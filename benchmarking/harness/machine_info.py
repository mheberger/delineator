"""
Collect the machine, software, and data-file details that a benchmark report
needs, so every results folder documents the conditions it was measured in.

Deliberately leaves out the host name and user name. Every lookup is
best-effort: a failure records None rather than stopping the benchmark.
"""

import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import winperf

PACKAGES = ["delineator", "pysheds", "numpy", "numba", "shapely", "geopandas", "pyogrio",
            "rasterio", "pyproj", "pandas", "scipy", "pyshp"]


def _powershell_json(command: str, pipe: bool = True):
    """Run *command* and parse its JSON output; pipe=False if the command calls ConvertTo-Json itself."""
    try:
        out = subprocess.run(["powershell", "-NoProfile", "-NonInteractive", "-Command",
                              command + (" | ConvertTo-Json -Depth 3" if pipe else "")],
                             capture_output=True, text=True, timeout=60).stdout
        return json.loads(out) if out.strip() else None
    except Exception:
        return None


def _is_admin() -> bool:
    try:
        import ctypes
        return bool(ctypes.windll.shell32.IsUserAnAdmin())
    except Exception:
        return False


def _defender(paths) -> dict:
    """
    Antivirus state. On-access scanning can slow reads of files just written
    or changed, so a benchmark report should say whether the data folders
    were excluded. Windows shows the exclusion list only to administrators;
    otherwise ``exclusions_readable`` is False and the report has to state
    the exclusions by hand.
    """
    status = _powershell_json("Get-MpComputerStatus | Select-Object AMRunningMode, "
                              "RealTimeProtectionEnabled, OnAccessProtectionEnabled, AMProductVersion") or {}
    products = _powershell_json("Get-CimInstance -Namespace root/SecurityCenter2 "
                                "-ClassName AntivirusProduct | Select-Object displayName")
    if isinstance(products, dict):
        products = [products]
    info = {"antivirus_products": [p.get("displayName") for p in products or []],
            "running_as_admin": _is_admin(), **status}

    exclusions = _powershell_json("ConvertTo-Json -InputObject @((Get-MpPreference).ExclusionPath)",
                                  pipe=False)
    if isinstance(exclusions, str):
        exclusions = [exclusions]
    readable = isinstance(exclusions, list) and not any(
        str(e).startswith("N/A") for e in exclusions)
    info["exclusions_readable"] = readable
    folders = sorted({str(Path(p) if Path(p).is_dir() else Path(p).parent) for p in paths})
    if readable:
        info["exclusion_paths"] = exclusions
        prefixes = [str(e).rstrip("\\").lower() for e in exclusions]
        # Plain prefix matching; wildcard or environment-variable exclusions are not expanded.
        info["data_folders_excluded"] = {
            f: any(f.lower() == e or f.lower().startswith(e + "\\") for e in prefixes) for f in folders}
    else:
        info["data_folders_excluded"] = {f: None for f in folders}
    return info


def _cpu() -> dict:
    info = {"logical_processors": os.cpu_count()}
    try:
        import winreg
        with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                            r"HARDWARE\DESCRIPTION\System\CentralProcessor\0") as key:
            info["name"] = winreg.QueryValueEx(key, "ProcessorNameString")[0].strip()
    except OSError:
        pass
    cim = _powershell_json("Get-CimInstance Win32_Processor | "
                           "Select-Object NumberOfCores, MaxClockSpeed")
    if isinstance(cim, list):
        cim = cim[0]
    if cim:
        info["physical_cores"] = cim.get("NumberOfCores")
        info["max_clock_mhz"] = cim.get("MaxClockSpeed")
    return info


def _disk_for(path: str):
    """The physical disk and file system behind *path*'s drive letter."""
    drive = Path(path).resolve().drive.rstrip(":")
    if not drive:
        return None
    return _powershell_json(
        f"$p = Get-Partition -DriveLetter {drive}; $v = Get-Volume -DriveLetter {drive}; "
        f"$d = Get-PhysicalDisk | Where-Object DeviceId -eq ([string]$p.DiskNumber); "
        f"[pscustomobject]@{{Drive='{drive}'; FileSystem=$v.FileSystem; "
        f"Disk=$d.FriendlyName; MediaType=[string]$d.MediaType; BusType=[string]$d.BusType; "
        f"SizeGB=[math]::Round($d.Size/1GB)}}")


def _file_version(exe: str | None):
    if not exe:
        return None
    return _powershell_json(f"(Get-Item '{exe}').VersionInfo | "
                            f"Select-Object FileVersion, ProductVersion, ProductName")


_POWER_MODES = {"00000000-0000-0000-0000-000000000000": "Balanced",
                "ded574b5-45a0-4f42-8737-46345c09c238": "Best performance",
                "961cc777-2547-4f9d-8174-7d86181b8a7a": "Best power efficiency"}


def _power_scheme() -> dict:
    """The power plan, plus the Windows 10/11 power-mode slider, which overrides it."""
    info = {}
    try:
        info["plan"] = subprocess.run(["powercfg", "/getactivescheme"], capture_output=True,
                                      text=True, timeout=30).stdout.strip()
    except Exception:
        info["plan"] = None
    try:
        import winreg
        with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE,
                            r"SYSTEM\CurrentControlSet\Control\Power\User\PowerSchemes") as key:
            for name in ("ActiveOverlayAcPowerScheme", "ActiveOverlayDcPowerScheme"):
                try:
                    guid = winreg.QueryValueEx(key, name)[0].lower()
                    info[name] = _POWER_MODES.get(guid, guid)
                except OSError:
                    pass
    except OSError:
        pass
    return info


def _git(repo: Path) -> dict:
    def run(*args):
        try:
            return subprocess.run(["git", "-C", str(repo), *args], capture_output=True,
                                  text=True, timeout=30).stdout.strip()
        except Exception:
            return None
    return {"commit": run("rev-parse", "HEAD"),
            "describe": run("describe", "--tags", "--always", "--dirty"),
            "dirty_files": run("status", "--porcelain", "--", "src")}


def _versions() -> dict:
    versions = {}
    for name in PACKAGES:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    try:
        import rasterio
        versions["GDAL"] = rasterio.__gdal_version__
    except ImportError:
        pass
    return versions


def file_inventory(paths) -> list[dict]:
    rows = []
    for p in sorted(set(map(str, paths))):
        try:
            st = os.stat(p)
            rows.append({"path": p, "bytes": st.st_size,
                         "modified": datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat()})
        except OSError:
            rows.append({"path": p, "bytes": None, "modified": None})
    return rows


def collect(data_paths, repo: Path, mpiexec: str = "mpiexec",
            gagewatershed: str = "GageWatershed") -> dict:
    """Everything except the data inventory, which the caller adds."""
    drives = {}
    for p in data_paths:
        drive = Path(p).resolve().drive
        if drive and drive not in drives:
            drives[drive] = _disk_for(p)
    mpi_exe, gw_exe = shutil.which(mpiexec), shutil.which(gagewatershed)
    return {
        "collected_utc": datetime.now(timezone.utc).isoformat(),
        "os": platform.platform(),
        "python": sys.version,
        "python_executable": sys.executable,
        "cpu": _cpu(),
        "memory_mb": winperf.memory_status_mb(),
        "on_ac_power": winperf.ac_power(),
        "power_scheme": _power_scheme(),
        "storage": drives,
        "antivirus": _defender(data_paths),
        "packages": _versions(),
        "taudem": {"gagewatershed": gw_exe, "version": _file_version(gw_exe)},
        "mpi": {"mpiexec": mpi_exe, "version": _file_version(mpi_exe)},
        "delineator_repo": _git(repo),
    }
