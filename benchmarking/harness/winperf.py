"""
Windows process accounting and file-cache control for benchmarking.

- JobAccounting puts the current process in a Windows job object, so every
  process it starts afterwards (e.g. mpiexec and the TauDEM ranks) is in the
  job too. The job's counters then give CPU time, bytes read, and peak
  committed memory for the whole process tree, measured the same way for
  every tool.
- purge_file_cache() evicts a file's pages from the Windows file cache
  without admin rights: opening a file with FILE_FLAG_NO_BUFFERING makes the
  cache manager flush and purge that file's cached pages. It fails silently
  if another process has the file memory-mapped, so verify_purged() checks
  the result with a few random reads.

Bytes read are *logical* reads (ReadFile calls, whether or not the data came
from the cache). Reads through memory-mapped files are not counted; SQLite
(mmap_size=0 by default) and GDAL do not use them here.
"""

import ctypes
import random
import sys
import time
from ctypes import wintypes

if sys.platform != "win32":
    raise ImportError("winperf works only on Windows")

_k32 = ctypes.WinDLL("kernel32", use_last_error=True)

_k32.CreateFileW.restype = wintypes.HANDLE
_k32.CreateFileW.argtypes = [wintypes.LPCWSTR, wintypes.DWORD, wintypes.DWORD, wintypes.LPVOID,
                             wintypes.DWORD, wintypes.DWORD, wintypes.HANDLE]
_k32.CloseHandle.argtypes = [wintypes.HANDLE]
_k32.CreateJobObjectW.restype = wintypes.HANDLE
_k32.CreateJobObjectW.argtypes = [wintypes.LPVOID, wintypes.LPCWSTR]
_k32.GetCurrentProcess.restype = wintypes.HANDLE
_k32.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
_k32.QueryInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, wintypes.LPVOID,
                                           wintypes.DWORD, wintypes.LPDWORD]
_k32.GetProcessIoCounters.argtypes = [wintypes.HANDLE, wintypes.LPVOID]
_k32.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.POINTER(wintypes.FILETIME)] * 4

_INVALID_HANDLE = wintypes.HANDLE(-1).value
_GENERIC_READ = 0x80000000
_SHARE_ALL = 0x7  # FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE
_OPEN_EXISTING = 3
_FILE_FLAG_NO_BUFFERING = 0x20000000


class _IoCounters(ctypes.Structure):
    _fields_ = [(name, ctypes.c_ulonglong) for name in (
        "ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
        "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]


class _BasicAccounting(ctypes.Structure):
    _fields_ = [("TotalUserTime", ctypes.c_longlong), ("TotalKernelTime", ctypes.c_longlong),
                ("ThisPeriodTotalUserTime", ctypes.c_longlong),
                ("ThisPeriodTotalKernelTime", ctypes.c_longlong),
                ("TotalPageFaultCount", wintypes.DWORD), ("TotalProcesses", wintypes.DWORD),
                ("ActiveProcesses", wintypes.DWORD), ("TotalTerminatedProcesses", wintypes.DWORD)]


class _BasicAndIoAccounting(ctypes.Structure):
    _fields_ = [("BasicInfo", _BasicAccounting), ("IoInfo", _IoCounters)]


class _BasicLimit(ctypes.Structure):
    _fields_ = [("PerProcessUserTimeLimit", ctypes.c_longlong), ("PerJobUserTimeLimit", ctypes.c_longlong),
                ("LimitFlags", wintypes.DWORD), ("MinimumWorkingSetSize", ctypes.c_size_t),
                ("MaximumWorkingSetSize", ctypes.c_size_t), ("ActiveProcessLimit", wintypes.DWORD),
                ("Affinity", ctypes.c_size_t), ("PriorityClass", wintypes.DWORD),
                ("SchedulingClass", wintypes.DWORD)]


class _ExtendedLimit(ctypes.Structure):
    _fields_ = [("BasicLimitInformation", _BasicLimit), ("IoInfo", _IoCounters),
                ("ProcessMemoryLimit", ctypes.c_size_t), ("JobMemoryLimit", ctypes.c_size_t),
                ("PeakProcessMemoryUsed", ctypes.c_size_t), ("PeakJobMemoryUsed", ctypes.c_size_t)]


_JobObjectBasicAndIoAccountingInformation = 8
_JobObjectExtendedLimitInformation = 9


def _filetime_s(ft: wintypes.FILETIME) -> float:
    return ((ft.dwHighDateTime << 32) | ft.dwLowDateTime) / 1e7


class JobAccounting:
    """
    CPU, I/O, and memory counters for this process and everything it starts.

    If the process cannot join a job, falls back to counters for this
    process alone and sets ``covers_children`` to False.
    """

    def __init__(self):
        self._job = _k32.CreateJobObjectW(None, None)
        self.covers_children = bool(self._job) and bool(
            _k32.AssignProcessToJobObject(self._job, _k32.GetCurrentProcess()))

    def snapshot(self) -> dict:
        """Cumulative counters: cpu_user_s, cpu_kernel_s, read_bytes, read_ops, write_bytes."""
        if self.covers_children:
            info = _BasicAndIoAccounting()
            if not _k32.QueryInformationJobObject(self._job, _JobObjectBasicAndIoAccountingInformation,
                                                  ctypes.byref(info), ctypes.sizeof(info), None):
                raise ctypes.WinError(ctypes.get_last_error())
            user, kernel = info.BasicInfo.TotalUserTime / 1e7, info.BasicInfo.TotalKernelTime / 1e7
            io = info.IoInfo
        else:
            io = _IoCounters()
            _k32.GetProcessIoCounters(_k32.GetCurrentProcess(), ctypes.byref(io))
            times = [wintypes.FILETIME() for _ in range(4)]
            _k32.GetProcessTimes(_k32.GetCurrentProcess(), *[ctypes.byref(t) for t in times])
            kernel, user = _filetime_s(times[2]), _filetime_s(times[3])
        return {"cpu_user_s": user, "cpu_kernel_s": kernel, "read_bytes": io.ReadTransferCount,
                "read_ops": io.ReadOperationCount, "write_bytes": io.WriteTransferCount}

    def peak_memory_mb(self) -> float | None:
        """Peak committed memory of the whole job since it was created (None without a job)."""
        if not self.covers_children:
            return None
        info = _ExtendedLimit()
        if not _k32.QueryInformationJobObject(self._job, _JobObjectExtendedLimitInformation,
                                              ctypes.byref(info), ctypes.sizeof(info), None):
            return None
        return info.PeakJobMemoryUsed / 2**20

    @staticmethod
    def delta(before: dict, after: dict) -> dict:
        return {k: after[k] - before[k] for k in before}


def purge_file_cache(path: str) -> None:
    """Evict *path* from the Windows file cache (see module docstring)."""
    handle = _k32.CreateFileW(str(path), _GENERIC_READ, _SHARE_ALL, None, _OPEN_EXISTING,
                              _FILE_FLAG_NO_BUFFERING, None)
    if handle == _INVALID_HANDLE:
        raise ctypes.WinError(ctypes.get_last_error())
    _k32.CloseHandle(handle)


def random_read_latency_us(path: str, n: int = 64, seed: int | None = None) -> float:
    """Mean latency (µs) of *n* random 4 KiB reads: microseconds if cached, far more if not."""
    import os
    pages = os.path.getsize(path) // 4096
    rng = random.Random(seed)
    offsets = [rng.randrange(pages) * 4096 for _ in range(n)]
    with open(path, "rb", buffering=0) as f:
        t = time.perf_counter()
        for off in offsets:
            f.seek(off)
            f.read(4096)
        return (time.perf_counter() - t) / n * 1e6


def calibrate(path: str) -> dict:
    """Measure warm and cold random-read latency on *path*; the threshold sits between them."""
    random_read_latency_us(path, n=256, seed=1)
    warm = random_read_latency_us(path, n=256, seed=1)
    purge_file_cache(path)
    cold = random_read_latency_us(path, n=256, seed=2)
    purge_file_cache(path)
    return {"file": str(path), "warm_us": warm, "cold_us": cold,
            "threshold_us": (warm * cold) ** 0.5, "ratio": cold / warm}


def purge_and_verify(path: str, threshold_us: float) -> dict:
    """Purge *path*, probe it, and purge it again. ok=False means it was probably still cached."""
    purge_file_cache(path)
    latency = random_read_latency_us(path, n=32)
    purge_file_cache(path)
    return {"file": str(path), "probe_us": latency, "ok": latency >= threshold_us}


def prevent_idle_sleep() -> bool:
    """
    Keep Windows from idle-sleeping while this process runs (the display may
    still turn off). Released automatically when the process exits; a
    closed lid or choosing Sleep still sleeps.
    """
    ES_CONTINUOUS, ES_SYSTEM_REQUIRED = 0x80000000, 0x00000001
    _k32.SetThreadExecutionState.restype = wintypes.DWORD
    _k32.SetThreadExecutionState.argtypes = [wintypes.DWORD]
    return bool(_k32.SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED))


def ac_power() -> bool | None:
    """True on mains power, False on battery, None if unknown."""
    class _Status(ctypes.Structure):
        _fields_ = [("ACLineStatus", ctypes.c_ubyte), ("BatteryFlag", ctypes.c_ubyte),
                    ("BatteryLifePercent", ctypes.c_ubyte), ("SystemStatusFlag", ctypes.c_ubyte),
                    ("BatteryLifeTime", wintypes.DWORD), ("BatteryFullLifeTime", wintypes.DWORD)]
    status = _Status()
    if not _k32.GetSystemPowerStatus(ctypes.byref(status)):
        return None
    return {0: False, 1: True}.get(status.ACLineStatus)


def memory_status_mb() -> dict:
    class _MemStatus(ctypes.Structure):
        _fields_ = [("dwLength", wintypes.DWORD), ("dwMemoryLoad", wintypes.DWORD),
                    ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
    status = _MemStatus()
    status.dwLength = ctypes.sizeof(status)
    _k32.GlobalMemoryStatusEx(ctypes.byref(status))
    return {"total_mb": status.ullTotalPhys / 2**20, "available_mb": status.ullAvailPhys / 2**20}
