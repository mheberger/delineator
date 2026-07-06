"""
Functions for managing delineator data files.
Specifically, this handles finding the correct data files for a given megabasin.

Basin 27 (Iceland) is bundled with the package, and other data files are downloaded
from the author's website as needed.

Data files are stored in a system default data directory, which is platform-dependent.

Users may also set a custom data directory via the environment variable
DELINEATOR_DATA_DIR.

"""

import os
import logging
from importlib.resources import files
from pathlib import Path
from platformdirs import user_data_dir
from delineator.download import _download_file, _local_path
from delineator.settings import DelineatorConfig
from delineator.constants import DATA_DIR_NAME

logger = logging.getLogger(__name__)

_STALE_WARNED = False


def _warn_stale_data(config: DelineatorConfig) -> None:
    """
    Warn once if data from an older package version is present.

    The current data lives in a version-named subdirectory of the data
    directory (see DATA_DIR_NAME, currently 'v2.2'). When the data files change, that name
    is bumped, so the corrected files download to a fresh folder and the
    buggy ones become unreachable. The old folders are not deleted
    automatically; this function logs a one-time warning listing them so
    the user can remove them to reclaim disk space.

    Detects both the legacy un-versioned layout ("vector/", "raster/")
    and any prior "delineator*" data folder. The warning fires at most
    once per process, guarded by the module-level _STALE_WARNED flag.

    Parameters
    ----------
    config : DelineatorConfig
        Configuration object; config.data_dir is scanned for stale folders.

    Returns
    -------
    None
    """

    global _STALE_WARNED
    if _STALE_WARNED:
        return

    # The data directory may not exist yet (e.g. a fresh install, before any
    # data has been downloaded). If so there are no stale folders to report.
    # Skip without consuming the one-time warning, so a later call can still
    # scan once the directory has been created.
    if not config.data_dir.is_dir():
        return

    _STALE_WARNED = True

    stale = [
        child for child in config.data_dir.iterdir()
        if child.is_dir()
        and child.name != DATA_DIR_NAME
        and (child.name in ("vector", "raster") or child.name.startswith("delineator"))
    ]
    if stale:
        listing = "\n".join(f"  {p}" for p in stale)
        txt = f"""NOTICE: Found delineator data from an older version. The current version
uses corrected data in: 
  {config.data_dir}

The old folders below are no longer used and can be deleted to reclaim disk space:
{listing}
"""
        print(txt)


def _find_data_file(relative_path: str, config: DelineatorConfig) -> Path | None:
    """
    Locate a data file in the data directory.
    Includes special handling for the data files for Iceland, megabasin 27,
    for which the data is bundled with the package.
    Otherwise, the script will have already set the data directory, either
    to a default directory or to the custom location specified by the user.

    Parameters
    ----------
    relative_path : str
        A path relative to the data directory, e.g. "rivers23.db" or "accum73.tif"
    config : DelineatorConfig dataclass object
        includes config.data_dir, the path to the data directory

    Returns
    --------
    strthe full path to the data file on the user's computer
    """
    if not config.custom_data_dir:

        _warn_stale_data(config)  # reclaims old-version disk; no-op after first call

    if '27' in relative_path:
        filepath = files('delineator').joinpath('data', relative_path)
    else:
        filepath = _local_path(relative_path, config)

    if filepath.is_file():
        return filepath

    # If the file is not found in the data directory, try to download it
    filepath = _download_file(relative_path, config)

    if filepath is not None:
        return filepath

    else:
        logger.warning(
            f"Data file not found: {relative_path}\n"
            f"and could not be downloaded. \n"
            f"Run `delineator_download --basin ##`'` to fetch required data files, \n"
            f"or visit https://mghydro.com/watersheds/data/ to download them manually, \n"
            f"or check your data directory: {config.data_dir}"
        )
        return None


def _get_data_dir() -> Path:
    """
    Return the *base* data directory for delineator data files.

    Note this is the base directory, not where the data files actually
    live. Data files are stored under a version-named subfolder of this
    directory (see DATA_DIR_NAME, e.g. "v2.2"); that subfolder is
    appended downstream in _local_path(), not here. The layout looks like:

        <base_data_dir>/            <- this function returns this
        ├── v2.2/          <- DATA_DIR_NAME: current data files
        │    ├── acccum73.tif      <- raster data files
        │    ├── basins73.db       <- vector data files
        ├── v2.1/               <- stale; flagged by _warn_stale_data()
        ├── vector/             <- legacy un-versioned; also flagged
        └── raster/

    Keeping this function version-agnostic is deliberate: _warn_stale_data()
    scans this base directory for the *siblings* of the current version
    folder to find data to clean up. If the version were baked in here, that
    scan would look inside the current folder and never see the old ones.

    The default location is the platform user-data dir (via platformdirs).
    Override it by setting the DELINEATOR_DATA_DIR environment variable:
        Windows:      set DELINEATOR_DATA_DIR=D:\\GIS\\delineator_data
        macOS/Linux:  export DELINEATOR_DATA_DIR=~/gis/delineator_data

    Returns
    -------
    Path
        The base data directory, created (including parents) if it does not
        already exist.
    """
    custom_path = os.environ.get("DELINEATOR_DATA_DIR", "").strip("")
    if custom_path:
        data_dir = Path(custom_path).expanduser()
    else:
        data_dir = Path(user_data_dir("delineator", appauthor=False))

    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir
