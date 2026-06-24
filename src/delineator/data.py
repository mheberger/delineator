from importlib.resources import files
import os
import logging
from pathlib import Path
from platformdirs import user_data_dir

from delineator.download import _download_file, _local_path
from delineator.settings import DelineatorConfig
import shutil
from delineator.constants import DATA_VERSION

logger = logging.getLogger(__name__)

_MIGRATION_DONE = False          # in-process guard
_LEGACY_DIRS = ("vector", "raster")
_MARKER_NAME = ".data_version"


def _migrate_data_dir(config) -> None:
    """Purge stale data when DATA_VERSION has changed. Idempotent & cheap."""
    global _MIGRATION_DONE
    if _MIGRATION_DONE:
        return

    base = config.data_dir       # already mkdir'd by _get_data_dir()
    marker = base / _MARKER_NAME

    try:
        current = marker.read_text().strip()
    except (FileNotFoundError, OSError):
        current = None

    if current == DATA_VERSION:
        _MIGRATION_DONE = True
        return

    for child in base.iterdir():
        if child.name == _MARKER_NAME:
            continue
        if child.name in _LEGACY_DIRS or (child.is_dir() and child.name != DATA_VERSION):
            logger.info("Removing stale delineator data: %s", child)
            shutil.rmtree(child, ignore_errors=True)

    try:
        marker.write_text(DATA_VERSION)
    except OSError as exc:
        logger.warning("Could not write data version marker %s: %s", marker, exc)

    _MIGRATION_DONE = True


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
        A path relative to the data directory, e.g. "vector/rivers23.db" or "raster/accum73.tif"
    config : DelineatorConfig dataclass object
        includes config.data_dir, the path to the data directory

    Returns
    -------
    the full path to the data file on the user's computer
    """
    #
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
            f"and could not be downloaded."
            f"Run 'delineator_download --basin ##' to fetch required data files, "
            f"or check your data directory: {config.data_dir}"
        )
        return None


def _get_data_dir() -> Path:
    """
    Return the data directory for delineator data files.
    
    Override the default location by setting the DELINEATOR_DATA environment variable:
        Windows:  set DELINEATOR_DATA=D:\\GIS\\delineator_data
        macOS/Linux: export DELINEATOR_DATA=~/gis/delineator_data
    """
    custom_path = os.environ.get("DELINEATOR_DATA")
    if custom_path:
        data_dir = Path(custom_path).expanduser()
    else:
        data_dir = Path(user_data_dir("delineator", appauthor=False))

    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir
