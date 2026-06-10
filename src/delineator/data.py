from importlib.resources import files
import os
import logging
from pathlib import Path
from platformdirs import user_data_dir

from delineator.download import _download_file
from delineator.settings import DelineatorConfig

logger = logging.getLogger(__name__)

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
        filepath = config.data_dir / relative_path

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
