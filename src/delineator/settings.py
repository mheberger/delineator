# src/delineator/settings.py
import os
from dataclasses import dataclass, field
from pathlib import Path
from platformdirs import user_data_dir


def _default_data_dir() -> Path:
    """Return the default data directory, respecting DELINEATOR_DATADIR env var."""
    custom = os.environ.get("DELINEATOR_DATADIR")
    if custom:
        return Path(custom).expanduser()
    return Path(user_data_dir("delineator", appauthor=False))


def _default_output_dir() -> Path:
    custom = os.environ.get("DELINEATOR_OUTDIR")
    if custom:
        return Path(custom).expanduser()
    output_dir = Path.cwd() / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


@dataclass
class DelineatorConfig:
    """
    Configuration settings for watershed delineation.

    Parameters
    ----------
    auto_download: bool
        If True, the script will attempt to download needed data files for
        watershed delineation. You can also manually download the data files,
        see documentation for download_data().

    clean : bool
        If True, the watershed boundary polygon will be "cleaned"
        which repairs artifacts such as seams that affect
        the appearance of the watershed polygon.

    data_dir : Path or str, optional
        Directory where the script will look for data files and store downloaded
        data files needed for watershed delineation. If not provided,
        the script will use the default cache directory. For example, on Windows:
        C:\\Users\\<username>\\.delineator\\data
        Or on macOS/Linux: ~/.delineator/data

    high_res : bool
        If True, the script will split the unit catchment around the outlet point
        for greater accuracy.
        If False, the script skips this step. The result will be a watershed that is
        slightly too large, and will include some area downstream of the outlet point.
        This error may be insignificant for very large watersheds.

    fill : bool
        MERIT-Hydro data tends to produce watersheds with small "donut holes" due
        to minor topological errors in the input data. If True, these holes will be
        filled in, resulting in a cleaner appearance and smaller output files.

    fill_threshold : int
        Only meaningful if fill=True. Holes smaller than this number of pixels
        (on the 3 arcsecond grid) will be filled. Larger holes are preserved.
        Set to 0 to fill all holes regardless of size. Default is 100 pixels.

    low_res_threshold : float
        Watershed area in km² above which the script will automatically switch to
        lower-resolution mode, regardless of the high_res setting.
        Default is 80,000 km², meaning that the script will automatically switch to
        lower-resolution mode for watersheds over 80,000 km². Set it to a large
        number (6e6 or higher) to disable lower-resolution mode. (The largest
        basin in the dataset is the Amazon at 5.9 × 10⁶ km².)

    num_stream_orders : int
        Number of stream orders to include in the river network. This will only
        have an effect if rivers=True.
        Default is 4. Higher values will result in more detailed stream network.
        Set this to a value of 9 or greater to get every river reach available.

    rivers: bool
        If True, return the river flow lines in the watershed. These will be saved
        as the "rivers" layer in the output GeoPackage, or to a separate file
        such as "rivers.shp" or "rivers.geojson", depending on  `output_format`.
    
    output_format: str
        Specify the output format for geodata files. Options are "gpkg" (default)
        or "geojson", "shp", "kml", and others supported by your installation of geopandas.
        See https://geopandas.org/en/stable/docs/user_guide/io.html#supported-drivers-file-formats
        for more details.

    search_dist : float
        If the outlet point does not fall inside any catchment, the script will
        search for the nearest catchment within this distance (in decimal degrees).
        Default is 0.025°, which is roughly 2-3 km near the equator.

    """

    auto_download: bool = True
    clean: bool = False
    data_dir: Path = field(default_factory=_default_data_dir)
    fill: bool = True
    fill_threshold: int = 100
    high_res: bool = True
    low_res_threshold: float = 6e6
    num_stream_orders: int = 4
    output_dir: Path = field(default_factory=_default_output_dir)
    output_format: str = "gpkg"
    outlets: bool = True
    rivers: bool = True
    search_dist: float = 0.025

