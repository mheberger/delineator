"""
Configuration settings for the `delineator` package

The `DelineatorConfig` is used by several functions in the `delineator` package
including:

  - `delineator.core.delineate_watershed`
  - `delineator.core.downloder`
  - `delineator.util.write_outputs`
  - `delineator.cli.delineate_dataframe`

Allows the user to specify various options for the watershed delineation process
such as the output format, whether to include rivers, directory locations, etc.
"""
import os
from dataclasses import dataclass, field
from pathlib import Path
from platformdirs import user_data_dir

from delineator.constants import SUPPORTED_OUTPUT_FORMATS


def _default_auto_download() -> bool:
    custom = os.environ.get("DELINEATOR_AUTO_DOWNLOAD")
    if custom:
        return custom.lower() == "true"
    return True


def _default_data_dir() -> Path:
    """Return the default data directory, respecting DELINEATOR_DATADIR env var."""
    custom = os.environ.get("DELINEATOR_DATA_DIR")
    if custom:
        return Path(custom).expanduser()
    return Path(user_data_dir("delineator", appauthor=False))


def _default_output_dir() -> Path:
    custom = os.environ.get("DELINEATOR_OUTPUT_DIR")
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

    calc_area : bool
        If True, the script will calculate the area of the watershed.
        Setting it to False will save a small amount of time.

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
        (on the 3 arcsecond grid) will be filled. For example, enter 100 to fill
        holes smaller than 100 pixels in size. Holes larger than this will
        be preserved. Set to 0 to fill all holes. Default is 0.

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

    simplify : bool
        Used by the `write_outputs` utility function to simplify the output
        before saving it to disk. If True, the output will be simplified using
        Douglas-Peucker algorithm. This will result in a smaller file size,
        and also less of a jagged "staircase" appearance in the output.

    simplify_tolerance: bool
        If `simplify=True`, this is the tolerance value for the simplification.
        With the Douglas-Peucker algorithm, the tolerance value is the maximum
        distance (in decimal degrees) that a vertex can be from the original.

    threshold_single : int
        Threshold for number of upstream pixels that defines a stream for
        purposes of snapping the outlet to a river centerline. The script
        uses the `threshold_single` value when the outlet falls within a
        unit catchment with no upstream contributing catchments. A default
        value of 300 - 500 worked well in my testing.

    threshold_multi : int
        Threshold for number of upstream pixels that defines a stream for
        snapping the outlet to a river centerline. The script uses the
        `threshold_multi` value when the outlet falls within a unit catchment
        that has upstream contributing catchments. A default value of 5000 or
        higher will typically force the outlet to be snapped to a larger river,
        i.e. one with a larger upstream area.

    verbose: bool
        Set verbose=True to see messages about the progress of the script.
        This is equivalent to setting the loggin level to "INFO" using
        the `logging` module.

    """

    auto_download: bool = field(default_factory=_default_auto_download)
    calc_area: bool = True
    clean: bool = False
    data_dir: Path | str = field(default_factory=_default_data_dir)
    fill: bool = True
    fill_threshold: int = 0
    high_res: bool = True
    low_res_threshold: float = 6e6
    num_stream_orders: int = 4
    output_dir: Path | str = field(default_factory=_default_output_dir)
    output_format: str = "gpkg"
    outlets: bool = True
    rivers: bool = True
    search_dist: float = 0.1
    simplify: bool = False
    simplify_tolerance: float = 0.0008
    threshold_single: int = 300
    threshold_multi: int = 5000
    user_id: bool = False
    verbose: bool = False

    def __post_init__(self) -> None:
        self.data_dir = self._coerce_path(self.data_dir, "data_dir")
        self.output_dir = self._coerce_path(self.output_dir, "output_dir")

        self.output_format = self.output_format.lower().lstrip(".")

        self._validate_bool("auto_download", self.auto_download)
        self._validate_bool("clean", self.clean)
        self._validate_bool("fill", self.fill)
        self._validate_bool("high_res", self.high_res)
        self._validate_bool("outlets", self.outlets)
        self._validate_bool("rivers", self.rivers)
        self._validate_bool("simplify", self.simplify)

        if not isinstance(self.fill_threshold, int):
            raise TypeError("fill_threshold must be an int.")
        if self.fill_threshold < 0:
            raise ValueError("fill_threshold must be greater than or equal to 0.")

        if not isinstance(self.num_stream_orders, int):
            raise TypeError("num_stream_orders must be an int.")
        if self.num_stream_orders < 1:
            raise ValueError("num_stream_orders must be greater than or equal to 1.")

        self.low_res_threshold = self._coerce_float(
            self.low_res_threshold,
            "low_res_threshold",
        )
        if self.low_res_threshold <= 0:
            raise ValueError("low_res_threshold must be greater than 0.")

        self.search_dist = self._coerce_float(self.search_dist, "search_dist")
        if self.search_dist <= 0:
            raise ValueError("search_dist must be greater than 0.")

        self.simplify_tolerance = self._coerce_float(
            self.simplify_tolerance,
            "simplify_tolerance",
        )
        if self.simplify_tolerance < 0:
            raise ValueError("simplify_tolerance must be greater than or equal to 0.")

        if not isinstance(self.output_format, str):
            raise TypeError("output_format must be a string.")
        if self.output_format not in SUPPORTED_OUTPUT_FORMATS:
            valid_formats = ", ".join(sorted(SUPPORTED_OUTPUT_FORMATS))
            raise ValueError(
                f"Unsupported output_format: {self.output_format!r}. "
                f"Supported formats are: {valid_formats}."
            )

    @staticmethod
    def _coerce_path(value: Path | str, name: str) -> Path:
        if isinstance(value, Path):
            return value
        if isinstance(value, str):
            return Path(value)
        raise TypeError(f"{name} must be a pathlib.Path or string.")

    @staticmethod
    def _coerce_float(value: float | int, name: str) -> float:
        if isinstance(value, bool):
            raise TypeError(f"{name} must be a number, not bool.")
        if isinstance(value, int | float):
            return float(value)
        raise TypeError(f"{name} must be a number.")

    @staticmethod
    def _validate_bool(name: str, value: bool) -> None:
        if not isinstance(value, bool):
            raise TypeError(f"{name} must be a bool.")
