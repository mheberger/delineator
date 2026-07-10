"""
Configuration settings for the `delineator` package

The `DelineatorConfig` is used by several functions in the `delineator` package
including:

  - `delineator.core.delineate`
  - `delineator.core.downloader`
  - `delineator.util.write_outputs`
  - `delineator.cli._delineate_dataframe`

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

    if custom is None:
        return True

    if custom.lower() == "false":
        return False
    if custom == "0":
        return False
    return True


def _resolve_default_data_dir() -> tuple[Path, bool]:
    """
    Return ``(data_dir, is_custom)``: the DELINEATOR_DATA_DIR environment
    variable if set to a non-blank value (custom, flat layout), otherwise
    the platform user-data directory (managed, versioned layout).

    This is the single source of truth for the environment-variable
    fallback; both DelineatorConfig and data._get_data_dir use it.
    """
    custom = os.environ.get("DELINEATOR_DATA_DIR", "").strip()
    if custom:
        return Path(custom).expanduser(), True
    return Path(user_data_dir("delineator", appauthor=False)), False


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

    Parameters:
    ----------
    auto_download: bool
        If True, the script will attempt to download needed data files for
        watershed delineation. You can also manually download the data files,
        see documentation for download_data().

    cache : bool
        If True, the script will cache the results of the dissolve operation
        for merging upstream unit catchment geometries. This will save time
        when delineating multiple watersheds whose outlets fall in the same
        unit catchment. The cache is in-memory, holds at most a few dozen
        watersheds, and can be emptied with `delineator.clear_cache()`.
        Default is False.

    calc_area : bool
        If True, the script will calculate the area of the watershed.
        Setting it to False will save a small amount of time.

    clean : bool
        If True, the watershed boundary polygon will be "cleaned"
        which repairs artifacts such as seams that affect
        the appearance of the watershed polygon. Default is True.

    data_dir : Path or str, optional
        Directory where the script will look for data files and store downloaded
        data files needed for watershed delineation. If not provided,
        the script will use the platform default data directory. For example:
        Windows: C:\\Users\\<username>\\AppData\\Local\\delineator
        Linux:   ~/.local/share/delineator
        macOS:   ~/Library/Application Support/delineator

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

    fill_area_max : float
        Only used when fill=True. Fill holes smaller than this value, in km².
        For example, enter 10 to fill holes  smaller than 10 km².
       **Set to 0 to fill all holes.**
       Default is 1.0 -- only holes smaller than 1.0 km² will be filled.

    num_stream_orders : int
        Number of stream orders to include in the river network. This will only
        have an effect if rivers=True.
        Default is 4. Higher values will result in more detailed stream network.
        Set this to a value of 9 or greater to get every river reach available.

    output_format: str
        Specify the output format for geodata files. Options are "gpkg" (default)
        or "geojson", "shp", "kml", and others supported by your installation of geopandas.
        See https://geopandas.org/en/stable/docs/user_guide/io.html#supported-drivers-file-formats
        for more details.

    rivers: bool
        If True, return the river flow lines in the watershed. These will be saved
        as the "rivers" layer in the output GeoPackage, or to a separate file
        such as "rivers.shp" or "rivers.geojson", depending on  `output_format`.

    round_coordinates: bool
        If True, the coordinates of the watershed and rivers will be rounded to
        6 decimal places. This will reduce the file size of the output files,
        with minimal impact on the quality of the output. Default is True.
        Rounds to 6 decimal places, which is about 10 cm precision near the equator.

    search_dist : float
        The distance (in km) to search for a catchment that contains the outlet point,
        or to snap the outlet point to a river centerline.
        If the outlet point does not fall inside any catchment. the script will
        search for the nearest catchment within this distance (in kilometers).
        Default is 5.0 km. Accepts values from 0 to 50 km.

    simplify : bool
        Used by the `write_outputs` utility function to simplify the output
        before saving it to disk. If True, the output will be simplified using
        Douglas-Peucker algorithm. This will result in a smaller file size,
        and also less of a jagged "staircase" appearance in the output.

    simplify_tolerance: float
        If `simplify=True`, this is the tolerance value for the simplification.
        With the Douglas-Peucker algorithm, the tolerance value is the maximum
        distance (in kilometers) that a vertex can be from the original.

    smooth: bool
        If True, smooth the output geometries for a more natural appearance:
        rivers with centripetal Catmull-Rom splines, and the watershed
        boundary with Chaikin corner cutting. Best used together with
        `simplify=True` -- smoothing is applied after simplification, and
        without it the smoothed line faithfully traces the raster
        "staircase" instead of removing it. Default is False.

    snapping: bool
        If snapping=False, the script will not attempt to snap the outlet point.
        Instead, it will use the coordinates provided by the user. This is
        provided as an alternative to the built-in pour-point snapping function
        for testing or for those who wish to implement a custom snapping routine.

    stream_threshold_km2: float
        Minimum upstream drainage area, in km², for a grid cell to be considered
        part of a stream when snapping the outlet. Cells below this are not valid
        snap targets. Applies to watersheds spanning multiple unit catchments;
        setting this too low causes incorrect results and topology problems.

    verbose: bool
        Set verbose=True to see messages about the progress of the script.
        This is equivalent to setting the loggin level to "INFO" using
        the `logging` module.

    """

    auto_download: bool = field(default_factory=_default_auto_download)
    cache: bool = False
    calc_area: bool = True
    clean: bool = True
    custom_data_dir: bool = field(init=False, default=False)
    data_dir: Path | str | None = None
    fill: bool = True
    fill_area_max: float = 1.0
    high_res: bool = True
    num_stream_orders: int = 4
    output_dir: Path | str = field(default_factory=_default_output_dir)
    output_format: str = "gpkg"
    outlets: bool = True
    rivers: bool = True
    round_coordinates: bool = True
    search_dist: float = 5.0
    simplify: bool = False
    simplify_tolerance: float = 0.10  # km; ≈ the old default of 0.0008 decimal degrees
    smooth: bool = False
    snapping: bool = True
    stream_threshold_km2: float = 25.0
    verbose: bool = False

    def __post_init__(self) -> None:
        self._resolve_data_dir()
        self.output_dir = self._coerce_path(self.output_dir, "output_dir")

        self.output_format = self.output_format.lower().lstrip(".")

        self._validate_bool("auto_download", self.auto_download)
        self._validate_bool("cache", self.cache)
        self._validate_bool("calc_area", self.calc_area)
        self._validate_bool("clean", self.clean)
        self._validate_bool("fill", self.fill)
        self._validate_bool("high_res", self.high_res)
        self._validate_bool("outlets", self.outlets)
        self._validate_bool("rivers", self.rivers)
        self._validate_bool("round_coordinates", self.round_coordinates)
        self._validate_bool("simplify", self.simplify)
        self._validate_bool("smooth", self.smooth)
        self._validate_bool("snapping", self.snapping)

        self.fill_area_max = self._coerce_float(
            self.fill_area_max,
            "fill_area_max",
        )
        if self.fill_area_max < 0:
            raise ValueError("fill_area_max must be greater than or equal to 0.")

        if not isinstance(self.num_stream_orders, int):
            raise TypeError("num_stream_orders must be an int.")
        if self.num_stream_orders < 1:
            raise ValueError("num_stream_orders must be greater than or equal to 1.")

        self.search_dist = self._coerce_float(self.search_dist, "search_dist")
        if self.search_dist <= 0:
            raise ValueError("search_dist (in km) must be greater than 0.")
        if self.search_dist > 50:
            raise ValueError("search_dist (in km) must be less than or equal to 50.")

        self.simplify_tolerance = self._coerce_float(
            self.simplify_tolerance,
            "simplify_tolerance",
        )
        if self.simplify_tolerance < 0:
            raise ValueError("simplify_tolerance must be greater than or equal to 0.")

        self.stream_threshold_km2 = self._coerce_float(
            self.stream_threshold_km2,
            "stream_threshold_km2",
        )
        if self.stream_threshold_km2 <= 0:
            raise ValueError("stream_threshold_km2 must be greater than 0.")
        if self.stream_threshold_km2 > 50000:
            raise ValueError("stream_threshold_km2 must be less than or equal to 50000.")

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

    def _resolve_data_dir(self) -> None:
        """
        Resolve data_dir and record whether it is user-supplied ("custom").

        Precedence:
          1. An explicit data_dir argument           -> custom (flat layout)
          2. The DELINEATOR_DATA_DIR environment var  -> custom (flat layout)
          3. The platform user-data directory         -> managed default (versioned)

        For a custom directory the data files are read and written directly in
        that directory. For the managed default they live in a version-named
        subfolder (DATA_DIR_NAME) so data updates land in a fresh folder; see
        download._local_path and data._warn_stale_data.
        """
        if self.data_dir is not None:
            self.custom_data_dir = True
            self.data_dir = self._coerce_path(self.data_dir, "data_dir")
            self.data_dir.mkdir(parents=True, exist_ok=True)
            return

        self.data_dir, self.custom_data_dir = _resolve_default_data_dir()
        self.data_dir.mkdir(parents=True, exist_ok=True)
        