"""
delineator 2.0: Global Watershed Delineation with Python and data from
MERIT-Hydro and MERIT-Basins

Matthew Heberger, major refactor in May 2026

Formerly a collection of Python scripts, rewritten as a Python package

"""

import logging
from collections import OrderedDict
from contextlib import closing
from pathlib import Path
from typing import Any
import geopandas as gpd
import shapely
from shapely.ops import unary_union
from shapely.geometry import Point
from shapely.geometry.base import BaseGeometry
import sqlite3
import warnings

from delineator.constants import (
    MEGABASINS_DB_FILE,
    VALID_MEGABASINS,
    ROUNDING_DECIMALS,
    KM_PER_DEGREE,
)
from delineator.data import _find_data_file
from delineator.util import _validate_outlet_coordinates, _close_holes, _get_polygon_area
from delineator.spatial import _point_in_polygon_analysis
from delineator.queries import get_upstream_comids, get_home_unit_catchment_geom_and_area, get_hiearchical_basins, \
    get_upstream_geometries, get_rivers, get_upstream_area
from delineator.merit_detailed import _split_catchment
from delineator.settings import DelineatorConfig
from delineator.smoothing import _smooth_river_geometry, _smooth_watershed_geometry

logger = logging.getLogger(__name__)

# Cache used when config.cache=True. Unit catchment comids are globally
# unique (the megabasin is encoded in the first two digits), and the home
# comid fully determines the set of upstream catchments, so the comid alone
# is a sufficient cache key. The cache is bounded because dissolved
# geometries can be large; it assumes one set of data files per process
# (call clear_cache() if you switch data directories mid-session).
_CACHE: OrderedDict[tuple[str, int], Any] = OrderedDict()
_CACHE_MAXSIZE = 64


def clear_cache() -> None:
    """Empty the upstream-catchment cache (see ``DelineatorConfig.cache``)."""
    _CACHE.clear()


def _cache_get(key: tuple[str, int]) -> Any | None:
    """Return the cached value for *key* and mark it most recently used, or None on a miss."""
    value = _CACHE.get(key)
    if value is not None:
        _CACHE.move_to_end(key)
    return value


def _cache_put(key: tuple[str, int], value: Any) -> None:
    """Store *value* under *key*, evicting the least recently used entries past ``_CACHE_MAXSIZE``."""
    _CACHE[key] = value
    while len(_CACHE) > _CACHE_MAXSIZE:
        _CACHE.popitem(last=False)


def _get_upstream_catchments(basins_db_conn: sqlite3.Connection, home_unit_catchment: int,
                             use_cache: bool) -> tuple[list[int], dict | None]:
    """
    Return ``(upstream_unit_catchments, basins_dict)`` for the home unit
    catchment: the recursive-CTE traversal of the drainage network and the
    hierarchical grouping of the result, optionally cached.
    """
    key = ('upstream', int(home_unit_catchment))
    if use_cache:
        cached = _cache_get(key)
        if cached is not None:
            logger.info(f"Using cached upstream catchments for unit catchment "
                        f"{home_unit_catchment}")
            return cached

    upstream_unit_catchments = get_upstream_comids(basins_db_conn, int(home_unit_catchment))
    basins_dict = None
    if len(upstream_unit_catchments) > 1:
        basins_dict = get_hiearchical_basins(basins_db_conn, upstream_unit_catchments[1:])

    if use_cache:
        _cache_put(key, (upstream_unit_catchments, basins_dict))

    return upstream_unit_catchments, basins_dict


def _get_upstream_union(basins_db_path: Path, home_unit_catchment: int, basins_dict: dict,
                        use_cache: bool) -> BaseGeometry:
    """
    Fetch the upstream catchment geometries and dissolve them into a single
    (Multi)Polygon, optionally caching the result. This is the most expensive
    vector operation in the pipeline, and it depends only on the home unit
    catchment, so repeated delineations with outlets in the same unit
    catchment can reuse it. The union never includes the home catchment
    itself; the caller merges in the home geometry (split or whole).
    """
    key = ('union', int(home_unit_catchment))
    if use_cache:
        cached = _cache_get(key)
        if cached is not None:
            logger.info(f"Using cached upstream geometry for unit catchment "
                        f"{home_unit_catchment}")
            return cached

    upstream_polygons = get_upstream_geometries(basins_db_path, basins_dict)
    union = unary_union(upstream_polygons)

    if use_cache:
        _cache_put(key, union)

    return union


def delineate(lat: float, lng: float, config: DelineatorConfig | None = None) -> tuple[
    gpd.GeoDataFrame | None,  # watershed
    gpd.GeoDataFrame | None,  # rivers
    gpd.GeoDataFrame | None,  # outlets
]:
    """
    Delineate the watershed upstream of an outlet point.

    Parameters
    ----------
    lat, lng : float
        Coordinates of the outlet point in decimal degrees (EPSG:4326).
        Note the order: latitude first.
    config : DelineatorConfig, optional
        Configuration object. If not provided, one is created with default
        settings. The object is not modified, so it can be reused across calls.

    Returns
    -------
    watershed : GeoDataFrame or None
        A single-row GeoDataFrame with the watershed boundary polygon, in
        EPSG:4326. In high-res mode it has ``outlet_lat`` and ``outlet_lng``
        columns with the snapped outlet. When ``config.calc_area`` is True it
        has an ``area_km2`` column with the watershed's area, and when
        ``config.fill`` is also True, a ``filled_area`` column with the part of
        that area (km²) added by filling holes.
        None if the watershed could not be delineated.
    rivers : GeoDataFrame or None
        The river reaches within the watershed. None if ``config.rivers`` is
        False, or if there are no rivers (e.g. in coastal catchments).
    outlets : GeoDataFrame or None
        Point features for the requested outlet and, in high-res mode, the
        snapped outlet, distinguished by the ``type`` column (``"requested"``
        or ``"snapped"``). None if ``config.outlets`` is False.

    Notes
    -----
    This function does not raise on failure. If the watershed cannot be
    delineated (invalid coordinates, an outlet over the ocean, missing data
    files, a failed snap), it logs a warning and returns ``(None, None, None)``.

    Missing data files for the outlet's megabasin are downloaded on first use
    unless ``config.auto_download`` is False.

    Examples
    --------
    >>> from delineator import delineate, DelineatorConfig
    >>> watershed, rivers, outlets = delineate(64.71072, -21.60337)
    >>> watershed, _, _ = delineate(64.71072, -21.60337,
    ...                             DelineatorConfig(rivers=False, outlets=False))
    """
    if config is None:
        config = DelineatorConfig()

    if config.verbose:
        pkg_logger = logging.getLogger("delineator")
        pkg_logger.setLevel(logging.INFO)
        if not pkg_logger.handlers:
            handler = logging.StreamHandler()
            handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
            pkg_logger.addHandler(handler)

    # Step 1: Validate the input coordinates (before opening any resources)
    try:
        lat, lng = _validate_outlet_coordinates(lat, lng)
    except Exception as e:
        logger.warning(f"ERROR: Invalid coordinates: {e}")
        return None, None, None

    print(f"Delineating watershed for outlet: {lat}, {lng}")

    # Step 2: Determine the megabasin
    requested_outlet = Point(lng, lat)
    with closing(sqlite3.connect(MEGABASINS_DB_FILE)) as megabasins_db_conn:
        megabasin = _point_in_polygon_analysis(megabasins_db_conn, requested_outlet, table='megabasins',
                                               geom_col='geometry', id_col='basin',
                                               search_dist_km=config.search_dist)

    if megabasin is None or megabasin not in VALID_MEGABASINS:
        logger.warning(
            "The requested outlet is not in any megabasin. Check that it is over land. Maybe you flipped the lat/lng?")
        return None, None, None

    logger.info(f"Your watershed is in megabasin {megabasin}")

    # Step 3: Determine the unit catchment of the outlet
    basins_db_file = f'basins{megabasin}.db'
    basins_db_path = _find_data_file(basins_db_file, config)
    if basins_db_path is None:
        logger.warning(f"Could not load unit catchment database: {basins_db_file}")
        return None, None, None
    # The connection is scoped to the queries that need it, so the database
    # file is closed (and unlocked) before the heavier geometry work begins.
    with closing(sqlite3.connect(basins_db_path)) as basins_db_conn:
        home_unit_catchment = _point_in_polygon_analysis(basins_db_conn, requested_outlet, table='l0_basins',
                                                         id_col='comid',
                                                         search_dist_km=config.search_dist)

        # If the outlet is not in any unit catchment, return None
        if home_unit_catchment is None:
            logger.warning("The requested outlet is not in any unit catchment. Check that it is over land.")
            return None, None, None

        logger.info(f"Your watershed is in unit catchment {home_unit_catchment}")

        # Step 4: Get the collection of upstream catchments while the
        # connection is open: the comid list, and for a non-singleton the
        # hierarchical nested set of catchments used in Step 6.
        upstream_unit_catchments, basins_dict = _get_upstream_catchments(
            basins_db_conn, home_unit_catchment, config.cache)
        is_singleton = len(upstream_unit_catchments) == 1

    # Get the upstream area of the home unit catchment based on the rivers database
    upstream_area = get_upstream_area(int(home_unit_catchment), config)

    # Step 5: Split the home unit catchment at the outlet
    # First, get the home unit catchment geometry as a shapely Polygon
    home_unit_catchment_polygon, home_catchment_area = get_home_unit_catchment_geom_and_area(basins_db_path,
                                                                                             home_unit_catchment)
    if config.high_res:
        # Perform the split - returns the portion of the home catchment upstream of the outlet using raster methods.
        # Also reports whether the outlet is on the home catchment's river reach. If it is not (it is on a
        # small tributary within the home catchment), none of the upstream unit catchments drain to it.
        split_catchment_polygon, lat_snapped, lon_snapped, on_mainstem = _split_catchment(
            megabasin, lat, lng, home_unit_catchment_polygon, config)
        if split_catchment_polygon is None:
            return None, None, None

        if lat_snapped is None or lon_snapped is None:
            return None, None, None

        if config.calc_area:
            split_catchment_area = _get_polygon_area(split_catchment_polygon)

    else:
        split_catchment_polygon = None

    # Step 6: Find the set of upstream catchments
    # For a singleton unit catchment (no upstream neighbors),
    # the watershed is just the split home catchment; otherwise
    # we assemble it from the set of upstream unit catchments.
    #
    # inflow_area is the drainage area (km²) entering the home catchment from the
    # upstream catchments that drain to the outlet, or None if unknown. By
    # default that is every upstream catchment: the upstream area at the home
    # catchment's outlet minus the home catchment itself.
    if upstream_area is None:
        inflow_area = None
    else:
        inflow_area = upstream_area - home_catchment_area

    if is_singleton:
        if split_catchment_polygon is None:
            return None, None, None

        watershed_polygon = split_catchment_polygon
    elif config.high_res and not on_mainstem:
        # Nothing upstream drains to the outlet: it lies on a small tributary
        # inside the home catchment (or off the river network, with snapping
        # disabled). The watershed is just the local split.
        logger.warning(
            "None of the upstream unit catchments drain to this outlet: it lies on a "
            "small tributary or off the river network within its unit catchment. "
            "Returning only the local catchment above the outlet. Raise "
            "stream_threshold_km2 (or move the outlet onto a larger river) to "
            "delineate the full upstream watershed."
        )
        upstream_unit_catchments = [home_unit_catchment]
        inflow_area = 0.0
        watershed_polygon = split_catchment_polygon

    else:
        # Step 7: Merge and dissolve the polygons. basins_dict, the
        # hierarchical nested set of catchments, was fetched in Step 4 while
        # the basins database connection was open. The union covers only the
        # upstream catchments; the home catchment geometry (split in
        # high-res mode, whole in low-res mode) is merged in here.
        upstream_union = _get_upstream_union(basins_db_path, home_unit_catchment,
                                             basins_dict, config.cache)
        home_geom = split_catchment_polygon if config.high_res else home_unit_catchment_polygon
        watershed_polygon = unary_union([upstream_union, home_geom])

    # Fill interior "donut holes". These appear both in raster-split home
    # catchments (including single-catchment watersheds) and in the gaps between
    # dissolved unit catchments, so filling applies to every watershed.
    if config.fill:
        watershed_polygon = _close_holes(watershed_polygon, config.fill_area_max)

    watershed_gdf = gpd.GeoDataFrame(geometry=[watershed_polygon], crs='epsg:4326')

    # The snapped outlet only exists in high-res mode (the low-res path does no
    # snapping and never defines lat_snapped/lon_snapped). Singletons are always
    # high-res, so this reports the snapped outlet for them too.
    if config.high_res:
        watershed_gdf["outlet_lat"] = round(lat_snapped, 4)
        watershed_gdf["outlet_lng"] = round(lon_snapped, 4)

    # Calculate the area of the watershed
    if config.calc_area:
        if inflow_area is None:
            # The rivers database (which stores the pre-computed upstream
            # area) was unavailable or had no record for this catchment;
            # return the watershed without area attributes rather than fail.
            logger.warning("Could not determine the upstream area; "
                           "the watershed is returned without area_km2.")
            watershed_area = None
        elif config.high_res:
            watershed_area = inflow_area + split_catchment_area
        else:
            # Low-res mode merges in the whole home catchment
            watershed_area = inflow_area + home_catchment_area

        if watershed_area is not None:
            watershed_gdf["area_km2"] = round(watershed_area, 1)

            # watershed_area is the drainage area summed from the source-dataset
            # unit catchments. When fill=True, _close_holes adds the area of the
            # filled interior holes to the polygon, so that sum no longer matches
            # the delivered polygon. area_km2 then reports the polygon's actual
            # area, and filled_area reports how much of it the filling added.
            if config.fill:
                filled_polygon_area = _get_polygon_area(watershed_polygon)
                watershed_gdf["area_km2"] = round(filled_polygon_area, 1)
                watershed_gdf["filled_area"] = round(filled_polygon_area - watershed_area, 1)

    # Step 7.5: optionally clean, simplify, and smooth the watershed
    if config.smooth and not config.simplify:
        logger.warning(
            "smooth=True works best together with simplify=True: without "
            "simplification, the smoothed geometry traces the raster "
            "'staircase' instead of removing it."
        )

    if config.simplify:
        tolerance_deg = config.simplify_tolerance / KM_PER_DEGREE
        watershed_gdf.geometry = watershed_gdf.geometry.simplify(
            tolerance=tolerance_deg, preserve_topology=True)

    if config.clean:
        # This is a bit hack-ish, but it seems to work to remove "seams"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS")
            watershed_gdf.geometry = watershed_gdf.geometry.buffer(0.0001, join_style="mitre").buffer(-0.0001, join_style="mitre")

    # Smoothing runs last, after simplify and clean: it is a cosmetic
    # refinement of whatever boundary the earlier steps settled on.
    if config.smooth:
        watershed_gdf.geometry = watershed_gdf.geometry.apply(_smooth_watershed_geometry)

    # Step 8: Get the rivers if requested by the user
    if config.rivers:
        logger.info("Getting rivers")
        # upstream_unit_catchments was narrowed above to the home catchment
        # alone when the outlet is off its reach, so other reaches are left out.
        rivers_gdf = get_rivers(upstream_unit_catchments, split_catchment_polygon, config)
        if rivers_gdf is not None:
            if config.simplify:
                rivers_gdf.geometry = rivers_gdf.geometry.simplify(
                    tolerance=tolerance_deg, preserve_topology=True)
            # Catmull-Rom must come after simplification (see smoothing.py).
            if config.smooth:
                rivers_gdf.geometry = rivers_gdf.geometry.apply(_smooth_river_geometry)
    else:
        rivers_gdf = None

    # Step 9: Create a GeoDataFrame of the requested and snapped outlet point
    if config.outlets:
        if config.high_res:
            snapped_outlet = Point(lon_snapped, lat_snapped)
            outlets_gdf = gpd.GeoDataFrame(
                [
                    {'type': 'requested', 'geometry': requested_outlet},
                    {'type': 'snapped', 'geometry': snapped_outlet},
                ],
                crs='EPSG:4326'
            )
        else:
            outlets_gdf = gpd.GeoDataFrame([{'type': 'requested', 'geometry': requested_outlet}], crs='epsg:4326')
    else:
        outlets_gdf = None

    # Round coordinates in the geodata files, if requested by the user
    if config.round_coordinates:
        grid_size = 10 ** -ROUNDING_DECIMALS
        watershed_gdf['geometry'] = shapely.set_precision(watershed_gdf['geometry'], grid_size=grid_size)
        # rivers_gdf is None when rivers were requested but none exist
        # (e.g. coastal catchments have no river geometries)
        if rivers_gdf is not None:
            rivers_gdf['geometry'] = shapely.set_precision(rivers_gdf['geometry'], grid_size=grid_size)
        
    return watershed_gdf, rivers_gdf, outlets_gdf


def downloader(basin: int, data_dir: str | None = None) -> dict | None:
    """
    Download the data files for one megabasin, if not already present.

    Fetches the four files that delineation needs for the megabasin: the
    vector data (unit catchments, rivers) and the raster data (flow
    direction, flow accumulation). Files already on disk are not downloaded
    again. Use this to prepare data ahead of time, for example before
    working offline; otherwise ``delineate`` downloads files on first use.

    Parameters
    ----------
    basin : int
        The megabasin ID, an integer from 11 to 86 (see
        ``constants.VALID_MEGABASINS``). For example, megabasin 56 covers
        Australia.
    data_dir : str, optional
        Directory to download the data files into. If not provided, uses the
        ``DELINEATOR_DATA_DIR`` environment variable if set, and otherwise your
        system's default data directory, for example on Windows:
        C:\\Users\\<username>\\AppData\\Local\\delineator

    Returns
    -------
    dict[str, Path | None] or None
        The local path of each data file, keyed by file type ("Unit catchments",
        "Rivers", "Flow Direction", "Flow Accumulation"). A path is None if that
        file could not be found or downloaded. Returns None if the megabasin ID
        is invalid.

    Examples
    --------
    >>> from delineator import downloader
    >>> paths = downloader(basin=56, data_dir="path/to/data/directory")
    """
    # Pass data_dir to the constructor, not set it afterwards, so the config
    # marks it as a custom (flat) directory, the same as delineate() does.
    # Otherwise the files would download into a v3/ subfolder of data_dir.
    config = DelineatorConfig(data_dir=data_dir)

    # Check that the megabasin is valid
    if basin not in VALID_MEGABASINS:
        logger.warning(f"Invalid megabasin ID. Must be one of: {VALID_MEGABASINS}")
        return None

    files = {
        "Unit catchments": f"basins{basin}.db",
        "Rivers": f"rivers{basin}.db",
        "Flow Direction": f"flowdir{basin}.tif",
        "Flow Accumulation": f"accum{basin}.tif",
    }

    paths = {}
    for filetype, file in files.items():
        path = _find_data_file(file, config)
        print(f"{filetype} file is at: {path}")
        paths[filetype] = path
    return paths
