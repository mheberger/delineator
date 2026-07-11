"""
delineator 2.0: Global Watershed Delineation with Python

Matthew Heberger, May 2026

Formerly a collection of Python scripts, rewritten as a Python package

"""

import logging
from collections import OrderedDict
from contextlib import closing
from pathlib import Path
import geopandas as gpd
import shapely
from shapely.ops import unary_union
from shapely.geometry import Point
import sqlite3
import warnings

from delineator.constants import (
    MEGABASINS_DB_FILE,
    VALID_MEGABASINS,
    ROUNDING_DECIMALS,
    KM_PER_DEGREE,
    DISCONNECT_TOL_DEG,
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
_CACHE: OrderedDict = OrderedDict()
_CACHE_MAXSIZE = 64


def clear_cache() -> None:
    """Empty the upstream-catchment cache (see ``DelineatorConfig.cache``)."""
    _CACHE.clear()


def _cache_get(key):
    value = _CACHE.get(key)
    if value is not None:
        _CACHE.move_to_end(key)
    return value


def _cache_put(key, value) -> None:
    _CACHE[key] = value
    while len(_CACHE) > _CACHE_MAXSIZE:
        _CACHE.popitem(last=False)


def _get_upstream_catchments(basins_db_conn, home_unit_catchment, use_cache: bool):
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


def _get_upstream_union(basins_db_path, home_unit_catchment, basins_dict,
                        use_cache: bool):
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
    Main watershed delineation routine.

    Parameters
    ----------
    lat, lng : float
        Coordinates of the outlet point in decimal degrees.
    config : DelineatorConfig, optional
        Configuration object. If not provided, one is created with default settings.

    Returns
    -------
    watershed : GeoDataFrame
        The watershed boundary polygon.
    rivers : GeoDataFrame
        The river reaches within the watershed.
    outlets : GeoDataFrame
        The requested and snapped watershed outlets.

    Notes
    -----
    Any of the returned objects will be None if it cannot be produced.
    This function logs a warning if the watershed cannot be created but
    does not raise an error.
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
        # Perform the split
        split_catchment_polygon, lat_snapped, lon_snapped = _split_catchment(megabasin, lat, lng,
                                                                             home_unit_catchment_polygon,
                                                                             home_catchment_area,
                                                                             config)
        if split_catchment_polygon is None:
            return None, None, None

        if lat_snapped is None or lon_snapped is None:
            return None, None, None

        if config.calc_area:
            split_catchment_area = _get_polygon_area(split_catchment_polygon)

    else:
        split_catchment_polygon = None

    # Step 6: Find the set of upstream catchments
    # For a singleton the watershed is just the split home catchment; otherwise
    # we assemble it from the set of upstream unit catchments.
    #
    # local_only records that we collapsed a non-singleton watershed down to just
    # the home catchment's local split because the outlet turned out to be off the
    # mainstem (see the disconnection check below). It changes how area and rivers
    # are handled downstream.
    local_only = False
    if is_singleton:
        if split_catchment_polygon is None:
            return None, None, None

        watershed_polygon = split_catchment_polygon
    else:
        # Step 7: Merge and dissolve the polygons. basins_dict, the
        # hierarchical nested set of catchments, was fetched in Step 4 while
        # the basins database connection was open. The union covers only the
        # upstream catchments; the home catchment geometry (split or whole)
        # is merged in below.
        upstream_union = _get_upstream_union(basins_db_path, home_unit_catchment,
                                             basins_dict, config.cache)

        if config.high_res:
            # The upstream neighbours are only genuinely upstream of the outlet
            # when the raster split of the home catchment reaches the shared edge
            # where their rivers enter -- which requires the outlet to lie on the
            # mainstem. Snapping usually guarantees that; with snapping disabled the user
            # may place the outlet off-stream, and then the split is a small local
            # catchment that does not connect to the upstream block. Splicing the
            # two together would produce a disconnected polygon and credit the
            # outlet with upstream drainage area that does not reach it. Detect
            # that and keep only the local catchment above the outlet.
            if (not config.snapping):  # Only do this expensive check when snapping is disabled.
                if split_catchment_polygon.distance(upstream_union) > DISCONNECT_TOL_DEG:
                    logger.warning(
                        "The outlet's local catchment is disconnected from its upstream "
                        "catchments: the outlet appears to lie off the river network "
                        "while snapping is disabled. Returning only the local catchment "
                        "above the outlet. Adjust snapping tolerance (or move the outlet onto a "
                        "river) to delineate the full upstream watershed."
                    )
                    local_only = True
                    watershed_polygon = split_catchment_polygon
                else:
                    pass
            else:
                # Merge the split home catchment into the upstream block.
                watershed_polygon = unary_union([upstream_union, split_catchment_polygon])
        else:
            # Low-res case: merge in the home catchment's original (unsplit)
            # geometry, fetched in Step 5.
            watershed_polygon = unary_union([upstream_union, home_unit_catchment_polygon])

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
        if local_only:
            # Collapsed to the home catchment's local split; the upstream
            # catchments are not actually upstream of this outlet.
            watershed_area = split_catchment_area
        elif upstream_area is None:
            # The rivers database (which stores the pre-computed upstream
            # area) was unavailable or had no record for this catchment;
            # return the watershed without area attributes rather than fail.
            logger.warning("Could not determine the upstream area; "
                           "the watershed is returned without area_km2.")
            watershed_area = None
        elif config.high_res:
            watershed_area = upstream_area - home_catchment_area + split_catchment_area
        else:
            watershed_area = upstream_area - home_catchment_area

        if watershed_area is not None:
            watershed_gdf["area_km2"] = round(watershed_area, 1)

            # area_km2 is the drainage area summed from the source-dataset unit
            # catchments. When fill=True, _close_holes adds the area of the filled
            # interior holes to the polygon, so its geometric area no longer matches
            # that figure. Rather than distort area_km2, report the added area
            # separately as area_filled (actual polygon area minus area_km2).
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
        # When collapsed to local-only, restrict rivers to the home catchment's
        # (split) reach; the upstream reaches do not drain to this outlet.
        rivers_catchments = [upstream_unit_catchments[0]] if local_only else upstream_unit_catchments
        rivers_gdf = get_rivers(rivers_catchments, split_catchment_polygon, config)
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


def downloader(basin: int, data_dir: str | None = None):
    """
    Utility function to download the data files for a given basin.
    This will download the data files for megabasin 11 to the specified directory,
    including the vector data (unit catchments, rivers) and raster data (flow
    direction, accumulation)

    Usage
    -----
    `downloader(basin=11, data_dir="path/to/data/directory")`

    Parameters
    ----------
    basin: int
        The megabasin ID, and integer from 11 to 86
        For example, megabasin 56 covers Australia

    data_dir: str, optional
        Location where the data files will be downloaded. If not provided,
        defaults to your system's default data directory, for example, on Windows:
        C:\\Users\\<username>\\AppData\\Local\\delineator

    Returns
    -------
    Nothing. Raises an error if something goes wrong.
    """
    config = DelineatorConfig()

    if data_dir is not None:
        config.data_dir = Path(data_dir)

    # Check that the megabasin is valid
    if basin not in VALID_MEGABASINS:
        raise ValueError(f"Invalid megabasin ID. Must be one of: {VALID_MEGABASINS}")

    files = {
        "Unit catchments": f"basins{basin}.db",
        "Rivers": f"rivers{basin}.db",
        "Flow Direction": f"flowdir{basin}.tif",
        "Flow Accumulation": f"accum{basin}.tif",
    }

    for filetype, file in files.items():
        path = _find_data_file(file, config)
        print(f"{filetype} file is at: {path}")
