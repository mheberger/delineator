"""
Version 2.0 of `delineator` Global Watershed Delineation with Python
  - refactored from a collection of scripts to a single Python package
  - uses sqlite .db files instead of shapefiles
  - implements "hierarchical nested" set of catchment boundaries to speed up vector operations.
"""

import logging
from pathlib import Path
import geopandas as gpd
from shapely.ops import unary_union
from shapely.geometry import Point
import sqlite3
import warnings

from delineator.constants import (
    MEGABASINS_DB_FILE,
    VALID_MEGABASINS,
    USE_SPATIALITE,
    PIXEL_AREA,
)
from delineator.data import _find_data_file
from delineator.util import _validate_outlet_coordinates, _close_holes, _get_polygon_area
from delineator.spatial import _point_in_polygon_analysis
from delineator.queries import get_upstream_comids, get_home_unit_catchment_geom_and_area, get_hiearchical_basins, \
    get_upstream_geometries, get_rivers, get_upstream_area
from delineator.merit_detailed import split_catchment
from delineator.settings import DelineatorConfig

logger = logging.getLogger(__name__)


def delineate(
    lat: float,
    lng: float,
    config: DelineatorConfig | None = None
) -> tuple[
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
        Configuration object. If not provided, one is created from kwargs.

    Returns
    -------
    3 GeoPandas GeoDataFrames:
        - watershed
        - rivers
        - outlets

    If any of these cannot be produced, they will be None.
    This function will log a warning if the watershed cannot be created
    but will not raise an error.
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

    # Set up database connection
    megabasins_db_conn = sqlite3.connect(MEGABASINS_DB_FILE)

    # Step 1: Validate the input coordinates
    try:
        lat, lng = _validate_outlet_coordinates(lat, lng)
    except Exception as e:
        logger.warning(f"ERROR: Invalid coordinates: {e}")
        return None, None, None

    print(f"Delineating watershed for outlet: {lat}, {lng}")

    # Step 2: Determine the megabasin
    requested_outlet = Point(lng, lat)
    megabasin = _point_in_polygon_analysis(megabasins_db_conn, requested_outlet, table='megabasins',
                                           geom_col='geometry', id_col='basin', use_spatialite=USE_SPATIALITE,
                                           search_dist=config.search_dist)
    if megabasin is None:
        logger.warning(
            "The requested outlet is not in any megabasin. Check that it is over land. Maybe you flipped the lat/lng?")
        return None, None, None

    logger.info(f"Your watershed is in megabasin {megabasin}")

    # Step 3: Determine the unit catchment of the outlet
    basins_db_file = f'basins{megabasin}.db'
    basins_db_path = _find_data_file(basins_db_file, config)
    basins_db_conn = sqlite3.connect(basins_db_path)
    home_unit_catchment = _point_in_polygon_analysis(basins_db_conn, requested_outlet, table='l0_basins',
                                                     id_col='comid', use_spatialite=USE_SPATIALITE,
                                                     search_dist=config.search_dist)

    # If the outlet is not in any unit catchment, return None
    if home_unit_catchment is None:
        logger.warning("The requested outlet is not in any unit catchment. Check that it is over land.")
        return None, None, None

    logger.info(f"Your watershed is in unit catchment {home_unit_catchment}")

    # Step 4: Get the collection of upstream catchment polygons
    upstream_unit_catchments = get_upstream_comids(basins_db_conn, int(home_unit_catchment))
    is_singleton = len(upstream_unit_catchments) == 1

    upstream_area = get_upstream_area(int(home_unit_catchment), config)

    # Step 4.5 Determine if the watershed is so huge that the user would prefer to use the low-res mode
    if config.high_res:
        # Get the area of the home unit catchment
        if upstream_area > config.low_res_threshold:
            logger.info(f"Switching to low-resolution mode; watershed area = {upstream_area:.2f} km²")
            config.high_res = False

    # Step 5: Split the home unit catchment at the outlet
    # First, we need to get the home unit catchment geometry as a shapely Polygon
    home_unit_catchment_polygon, home_catchment_area = get_home_unit_catchment_geom_and_area(basins_db_path,
                                                                                             home_unit_catchment)
    if config.high_res:
        # Perform the split
        split_catchment_polygon, lat_snapped, lon_snapped = split_catchment(megabasin, lat, lng,
                                                                            home_unit_catchment_polygon, is_singleton,
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
    # If the watershed is a singleton, we don't need to do anything
    if is_singleton:
        if split_catchment_polygon is None:
            return None, None, None

        watershed_gdf = gpd.GeoDataFrame(geometry=[split_catchment_polygon], crs='epsg:4326')
    else:
        # Otherwise, we need to find the set of upstream catchments
        # First, we need to get the hierarchical nested set of catchments
        basins_dict = get_hiearchical_basins(basins_db_conn, upstream_unit_catchments[1:])

        # In low-res mode, we do need to include the home unit catchment, as we will use the original geometry
        if not config.high_res:
            basins_dict['L0'].append(home_unit_catchment)

        upstream_polygons = get_upstream_geometries(basins_db_path, basins_dict)

        # In high-res mode, do not include the geoometries of the home unit catchment, as we will use the
        # split unit catchment polygon obtained in step 5
        if config.high_res:
            upstream_polygons.append(split_catchment_polygon)

        # Step 7: Merge and dissolve the polygons
        watershed_polygon = unary_union(upstream_polygons)

        if config.fill:
            watershed_polygon = _close_holes(watershed_polygon, config.fill_threshold * PIXEL_AREA)

        watershed_gdf = gpd.GeoDataFrame(geometry=[watershed_polygon], crs='epsg:4326')
        
        watershed_gdf["outlet_lat"] = round(lat_snapped, 4)
        watershed_gdf["outlet_lng"] = round(lon_snapped, 4)

    # Calculate the area of the watershed
    if config.calc_area:
        if config.high_res:
            watershed_area = upstream_area - home_catchment_area + split_catchment_area
        else:
            watershed_area = upstream_area - home_catchment_area  # Check this!!!

        watershed_gdf["area_km2"] = round(watershed_area, 1)

    # Step 7.5: optionally clean and simplify the watershed
    if config.simplify:
        watershed_gdf.geometry = watershed_gdf.geometry.simplify(
            tolerance=config.simplify_tolerance, preserve_topology=True)

    if config.clean:
        # This is a bit hack-ish, but it seems to work to remove "seams"
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Geometry is in a geographic CRS")
            watershed_gdf.geometry = watershed_gdf.geometry.buffer(0.0001).buffer(-0.0001)

    # Step 8: Get the rivers if requested by the user
    if config.rivers:
        logger.info("Getting rivers")
        rivers_gdf = get_rivers(upstream_unit_catchments, split_catchment_polygon, config)
        if config.simplify:
            rivers_gdf.geometry = rivers_gdf.geometry.simplify(
                tolerance=config.simplify_tolerance, preserve_topology=True)
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

    return watershed_gdf, rivers_gdf, outlets_gdf


def downloader(basin: int, data_dir: str | None = None):
    """
    Utility function to download the data files for a given basin.

    Parameters:
    -----------
    basin: int
        The megabasin ID, and integer from 11 to 86
        For example, megabasin 56 covers Australia

    data_dir: str, optional
        Location where the data files will be downloaded. If not provided,
        defaults to your system's default data directory, for example, on Windows:
        C:\\Users\\<username>\\AppData\\Local\\delineator

    Usage:
    ------
        `downloader(basin=11, data_dir="path/to/data/directory")`

        This will download the data files for megabasin 11 to the specified directory,
        including the vector data (unit catchments, rivers) and raster data (flow direction, accumulation)

    Returns:
    --------
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
