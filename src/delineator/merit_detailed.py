"""
Performs a detailed, raster-based watershed delineation with `pysheds`,
but only inside of a *single* unit catchment.
This is my implementation of the *hybrid* method that I "invented" -- but which was actually
first described in a a paper by Djokic and Ye at the 1999 ESRI User Conference
and is used by the USGS as part of their watershed delineation API in the NLDI.
Raster-based delineation is slow and requires a lot of memory. So we only do the bare minimum,
and use vector data for the rest of the upstream watershed.
"""
import os
import logging
import numpy as np
from numpy import floor, ceil
from pysheds.grid import Grid
from shapely.geometry import Polygon, MultiPolygon
from shapely import wkb, ops

from delineator.constants import THRESHOLD_SINGLE, THRESHOLD_MULTIPLE, HALF_PIXEL, DIRMAP
from delineator.data import _find_data_file
from delineator.settings import DelineatorConfig
from delineator.util import _close_holes

logger = logging.getLogger(__name__)


def split_catchment(
    basin: int,
    lat: float,
    lng: float,
    catchment_poly: Polygon | MultiPolygon,
    bSingleCatchment: bool,
    config: DelineatorConfig,
) -> tuple[Polygon | MultiPolygon | None, float | None, float | None]:
    """
    Performs the detailed pixel-scale raster-based delineation for a watershed.

    To efficiently delineate large watersheds, we only use raster-based methods in a small area,
    the size of a single unit catchment that is most downstream. This results in big
    savings in processing time and memory use, making it possible to delineate even large watersheds
    on a laptop computer.

    Args:
        wid: the watershed id, a string
        basin: 2-digit Pfafstetter code for the level 2 basin we're in (tells us which raster files to open)
        lat: latitude
        lng: longitude
        catchment_poly: a Shapely polygon; we'll use it to clip the flow accumulation raster to get an accurate snap
        bSingleCatchment: is the watershed small, i.e. there is only one unit catchment in it?
            If so, we'll use a lower snap threshold for the outlet.
        config: a DelineatorConfig dataclass object

    Returns:
        poly: a shapely Polygon or MultiPolygon representing the part of the terminal unit catchment
            that is upstream of the outlet point
        lat_snap:  latitude of the outlet, snapped to the river centerline in the accumulation raster
        lng_snap: longitude of the outlet, snapped to the river centerline in the accumulation raster

    """

    # Get a bounding box for the unit catchment
    bounds = catchment_poly.bounds
    bounds_list = [float(i) for i in bounds]

    # The coordinates of the bounding box edges that we get from the above query
    # do not correspond well with the edges of the grid pixels.
    # We need to round them to the nearest whole pixel and then
    # adjust them by a half-pixel width to get good results in pysheds.

    # Bounding box is xmin, ymin, xmax, ymax
    # round the elements DOWN, DOWN, UP, UP
    # The number 1200 is because the MERIT-Hydro rasters have 3 arsecond resolution, or 1/1200 of a decimal degree.
    # So we just multiply it by 1200, round up or down to the nearest whole number, then divide by 1200
    # to put it back in its regular units of decimal degrees. Then, since pysheds wants the *center*
    # of the pixel, not its edge, add or subtract a half-pixel width as appropriate.
    # This took me a while to figure out but was essential to getting results that look correct
    bounds_list[0] = floor(bounds_list[0] * 1200) / 1200 - HALF_PIXEL
    bounds_list[1] = floor(bounds_list[1] * 1200) / 1200 - HALF_PIXEL
    bounds_list[2] = ceil(bounds_list[2] * 1200) / 1200 + HALF_PIXEL
    bounds_list[3] = ceil(bounds_list[3] * 1200) / 1200 + HALF_PIXEL
    # The bounding box needs to be a tuple for pysheds.
    bounding_box = tuple(bounds_list)

    # Open the flow direction raster *using windowed reading mode*
    fdir_file = f"raster/flowdir{basin}.tif"
    fdir_path = _find_data_file(fdir_file, config)

    logger.info("Loading flow direction raster from: {}".format(fdir_path))

    if not os.path.isfile(fdir_path):
        raise Exception("Could not find flow flow direction raster: {}".format(fdir_path))

    # Get the Grid object that matches the flow direction grid, which we
    # can use to rasterize the unit catchment polygon
    grid = Grid.from_raster(str(fdir_path), window=bounding_box, nodata=0)

    # Now "clip" the rectangular flow direction grid even further so that it ONLY contains data
    # inside the boundaries of the terminal unit catchment.
    hexpoly = catchment_poly.wkb_hex
    poly = wkb.loads(hexpoly, hex=True)

    # Fix any holes in the polygon by taking the exterior coordinates.
    filled_poly = _close_holes(poly, area_max=0)

    # It needs to be of type MultiPolygon to work with rasterio apparently
    if isinstance(filled_poly, Polygon):
        filled_poly = MultiPolygon([filled_poly])

    polygon_list = list(filled_poly.geoms)

    # Convert the polygon into a pixelized raster "mask".
    mymask = grid.rasterize(polygon_list)

    # LOAD Flow Direction Grid
    fdir = grid.read_raster(fdir_path, window=bounding_box, nodata=0)

    # Not clear if this this step was unnecessary, but it makes the plots look nicer
    fdir[mymask == 0] = 0

    logger.info("Snapping pour point")

    # Open the accumulation raster, again using windowed reading mode.
    accum_file = f'raster/accum{basin}.tif'
    accum_path = _find_data_file(accum_file, config)

    if not os.path.isfile(accum_path):
        raise Exception("Could not find accumulation raster: {}".format(accum_path))

    # pysheds grids expects a string, not a Path object
    acc = grid.read_raster(str(accum_path), window=bounding_box, window_crs=grid.crs, nodata=0)

    # MASK the accumulation raster to the unit catchment POLYGON.
    acc[mymask == 0] = 0

    # Snap the outlet to the nearest stream.
    if bSingleCatchment:
        numpixels = THRESHOLD_SINGLE
    else:
        # Case where there are 2 or more unit catchments in the watershed
        # setting this value too low causes incorrect results and weird topology problems in the output
        numpixels = THRESHOLD_MULTIPLE

    logger.info("Using threshold of {} for number of upstream pixels.".format(numpixels))

    # Snap the pour point to a point on the accumulation grid where accum (# of upstream pixels)
    # is greater than our threshold
    streams = acc > numpixels
    xy = (lng, lat)
    try:
        [lng_snap, lat_snap] = grid.snap_to_mask(streams, xy)  # New version does not give you the snap distance.
    except Exception as e:
        logger.warning(f"Could not snap the pour point. Error: {e}")
        return None, None, None

    # Finally, here is the raster based watershed delineation with pysheds!
    logger.info("Splitting home unit catchment")
    try:
        catch = grid.catchment(fdir=fdir, x=lng_snap, y=lat_snap, dirmap=DIRMAP,
                               xytype='coordinate', recursionlimit=15000)

        # Clip the bounding box to the catchment
        # Seems optional, but turns out this line is essential.
        grid.clip_to(catch)
        clipped_catch = grid.view(catch, dtype=np.uint8)
    except Exception as e:
        logger.error(f"ERROR: something went wrong during pysheds grid.catchment(). Error: {e}")
        return None, lng_snap, lat_snap

    # Convert high-precision raster subcatchment to a polygon using pysheds method .polygonize()
    logger.info("Convert split catchment raster to polygon")
    shapes = grid.polygonize(clipped_catch)

    # The output from pysheds is creating MANY shapes.
    # Dissolve them together with the unary union operation in shapely
    shapely_polygons = []

    shape_count = 0

    # The snapped vertices need to be nudged down and right by one half pixel
    lng_snap += HALF_PIXEL
    lat_snap -= HALF_PIXEL

    # Convert the result from pysheds into a list of shapely polygons
    for shape, value in shapes:
        pysheds_polygon = shape
        shape_count += 1
        # The pyshseds polygon can be converted to a shapely Polygon in this one-liner
        shapely_polygon = Polygon([[p[0], p[1]] for p in pysheds_polygon['coordinates'][0]])
        shapely_polygons.append(shapely_polygon)

    if shape_count > 1:
        # If pysheds returned multiple polygons, dissolve them using shapely's unary_union() function
        # Note that this can sometimes return a MultiPolygon, which we'll need to fix later
        result_polygon = ops.unary_union(shapely_polygons)

    else:
        # If pysheds generated a single polygon, that is our answer
        result_polygon = shapely_polygons[0]

    return result_polygon, lat_snap, lng_snap
