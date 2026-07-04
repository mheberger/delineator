"""
Performs detailed raster-based watershed delineation with `pysheds`,
but only inside of a *single* unit catchment.
"""
import os
import logging
import numpy as np
from numpy import floor, ceil
from pysheds.grid import Grid
from shapely.geometry import Polygon, MultiPolygon, shape
from shapely import wkb

from delineator.constants import HALF_PIXEL, DIRMAP
from delineator.data import _find_data_file
from delineator.settings import DelineatorConfig

logger = logging.getLogger(__name__)

# Suppress a harmless warning from rasterio:
# "DeprecationWarning: 'Memory' driver is deprecated since GDAL 3.11. Use 'MEM' onwards. Further messages of this type will be suppressed."
logging.getLogger("rasterio._env").setLevel(logging.ERROR)

def _split_catchment(basin: int, lat: float, lng: float, catchment_poly: Polygon | MultiPolygon,
                    bSingleCatchment: bool, config: DelineatorConfig,) -> \
        tuple[Polygon | MultiPolygon | None, float | None, float | None]:
    """
    Performs the detailed pixel-scale raster-based delineation for a watershed.

    To efficiently delineate large watersheds, we only use raster-based methods in a small area,
    the size of a single unit catchment that is most downstream. This results in big
    savings in processing time and memory use, making it possible to delineate even large watersheds
    on a laptop computer.

    Parameters
    ----------
    basin : int
        2-digit Pfafstetter code for the level 2 basin we're in (tells us which raster files to open)
    lat : float
        latitude
    lng : float
        longitude
    catchment_poly : Shapely Polygon or MultiPolygon
        a Shapely polygon; we'll use it to clip the flow accumulation raster to get an accurate snap
    bSingleCatchment : bool
        is the watershed small, i.e. there is only one unit catchment in it?
            If so, we'll use a lower snap threshold for the outlet.
    config : DelineatorConfig
        a DelineatorConfig dataclass object

    Returns
    -------
    poly : shapely Polygon or MultiPolygon
        the part of the terminal unit catchment that is upstream of the outlet point
    lat_snap : float
        latitude of the outlet, snapped to the river centerline in the accumulation raster
    lng_snap : float
        longitude of the outlet, snapped to the river centerline in the accumulation raster

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
    fdir_file = f"flowdir{basin}.tif"
    fdir_path = _find_data_file(fdir_file, config)
    if fdir_path is None:
        logging.warning(f"Could not load flow direction raster: {fdir_file}")
        return None, None, None
    logger.info("Loading flow direction raster from: {}".format(fdir_path))

    if not os.path.isfile(fdir_path):
        logging.warning("Could not find flow flow direction raster: {}".format(fdir_path))
        return None, None, None

    # Get the Grid object that matches the flow direction grid, which we
    # can use to rasterize the unit catchment polygon
    grid = Grid.from_raster(str(fdir_path), window=bounding_box, nodata=0)

    # Now "clip" the rectangular flow direction grid even further so that it ONLY contains data
    # inside the boundaries of the terminal unit catchment.
    hexpoly = catchment_poly.wkb_hex
    poly = wkb.loads(hexpoly, hex=True)

    # Fix any holes in the polygon by taking the exterior coordinates.
    # Updated 2026-07-02 to remove. Now we are handling MultiPolygons properly. No longer need this hack!
    #filled_poly = _close_holes(poly, area_max=0)

    # It needs to be of type MultiPolygon to work with rasterio apparently
    if isinstance(poly, Polygon):
        poly = MultiPolygon([poly])

    polygon_list = list(poly.geoms)

    # Convert the polygon into a pixelized raster "mask".
    mymask = grid.rasterize(polygon_list)

    # LOAD Flow Direction Grid
    fdir = grid.read_raster(fdir_path, window=bounding_box, nodata=0)

    # Not clear if this this step was unnecessary, but it makes the plots look nicer
    fdir[mymask == 0] = 0

    # Pour point snapping
    def snap_outlet() -> tuple[float | None, float | None]:
        """Snap the user-provided lat/lng coordinates to a valid point on the grid."""
        logger.info("Snapping pour point")

        # Open the accumulation raster, again using windowed reading mode.
        accum_file = f"accum{basin}.tif"
        accum_path = _find_data_file(accum_file, config)
        if accum_path is None:
            logging.warning(f"Could not load accumulation raster: {accum_file}")
            return None, None

        if not os.path.isfile(accum_path):
            raise Exception("Could not find accumulation raster: {}".format(accum_path))

        # pysheds grids expects a string, not a Path object
        acc = grid.read_raster(
            str(accum_path), window=bounding_box, window_crs=grid.crs, nodata=0
        )

        # MASK the accumulation raster to the unit catchment POLYGON.
        # I found that this step is essential for good pour point snapping
        acc[mymask == 0] = 0

        # Snap the outlet to the nearest stream.
        if bSingleCatchment:
            numpixels = config.threshold_single
        else:
            # Case where there are 2 or more unit catchments in the watershed
            # setting this value too low causes incorrect results and weird topology problems in the output
            numpixels = config.threshold_multi

        logger.info(
            "Using threshold of {} for number of upstream pixels.".format(numpixels)
        )

        # Snap the pour point to a point on the accumulation grid where accum (# of upstream pixels)
        # is greater than our threshold. This is conventionally called the streams layer.
        streams = acc > numpixels
        xy = (lng, lat)
        try:
            [lng_snap, lat_snap] = grid.snap_to_mask(
                streams, xy
            )  # New version of pysheds does not give you the snap distance.
            return lat_snap, lng_snap

        except Exception as e:
            logger.warning(f"Could not snap the pour point. Error: {e}")
            return None, None

    # Snap the pour point
    if config.snapping:
        lat_snap, lng_snap = snap_outlet()
        if lat_snap is None:
            return None, None, None
        snap_setting = "corner"
        #snap_setting = "center"
    else:
        # Snapping to center is consistent with TauDEM
        lat_snap, lng_snap = lat, lng
        snap_setting = "center"

    # Finally, here is the raster based watershed delineation with pysheds!
    logger.info("Splitting home unit catchment")
    try:
        catch = grid.catchment(fdir=fdir, x=lng_snap, y=lat_snap, dirmap=DIRMAP,
                               xytype='coordinate', recursionlimit=15000, snap=snap_setting)
        # Clip the bounding box to the catchment
        grid.clip_to(catch)
        clipped_catch = grid.view(catch, dtype=np.uint8)
    except Exception as e:
        logger.error(f"ERROR: something went wrong during pysheds grid.catchment(). Error: {e}")
        return None, lng_snap, lat_snap

    # Convert high-precision raster subcatchment to a polygon using pysheds method .polygonize()
    logger.info("Convert split catchment raster to polygon")
    shapes = grid.polygonize(clipped_catch)

    polygons = [shape(geom) for geom, value in shapes]
    multi = MultiPolygon(polygons)

    # The snapped vertices need to be nudged down and right by one half pixel
    lng_snap += HALF_PIXEL
    lat_snap -= HALF_PIXEL

    return multi, lat_snap, lng_snap
