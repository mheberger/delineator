"""
Performs detailed raster-based watershed delineation with `pysheds`,
but only inside of a *single* unit catchment.
"""
import os
import logging
import numpy as np
import rasterio
from numpy import floor, ceil
from pysheds.grid import Grid
from shapely.geometry import Polygon, MultiPolygon, shape
from shapely import wkb

from delineator.constants import DIRMAP
from delineator.data import _find_data_file
from delineator.settings import DelineatorConfig

logger = logging.getLogger(__name__)

# Radius of the authalic sphere in km, used for haversine distances.
EARTH_RADIUS_KM = 6371.0088


def _snap_to_stream(streams: np.ndarray, affine, lat: float, lng: float,
                    max_dist_km: float) -> tuple[float | None, float | None, float | None]:
    """
    Snap a pour point to the nearest stream cell by *geodesic* distance,
    subject to a maximum snap distance in kilometers.

    This replaces pysheds' `snap_to_mask()`, which finds the nearest cell in
    raw coordinate (degree) space. Degree-space distance is anisotropic: away
    from the equator, one degree of longitude is much shorter than one degree
    of latitude, so degree-nearest snapping is biased toward north-south
    neighbors by up to a factor of 1/cos(latitude). This function measures
    true (haversine) distances instead, and also returns the snap distance,
    which newer versions of pysheds no longer provide.

    Parameters
    ----------
    streams : array-like of bool
        Boolean raster; True where flow accumulation exceeds the stream threshold.
    affine : affine.Affine
        The affine transform of the streams raster.
    lat, lng : float
        Coordinates of the pour point provided by the user.
    max_dist_km : float
        Maximum allowed snap distance in kilometers. If the nearest stream
        cell is farther than this, the snap fails.

    Returns
    -------
    (lat_snap, lng_snap, dist_km), or (None, None, None) if there are no
    stream cells or the nearest one is beyond max_dist_km.

    Coordinate convention: the returned coordinates are the *center* of the
    snapped cell, i.e. affine * (col + 0.5, row + 0.5). Distances are also
    measured to cell centers. Callers passing the result to pysheds'
    grid.catchment() should use snap='center', which resolves a cell-center
    coordinate to its cell unambiguously.
    """
    rows, cols = np.nonzero(streams)
    if rows.size == 0:
        logger.warning("No stream cells above the accumulation threshold in the unit catchment.")
        return None, None, None

    # Cell centers for the distance calculation
    x_ctr = affine.c + (cols + 0.5) * affine.a
    y_ctr = affine.f + (rows + 0.5) * affine.e

    # Haversine distance from the pour point to every candidate cell.
    # Unit catchments are small, so this is at most a few thousand cells.
    phi1 = np.deg2rad(lat)
    phi2 = np.deg2rad(y_ctr)
    dphi = phi2 - phi1
    dlam = np.deg2rad(x_ctr - lng)
    h = np.sin(dphi / 2.0) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlam / 2.0) ** 2
    dist_km = 2.0 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(h))

    i = int(np.argmin(dist_km))
    d = float(dist_km[i])
    if d > max_dist_km:
        logger.warning(
            f"Nearest stream cell is {d:.2f} km from the pour point, "
            f"beyond the maximum snap distance of {max_dist_km:.2f} km."
        )
        return None, None, None

    # Return the cell center (we already computed centers for the distances)
    return float(y_ctr[i]), float(x_ctr[i]), d

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

    # Locate the flow direction raster (we will open it *using windowed
    # reading mode*, but first we need its transform to align the window).
    fdir_file = f"flowdir{basin}.tif"
    fdir_path = _find_data_file(fdir_file, config)
    if fdir_path is None:
        logger.warning(f"Could not load flow direction raster: {fdir_file}")
        return None, None, None
    logger.info("Loading flow direction raster from: {}".format(fdir_path))

    if not os.path.isfile(fdir_path):
        logger.warning("Could not find flow flow direction raster: {}".format(fdir_path))
        return None, None, None

    # Align the bounding box to the source raster's own pixel grid.
    #
    # pysheds uses the window bbox coordinates VERBATIM as the affine origin
    # of the windowed grid, pulling source pixels by nearest-index rounding.
    # If the bbox does not fall exactly on the source raster's pixel edges,
    # the windowed grid is georeferenced up to half a pixel away from where
    # the data truly sits, and everything downstream (snapping, delineation,
    # polygonization) comes out misregistered by that amount.
    #
    # Crucially, we make NO assumption about the raster's registration
    # (pixel edges vs. pixel centers on round-number coordinates). We read
    # the raster's affine transform and round the bounds outward to whole
    # pixel indices of THAT grid, padded by one pixel so that boundary-pixel
    # inclusion is robust to floating-point jitter when a polygon bound
    # lands exactly on a grid line.
    with rasterio.open(fdir_path) as src:
        T = src.transform
        src_height, src_width = src.height, src.width

    # Fractional pixel indices of the polygon bounds. T.e is negative for
    # north-up rasters, so ymax (north) maps to the smallest row index.
    col_min = int(floor((bounds[0] - T.c) / T.a)) - 1
    col_max = int(ceil((bounds[2] - T.c) / T.a)) + 1
    row_min = int(floor((bounds[3] - T.f) / T.e)) - 1
    row_max = int(ceil((bounds[1] - T.f) / T.e)) + 1

    # Clamp to the raster extent
    col_min = max(col_min, 0)
    row_min = max(row_min, 0)
    col_max = min(col_max, src_width)
    row_max = min(row_max, src_height)

    # Convert back to coordinates: these lie exactly on the source pixel grid
    xmin_aligned = T.c + col_min * T.a
    xmax_aligned = T.c + col_max * T.a
    ymax_aligned = T.f + row_min * T.e
    ymin_aligned = T.f + row_max * T.e

    # Bounding box is xmin, ymin, xmax, ymax; it needs to be a tuple for pysheds.
    bounding_box = (xmin_aligned, ymin_aligned, xmax_aligned, ymax_aligned)

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

        # Snap to the nearest stream cell by true (geodesic) distance, limited
        # to a maximum snap distance in km. See _snap_to_stream for details.
        max_dist_km = config.search_dist
        try:
            lat_snap, lng_snap, snap_dist_km = _snap_to_stream(
                streams, acc.affine, lat, lng, max_dist_km
            )
            if lat_snap is None:
                return None, None
            logger.info(f"Snapped pour point moved {snap_dist_km:.3f} km.")
            return lat_snap, lng_snap

        except Exception as e:
            logger.warning(f"Could not snap the pour point. Error: {e}")
            return None, None

    # Snap the pour point.
    # Both paths use the cell-center coordinate convention: _snap_to_stream
    # returns cell centers, and raw user coordinates resolved with
    # snap='center' is consistent with TauDEM's behavior.
    if config.snapping:
        lat_snap, lng_snap = snap_outlet()
        if lat_snap is None:
            return None, None, None
    else:
        lat_snap, lng_snap = lat, lng

    # Finally, here is the raster based watershed delineation with pysheds!
    logger.info("Splitting home unit catchment")
    try:
        catch = grid.catchment(fdir=fdir, x=lng_snap, y=lat_snap, dirmap=DIRMAP,
                               xytype='coordinate', recursionlimit=15000, snap='center')
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

    # No coordinate adjustment needed: the window grid is aligned to the
    # source raster, snapped coordinates are true cell centers, and the
    # polygonized boundary vertices fall on true pixel edges.
    return multi, lat_snap, lng_snap
