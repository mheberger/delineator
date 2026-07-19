"""
Performs detailed raster-based watershed delineation with `pysheds`,
but only inside of a *single* unit catchment, using GRIT rasters.

This is the GRIT counterpart of `merit_detailed.py`. GRIT distributes its
rasters differently from MERIT-Hydro, and several assumptions change:

- CRS: GRIT rasters are in EPSG:8857 (Equal Earth), 30 m pixels, while the
  rest of the pipeline works in EPSG:4326. This module projects the outlet
  point and the unit catchment polygon into EPSG:8857, does all raster work
  there, and projects the results back.
- Tiling: GRIT ships 300 km x 300 km tiles on a common aligned grid. Callers
  should pass a VRT mosaic over the tiles (see data_mgmt/build_grit_vrt.py);
  windowed reads through the VRT are seamless across tile boundaries, so
  this module never needs to know about tiles.
- Flow direction encoding: GRIT uses values 1-8 (see GRIT_DIRMAP below),
  not MERIT's ESRI-style powers of two.
- Pixel area: Equal Earth is an equal-area projection, so every pixel is
  exactly 900 m2 regardless of latitude. No geodesic pixel-area correction
  is needed when converting a drainage-area threshold to a pixel count.

Like `merit_detailed.py`, functions here log warnings and return None on
failure rather than raising.
"""
import os
import logging
import numpy as np
import rasterio
from numpy import floor, ceil
from pyproj import Transformer
from pysheds.grid import Grid
from shapely.geometry import Polygon, MultiPolygon, shape
from shapely.ops import transform as _shapely_transform

from delineator.constants import ADAPTIVE_THRESHOLD_FRACTION
from delineator.settings import DelineatorConfig

logger = logging.getLogger(__name__)

# The CRS of all GRIT raster products: Equal Earth (equal-area).
GRIT_CRS = "EPSG:8857"

# Nodata value in the GRIT drainage direction rasters. The value 0 also
# occurs and means "no flow" (sinks, ocean); both are treated as no-flow.
GRIT_NODATA = 255

# Area of one 30 m x 30 m pixel in km2. Exact everywhere, because Equal
# Earth preserves area.
GRIT_PIXEL_AREA_KM2 = 0.0009

# GRIT drainage direction encoding, from the GRIT documentation:
#
#     3 2 1
#     4 . 8
#     5 6 7
#
# i.e. 1=NE, 2=N, 3=NW, 4=W, 5=SW, 6=S, 7=SE, 8=E. Expressed in pysheds'
# dirmap order (N, NE, E, SE, S, SW, W, NW). Confirmed empirically: of all
# 16 cyclic 1-8 encodings, only this one concentrates flow accumulation
# into the Mississippi channel (see data_mgmt/build_grit_vrt.py notes).
GRIT_DIRMAP = (2, 1, 8, 7, 6, 5, 4, 3)

# Coordinate transformers, created once. always_xy gives (lng, lat) order.
_TO_GRIT = Transformer.from_crs("EPSG:4326", GRIT_CRS, always_xy=True)
_TO_WGS84 = Transformer.from_crs(GRIT_CRS, "EPSG:4326", always_xy=True)


def _snap_to_stream(streams: np.ndarray, affine, x: float, y: float,
                    max_dist_km: float) -> tuple[float | None, float | None, float | None]:
    """
    Snap a pour point to the nearest stream cell, in projected coordinates.

    The Equal Earth analog of merit_detailed._snap_to_stream. Coordinates
    here are EPSG:8857 meters, so plain Euclidean distance is used instead
    of haversine. Equal Earth preserves area, not distance, so the measured
    distance carries some scale distortion (growing toward high latitudes);
    over a snap radius of a few km this is an acceptable approximation.

    Parameters
    ----------
    streams : array-like of bool
        Boolean raster; True where flow accumulation exceeds the stream threshold.
    affine : affine.Affine
        The affine transform of the streams raster (EPSG:8857).
    x, y : float
        Pour point in EPSG:8857 coordinates.
    max_dist_km : float
        Maximum allowed snap distance in kilometers.

    Returns
    -------
    (x_snap, y_snap, dist_km) in EPSG:8857, or (None, None, None) if there
    are no stream cells or the nearest one is beyond max_dist_km. As in
    merit_detailed, the returned coordinates are the *center* of the snapped
    cell, for use with pysheds' snap='center'.
    """
    rows, cols = np.nonzero(streams)
    if rows.size == 0:
        logger.warning("No stream cells above the accumulation threshold in the unit catchment.")
        return None, None, None

    # Cell centers for the distance calculation
    x_ctr = affine.c + (cols + 0.5) * affine.a
    y_ctr = affine.f + (rows + 0.5) * affine.e

    dist_km = np.hypot(x_ctr - x, y_ctr - y) / 1000.0

    i = int(np.argmin(dist_km))
    d = float(dist_km[i])
    if d > max_dist_km:
        logger.warning(
            f"Nearest stream cell is {d:.2f} km from the pour point, "
            f"beyond the maximum snap distance of {max_dist_km:.2f} km."
        )
        return None, None, None

    return float(x_ctr[i]), float(y_ctr[i]), d


def _split_catchment(lat: float, lng: float, catchment_poly: Polygon | MultiPolygon,
                     home_catchment_area_km2: float, config: DelineatorConfig,
                     fdir_path: str | os.PathLike,
                     accum_path: str | os.PathLike | None = None) -> \
        tuple[Polygon | MultiPolygon | None, float | None, float | None]:
    """
    Performs the detailed pixel-scale raster-based delineation for a watershed,
    using GRIT drainage direction (and optionally upstream area) rasters.

    As in merit_detailed, raster work is confined to the single most
    downstream unit catchment, so even large watersheds stay fast.

    Parameters
    ----------
    lat : float
        latitude of the outlet (EPSG:4326)
    lng : float
        longitude of the outlet (EPSG:4326)
    catchment_poly : Shapely Polygon or MultiPolygon
        the home unit catchment, in EPSG:4326; used to clip the rasters
    home_catchment_area_km2 : float
        area of the home unit catchment in km2, used to adapt the stream
        threshold for small watersheds (see merit_detailed for details)
    config : DelineatorConfig
        a DelineatorConfig dataclass object
    fdir_path : str or PathLike
        path to the GRIT drainage direction raster: either a single tile or,
        normally, a VRT mosaic over the tiles (EPSG:8857)
    accum_path : str or PathLike, optional
        path to the GRIT drainage area ("mainstem") raster or, normally, a
        VRT mosaic over its tiles. Required when config.snapping is True.
        The product is float32 in km2 and sparse: values exist only on
        channel cells, NaN everywhere else. Build the mosaic with
        --extent-from the flow direction mosaic (see
        data_mgmt/build_grit_vrt.py) so both cover the same extent.

    Returns
    -------
    poly : shapely Polygon or MultiPolygon
        the part of the home unit catchment upstream of the outlet, EPSG:4326
    lat_snap : float
        latitude of the outlet after snapping (EPSG:4326)
    lng_snap : float
        longitude of the outlet after snapping (EPSG:4326)
    """

    fdir_path = str(fdir_path)
    if not os.path.isfile(fdir_path):
        logger.warning(f"Could not find flow direction raster: {fdir_path}")
        return None, None, None
    logger.info(f"Loading flow direction raster from: {fdir_path}")

    # Project the outlet and the unit catchment polygon into the raster CRS.
    x, y = _TO_GRIT.transform(lng, lat)
    poly_proj = _shapely_transform(_TO_GRIT.transform, catchment_poly)

    # Get a bounding box for the unit catchment, in EPSG:8857
    bounds = poly_proj.bounds

    # Align the bounding box to the source raster's own pixel grid, padded
    # by one pixel. Same reasoning as in merit_detailed: pysheds uses the
    # window bbox verbatim as the affine origin of the windowed grid, so the
    # bbox must fall exactly on source pixel edges or everything downstream
    # is misregistered by up to half a pixel. The arithmetic is pure affine
    # algebra and works the same in a projected CRS.
    with rasterio.open(fdir_path) as src:
        T = src.transform
        src_height, src_width = src.height, src.width

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

    bounding_box = (xmin_aligned, ymin_aligned, xmax_aligned, ymax_aligned)

    # Get the Grid object that matches the flow direction grid
    grid = Grid.from_raster(fdir_path, window=bounding_box, nodata=0)

    # Rasterize the unit catchment polygon into a mask on the windowed grid
    if isinstance(poly_proj, Polygon):
        poly_proj = MultiPolygon([poly_proj])
    mymask = grid.rasterize(list(poly_proj.geoms))

    # LOAD Flow Direction Grid
    fdir = grid.read_raster(fdir_path, window=bounding_box, nodata=0)

    # GRIT uses 255 for nodata and 0 for no-flow; fold both into no-flow,
    # and blank out everything outside the unit catchment polygon.
    fdir[fdir == GRIT_NODATA] = 0
    fdir[mymask == 0] = 0

    # Pour point snapping
    def snap_outlet() -> tuple[float | None, float | None]:
        """Snap the outlet to a stream cell. Works in EPSG:8857 coordinates."""
        logger.info("Snapping pour point")

        if accum_path is None:
            logger.warning("Snapping requires an upstream area raster, but accum_path is None.")
            return None, None
        if not os.path.isfile(str(accum_path)):
            logger.warning(f"Could not find upstream area raster: {accum_path}")
            return None, None

        acc = grid.read_raster(
            str(accum_path), window=bounding_box, window_crs=grid.crs, nodata=0
        )

        # GRIT's drainage-area raster holds km2 on channel cells and NaN
        # everywhere else; fold NaN to 0 so threshold comparisons behave.
        acc[~np.isfinite(acc)] = 0

        # MASK the accumulation raster to the unit catchment polygon
        acc[mymask == 0] = 0

        # Cap the stream threshold for small catchments, as in merit_detailed
        threshold_km2 = min(
            config.stream_threshold_km2,
            ADAPTIVE_THRESHOLD_FRACTION * home_catchment_area_km2,
        )
        logger.info(f"Using stream threshold of {threshold_km2:.4g} km2")

        # The raster holds drainage area in km2 directly (unlike MERIT's
        # pixel counts), so no pixel-area conversion is needed.
        streams = acc > threshold_km2

        max_dist_km = config.search_dist
        try:
            x_snap, y_snap, snap_dist_km = _snap_to_stream(
                streams, acc.affine, x, y, max_dist_km
            )
            if x_snap is None:
                return None, None
            logger.info(f"Snapped pour point moved {snap_dist_km:.3f} km.")
            return x_snap, y_snap

        except Exception as e:
            logger.warning(f"Could not snap the pour point. Error: {e}")
            return None, None

    # Snap the pour point. Both paths use cell-center coordinates with
    # pysheds' snap='center', as in merit_detailed.
    if config.snapping:
        x_snap, y_snap = snap_outlet()
        if x_snap is None:
            return None, None, None
    else:
        x_snap, y_snap = x, y

    # Raster-based watershed delineation with pysheds
    logger.info("Splitting home unit catchment")
    try:
        catch = grid.catchment(fdir=fdir, x=x_snap, y=y_snap, dirmap=GRIT_DIRMAP,
                               xytype='coordinate', recursionlimit=15000, snap='center')
        # Clip the bounding box to the catchment
        grid.clip_to(catch)
        clipped_catch = grid.view(catch, dtype=np.uint8)
    except Exception as e:
        logger.error(f"ERROR: something went wrong during pysheds grid.catchment(). Error: {e}")
        if config.snapping:
            lng_snap, lat_snap = _TO_WGS84.transform(x_snap, y_snap)
        else:
            lat_snap, lng_snap = lat, lng
        return None, lng_snap, lat_snap

    # Convert the raster subcatchment to polygons
    logger.info("Convert split catchment raster to polygon")
    shapes = grid.polygonize(clipped_catch)
    polygons = [shape(geom) for geom, value in shapes]
    multi = MultiPolygon(polygons)

    # Project the result and the snapped outlet back to EPSG:4326. The
    # polygon vertices fall on true pixel edges of the EPSG:8857 grid; the
    # reprojection preserves that georeferencing exactly (pyproj transforms
    # each vertex independently).
    multi_wgs84 = _shapely_transform(_TO_WGS84.transform, multi)
    if config.snapping:
        lng_snap, lat_snap = _TO_WGS84.transform(x_snap, y_snap)
    else:
        # Echo the user's coordinates back exactly, rather than round-tripping
        # them through EPSG:8857 (which introduces float noise).
        lat_snap, lng_snap = lat, lng

    return multi_wgs84, lat_snap, lng_snap
