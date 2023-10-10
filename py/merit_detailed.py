"""
Performs a detailed, raster-based delineation, but only inside of a *single* unit catchment.
This is my implementation of the hybrid method, which I "invented" but which was actually
first described by Djokic and Ye (1999). Doing the raster math is slow
and takes a lot of memory. So we only do the bare minimum that we have to, and use the
vector data that has already been processed for most of the upstream watershed.
"""
import copy  # TODO: undo this modification?

from pysheds.pgrid import Grid as Grid
from shapely.geometry import Polygon, MultiPolygon
from shapely import wkb, ops
from config import *
import numpy as np
import os
from config import PLOTS
from py.raster_plots import *


def split_catchment(wid: str, basin: int, lat: float, lng: float, catchment_poly: Polygon,
                                 bSingleCatchment: bool) -> (object or None, float, float):
    """
    Performs the detailed pixel-scale raster-based delineation for a watershed.

    To efficiently delineate large watersheds, we only use raster-based methods in a small area,
    the size of a single unit catchment that is most downstream. This results in big
    savings in processing time and memory use, making it possible to delineate even large watersheds
    on a laptop computer.

    Args:
        wid: the watershed id, a string
        basin: 2-digit Pfafstetter code for the level 2 basin we're in (tells us what files to open)
        lat: latitude
        lng: longitude
        catchment_poly: a Shapely polygon; we'll use it to clip the flow accumulation raster to get an accurate snap
        bSingleCatchment: is the watershed small, i.e. there is only one unit catchment in it?
            If so, we'll use a lower snap threshold for the outlet.

    Returns:
        poly: a shapely polygon representing the part of the terminal unit catchment that is upstream of the
            outlet point
        lat_snap:  latitude of the outlet, snapped to the river centerline in the accumulation raster
        lng_snap: longitude of the outlet, snapped to the river centerline in the accumulation raster

    For HIGH PRECISION, we will discard the most downstream catchment
    polygon, and replace it with a more detailed delineation.
    For this, we will use a small piece of the flow direction raster
    that fully contains the terminal unit catchment's boundaries,
    and do a raster-based delineation.
    This finds all of the pixels in the unit catchments that contribute flow to our
    pour point. Then, we will convert that collection of pixels to a polygon,
    and merge it with the upstream unit catchment polygons.

    Read the flow direction raster, but use "Windowed Reading",
    where we only read in the portion
    of interest, surrounding our catchment. See:
    https://mattbartos.com/pysheds/file-io.html
    There is no need to read the whole file into memory, because we are only
    interested in delineating a watershed within our little most-downstream unit catchment.
    Upstream of that, we used vector-based data
    For the window, get the bounding box for our catchment
    a tuple with 4 floats: (Left, Bottom, Right, Top)
    Note pysheds lets you read in a rectangular portion (and not a portion based
    on polygon geometry. We will clip the accumulation raster with the unit catchment
    polygon in a separate step below. (Not a clip per se, but we will replace the values
    in cells that are outside the unit catchment with NaN. This way, these cells will be
    ignored during the "snap pour point" routine, and we will only find rivers that are
    inside our unit catchment. It took me a while to figure out this workflow, but it is
    the key to getting accurate results!
    """

    # Get a bounding box for the unit catchment
    bounds = catchment_poly.bounds
    bounds_list = [float(i) for i in bounds]

    # The coordinates of the bounding box edges that we get from the above query
    # do not correspond well with the edges of the grid pixels.
    # We need to round them to the nearest whole pixel and then
    # adjust them by a half-pixel width to get good results in pysheds.

    # Distance of a half-pixel
    halfpix = 0.000416667

    # Bounding box is xmin, ymin, xmax, ymax
    # round the elements DOWN, DOWN, UP, UP
    # The number 1200 is because the MERIT-Hydro rasters have 3 arsecond resolution, or 1/1200 of a decimal degree.
    # So we just multiply it by 1200, round up or down to the nearest whole number, then divide by 1200
    # to put it back in its regular units of decimal degrees. Then, since pysheds wants the *center*
    # of the pixel, not its edge, add or subtract a half-pixel width as appropriate.
    # This took me a while to figure out but was essential to getting results that look correct
    bounds_list[0] = np.floor(bounds_list[0] * 1200) / 1200 - halfpix
    bounds_list[1] = np.floor(bounds_list[1] * 1200) / 1200 - halfpix
    bounds_list[2] = np.ceil( bounds_list[2] * 1200) / 1200 + halfpix
    bounds_list[3] = np.ceil( bounds_list[3] * 1200) / 1200 + halfpix

    # The bounding box needs to be a tuple for pysheds.
    bounding_box = tuple(bounds_list)

    # Open the flow direction raster *using windowed reading mode*
    fdir_fname = "{}/flowdir{}.tif".format(MERIT_FDIR_DIR, basin)
    if VERBOSE: print("Loading flow direction raster from: {}".format(fdir_fname))
    if VERBOSE: print(" using windowed reading mode with bounding_box = {}".format(repr(bounding_box)))

    if not os.path.isfile(fdir_fname):
        raise Exception("Could not find flow flow direction raster: {}".format(fdir_fname))

    # The pysheds documentation was not up-to-date. Seems there were some changes in the API
    # for the versions with and without the numba library (sgrid and pgrid)
    # The first line did not work, but the following does. Took ages to figure this out! :(
    # I think it had to do with when the developer added the ability to use numba, the code forked.
    # You can still use it without numba, but the code is older and has not evolved with the new stuff (?)
    # Anyhow, the old version worked better for me in my testing.
    # grid = Grid.from_raster(path=fdir_fname, data=fdir_fname, data_name="myflowdir", window=bounding_box,nodata=0)
    grid = Grid.from_raster(fdir_fname, data_name="fdir", window=bounding_box, nodata=0)

    # Now "clip" the rectangular flow direction grid even further so that it ONLY contains data
    # inside the bounaries of the terminal unit catchment.
    # This prevents us from accidentally snapping the pour point to a neighboring watershed.
    # This was especially a problem around confluences, but this step seems to fix it.
    # (Seems I had to first convert it to hex format to get this to work...)
    hexpoly = catchment_poly.wkb_hex
    poly = wkb.loads(hexpoly, hex=True)
    # coerce this into a single-part polygon, in case the geometry is a MultiPolygon
    poly = get_largest(poly)

    # Fix any holes in the polygon by taking the exterior coordinates.
    # One of the annoyances of working with GeoPandas and pysheds is that you have
    # to constantly switch back and forth between Polygons and MultiPolygons...
    filled_poly = Polygon(poly.exterior.coords)

    # It needs to be of type MultiPolygon to work with rasterio apparently
    multi_poly = MultiPolygon([filled_poly])

    # Convert the polygon into a pixelized raster "mask".
    mymask = grid.rasterize(multi_poly)
    grid.add_gridded_data(mymask, data_name="mymask", affine=grid.affine, crs=grid.crs, shape=grid.shape)

    # I believe this step was unnecessary, but it makes the plots look a little nicer
    m, n = grid.shape
    for i in range(0, m):
        for j in range(0, n):
            if int(grid.mymask[i, j]) == 0:
                grid.fdir[i, j] = 0

    # Plot the mask that I created from rasterized vector polygon
    if PLOTS:
        plot_mask(grid, catchment_poly, lat, lng, wid)

    # MERIT-Hydro flow direction uses the old ESRI standard for flow direction...
    dirmap = (64, 128, 1, 2, 4, 8, 16, 32)

    # Plot the flow-direction raster, for debugging
    if PLOTS:
        plot_flowdir(grid, lat, lng, wid, dirmap, catchment_poly)

    # Making a copy of the grid before clipping to split catchment boundary, for plots
    if PLOTS:
        grid2 = copy.deepcopy(grid)

    if VERBOSE: print("Snapping pour point")

    # Open the accumulation raster, again using windowed reading mode.
    accum_fname = '{}/accum{}.tif'.format(MERIT_ACCUM_DIR, basin)
    if not os.path.isfile(accum_fname):
        raise Exception("Could not find accumulation raster: {}".format(accum_fname))

    grid.read_raster(accum_fname, data_name="acc", window=bounding_box, window_crs=grid.crs)

    # Clips the flow direction grid to a new rectangular bounding box.
    # that corresponds to the mask of the unit catchment.
    grid.clip_to("mymask", inplace=True, apply_mask=True)

    # MASK the accumulation raster to the unit catchment POLYGON. Set any pixel that is not
    # in 'mymask' to zero. That way, the pour point will always snap to a grid cell that is
    # inside our polygon for the unit catchment, and will not accidentally snap
    # to a neighboring watershed. It took me a bunch of experimenting to realize
    # that this is the key to getting good results in small watersheds, especially
    # when there are other streams nearby.
    # The approach I used (looping over every pixel in the grid) is simple but a little hackish.
    # Would be better implemented as a method in pysheds.
    m, n = grid.shape
    for i in range(0, m):
        for j in range(0, n):
            if int(grid.mymask[i, j]) == 0:
                grid.acc[i, j] = 0

    # Snap the outlet to the nearest stream. This function depends entirely on the threshold
    # that you set for how minimum number of upstream pixels to define a waterway.
    # If the user is looking for a small headwater stream, we can use a small number.
    # In this case, there will be only one unit catchment in the watershed.
    # In most other circumstances (num unit catchments > 1), a much larger value gives better results.
    # The values here work OK, but I did not test very extensively...
    # Using a minimum value like 500 prevents the script from finding little tiny watersheds.
    if bSingleCatchment:
        numpixels = 500
    else:
        # Case where there are 2 or more unit catchments in the watershed
        # setting this value too low causes incorrect results and weird topology problems in the output
        numpixels = 5000

    if VERBOSE: print("Using threshold of {} for number of upstream pixels.".format(numpixels))

    # Snap the pour point to a point on the accumulation grid where accum (# of upstream pixels)
    # is greater than our threshold
    streams = grid.acc > numpixels
    xy = (lng, lat)
    try:
        [lng_snap, lat_snap], dist = grid.snap_to_mask(streams, xy)
    except Exception as e:
        if VERBOSE: print(f"Could not snap the pour point. Error: {e}")
        return None, None, None

    # Plot the accumulation grid, for debugging
    if PLOTS:
        plot_accum(grid, lat, lng, lat_snap, lng_snap, wid, catchment_poly)

    # Plot the streams!
    if PLOTS:
        plot_streams(grid, streams, catchment_poly, lat, lng, lat_snap, lng_snap, wid, numpixels)

    # Finally, here is the raster based watershed delineation with pysheds!
    if VERBOSE: print("Delineating catchment")
    try:
        grid.catchment(data='fdir',
                       x=lng_snap,
                       y=lat_snap,
                       dirmap=dirmap,
                       xytype='label',
                       out_name='catch',
                       recursionlimit=15000,
                       nodata_out=0)

        # Clip the bounding box to the catchment
        # Seems optional, but turns out this line is essential.
        grid.clip_to('catch')
        clipped_catch = grid.view('catch', dtype=np.uint8)
    except Exception as e:
        if VERBOSE: print(f"ERROR: something went wrong during pysheds grid.catchment(). Error: {e}")
        return None, lng_snap, lat_snap

    # Convert high-precision raster subcatchment to a polygon using pysheds method .polygonize()
    if VERBOSE: print("Converting to polygon")
    shapes = grid.polygonize()


    # The output (from pysheds without numba) is creating MANY shapes.
    # Dissolve them together with the unary union operation in shapely
    # (Could not install numba on my web server, but performance does
    # not seem to suffer too much, presumably because I'm using small
    # little clipped sections of the rasters.)
    # I believe that the reason that it sometimes prodcuces many polygons is
    # because MERIT-Hydro flow-direction grids, while being a very nice dataset,
    # can produce polygons with dangles, and with donut holes.
    # The solution "discard all but the largest polygon" is not ideal from a
    # theoretical standpoint, because we could discard a piece of the drainage
    # area, but in my testing, it usually only results in the loss of a few pixels
    # here and there, in other words, trivial differences. The tradeoff seemed
    # worthwhile -- we lose a bit of accuracy, but in exchange, we gain the simplicity
    # of working with Polygon geometries, rather than MultiPolygons.

    shapely_polygons = []

    shape_count = 0

    # The snapped vertices look better if we nudge them one half pixels
    lng_snap += halfpix
    lat_snap -= halfpix

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

        if result_polygon.geom_type == "MultiPolygon":
            if PLOTS:
                polygons = list(result_polygon)
                plot_polys(polygons, wid)

            result_polygon = get_largest(result_polygon)
    else:
        # If pysheds generated a single polygon, that is our answer
        result_polygon = shapely_polygons[0]

    if PLOTS:
        plot_catchment(grid, grid2, catchment_poly, result_polygon, lat, lng, lat_snap, lng_snap, wid, dirmap)
        plot_clipped(grid, clipped_catch, catchment_poly, lat, lng, lat_snap, lng_snap, wid, result_polygon)


    return result_polygon, lat_snap, lng_snap


def get_largest(input_poly: MultiPolygon or Polygon) -> Polygon:
    """
    Converts a Shapely MultiPolygon to a Shapely Polygon
    For multipart polygons, will only keep the largest polygon
    in terms of area. In my testing, this was usually good enough

    Note: can also do this via PostGIS query... see myqueries_merit.py, query19a and query19b
          Not sure one approach is better than the other. They both seem to work well.
    Args:
        input_poly: A Shapely Polygon or MultiPolygon

    Returns:
        a shapely Polygon
    """
    if input_poly.geom_type == "MultiPolygon":
        areas = []
        polygons = list(input_poly)

        for poly in polygons:
            areas.append(poly.area)

        max_index = areas.index(max(areas))

        return polygons[max_index]
    else:
        return input_poly


