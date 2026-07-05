"""
Benchmarking experiment with `pysheds`
to compare how long it takes to do watershed delineation
with different methods, in order to compare them.

 - Reads a CSV file with a list of outlets (gaging stations)
 - Does the raster-based delineation with the Python pysheds library
 - Exports the results to a GeoPackage.
 - Calculates the area of the result for comparison with other methods.
 - Logs the time that it takes for each step and writes a CSV file with the results.

"""

import pyproj
from pysheds.grid import Grid
import matplotlib.pyplot as plt
import fiona
import numpy as np
import pandas
import time
from shapely.geometry import Polygon
import shapely.ops
import shapely.geometry


def timer(func):
    """
    Decorator function to time another Python function
    Usage:
       @timer
       def myfunction()

    timer(myfunction)
    """
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = round(end_time - start_time, 3)
        print(f"Execution time of '{func.__name__}': {execution_time} seconds")
        return execution_time, result
    return wrapper


def main():
    # This is the set of sample outlet locations
    megabasin = 23
    csv_fname = f'outlets{megabasin}.csv'

    # Get the data for the set of gages.
    gages_df = pandas.read_csv(csv_fname, header=0,
                               dtype={'id': 'str', 'lat': 'float', 'lng': 'float', 'area': 'float'})

    # We'll record the time in this field
    gages_df['time'] = np.nan
    gages_df['delineated_area'] = np.nan

    n = len(gages_df)

    # Read the Raster data files. Here, all of the watersheds are in the same Level 2 megabasin
    accum_fname = f'C:/Data/GIS/MERITHydro/accum_basins/accum{megabasin}.tif'
    fdir_fname = f"C:/Data/GIS/MERITHydro/flow_dir_basins/flowdir{megabasin}.tif"
    # Get the stream centerlines from the flow accumulation grid
    print("Reading flow accumulation raster")
    grid = Grid.from_raster(accum_fname)
    elapsed_time, acc = read_accum(accum_fname, grid)
    NUMPIXELS = 5000
    streams = acc > NUMPIXELS
    # Read the flow direction grid
    print("Reading flow direction raster")
    elapsed_time, fdir = read_fdir(fdir_fname, grid)

    for i in range(0, n):
        print(f'***Delineating {i + 1} of {n}')
        wid = gages_df['id'].iloc[i]
        y = gages_df['lat'].iloc[i]
        x = gages_df['lng'].iloc[i]

        # Begin timer
        start_time = time.time()

        # Snap pour point
        x_snap, y_snap = grid.snap_to_mask(streams, (x, y))

        # Delineation with pysheds catchment() command
        catch = grid.catchment(x=x_snap, y=y_snap, fdir=fdir, xytype='coordinate')

        # This step seems critical to getting a clean export
        grid.clip_to(catch)
        catch_view = grid.view(catch, dtype=np.uint8)

        # Draw a plot of the delineated watershed (uncomment for debugging)
        # plot_catchment(catch_view, extent=grid.extent, wid=wid)

        # Create a vector representation of the catchment
        shapes = grid.polygonize(catch_view)
        # shapes is a "generator" ... iterable, but not exactly a Python list.

        # Log the elapsed time
        end_time = time.time()
        duration = end_time - start_time
        gages_df.loc[i, 'delineation_time'] = round(duration, 2)
        print(f'pysheds delineation time: {duration}')

        # Save the catchment as a GeoJSON file.
        print('  Exporting')
        elapsed_time, result = export(shapes, wid, crs=grid.crs.srs)
        gages_df.loc[i, 'export_time'] = elapsed_time

        # Convert to a shapely polygon so we can calculate area
        shapes = grid.polygonize(catch_view)  # bizarrely, generator works one time only! have to recreate
        catchment_polygon = shapely.ops.unary_union([shapely.geometry.shape(shape)
                                                     for shape, value in shapes])

        delineated_area = get_area(catchment_polygon)
        gages_df.loc[i, 'delineated_area'] = delineated_area

        # Restore the old grid dimensions
        grid.viewfinder = fdir.viewfinder

    # Write the timer data to a csv file
    gages_df.to_csv(f'pysheds_times{megabasin}.csv')


@timer
def read_fdir(fdir_fname, grid):
    return grid.read_raster(fdir_fname, nodata=0)


@timer
def read_accum(accum_fname, grid):
    return grid.read_raster(accum_fname)


def get_area(poly: Polygon) -> float:
    """
    Projects a Shapely polygon in raw lat, lng coordinates, and calculates its
    (approximate) area.

    Args:
        poly: Shapely polygon
    :return: area in km²
    """
    wgs84 = pyproj.CRS('EPSG:4326')
    albers = pyproj.Proj(proj='aea', lat_1=poly.bounds[1], lat_2=poly.bounds[3])
    project = pyproj.Transformer.from_proj(wgs84, albers, always_xy=True).transform

    projected_poly = shapely.ops.transform(project, poly)

    # Get the area in km^2
    return projected_poly.area / 1e6


@timer
def export(shapes, wid: str, crs: str):
    """
    Export a delineated watershed as a GeoJSON file
    :param shapes: pysheds generator
    :param wid: the ID of the watershed, an string
    :param crs: string like 'epsg4326'
    :return: nothing
    """
    schema = {'geometry': 'Polygon', 'properties': {'LABEL': 'float:16'}}

    # Write geopackage
    fname = f"{wid}.gpkg"
    with fiona.open(fname, 'w', crs=crs, driver="GPKG",
                    schema=schema) as c:
        i = 0
        for shape, value in shapes:
            rec = {'geometry': shape, 'properties': {'LABEL': str(value)}, 'id': str(i)}
            c.write(rec)
            i += 1

    return True


def plot_catchment(catch_view, extent, wid: str):
    """
    Draws a map of the delineated pysheds catchment
    :param catch_view: a pysheds object.
    :param extent:
    :param wid: the ID of the watershed
    :return: nothing
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_alpha(0)

    ax.imshow(np.where(catch_view, catch_view, np.nan), extent=extent,
              zorder=1, cmap='Greys_r')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.grid(visible=True, zorder=0)
    plt.title(label='Delineated Catchment', size=14)
    # plt.show()
    fname = f"plots/{wid}.jpg"
    plt.savefig(fname)
    plt.close(fig)


if __name__ == "__main__":
    main()
