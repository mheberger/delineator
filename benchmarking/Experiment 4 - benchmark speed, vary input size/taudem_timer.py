"""
TauDEM timer for Benchmarking Experiment #4

Here we are delienating at one point in each megabasin, with the same upstream accum.

"""
import subprocess
import pandas
import geopandas
import pyproj
from shapely.geometry import Point, shape
import rasterio
from rasterio.features import shapes
import fiona
import os
import time
from sigfig import round
import numpy as np


wgs84 = pyproj.CRS('EPSG:4326')

INFILE = "expt4_outlets.csv"
OUTFILE= "taudem_results.csv"


def timer(func):
    """
    Decorator function to time another Python function
    Usage:
       @timer
       def myfunction()

    Returns:
        - time (in seconds), rounded to 3 sig figs
        - the result of the function it decorates.
    """

    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        execution_time = round(end_time - start_time, 3)
        print(f"Execution time of '{func.__name__}': {execution_time} seconds")
        return execution_time, result
    return wrapper


@timer
def make_gage_shapefile():
    """"
    Create a shapefile of of the watershed outlets from the .csv file
    """

    # My gage data is in a CSV file, but TauDEM wants a shapefile. 
    # Use GeoPandas to convert
    data = pandas.read_csv(INFILE)

    # Create a GeoDataFrame
    geometry = [Point(xy) for xy in zip(data['lng'], data['lat'])]
    geo_data = geopandas.GeoDataFrame(data, crs=wgs84, geometry=geometry)

    outfile = f"shp/outlets.shp"
    # Save the GeoDataFrame as a shapefile
    geo_data.to_file(outfile, driver='ESRI Shapefile')


def delineate() -> pd.DataFrame:
    """
    Delineates a set of watersheds using TauDEM
    expects there to be a shapefile with the outlets that have
    been snapped to streams.
    - Runs TauDEM `GageWatershed`, this outputs a .tif file.
    - Python function to convert this to vector data
    - includes a timer for benchmarking
    """

    infile = f'shp/outlets.shp'

    # Get the data for the set of gages.
    gages_gdf = geopandas.read_file(infile)
    n = len(gages_gdf)

    # Make a regular Pandas DataFrame to hold the outlets and benchmarking results
    output_df = gages_gdf.drop(columns='geometry').copy()
    output_df['time_delineate'] = np.nan
    output_df['time_vectorize'] = np.nan
    output_df['delineated_area'] = np.nan

    # Iterate over the gages
    for i in range(0, n):
        print(f"Delineating {i + 1} of {n}")

        # Create a shapefile with ONE point in it (what TauDEM wants)
        selected_row = gages_gdf.iloc[i]
        point_geometry = selected_row['geometry']
        single_point_gdf = geopandas.GeoDataFrame(geometry=[point_geometry])
        single_point_gdf.crs = gages_gdf.crs  # Assign the CRS from the original GeoDataFrame
        single_point_gdf.to_file("outlet.shp")  # Export the single point as a shapefile

        # Run TauDEM
        watershed_id = selected_row['id']
        megabasin = selected_row['basin']

        print("Running TauDEM")
        elapsed_time, result = run_taudem(megabasin)
        output_df.at[i, 'time_delineate'] = elapsed_time
        print("Vectorizing")
        elapsed_time, result = vectorize('watershed.tif', f'taudem_out/{watershed_id}.gpkg')
        output_df.at[i, 'time_vectorize'] = elapsed_time
        os.remove('watershed.tif')
        delete_shapefile('outlet.shp')

        # Finally, get the area of the delineated watershed, to add to results table
        watershed_gdf = geopandas.read_file(f'taudem_out/{watershed_id}.gpkg')
        area = get_area(watershed_gdf)
        print(area)
        output_df.at[i, 'delineated_area'] = area

    return output_df


def get_area(gdf: geopandas.GeoDataFrame) -> float:
    """
    Calculates the approximate area of all polygon features in a GeoDataFrame
    :param gdf: geopandas GeoDataFrame in WGS84 projection, in square kilometers
    :return: total area, km²
    """

    # Get the bounding box coordinates of the GeoDataFrame
    bbox = gdf.total_bounds  # Returns [minx, miny, maxx, maxy]

    # Define parameters for the custom Albers projection
    # Adjust these values based on the bounding box of your data
    # You may need to fine-tune these parameters for your specific area
    lon_0 = (bbox[0] + bbox[2]) / 2  # Center longitude
    lat_0 = (bbox[1] + bbox[3]) / 2  # Center latitude

    # Define the custom projection using Proj string
    custom_proj = (
        "+proj=aea "  # Albers Equal Area projection
        f"+lat_1={bbox[1]} "  # First standard parallel (minimum latitude)
        f"+lat_2={bbox[3]} "  # Second standard parallel (maximum latitude)
        f"+lat_0={lat_0} "  # Latitude of origin (center latitude)
        f"+lon_0={lon_0} "  # Central meridian (center longitude)
        "+x_0=0 +y_0=0"  # False Easting and Northing (usually 0)
    )

    # Create a custom CRS (Coordinate Reference System) using the Proj string
    custom_crs = pyproj.CRS.from_proj4(custom_proj)

    # Reproject the GeoDataFrame to the custom Albers projection
    gdf_custom_albers = gdf.to_crs(custom_crs)

    # Calculate the area in square kilometers
    total_area_km2 = gdf_custom_albers.geometry.area.sum() / 10 ** 6  # Convert square meters to square kilometers
    return round(total_area_km2, 6)


@timer
def run_taudem(megabasin) -> bool:
    """
    This will run TauDEM delineation for a single outlet point, in file `outlet.shp`
    and create the file `watershed.tif`.
    :return: True if successful, False if some kind of error
    """

    fdir = fr"C:\Data\GIS\MERITHydro\fdir_taudem\flowdir{megabasin}.tif"

    # TauDEM command to execute in the terminal; my laptop has 4 physical cores, so use max
    command = f"mpiexec -n 4 GageWatershed -p {fdir} -o outlet.shp -gw watershed.tif"
    return run_command(command)


def run_command(command: str) -> bool:
    """
    Runs a shell (terminal) command
    :param command: the command to run
    :return: True if it ran OK, False if some kind of error
    """
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)
    # Read the output and errors (if any)
    output, errors = process.communicate()
    # Print output and errors
    print("Output:")
    print(output)
    # Get the return code
    return_code = process.returncode
    if return_code != 0:
        print("Errors:")
        print(errors)
        print("Return Code:", return_code)
        return False
    else:
        return True


@timer
def vectorize(input_raster: str, output_geopackage: str) -> bool:
    """
    Converts a raster .tif file into a vector Geopackage
    :param input_raster:
    :param output_geopackage:
    :return: Boolean, indicating if it succeeded or not!
    """

    # Open the raster file
    with rasterio.open(input_raster) as src:
        # Read the raster data as numpy array
        image = src.read(1)  # Assuming you're working with the first band

        # Get metadata for the output GeoPackage file
        transform = src.transform

        # Generate vector shapes from raster data
        out_shapes = list(shapes(image, transform=transform))
        # out_shapes are tuples (geom, val).
        # Where val = -2147483647.0, that is the rest of the raster extent that is NOT the watershed.

    # Create output GeoPackage file
    schema = {'geometry': 'Polygon', 'properties': {'id': 'int'}}
    with fiona.open(output_geopackage, 'w', driver='GPKG', crs=wgs84, schema=schema) as dst:
        for geom, val in out_shapes:

            if geom is not None and val > 0:
                # Convert raster shapes to Shapely geometries
                geom = shape(geom)
                # Create a feature and write it to the output GeoPackage file
                feature = {'geometry': geom.__geo_interface__, 'properties': {'id': val}}
                dst.write(feature)

    return True


def delete_shapefile(shapefile_path):
    """
    Utility function to delete a shapefile, removing all of the various files
    :param shapefile_path:
    :return:
    """
    # Get the base filename without extension
    file_prefix = os.path.splitext(shapefile_path)[0]

    # List all files associated with the shapefile
    shapefile_extensions = ['.shp', '.shx', '.dbf', '.prj', '.cpg', '.sbn', '.sbx', '.shp.xml', '.qpj', '.shp.lock']

    # Delete all associated files
    for extension in shapefile_extensions:
        file_to_delete = file_prefix + extension
        if os.path.exists(file_to_delete):
            os.remove(file_to_delete)
            # print(f"Deleted {file_to_delete}")
        # else:
            # print(f"File {file_to_delete} does not exist.")


if __name__ == "__main__":

    make_gage_shapefile()
    results_df = delineate()
    results_df.to_csv(OUTFILE, index=False)
