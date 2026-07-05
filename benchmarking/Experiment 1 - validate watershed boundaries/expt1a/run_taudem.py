"""
Created by Matthew Heberger, November 2023
Updated June 2026

This script does raster-based watershed delineation using TauDEM, which has been
the industry-standard software for delineation for the last decade, as far as I can tell.

Here I'm using MERIT-Hydro raster data for flow direction and flow accumulation. These
have already been created by the publisher of this dataset, no need to rebuild from the
DEM.

Before I ran this, I had to prepare the MERIT-Hydro data, because the flow direction
rasters used the ESRI standard, and needed to be converted to TauDEM standard.
I also needed to create a set of "stream centerline" rasters in order to use the
TauDEM "Move outlet to stream" routine.

 - Input is a CSV file with a list of gages.
 - First, convert this to a shapefile, required by TauDEM
 - Next, "snap pour point", or in TauDEM "Move outlet to stream"
 - Script outputs a shapefile with a single point in it.
 - We then run TauDEM with a terminal command.
 - The output from TauDEM is a raster (.tif) file.
 - TauDEM doesn't have any built-in methods for converting this to a vector format,
     which is what we are interested in.
 - Convert the raster to vector (GeoPackage instead of shapefile because it's 2023)
     using Python libraries rasterio, fiona, and shapely.
 - The result is "blocky" from vectorizing the raster data. If desired, one could
    run this through a simplify routine to reduce the number of edge vertices for
    faster display and smaller file size.

"""
import subprocess
import pandas
import geopandas
import pyproj
from shapely.geometry import Point, shape
import os
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.features import shapes

wgs84 = pyproj.CRS('EPSG:4326')
# I set up files with outlets in 3 basins: 27 (Iceland), 23 (Western Europe), and 62 (Amazon)
megabasin = 27


def make_gage_shapefile():
    """"
    Create a shapefile of my gage data in .CSV format.
    A shapefile is required by TauDEM's `MoveOutletstoStreams` routine
    """

    # My gage data is in a CSV file, but TauDEM wants a shapefile. Let's use GeoPandas to convert
    infile = f"outlets.csv"
    outfile = f"shp/outlets.shp"

    # Read the data into a Pandas DataFrame
    df = pandas.read_csv(infile)

    # Create a GeoDataFrame
    geometry = [Point(xy) for xy in zip(df['lng'], df['lat'])]
    gdf = geopandas.GeoDataFrame(df, crs=wgs84, geometry=geometry)

    # Save the GeoDataFrame as a shapefile
    gdf.to_file(outfile, driver='ESRI Shapefile')


def snap_gages():
    """
     Runs the TauDEM 'MoveOutletstoStreams' routine to snap the pour point
    to the stream centerline in the raster data
     - Expects there to be a shapefile with points in it.
     - Outputs a new shapefile with the *snapped* points

    :return:
    """
    outlets = f"shp/outlets.shp"
    outlets_snapped = f"shp/outlets_snapped.shp"
    fdir = fr"C:\Data\GIS\MERITHydro\fdir_taudem\flowdir{megabasin}.tif"
    streams = fr"C:\Data\GIS\MERITHydro\streams\streams{megabasin}.tif"

    # My laptop has 4 physical cores so use `mpiexec -4`
    command = f"mpiexec -n 4 MoveOutletstoStreams -p {fdir} -src {streams} -md 50 -o {outlets} -om {outlets_snapped}"
    print(command)
    run_command(command)


def delineate_taudem() -> pandas.DataFrame:
    """
    Delineates a set of watersheds using TauDEM
    expects there to be a shapefile with the outlets that have
    been snapped to streams.
    - Runs TauDEM GageWatershed, this outputs a .tif file.
    - Python function to convert this to vector data
    - includes a timer for benchmarking
    """

    infile = f'shp/outlets_snapped.shp'

    # Get the data for the set of gages.
    gages_gdf = geopandas.read_file(infile)

    # Add the fields `lat_snapped` and `lng_snapped` to the GeoDataFrame
    lat_snapped = gages_gdf['geometry'].apply(lambda geom: geom.y)
    lng_snapped = gages_gdf['geometry'].apply(lambda geom: geom.x)
    gages_gdf['lat_snapped'] = lat_snapped
    gages_gdf['lng_snapped'] = lng_snapped

    # TauDEM's `MoveOutletstoStreams` routine reports the distance the outlet was moved,
    # in terms of the number of grid cells (Manhattan distance), but I want to report
    # the distance in kilometers.
    p1 = gpd.GeoSeries(gpd.points_from_xy(gages_gdf['lng'], gages_gdf['lat']), crs='EPSG:4326').to_crs(3057)
    p2 = gpd.GeoSeries(gpd.points_from_xy(gages_gdf['lng_snapped'], gages_gdf['lat_snapped']), crs='EPSG:4326').to_crs(3057)

    gages_gdf['dist_m'] = round(p1.distance(p2, align=False))

    n = len(gages_gdf)

    # Make a regular Pandas DataFrame to hold the benchmarking data
    output_df = gages_gdf.drop(columns='geometry').copy()
    output_df['area_taudem'] = np.nan

    for i in range(0, n):
        print(f"Delineating {i + 1} of {n}")
        selected_row = gages_gdf.iloc[i]
        point_geometry = selected_row['geometry']
        single_point_gdf = geopandas.GeoDataFrame(geometry=[point_geometry])
        single_point_gdf.crs = gages_gdf.crs  # Assign the CRS from the original GeoDataFrame
        single_point_gdf.to_file("tmp/outlet.shp")  # Export the single point as a shapefile

        watershed_id = selected_row['id']
        print("Running TauDEM")
        run_taudem()

        print("Vectorizing")
        watershed_gdf = vectorize('tmp/watershed.tif')
        os.rename('tmp/watershed.tif', f'taudem_out/{watershed_id}.tif')

        delete_shapefile('tmp/outlet.shp')

        # Finally, get the area of the delineated watershed, to add to results table
        area = get_area(watershed_gdf)
        print(f"Watershed area: {area}")
        output_df.at[i, 'area_taudem'] = round(area)

        geojson_file = f'taudem_out/{watershed_id}.geojson'
        watershed_gdf.to_file(geojson_file, driver='GeoJSON')

    return output_df


def get_area(gdf: geopandas.GeoDataFrame) -> float:
    """
    Calculates the approximate area of all polygon features
    in a GeoDataFrame.

    Projects to a coordinate system appropriate for Iceland

    :param gdf: geopandas GeoDataFrame in WGS84 projection,
    in square kilometers
    :return: total area, km²
    """

    # Create a custom CRS (Coordinate Reference System) using the Proj string
    iceland_crs = pyproj.crs.CRS.from_string('EPSG:3057')

    # Reproject the GeoDataFrame to the custom Albers projection
    gdf_projected = gdf.to_crs(iceland_crs)

    # Calculate the area in square kilometers
    total_area_km2 = gdf_projected.geometry.area.sum() / 10 ** 6  # Convert square meters to square kilometers
    return round(total_area_km2, 6)


def run_taudem() -> bool:
    """
    This will run TauDEM delineation for a single outlet point, in file `outlet.shp`
    and create the file `watershed.tif`.
    :return: True if successful, False if some kind of error
    """

    fdir = fr"C:\Data\GIS\MERITHydro\fdir_taudem\flowdir{megabasin}.tif"

    # TauDEM command to execute in the terminal
    command = fr"mpiexec -n 4 GageWatershed -p {fdir} -o tmp\outlet.shp -gw tmp\watershed.tif"
    result = run_command(command)
    return result


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


def vectorize(input_raster: str) -> gpd.GeoDataFrame:
    """
    Converts a raster .tif file into a vector GeoDataFrame

    :param input_raster: Path to the input raster (.tif).
    :return: a GeoDataFrame containing the vectorized data.
    """
    # Read the raster and extract polygon shapes
    with rasterio.open(input_raster) as src:
        image = src.read(1)            # first band
        transform = src.transform

    # shapes() yields (geom, val) tuples. The background (everything outside
    # the watershed) has val == -2147483647.0, so keep only val > 0.
    records = [
        {"id": int(val), "geometry": shape(geom)}
        for geom, val in shapes(image, transform=transform)
        if geom is not None and val > 0
    ]

    if not records:
        return False  # nothing to write

    gdf = gpd.GeoDataFrame(records, geometry="geometry", crs="EPSG:4326")
    return gdf


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
    snap_gages()
    results_df = delineate_taudem()
    results_df.to_csv(f"taudem_results.csv", index=False)
