"""
Watershed delineation validation experiment with `pysheds`

Matthew Heberger, June 2026

For the first experiment to validate `delineator`, this script performs watershed delineation
using the `pysheds` library, and compares the results to the results of the `delineator` library.

 - Reads a CSV file with a list of outlets (gaging stations in Iceland from the GRDC catalog)
 - Does the raster-based delineation with the Python pysheds library
 - Exports the results to a GeoPackage.
 - Calculates the area of the result for comparison with other methods.
 - Logs the time that it takes for each step and writes a CSV file with the results.

"""

import pyproj
from pysheds.grid import Grid
import geopandas as gpd
from shapely.geometry import shape
import numpy as np
import pandas
import geopandas


def get_area(gdf: geopandas.GeoDataFrame) -> float:
    """
    Calculates the approximate area of all polygon features
    in a GeoDataFrame.

    Projects to a coordinate system appropriate for Iceland

    Args:
        gdf: geopandas GeoDataFrame in WGS84 projection

    Returns:
        total area, km²
    """

    # Create a custom CRS (Coordinate Reference System) using the Proj string
    iceland_crs = pyproj.crs.CRS.from_string('EPSG:3057')

    # Reproject the GeoDataFrame to the custom Albers projection
    gdf_projected = gdf.to_crs(iceland_crs)

    # Calculate the area in square kilometers
    total_area_km2 = gdf_projected.geometry.area.sum() / 10 ** 6  # Convert square meters to square kilometers
    return round(total_area_km2, 3)


def pysheds_to_geodataframe(shapes, wid: str) -> gpd.GeoDataFrame:
    """
    Convert a pysheds delineation into a GeoDataFrame in WGS84 (EPSG:4326).

    :param shapes: pysheds generator yielding (geometry, value) tuples
    :param wid: the ID of the watershed, a string
    :return: a GeoDataFrame with LABEL and geometry columns, CRS EPSG:4326
    """
    records = [
        {"LABEL": value, "geometry": shape(geom)}
        for geom, value in shapes
    ]

    gdf = gpd.GeoDataFrame(records, geometry="geometry", crs="EPSG:4326")

    return gdf


def main():
    # This is the set of sample outlet locations
    megabasin = 27
    csv_fname = 'outlets_snapped.csv'

    # Get the data for the set of gages.
    gages_df = pandas.read_csv(csv_fname, header=0,
                               dtype={'id': 'str', 'lat': 'float', 'lng': 'float', 'area': 'float'})

    gages_df['delineated_area'] = np.nan

    n = len(gages_df)

    # Read the Raster data files.
    fdir_fname = f"C:/Data/GIS/MERITHydro/flow_dir_basins/flowdir{megabasin}.tif"

    # Get the stream centerlines from the flow accumulation grid
    print("Reading flow accumulation raster")
    grid = Grid.from_raster(fdir_fname)

    # Read the flow direction grid
    print("Reading flow direction raster")
    fdir = grid.read_raster(fdir_fname, nodata=0)

    for i in range(0, n):
        print(f'Delineating {i + 1} of {n}')
        wid = gages_df['id'].iloc[i]
        lat = gages_df['lat'].iloc[i]
        lng = gages_df['lng'].iloc[i]

        # Skip pour point snapping
        lat_snap = lat
        lng_snap = lng

        gages_df.at[i, 'lat_snapped'] = lat_snap
        gages_df.at[i, 'lng_snapped'] = lng_snap

        # Delineation with pysheds catchment() command
        catch = grid.catchment(x=lng_snap, y=lat_snap, fdir=fdir, xytype='coordinate', snap='center')

        # This step seems critical to getting a clean export
        grid.clip_to(catch)
        catch_view = grid.view(catch, dtype=np.uint8)

        # Create a vector representation of the catchment
        shapes = grid.polygonize(catch_view)

        # Export the catchment mask to a GeoTIFF file
        tiff_name = f'pysheds_out/{wid}_py.tif'
        grid.to_raster(catch_view, tiff_name, dtype=np.uint8, blockxsize=16, blockysize=16, compress='deflate')

        # Save the catchment as a GeoJSON file.
        gdf = pysheds_to_geodataframe(shapes, wid)
        gdf.to_file(f'pysheds_out/{wid}_py.geojson', driver='GeoJSON')

        delineated_area = get_area(gdf)
        gages_df.loc[i, 'area_pysheds'] = round(delineated_area)

        # Restore the old grid dimensions
        grid.viewfinder = fdir.viewfinder

    # Write the summary to a csv file
    gages_df.to_csv('pysheds_results.csv')


if __name__ == "__main__":
    main()
