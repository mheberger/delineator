"""
Experiment 2 comarisons

Script to compare the results of 2 competing methods of watershed delineation
using the Intersection over Union (IOU) method to compare 2 different polygons.
Here I am comparing watersheds created with `delineator` and with the 
ESRI ArcGIS Online "Create Watersheds" tool. 

"""

import os
import geopandas as gpd
import pandas as pd
from shapely.ops import transform
import pyproj


def calculate_iou(feature, gpkgfile):

    ECKERT_IV = "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

    gdf1 = gpd.GeoDataFrame([feature], geometry='geometry', crs='epsg:4326').to_crs(ECKERT_IV)
    gdf2 = gpd.read_file(gpkgfile).to_crs(ECKERT_IV)

    intersection = gpd.overlay(gdf1, gdf2, how='intersection')
    intersection_area = intersection.geometry.area.sum()

    union = gpd.overlay(gdf1, gdf2, how='union')
    union_area = union.geometry.area.sum()

    #iou = intersection_area / union_area
    total_area1 = gdf1.geometry.area.sum()
    total_area2 = gdf2.geometry.area.sum()

    return total_area1 / 1e6, total_area2 / 1e6, intersection_area / 1e6, union_area / 1e6  # Convert to square kilometers


def main():

    results = []

    shapefile = 'C:/Data/Discharge/Observed/Caravan/shapefiles/lamah/lamah_basin_shapes.shp'
    gpkg_folder = 'C:/Users/mheberger/Documents/delineator/output'

    results = []

    shapefile_gdf = gpd.read_file(shapefile)

    for index, feature in shapefile_gdf.iterrows():
        gauge_id = feature['gauge_id']
        gpkg_filename = f"{gauge_id}.gpkg"
        gpkg_path = os.path.join(gpkg_folder, gpkg_filename)

        if os.path.exists(gpkg_path):
            area1, area2, intersection, union,= calculate_iou(feature, gpkg_path)
            print(gauge_id, area1, area2, intersection, union)
            results.append({'gauge_id': gauge_id, 'area1': area1,
                            'area2': area2,  'intersection': intersection, 'union': union})
        else:
            print(f"Skipped: {gauge_id} - Corresponding GeoPackage not found.")



    results_df = pd.DataFrame(results)

    print(results_df)

    results_df.to_csv("comparison.csv")


if __name__ == '__main__':
    main()