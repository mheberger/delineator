"""
Watershed Delineation Validation Experiment with Delineator

Matthew Heberger, June 2026

This is for experiment 1-b, where I am trying to see whether `delineator`
is producing the same results as purely raster-based methods.

"""

from pathlib import Path

from delineator import delineate, DelineatorConfig
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize

example_raster = r'C:\Users\mheberger\Dropbox\RESEARCH\Watershed Article\benchmarking\Experiment 1 - validate watershed boundaries\expt1b\taudem_out\6401070_tau.tif'

def rasterize_to_template(gdf:gpd.GeoDataFrame, template_fname: str | Path,
                          out_fname: str | Path, burn_value: int | float = 1,
                          all_touched: bool = False):
    """
    Create a raster from vector data in a geodataframe.
    and give it the same dimensions and resolution as an existing raster.

    :param gdf: A geodataframe containing polygons or multipolygons.
    :param template_fname: str or Path to the example raster to use as a template
    :param out_fname: str or Path to the output raster
    :param burn_value: int or float, the value to burn into the raster
    :param all_touched: If True, all pixels touched by geometries will be burned in.
      If False, only pixels whose center is within the polygon will be burned in.
    :return: Nothing, but will raise an error if something goes wrong.
    """

    # Read the grid definition from the existing raster
    with rasterio.open(template_fname) as template:
        meta = template.meta.copy()
        transform = template.transform
        out_shape = (template.height, template.width)
        template_crs = template.crs

    # Make sure the vector is in the same CRS as the raster grid
    if gdf.crs != template_crs:
        gdf = gdf.to_crs(template_crs)

    # rasterize() wants an iterable of (geometry, value) pairs
    shapes = ((geom, burn_value) for geom in gdf.geometry)

    burned = rasterize(
        shapes=shapes,
        out_shape=out_shape,
        transform=transform,
        fill=0,            # value for cells outside any geometry
        all_touched=all_touched,
        dtype='uint8',
    )

    meta.update(dtype='uint8', count=1, nodata=0, compress='deflate')
    with rasterio.open(out_fname, 'w', **meta) as dst:
        dst.write(burned, 1)

def main():
    config = DelineatorConfig(
        rivers=False,
        threshold_multi=3000,
        calc_area=True,
        data_dir=r'C:\Users\mheberger\Documents\watershed_app\static\data',
        fill=False,
        snapping=False
    )

    df = pd.read_csv('outlets_snapped.csv', index_col=0)

    for id, row in df.iterrows():
        lat = row['lat']
        lng = row['lng']

        print(f'Delineating {id}')
        watershed_gdf, _, outlets_gdf = delineate(lat, lng, config)

        fname = f'output/{id}_del.geojson'
        outlets_fname = f'output/{id}_outlets.geojson'

        watershed_gdf.to_file(fname, driver='GeoJSON')
        outlets_gdf.to_file(outlets_fname, driver='GeoJSON')

        # Add the area of the delineated watershed to the results table
        area = float(watershed_gdf['area_km2'][0])
        df.at[id, 'area_delineator'] = area

        # Add the snapped coordinates to the results table (row 0 is requested point, row 1 is snapped point)
        lat_snapped = outlets_gdf.iloc[1].geometry.y
        lng_snapped = outlets_gdf.iloc[1].geometry.x

        df.at[id, 'lat_snapped'] = lat_snapped
        df.at[id, 'lng_snapped'] = lng_snapped

        # Export a GeoTiff, and make it the same dimensions and resolution as the example raster
        rasterize_to_template(watershed_gdf, example_raster, f'output/{id}_del.tif')

    # Add the snap distance in meters to the dataframe
    p1 = gpd.GeoSeries(gpd.points_from_xy(df['lng'], df['lat']),
                       index=df.index, crs='EPSG:4326').to_crs(3057)
    p2 = gpd.GeoSeries(gpd.points_from_xy(df['lng_snapped'], df['lat_snapped']),
                       index=df.index, crs='EPSG:4326').to_crs(3057)

    # This line is not adding the distance to my dataframe, but I don't know why
    df['dist_m'] = round(p1.distance(p2, align=False))

    df.to_csv('delineator_results.csv')


if __name__ == '__main__':
    main()
