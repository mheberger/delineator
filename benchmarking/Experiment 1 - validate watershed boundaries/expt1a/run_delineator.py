from pathlib import Path

from delineator import delineate, DelineatorConfig
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize

example_raster = 'taudem_out/6401070_tau.tif'
OUTPUT_DIR = "delineator_out"


def rasterize_to_template(gdf:gpd.GeoDataFrame, template_fname: str | Path,
                          out_fname: str | Path, burn_value: int | float = 1,
                          all_touched: bool = False):
    """
    Create a raster from vector data in a geodataframe and give it the same
    dimensions and resolution as an existing raster.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        A geodataframe containing polygons or multipolygons.
    template_fname : str or Path
        Path to the example raster to use as a template.
    out_fname : str or Path
        Path to the output raster.
    burn_value : int or float, default 1
        The value to burn into the raster.
    all_touched : bool, default False
        If True, all pixels touched by geometries will be burned in.
        If False, only pixels whose center is within the polygon will be burned in.

    Returns
    -------
    None
        Nothing, but will raise an error if something goes wrong.
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
    """
    Validation experiment 1a.

    Uses delineator Python package to delineate watersheds for 20 outlets
    in Iceland. Captures snapped coordinates from the delineation process,
     and calculates the distance between the input and snapped coordinates.
    Additionally, it generates two outputs: GeoJSON files (vector)
    and GeoTiff files (raster), each containing the delineated watersheds,
    and a final CSV file summarizing results.
    """
    config = DelineatorConfig(
        rivers=False,
        fill=False,
        round_coordinates=False,
        threshold_multi=3000,
        calc_area=True,
        data_dir=r'C:\Users\mheberger\Documents\watershed_app\static\data'
    )

    # Read the outlets from the CSV file to a Pandas DataFrame
    # We'll also use this DataFrame to store the results of the delineation process
    df = pd.read_csv('outlets.csv', index_col=0)

    for id, row in df.iterrows():
        lat = row['lat']
        lng = row['lng']

        print(f'Delineating {id}')
        watershed_gdf, _, outlets_gdf = delineate(lat, lng, config)

        fname = f'{OUTPUT_DIR}/{id}_del.geojson'

        watershed_gdf.to_file(fname, driver='GeoJSON')

        # Add the area of the delineated watershed to the results table
        area = float(watershed_gdf['area_km2'][0])
        df.at[id, 'area_delineator'] = area

        # Add the snapped coordinates to the DataFrame
        lat_snapped = outlets_gdf.iloc[1].geometry.y
        lng_snapped = outlets_gdf.iloc[1].geometry.x

        df.at[id, 'lat_snapped'] = lat_snapped
        df.at[id, 'lng_snapped'] = lng_snapped

        # Export a GeoTiff, and make it the same dimensions and resolution as the example raster
        rasterize_to_template(watershed_gdf, example_raster, f"{OUTPUT_DIR}/{id}_del.tif")

    # Add the snap distance in meters to the dataframe
    requested_point = gpd.GeoSeries(gpd.points_from_xy(df['lng'], df['lat']),
                                    index=df.index, crs='EPSG:4326').to_crs(3057)
    snapped_point = gpd.GeoSeries(gpd.points_from_xy(df['lng_snapped'], df['lat_snapped']),
                       index=df.index, crs='EPSG:4326').to_crs(3057)

    # Adding the distance to the DataFrame
    df['dist_m'] = round(requested_point.distance(snapped_point, align=False))

    # Save the results to a CSV file
    df.to_csv('delineator_results.csv')


if __name__ == '__main__':
    main()
