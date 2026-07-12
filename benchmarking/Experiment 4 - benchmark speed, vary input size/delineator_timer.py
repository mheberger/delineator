"""
Benchmarking delineator

Matthew Heberger, June 2026

For experiment 4, where we vary the size of the raster, while holding
the watershed size roughly equal at 10 million pixels. 
"""
from pathlib import Path
import pandas
import time
from delineator import delineate, write_outputs, DelineatorConfig
import rasterio
import numpy as np
import geopandas as gpd
from rasterio.features import rasterize
from rasterio.transform import from_origin

OUTPUT_DIR = "delineator_out"


def rasterize_to_geotiff(gdf: gpd.GeoDataFrame, path, resolution: float =3 / 3600,
                         all_touched: bool=False, snap: bool=True):
    """
    Rasterize a GeoDataFrame to a uint8 0/1 mask and write it as a GeoTIFF
    in EPSG:4326.

    Parameters
    ----------
    gdf : GeoDataFrame
        Input geometries (any CRS; reprojected to 4326 if needed).
    path: Path or str
        Path of the output GeoTiff file.
    resolution : float
        Pixel size in degrees. Default 3 arcsec = 3/3600 = 1/1200.
    all_touched : bool
        False (default): burn a pixel only if its center falls inside a geometry.
        True: burn every pixel the geometry touches.
    snap : bool
        Snap bounds outward to a global grid aligned to multiples of
        `resolution`, so masks stay coincident with MERIT/HydroSHEDS tiles.

    Returns
    -------
    The path of the file that was successfully written.
    """
    if gdf.crs is None:
        raise ValueError("Input GeoDataFrame has no CRS.")
    if gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(4326)

    minx, miny, maxx, maxy = gdf.total_bounds
    if snap:
        minx = np.floor(minx / resolution) * resolution
        miny = np.floor(miny / resolution) * resolution
        maxx = np.ceil(maxx / resolution) * resolution
        maxy = np.ceil(maxy / resolution) * resolution

    width = int(round((maxx - minx) / resolution))
    height = int(round((maxy - miny) / resolution))
    transform = from_origin(minx, maxy, resolution, resolution)

    mask = rasterize(
        ((geom, 1) for geom in gdf.geometry if geom is not None),
        out_shape=(height, width),
        transform=transform,
        fill=0,
        all_touched=all_touched,
        dtype="uint8",
    )

    with rasterio.open(
        path, "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype="uint8",
        crs="EPSG:4326",
        transform=transform,
        nodata=None,        # 0 is a real class here, not nodata
        compress="deflate",
        predictor=1,        # keep at 1 for masks; 2 doesn't help binary data
        tiled=True,
        blockxsize=256,
        blockysize=256,
    ) as dst:
        dst.write(mask, 1)

    return path


# Read the watershed outlet data, corresponding to pixels where accum ~ 10E6 pixels
csv_file = "expt4_outlets.csv"
outlets_df = pandas.read_csv(csv_file, sep=",")

# Add columns for the benchmarked time (in seconds) and the watershed area
outlets_df['time'] = 0.0
outlets_df['area_calc'] = 0.0

# Set delineator options
config = DelineatorConfig(
    rivers=False,
    outlets=False,
    simplify=False,
    fill=False,
    clean=False,
    calc_area=True,
    round_coordinates=False,
    output_format="geojson",
    data_dir=r"C:\Users\mheberger\Documents\watershed_app\static\data"
)

# Warm up. First one is always slower due to cold start
delineate(63.938, -21.004, config)

for i in range(0, len(outlets_df)):
    # Use basin_id as variable name to avoid confusion, since id() is a Python function 
    basin_id = outlets_df['id'][i]
    lat = outlets_df['lat'][i]
    lng = outlets_df['lng'][i]

    print(f"Processing watershed {basin_id}")

    # Start the timer
    start_time = time.time()

    # Perform watershed delineation
    watershed_gdf, _, _ = delineate(lat, lng, config)

    # Save the watershed as a geopackage
    #write_outputs(watershed_gdf, id=basin_id)

    # End the timer and log the elapsed time in seconds
    end_time = time.time()
    duration = end_time - start_time
    outlets_df.loc[i, 'time'] = round(float(duration), 2)

    # Add the watershed area to the output table
    # Helps ensure we did not accidentally create spurious small watersheds
    # which are incorrect and would skew the results.
    area_calc = watershed_gdf['area_km2'].loc[0]
    outlets_df.loc[i, 'area_calc'] = area_calc

    output_fname = f"{OUTPUT_DIR}/{basin_id}.tif"
    rasterize_to_geotiff(watershed_gdf, output_fname)

# Save the results
#outlets_df.to_csv('delineator_results.csv')
