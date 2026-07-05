"""
Benchmark how long it takes to generate watersheds using delineator

Matthew Heberger
created Nov 2023
updated June 2026
"""

import pandas
import time
import logging

from delineator import delineate, write_outputs, DelineatorConfig

# Read the watershed outlet data corresponding to GRDC gages
csv_file = "expt3_outlets.csv"
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
    output_format="geojson",
    data_dir=r"C:\Users\mheberger\Documents\watershed_app\static\data"
)

# Warm up. First one is always slower due to cold start
delineate(63.938, -21.004, config)

for i in range(0, len(outlets_df)):
    # Use basin_id as variable name to avoid confusion, since id() is a Python function 
    basin_id = outlets_df['id'][i]
    #name = outlets_df['name'][i]
    lat = outlets_df['lat'][i]
    lng = outlets_df['lng'][i]

    print(f"Processing watershed {basin_id}")

    # Start the timer
    start_time = time.time()

    # Perform watershed delineation
    watershed_gdf, _, _ = delineate(lat, lng, config)
    write_outputs(watershed_gdf, id=basin_id, config=config)

    # End the timer and log the elapsed time in seconds
    end_time = time.time()
    duration = end_time - start_time
    outlets_df.loc[i, 'time'] = round(float(duration), 2)

    # Add the watershed area to the output table
    # Helps ensure we did not accidentally create a spurious small watershed
    area_calc = watershed_gdf['area_km2'].loc[0]
    outlets_df.loc[i, 'area_calc'] = area_calc

# Save the results of the benchmarking experiment
outlets_df.to_csv('delineator_times.csv')
