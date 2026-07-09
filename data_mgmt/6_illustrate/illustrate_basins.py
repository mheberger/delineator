"""
Map the catchments for each of the aggregation levels, L0 through L4.

Just for visualizing the data and verifying that it looks right.
Works OK but I found that QGIS is better for visualizing the geodata.
"""

import geopandas as gpd
from matplotlib import pyplot as plt


basins = [11,12,13,14,15,16,17,18,21,22,23,24,25,26,28,29,31,32,33,34,35,36,41,42,43,44,45,46,47,48,49,51,52,53,54,
          55,56,57,61,62,63,64,65,66,67,71,72,73,74,75,76,77,78,81,82,83,84,85,86]

# Create a figure with 5 subplots
fig, axes = plt.subplots(2, 3, figsize=(50, 30), layout="constrained")

for basin in [62]:
    print(basin)
    fname = rf"C:\Users\mheberger\Documents\watershed_app\static\data\basins{basin}.db"

    fig.suptitle(f"Basin {basin}")

    # Iterate over the tables in the database
    for i, ax in enumerate(axes.flatten()):
        ax.clear()
        if i == 5:
            continue
        table = f"l{i}_basins"
        try:
            gdf = gpd.read_file(fname, layer=table)
            gdf.plot(ax=ax, edgecolor='black', facecolor='lightblue', linewidth=0.2)
            ax.set_title(table)
        except Exception as e:
            print(e)
            continue

    plt.savefig(f"illustrate_basins/basin{basin}.jpg")

