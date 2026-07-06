"""
Merge MERIT-Hydro flow direction rasters into larger files covering
megabasins, or Level2 basins. 
"""

import csv

path = "C:/Data/GIS/MERITHydro/upstream_grid"

#First step is to import the "mapping" data of Level2 basins to MERIT tiles

# Let D be a dictionary mapping basin_id to a list of tiles
# e.g.: D[11] = ['n00e030_upg.tif', ..., 's20e040_upg.tif']

D = {}

# Helpers for constructing the dictionary from the CSV data
current_basin = ""
tile_list = []

fname = "C:/Data/GIS/GRASS/basins_tiles.csv"
with open(fname, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=",")
    for [basin, tile] in reader:
        if basin != current_basin:
            # start a new dictionary entry
            if current_basin != "":
                D[current_basin] = tile_list
            
            current_basin = basin
            tile_list = [ "{}/{}_upg.tif".format(path, tile) ]
        else:
            # expand the list for the current dictionary entry
            tile_list.append( "{}/{}_upg.tif".format(path, tile) )

# add the last row
D[current_basin] = tile_list
        
# Now, iterate through the dictionary, and use the information to build our merged rasters!

for basin, tiles in D.items():
    
    print("Now processing basin {}".format(basin))
    
    merged_raster_filename = "C:/Data/GIS/MERITHydro/accum_basins/accum{}.tif".format(basin)
 
    processing.run("gdal:merge", {
        'INPUT': tiles,
        'PCT':False,
        'SEPARATE':False,
        'NODATA_OUTPUT':0,
        'OPTIONS':'COMPRESS=DEFLATE',
        'DATA_TYPE':4,
        'OUTPUT': merged_raster_filename})

print("Finished!")