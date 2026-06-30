"""
Import the list of rasters as externally-linked files.
"""

import csv

path = "C:/Data/GIS/MERITHydro/flow_direction"

#First step is to import the "mapping" data of Level2 basins to MERIT tiles

D = {}

current_basin = ""
tile_list = []

fname = "C:/Data/GIS/GRASS/basins_tiles.csv"
with open(fname, newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=",")
    for [basin, tile] in reader:
        if basin != current_basin:
            if current_basin != "":
                D[current_basin] = tile_list
            
            current_basin = basin
            tile_list = [ "{}/{}_dir.tif".format(path, tile) ]
        else:
            tile_list.append( "{}/{}_dir.tif".format(path, tile) )

# add the last row
D[current_basin] = tile_list
        
# Now, iterate through the dictionary, and use the information to build our merged rasters!

for basin, tiles in D.items():

     
    print("Now processing basin {}".format(basin))
    
    merged_raster_filename = "C:/Data/GIS/MERITHydro/flow_dir_basins/flowdir{}.tif".format(basin)
 
    processing.run("gdal:merge", {
        'INPUT': tiles,
        'PCT':False,
        'SEPARATE':False,
        'NODATA_OUTPUT':247,
        'OPTIONS':'COMPRESS=DEFLATE',
        'DATA_TYPE':0,
        'OUTPUT': merged_raster_filename})

print("Finished!")