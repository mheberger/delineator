"""
Converts the MERIT-Hydro flow direction rasters (basin scale)
to the TauDEM standard, using Raster Reclassify by Table

Matthew Heberger, Dec. 2023

This script is meant to be run in the Python window in QGIS.

"""

import os

os.chdir("C:/Data/GIS/MERITHydro")

basins = [11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 41, 42, 43, 44, 
          45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 56, 57, 61, 62, 63, 64, 65, 66, 67, 71, 72, 73, 74, 75, 76, 77, 78, 
          81, 82, 83, 84, 85, 86, 91]

n = len(basins)

for basin in basins:
    infile  = "flow_dir_basins/flowdir%s.tif" % basin
    tmp_file = "fdir_taudem/temp.tif"
    print(infile)
    
    # Run the RECLASSIFY routine. Even if I provide a COMPRESS argument, it does not do anything!
    # As a result, the output from this routine is huge, often 3GB+
    processing.run("native:reclassifybytable", {
        'INPUT_RASTER': infile,
        'RASTER_BAND':1,
        'TABLE':['1','1','1','2','2','8','4','4','7','8','8','6','16','16','5','32','32','4','64','64','3','128','128','2','0','0','0','247','247','-32768'],
        'NO_DATA':-9999,
        'RANGE_BOUNDARIES':2,
        'NODATA_FOR_MISSING':False,
        'DATA_TYPE':1,
        'OUTPUT': tmp_file,
    })
    
    outfile = "fdir_taudem/flowdir%s.tif" % basin
    
    # So, write to a temporary file and use GDAL translate to do compression. 
    # The resulting files are much smaller.
    processing.run("gdal:translate", {
        'INPUT':tmp_file,
        'TARGET_CRS':None,
        'NODATA':None,
        'COPY_SUBDATASETS':False,
        'OPTIONS':'COMPRESS=DEFLATE',
        'EXTRA':'',
        'DATA_TYPE': 1,
        'OUTPUT': outfile
    })
        
    os.remove(tmp_file)

print("Finished!")
