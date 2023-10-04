"""
Configuration / Settings for delineator.py

Edit this file carefully before running delineator.py

See README for more information about the options and for instructions
on how to download the input data. 

"""

# Path to your CSV file with the watershed outlet data
OUTLETS_CSV = 'outlets.csv'

# Set to True for "higher resolution" mode or False for "lower resolution."
HIGH_RES = True

# Directory containing the MERIT basin-scale flow direction rasters (.tif).
# Download from 
# For all paths, do not include a trailing slash.
MERIT_FDIR_DIR = "C:/Data/GIS/MERITHydro/flow_dir_basins"
#MERIT_FDIR_DIR = "data/raster/flowdir_basins"

# Directory containing the MERIT the flow accumulation rasters (.tif files).
MERIT_ACCUM_DIR = "C:/Data/GIS/MERITHydro/accum_basins"
#MERIT_ACCUM_DIR = "data/raster/accum_basins"

# Set to True if you want the script to write status messages to the console
VERBOSE = True

# Set to True to make a bunch of plots of each watershed.
# (Just for debugging.)
PLOTS = True

# Folder where you have stored the Merit-BASINS catchment shapefiles.
# These files need to be downloaded from: https://www.reachhydro.org/home/params/merit-basins
HIGHRES_CATCHMENTS_DIR = "C:/Data/GIS/MERITBasins/catchments"
#HIGHRES_CATCHMENTS_DIR = "data/shp/merit_catchments"

# Location of simplified catchment boundaries. Download the files from
# https://mghydro.org/watersheds/share/catchments_simplified.zip
LOWRES_CATCHMENTS_DIR = "C:/Data/GIS/MERITBasins/catchments"
#LOWRES_CATCHMENTS_DIR = "data/shp/catchments_simplified"

# Folder where you have stored the MERIT-Basins rivers shapefiles
RIVERS_DIR = "C:/Data/GIS/MERITBasins/rivers"

# Folder where the script will write the output GeoJSON files or shapefiles
OUTPUT_DIR = "output"

# The file extension will determine the types of files the script creates.
#   "geojson" for GeoJSON files
#   "shp" for shapefile
#   "gpkg" for GeoPackage (recommended)
# Use a blank string "" if you DO NOT want any output (for example,
# you are only making the interactive map).  Other file formats are available;
# see: https://geopandas.org/en/stable/docs/user_guide/io.html#writing-spatial-data
OUTPUT_EXT = "gpkg"

# Set to True to ouput a summary of the delineation in OUTPUT.CSV
OUTPUT_CSV = True

# Directory to store Python pickle files. Because it can be slow for Python to
# read shapefiles and create a GeoDataFrame. Once you have done this once, you
# can save time in the future by storing the GeoDataFrame as a .pkl file.
# The script will not search for pickle files if you leave this as a blank string, ''
PICKLE_DIR = 'pkl'

# Threshold for watershed size in km² above which the script will revert to
# low-resolution mode 
LOW_RES_THRESHOLD = 50000

# If the requested watershed outlet is not inside a catchment, how far away 
# from the point should we look for the nearest catchment (in degrees)
SEARCH_DIST = 0.025

# Watersheds created with Merit-Hydro data tend to have many "donut holes"
# ranging from one or two pixels to much larger.
FILL = True

# If FILL = True, you many choose to to fill donut holes that are below a
# certain size. This is the number of pixels, on the 3 arcsecond grid. 
# Set to 0 to fill ALL holes.
FILL_THRESHOLD = 100

# Simplify the watershed boundary? This will remove some vertices 
# from the watershed boundary and output smaller files.
SIMPLIFY = False

# If SIMPLIFY is True, set SIMPLIFY_TOLERANCE to a value in decimal degrees. 
# Note that the vector polygons
SIMPLIFY_TOLERANCE = 0.0008

# Set to TRUE if you want the script to create a local web page where you 
# can review the results
MAKE_MAP = True

# Folder where the script should put the map files. (MAKE sure it exists!)
# The mapping routine will make _viewer.html and .js files for every watershed
MAP_FOLDER = "map"

# On the map, do you also want to include the rivers?
MAP_RIVERS = True

# If you mapped the rivers, how many stream orders to include?
# I recommend 4 or 5. More than this and the browser may not display all the rivers in a large watershed.
NUM_STREAM_ORDERS = 3

# Set to True to use the experimental match areas feature. 
# You must include watershed areas in your outlets CSV file to use this feature. 
MATCH_AREAS = False

# If you set MATCH_AREAS = True, how close of a match should the script look for?
# Enter 0.25 for 25%. If you have not entered areas in your CSV file, you can ignore this parameter.
AREA_MATCHING_THRESHOLD = 0.25

# If you set MATCH_AREAS = True, how far away from the original outlet point should the script look 
# for a river reach that is a better match in terms of upstream area?
# Units are decimal degrees (not a proper distance measurement!)
MAX_DIST = 0.075
