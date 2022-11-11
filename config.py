"""
Configuration / Settings for delineator.py

Edit this file carefully before running delineator.py

See README for more information about the options and for instructions
on how to download the input data. 

"""

# Path to your CSV file with the watershed outlet data
OUTLETS_CSV = 'outlets_sample.csv'

# Set to True for "high-resolution" mode or False for "low-resolution."
HIGH_RES = True

# Directory containing the MERIT basin-scale flow direction rasters (.tif).
# Download from 
# For all paths, do not include a trailing slash.
MERIT_FDIR_DIR =  "data/raster/flowdir_basins"

# Directory containing the MERIT the flow accumulation rasters (.tif files).
MERIT_ACCUM_DIR = "data/raster/accum_basins"

# Set to True if you want the script to write status messages to the console
VERBOSE = True

# Set to True to make a bunch of plots of each watershed.
# (Just for debugging.)
PLOTS = False

# Folder where you have stored the Merit-BASINS catchment shapefiles.
HIGHRES_CATCHMENTS_DIR = "data/shp/merit_catchments"

# Folder where you have stored the MERIT-Basins rivers shapefiles
RIVERS_DIR = "data/shp/merit_rivers"

# Location of simplified catchment boundaries. Download the files from
# https://mghydro.org/watersheds/share/catchments_simplified.zip
LOWRES_CATCHMENTS_DIR = "data/shp/catchments_simplified"

# Folder where the script will write the output GeoJSON files or shapefiles
OUTPUT_DIR = "output"

# The file extension will determine the types of files the script creates.
# "geojson" for GeoJSON files (recommended) or "shp" for shapefiles
# Use a blank string "" if you don't want any output 
# (for example, you are only making the map)
OUTPUT_EXT = "shp"

# Set to True to ouput a summary of the delineation in OUTPUT.CSV
OUTPUT_CSV = True

# Threshold for watershed size in km² above which the script will revert to
# low-resolution mode 
LOW_RES_THRESHOLD = 50000

# If the requested watershed outlet is not inside a catchment, how far away 
# from the point should we look for the nearest catchment (in degrees)
SEARCH_DIST = 0.01

# Watersheds created with Merit-Hydro data tend to have many "donut holes"
# ranging from one or two pixels to much larger.
FILL = True

# If FILL = True, you many choose to to fill donut holes that are below a
# certain size. This is the number of pixels, on the 3 arcsecond grid.
FILL_THRESHOLD = 0

# Simplify the watershed boundary? This will remove some vertices from the watershed boundary and output smaller files. 
SIMPLIFY = False

# If SIMPLIFY is True, set SIMPLIFY_TOLERANCE to a value in decimal degrees. Note that the vector polygons
SIMPLIFY_TOLERANCE = 0.0008

# Set to TRUE if you want the script to create a local web page where you 
# can review the results
MAKE_MAP = True

# Folder where the script should put the map files.
# The mapping routine will make _viewer.html and .js files for every watershed
MAP_FOLDER = "map"

# On the map, do you also want to include the rivers?
MAP_RIVERS = True

# If you mapped the rivers, how many stream orders to include?
# I recommend 4. More than this and the browser may not display all the rivers in a large watershed.
NUM_STREAM_ORDERS = 4

# Set to True to use the experimental match areas feature. 
# You must include watershed areas in your outlets CSV file to use this feature. 
MATCH_AREAS = False

# If you set MATCH_AREAS = True, how close of a match should the script look for?
# Enter 0.25 for 25%. If you have not entered areas in your CSV file, you can ignore this parameter.
AREA_MATCHING_THRESHOLD = 0.35

# If you set MATCH_AREAS = True, how far away from the original outlet point should the script look 
# for a river reach that is a better match in terms of upstream area?
MAX_DIST = 0.15
