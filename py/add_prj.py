# Creates a .prj file defining the projection (WGs84 = CRS 4326) for all the shapefiles
# in a designated directory.

import os, glob

# This is the entire contents that needs to go in the prj file
prj_string = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]'

# Set the working directory
path = 'C:/Data/GIS/MERITBasins/catchments'

os.chdir(path)

# List all of the shapefiles.

file_pattern = '*.shp'
file_list = glob.glob(os.path.join(path, file_pattern))

for file in file_list:
    prj_file = file[:-3] + "prj"
    if not os.path.isfile(prj_file):
        with open(prj_file, 'w') as file:
            file.write(prj_string)