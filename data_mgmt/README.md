# Vector Data Processing for `delineator`

Here is a description of the data management pipeline I used for the vector 
data for `delineator`. This is for the MERIT-Basins unit catchments and rivers

1. Download 

Get the shapefiles from: https://www.reachhydro.org/home/params/merit-basins
   There is a link to a Google Drive folder. Get the "bugfix" files and the 
   coastal catchments:
   
   - MERIT-Basins hydrography (based on MERIT-Hydro v0.7/v1.0_bugfix1
   
   - Coastal Hillslopes (below 25 km2 channelization threshold)

The files are divided into "Pfafstetter Level 2" basins. I refere to these as
"megabasins." There are 61 megabasins, numbered 11 to 91. 

2. Add projection files

Some of the downloaded shapefiles from MERIT-Basins are missing `.prj` files. 
Run the sript `add_prj.py` to add these where needed. All the files are in
CRS 4326, the standard geographic coordinate reference system using latitude
and longitude coordinates.

3. Fix Geometries

Based on my previous experience working with these datasets, I recall that they
have some issues, the main one being with "bow-tie" polygons. I verified that 
this is still an issue with the 'bugfix' shapefiles with this script, created
with the help of Claude.ai: `shapefile_check.py` 

A reliable way to fix these is with QGIS's "Fix Geometries" tool. I used a 
Python script to batch process the shapefiles: `qgis_fix_shapefiles.py`

After this step, the files still have some issues with ring orientation, 
where the polyon boundary is counter-clockwise when it's supposed to be clockwise,
but this is an ESRI convention, and does not matter for our use case, which 
is to put the data into a spatialite-enabled sql database. 

4. Load the unit catchments to a spatialite database


5. Load the rivers to a spatialite database

The rivers geodata is where we have the

