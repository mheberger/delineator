# Vector Data Processing for `delineator`

Here is a description of the data management pipeline I used for the vector 
data for `delineator`. This is for the MERIT-Basins unit catchments and rivers

1. Download MERIT-Basins data

Get the shapefiles from: https://www.reachhydro.org/home/params/merit-basins
   There is a link on the site to a Google Drive folder. Get the "bugfix" files and the 
   coastal catchments:
   
   - MERIT-Basins hydrography (based on MERIT-Hydro v0.7/v1.0_bugfix1
     - rivers
     - unit catchments
   
   - Coastal Hillslopes (below 25 km2 channelization threshold)

The files are divided into "Pfafstetter Level 2" basins. I refere to these as
"megabasins." There are 61 megabasins, numbered 11 to 91. I am not carrying forward
basin 91, Greenland, as the file is huge and it does not describe realistic 
watersheds as Greenland is mostly ice-covered.

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

4. Create catchment boundaries and load to spatialite databases

The script `aggregate_basins.py` does a lot of heavy lifting here. For each
set of megabasin files, it creates a spatialite-enabled sqlite database 
`basins##.db`. Then we load in the unit catchment geometries into this database
as `l0_basins` for Level 0. Next, we create the larger catchment boundaries, 
Levels 1 to 4, calling the functions in `consolidate.py` which I developed for 
a previous project.

5. Load the rivers to a spatialite database with `make_sqlite_rivers.py`

6. Load the megabasins to a spatialite database with `make_megabasins_db.py`

7. Create extra database indices with `create_merit_sqlite_indexes.py`

This step was kind of a workaround to overcome a limitation with geopandas.
We can't set the field `comid` in the table `l0_basins` to be a primary key,
because then it won't get read in properly by geopandas. So we create a 
secondary index on the field for faster lookups.

8. Make the SHA-256 hashes for the files

These are needed for two reasons. First, they are used by the `pooch` package
to download the files and verify their integrity. Second, any user can use 
this information to verify the integrity of the files, i.e. there were no 
problems with the download and they have not been tampered with.

9. Create some basic visualizations with `illustrate_basins.py`

This was just to get a sense of the data and how the different aggregation 
levels look. Ultimately, not very satisfying; I preferred to make a map in 
QGIS, where I could adjust the lines, colors, etc.