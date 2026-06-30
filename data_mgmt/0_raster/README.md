# Merge MERIT-Hydro tiles 

Matthew Heberger
June 2022

To use MERIT-Hydro raster data for watershed delineation, we need to first
merge the 5° × 5° tiles to cover larger areas, or what I call "megabasins."

I previously experimented with the GRASS GIS "virtual raster" feature, 
but could not get it to work. 

I used QGIS to do some sample merges, then scripted the remaining
ones with Python. QGIS has a very helpful feature where you can 
get the Python equivalent of any processing step you've set up using
the user interface. In the command window, at the bottom left, click:
Advanced > Copy as Python command. You can then take this command, 
copy it to the Python console, and customize it, for example, insert
it in a loop to run repeatedly. 

The lookup table that shows which tiles belong to which basin are in
`basins_tiles.csv`. 

I created this by doing a Union in QGIS between two data layers:

- `wmobb_basins.shp` - World Meteorological Organization (WMO) Basins - 
MERIT follows the same organization. (Source)[https://www.bafg.de/GRDC/EN/02_srvcs/22_gslrs/223_WMO/wmo_regions_node.html]

- `merit_tiles_5_degree.shp` - I created this polygon layer in QGIS with the 
command Vector > Research Tools > Create Grid. 

I output the merged, megabasin-scale flow direction and flow accumulation
raster data as GeoTiff files, using two scripts. Note that these are
designed to be run from the Python console in QGIS. Access with: 
Plugins > Python console.

- qgis_merge_accum_rasters.py
- qgis_merge_flowdir_rasters.py

