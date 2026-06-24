from shapefile_check import check_shapefile

shp = r"C:\Data\GIS\MERITBasins\rivers\bugfix\riv_pfaf_11_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"

report = check_shapefile(shp)
print(report)
