"""
Repairs shapefiles

Matthew Heberger
2021-09-08

I got error messages when trying to do any kind of geoprocessing with them, like spatial joins
Something about invalid geometries or self-intersections. 
The advice on stackoverflow was that to fix it, you can do a zero-distance buffer. 
So I went ahead and did that for every shapefile.

This script is meant to be run in QGIS.


"""

in_path = "C:/Data/GIS/Princeton/Merit_Basins/pfaf_level_02/catchments"
out_path= "C:/Data/GIS/Princeton/Merit_Basins/pfaf_level_02/buffered"

i = 0
files = os.listdir(in_path)
n = len(files)/5

for file in files:
    if file.endswith('.shp'):
        i = i + 1
        
        infile = in_path + '/' + file
        outfile = out_path + '/' + file
        
        print("Processing %s, %s of n" % (i,n))
        
        processing.run("native:buffer", {'INPUT': infile,
            'DISTANCE':0,
            'SEGMENTS':5,
            'END_CAP_STYLE':0,
            'JOIN_STYLE':0,
            'MITER_LIMIT':2,
            'DISSOLVE':False,
            'OUTPUT': outfile})
        
print("Finished!")