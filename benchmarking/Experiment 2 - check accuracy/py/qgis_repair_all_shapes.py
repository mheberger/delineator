"""
Iterates over the shapefiles in a directory and repairs them. 
"""
    
inpath =  "C:/Users/mheberger/Documents/delineator/paper/watersheds/esri"
outpath = "C:/Users/mheberger/Documents/delineator/paper/watersheds/temp"


files = os.listdir(inpath)

for file in files:
    if file.endswith('.shp'):
     
        infile = inpath + '/' + file
        outfile = outpath + '/' + file
        
        
        # REPAIR
        processing.run("native:fixgeometries", {
            'INPUT': infile,
            'OUTPUT': outfile
        })
        
print("Finished!")