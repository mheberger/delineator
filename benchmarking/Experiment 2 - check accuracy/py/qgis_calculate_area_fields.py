"""
Iterates over a bunch of shapefiles and 
(1) adds a new 'shp_area' field to the shapefile's attribute table, and 
(2) calculates the features' area and updates the new shp_area filed

also

(3) collects the areas for each shapefile 

Updated the script, because each shapefile may contain multiple features.
In previous versions, I guaranteed that shapefiles had only one single feature
by running a "dissolve all" routine. Here, it does not seem necessary however,
as there are no donut holes to take care of.

"""

csv_file = r'areas_usgs.csv'

inpath = r"C:\Users\mheberger\Documents\delineator\paper\watersheds\usgs"
inpath = inpath.replace("\\", "/")

# variables i and n just to keep track of our progress. 
# (but does not work well, since the QGIS script window does not update after the first few lines. )
#i = 0
n = len(files) / 4  #each shapefile has 4 individual files. 

files = os.listdir(inpath)

A = {}
i = 0

for file in files:
    if file.endswith('.shp'):
        i = i + 1
        
        #gage is the gage id
        gage = file[0:-4]
        
        infile = inpath + '/' + file
        
        print('Running %s of %s' % (i, n))
        
        #Open the layer. 
        layer = QgsVectorLayer(infile, '', 'ogr')
        
        #Check whether the field shp_area already exists. 
        idx = layer.fields().indexFromName('shp_area')
        #print(gage, idx)
        
        #If the field is missing, idx will be -1, and we will add a new field
        if idx == -1:
        
            # Add the field shp_area to the attribute table
            pv = layer.dataProvider()
            pv.addAttributes([QgsField('shp_area', QVariant.Double, len=10, prec=0)])
        
            # Calculate the area of 
            expression1 = QgsExpression('$area/1e6')
            context = QgsExpressionContext()
            context.appendScopes(QgsExpressionContextUtils.globalProjectLayerScopes(layer))
            
            area = 0
            
            with edit(layer):
                for f in layer.getFeatures():
                    context.setFeature(f)
                    val = expression1.evaluate(context)
                    area += val
                    f['shp_area'] = val
                    layer.updateFeature(f)
                    #Also, add it to our dictionary of watershed areas (convert m2 to km2)
                    A[gage] = area 
        
        # Case where the field shp_area already exists
        else:
            #Get the value from the field and save it
            idx = layer.fields().indexFromName('shp_area')
            area = 0
            for feature in layer.getFeatures():
                shp_area = feature.attributes()[idx]
                #Also, add it to our dictionary of watershed areas (convert m2 to km2)
                area += shp_area
                
            A[gage] = area  
            #print(gage, shp_area/1e6)
                
# Write the areas to a CSV file

csv_columns = ('gage','shp_area')
f = open(csv_file, 'w')
f.write('%s, %s\n' % csv_columns)
for k, v in A.items(): 
    f.write('%s, %s\n' % (k,v))
    
f.close()

print("Done!")
