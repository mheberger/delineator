"""
Script meant to be run in the QGIS Python console.

UNIONS a series of shapefiles to determine the percent overlap. 

Was meant to help me to compare the overlap of watersheds from different sources.

"""
import os
import csv
from time import sleep

def getGageIDs():
    #Import the data on gages and their coordinates
    fname = r"C:\Users\mheberger\Documents\delineator\paper\gages_sampled.csv"

    gage_ids =[]
    #gage_nums = []

    with open(fname, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=",")
        next(reader, None)  # skip the headers
        
        for row in reader:
            #gage_nums.append(int(row[0]))
            gage_ids.append(row[0])
        
    return gage_ids


os.chdir("C:/Users/mheberger/Documents/delineator/paper/temp")

gage_ids = getGageIDs()
n = len(gage_ids)
print(n)

# Let O be the area in square kilometers of the overlap
O = []

for i in range(0, n):
    gage_id = gage_ids[i]
    print(gage_id)
    fname1 = "C:/Users/mheberger/Documents/delineator/paper/watersheds/mine/gid_{}.shp".format(gage_id)
    fname2 = "C:/Users/mheberger/Documents/delineator/paper/watersheds/esri/gid_{}.shp".format(gage_id)

    outfile = "{}.shp".format(gage_id)

    # Intersect, create temp file, add area field, calculate area, save the value,
    processing.run("native:union", {
        'INPUT': fname1,
        'OVERLAY': fname2,
        'INPUT_FIELDS':['LABEL'],
        'OVERLAY_FIELDS':['Id'],
        'OVERLAY_FIELDS_PREFIX':'',
        'OUTPUT': outfile
    })

    sleep(2)

    #Open the layer with QGIS
    layer = QgsVectorLayer(outfile, '', 'ogr')

    fields = layer.fields()
    fields_names = fields.names()


    # Add the field shp_area to the attribute table
    if 'shp_area' not in fields_names:
        pv = layer.dataProvider()
        pv.addAttributes([QgsField('shp_area', QVariant.Double, len=10, prec=0)])
        layer.updateFields()

    layer.commitChanges()

    # Calculate the area
    expression1 = QgsExpression('$area/1e6')
    context = QgsExpressionContext()
    context.appendScopes(QgsExpressionContextUtils.globalProjectLayerScopes(layer))

    area = 0


    with edit(layer):
        for feature in layer.getFeatures():
            context.setFeature(feature)
            val = expression1.evaluate(context)
            area += val
            feature['shp_area'] = val
            layer.updateFeature(feature)
            #Also, add it to our dictionary of watershed areas (convert m2 to km2)

    O.append(area)

    # Delete the temp shapefile
    #files = ['intersection.dbf', 'intersection.shp', 'intersection.shx', 'intersection.prj']
    #for file in files:
    #    os.remove(file)

csv_file = "../union_areas.csv"

csv_columns = ('gage','union_area')
f = open(csv_file, 'w')

f.write('%s, %s\n' % csv_columns)
for i in range(0, n):
    f.write('{}, {}\n'.format(i+1, O[i]))

f.close()