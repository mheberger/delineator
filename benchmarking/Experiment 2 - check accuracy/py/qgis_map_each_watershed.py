"""
Routine for making maps.
"""

# Import necessary modules
from qgis.core import QgsVectorLayer
import time

mypath = r'C:\Users\mheberger\Dropbox\RESEARCH\Watershed Article\benchmarking\Experiment 2 - check accuracy\watersheds'

canvas = iface.mapCanvas()

layer_names = ['usgs', 'esri', 'mine']

gages = ["06906500","09519800","07331600","09379500","06856600","03086000","08407500","06478000","07242000","07238000","05122000","02131000","01170500","08026000","07331000","07156900","06336000","06876000","06386500","07328500","03072655","07326500","01567000","02365200","01092000","06481000","08449000","03171000","08287000","05053000","04099000","07026040","01318500","03155000","06620000","05060000","06904500","05127000","05479000","07063000","02226500","02225500","02110500","06120500","08390800","03410500","01015800","06600500","05481000","03335000","01403060","01520500","07288500","07312200","03252500","04100500","09512500","05275000","05515500","07195400","07086000","05495000","06719505","08111700","06178000","08176900","05463000","03588500","08044000","04172000","08279000","09050700","03208500","01431500","02070500","07375300","07274000","05433000","03424730","03050000","02374500","07282000","01094000","01555500","09115500","03102500","03180500","09217900","06870300","03479000","01480870","06715000","02299950","01488500","06050000","01209901","09306242","01464907","05534500","06278300",
    "10353800",
    "12112600",
    "11208730",
    "14162200",
    "11414210",
    "11314500",
    "11118500",
    "11439500",
    "10297500",
    "13049500",
    "11451000",
    "12149000",
    "13186000",
    "13055000",
    "12452500",
    "11303000",
    "13152500",
    "12199000",
    "10068500",
    "13087995"
]

for gage in gages:

    for layer_name in layer_names:

        layer = QgsProject.instance().mapLayersByName(layer_name)[0]

        if layer:
            # Provide the new data source path
            new_data_source = f'{mypath}/{layer_name}/gid_{gage}.shp'

             # Get layer's symbology and transparency
            renderer = layer.renderer().clone()

            # Remove the current layer from the map
            QgsProject.instance().removeMapLayer(layer.id())

            # Load the layer again with the new data source
            new_layer = QgsVectorLayer(new_data_source, layer_name, 'ogr')

            if new_layer.isValid():
                # Apply the saved symbology and transparency to the new layer
                new_layer.setRenderer(renderer)

                # Add the updated layer to the map
                QgsProject.instance().addMapLayer(new_layer)
                
                # Get the current extent
                extent = new_layer.extent()
                # Add a buffer around the extent (e.g., 10% of the width/height)
                buffer_x = extent.width() * 0.1
                buffer_y = extent.height() * 0.1
                xmin = extent.xMinimum() - buffer_x
                xmax = extent.xMaximum() + buffer_x
                ymin = extent.yMinimum() - buffer_y
                ymax = extent.yMaximum() + buffer_y
                
                # Zoom to the extent of the USGS watershed layer
                if layer_name == 'usgs':
                    extent.setXMinimum(xmin)
                    extent.setXMaximum(xmax)
                    extent.setYMinimum(ymin)
                    extent.setYMaximum(ymax)
                    # Set the canvas extent with the buffered extent
                    canvas.setExtent(extent)
                    canvas.refresh()
                else:
                    # Subsequently, if the new layers are bigger, zoom out more
                    if xmin < extent.xMinimum():
                        extent.setXMinimum(xmin)
                    if xmax > extent.xMaximum():
                        extent.setXMaximum(xmax)
                    if ymin < extent.yMinimum():
                        extent.setYMinimum(ymin)
                    if ymax > extent.yMaximum():
                        extent.setYMaximum(xmax)
                
                print(f"Updated {layer_name}")
                
    canvas.refresh()
    
    canvas.waitWhileRendering()
    
                    
    # Save a screenshot!
    fname = f'C:/Users/mheberger/Desktop/{gage}.png'
    iface.mapCanvas().saveAsImage(fname)

