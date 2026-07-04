import rasterio
import numpy as np

data_dir = r'C:\Users\mheberger\Documents\watershed_app\static\data'

output_file = 'fdir_metadata.csv'
o = open(output_file, 'w')
o.write("basin, width, height, total_pixels\n")
basins = [11,12,13,14,15,16,17,18,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,41,42,43,44,45,46,47,48,49,51,52,53,54,
          55,56,57,61,62,63,64,65,66,67,71,72,73,74,75,76,77,78,81,82,83,84,85,86]

for basin in basins:
    file_path = f'{data_dir}/flowdir{basin}.tif'
    with rasterio.open(file_path) as src:
        width = src.width
        height = src.height
        total_pixels = width * height

        o.write(f"{basin}, {width}, {height}, {total_pixels}\n")

o.close()
