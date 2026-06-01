# Usage

For more control and flexibility, you can use the `delineator` package in your own 
Python scripts or notebooks.

```python
from delineator import delineate, write_outputs

# The basic delineate function returns three GeoDataFrames
watershed_gdf, rivers_gdf, outlets_gdf = delineate(63.938, -21.004)

# You can do whatever you want with the resulting GeoDataFrames.
# This utility function will write them to disk in one line. 
write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, id="test")

```