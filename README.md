# delineator: a Python Package for Global Watershed Delineation

Fast, accurate watershed delineation nearly anywhere on Earth, using 
state-of-the-art datasets and a combination of vector and raster methods.

  - uses state-of-the-art datasets from MERIT-Hydro and MERIT-Basins
  - near global coverage (excludes Greenland, Antarctica, and some small islands)
  - uses a hybrid of vector and raster methods for speed and accuracy

## Installation

```
pip install delineator
```


## Quick Start

The `delineator` package comes bundled with sample data for Iceland. To test the package,
run the following to delineate a watershed in Iceldand (the Ölfusá River where it crosses
Route 1). You can use the `delineate` command from the command line and pass the 
latitude and longitude of a single *point*, the watershed outlet.

```bash
delineate --point 63.938 -21.004
```

For more control and flexibility, you can use the `delineator` package in your own 
Python scripts or notebooks.

```python
from delineator import delineate, write_outputs

# The basic delineate function returns three GeoDataFrames
watershed_gdf, rivers_gdf, outlets_gdf = delineate(63.938, -21.004)

# You can do whatever you want with the resulting geodata.
# This utility function will write them to disk in one line. 
write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, id="test")

```

## Read the Docs

For more options and information, see:

https://delineator.readthedocs.io
