# delineator: Global Watershed Delineation with Python

Performs fast, accurate watershed delineation nearly anywhere on Earth.

  - packages state-of-the-art datasets from MERIT-Hydro and MERIT-Basins
  - near global coverage (excludes Greenland, Antarctica, and some small islands)
  - uses a hybrid of vector and raster methods for speed and accuracy

## Installation

Python 3.12 or higher is required. Quick install with `pip`:

```bash
pip install delineator
```

I highly recommended installing `delineator` in a fresh virtual environment.

Mac/Linux:
```
python3 -m venv venv
source venv/bin/activate
pip install delineator
```

Windows:
```aiignore
python -m venv venv
venv\Scripts\activate
pip install delineator
```

## Quick Start

After installing the package, use the `delineate` command-line utility to 
delineate a watershed. Open a terminal (or command prompt on Windows) and run:

```bash
delineate --point 63.938 -21.004  # Iceland's Ölfusá River at Route 1
```

This will create a single watershed for an outlet point with a given latitude
and longitude (note the order!) and saves it to `outlet/watershed.gpkg` in your 
current working directory.

To also produce rivers and outlet points, run:

```bash
delineate --point 63.938 -21.004 --rivers --outlets
```

Alternatively, you can use the `delineate()` function in your own 
Python scripts or notebooks.

```python
from delineator import delineate, write_outputs

# The basic delineate function returns three GeoDataFrames
# Note the order of latitude, longitude!
watershed_gdf, rivers_gdf, outlets_gdf = delineate(63.938, -21.004)

# Do whatever you wish with the resulting GeoDataFrames.
# This utility function will write them to disk in one line. 
write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, id="test")

```
You should get geodata for a watershed (polygon), rivers (poylines), and the 
watershed outlet (point). Here these data are displayed in QGIS:

![](docs/example_output.png)

## Datasets

The `delineator` package uses state-of-the art global hydrogrphy datasets from
[MERIT-Hydro](https://global-hydrodynamics.github.io/MERIT_Hydro/) and 
[MERIT-Basins](https://www.reachhydro.org/home/params/merit-basins) , 
and comes bundled with data files for Iceland. 

For allother regions, `delineator` will attempt to automaticallydownload the datasets
when they are needed. In some regions, the datasets are up to 3 GB, so this may take a while.

You may also download the datasets in advance so that they are staged and ready. 
These datasets have specially processed to work with the 
`delineator` package, and are available from the author's website.
Use the command-line utility `delineator download` to download the datasets for a 
given region. For example, to download the datasets for the Amazon Basin:

```bash
download_data --basin 62
```

The world's "megabasins" have IDs ranging from 11 to 86, as shown here (Greenland and Antarctica are 
excluded):

![](docs/megabasins.jpg)

By default, the datasets will be downloaded to your system's default data directory.
The usual locations are: 

  - Windows: `C:\Users\<username>\AppData\Local\delineator`.
  - Linux: `~/.local/share/delineator`.
  - Mac: `~/Library/Application Support/delineator`.

You can change this by setting the `data_dir` option. See below for more configuration
options.

# Command Line Examples

Delineate a watershed at a single outlet point:

```bash
delineate --point 63.938 -21.004
```

List the command line options:

```bash
delineate --help
```

Delineate watersheds for a set of outlet points in a CSV file:

```bash
delineate --csv outlets.csv
```

Create rivers and outlet points:
```bash
delineate --point 63.938 -21.004 --rivers --outlets
```

Output a different file format:
```bash
delineate --point 63.938 -21.004 --output-format geojson
delineate --point 63.938 -21.004 --output-format shp
delineate --point 63.938 -21.004 --output-format kml
delineate --point 63.938 -21.004 --output-format parquet
```

The `outlets.csv` file should contain, at minimum, the following columns:

- `id`
- `latitude`
- `longitude`

In this example, the fields `name` and `area` are not required and will be ignored.

```
id,lat,lng,name,area
6401070,64.71072,-21.60337,"Nordhura River at Stekkur",507
6401080,64.69229,-21.41046,"Hvita River at Kljafoss",1574
6401090,63.93796,-21.00666,"Olfusa River at Selfoss",5662
```

Change where the output files are written:

```bash
delineate --point 63.938 -21.004 --output-dir /path/to/output/dir`
delineate --csv outlets.csv --output-dir /path/to/output/dir`
```

## Configuration Options

The `delineate` function accepts a number of options to customize the output. A
limited set of options are available for the command line interface, and a more
complete set of options are available for the Python API.

To see a list of all options for the command line interface, run:

```bash
delineate --help
```

When using the Python API, options are set via the `DelineatorConfig` object, 
which is passed to the `delineate` function, as shown below.

```python
from delineator import DelineatorConfig, delineate

config = DelineatorConfig(
    data_dir="/path/to/data",
    rivers=False,  # the delineate function will not return rivers
    fill=False,    # do not fill in the gaps in the watersheds
    output_dir=r"C:\Users\user\Desktop\output"  # Set a custom output directory
    # Default output directory is the current working directory + "output"
)

watershed_gdf, rivers_gdf, outlets_gdf = delineate(63.938, -21.004, config)

# Update the config object to change options
config.rivers = True
watershed_gdf, _, outlets_gdf = delineate(63.938, -21.59, config)
```

The `DelineatorConfig` object has the following options:

```
auto_download: bool
    If True, the script will attempt to download needed data files for
    watershed delineation. You can also manually download the data files,
    see documentation for download_data().

clean : bool
    If True, the watershed boundary polygon will be "cleaned"
    which repairs artifacts such as seams that affect
    the appearance of the watershed polygon.

data_dir : Path or str, optional
    Directory where the script will look for data files and store downloaded
    data files needed for watershed delineation. If not provided,
    the script will use the default cache directory. For example, on Windows:
    C:\\Users\\<username>\\.delineator\\data
    Or on macOS/Linux: ~/.delineator/data

high_res : bool
    If True, the script will split the unit catchment around the outlet point
    for greater accuracy.
    If False, the script skips this step. The result will be a watershed that is
    slightly too large, and will include some area downstream of the outlet point.
    This error may be insignificant for very large watersheds.

fill : bool
    MERIT-Hydro data tends to produce watersheds with small "donut holes" due
    to minor topological errors in the input data. If True, these holes will be
    filled in, resulting in a cleaner appearance and smaller output files.

fill_threshold : int
    Only meaningful if fill=True. Holes smaller than this number of pixels
    (on the 3 arcsecond grid) will be filled. Larger holes are preserved.
    Set to 0 to fill all holes regardless of size. Default is 100 pixels.

low_res_threshold : float
    Watershed area in km² above which the script will automatically switch to
    lower-resolution mode, regardless of the high_res setting.
    Default is 80,000 km², meaning that the script will automatically switch to
    lower-resolution mode for watersheds over 80,000 km². Set it to a large
    number (6e6 or higher) to disable lower-resolution mode. (The largest
    basin in the dataset is the Amazon at 5.9 × 10⁶ km².)

num_stream_orders : int
    Number of stream orders to include in the river network. This will only
    have an effect if rivers=True.
    Default is 4. Higher values will result in more detailed stream network.
    Set this to a value of 9 or greater to get every river reach available.

rivers: bool
    If True, return the river flow lines in the watershed. These will be saved
    as the "rivers" layer in the output GeoPackage, or to a separate file
    such as "rivers.shp" or "rivers.geojson", depending on  `output_format`.

output_format: str
    Specify the output format for geodata files. Options are "gpkg" (default)
    or "geojson", "shp", "kml", and others supported by your installation of geopandas.
    See https://geopandas.org/en/stable/docs/user_guide/io.html#supported-drivers-file-formats
    for more details.

search_dist : float
    If the outlet point does not fall inside any catchment, the script will
    search for the nearest catchment within this distance (in decimal degrees).
    Default is 0.025°, which is roughly 2-3 km near the equator.

simplify : bool
    Used by the `write_outputs` utility function to simplify the output
    before saving it to disk. If True, the output will be simplified using
    Douglas-Peucker algorithm. This will result in a smaller file size,
    and also less of a jagged "staircase" appearance in the output.

simplify_tolerance: bool
    If `simplify=True`, this is the tolerance value for the simplification.
    With the Douglas-Peucker algorithm, the tolerance value is the maximum
    distance (in decimal degrees) that a vertex can be from the original.
```

Here is a more detailed description of a few of the options.

## Search Distance

Sometimes, your watershed outlet point will not fall inside one of the 
MERIT-Basins unit catchments. Your point could be just offshore in the ocean 
or an estuary, or it may happen to fall into one of the many small gaps and 
slivers in the source dataset. 

The variable `search_dist` controls how far away the script will look for the 
nearest catchment (in decimal degrees). Setting `search_dist=0` means that 
the point MUST fall inside a unit catchment to return results. I recommend 
setting it to at least 0.005 for reliable results, especially if your 
watersheds are near the coast.

## Filling holes in watershed polygons

Set `config(fill=True)` to fill in small gaps or "donut holes" in the watershed. 

Oftentimes, watersheds created by computer will contain internal gaps in the 
watershed polygon. Many of these come from the source data. MERIT-Basins has 
many small gaps and slivers between unit catchments, often only a few pixels wide. 
I recommend allowing the program to fill these in. To do this, in `config.py`, 
set `FILL = True`, and enter a non-zero value for `FILL_THRESHOLD`.

However, there are also larger donut holes inside a watershed, representing 
internal sinks, out of which water cannot flow. 
 How to handle these holes is somewhat of an outstanding question in hydrology, 
but my view is that most of them are unwanted, especially the smaller ones.

Larger gaps may be important, and you may not want to merge these with your 
watershed if they do not contribute to surface flow. As an example, here is 
a map of the Rio Grande watershed, with outlet coordinates at 26.05, -97.2.

![Rio Grande Watershed](docs/rio_grande.jpg)

Between the two main branches, the Rio Grande in the west and the Pecos River 
in the east, there is an endorheic basin that runs north-south for 560 km from 
New Mexico to Texas. Within this basin, there are several alkaline lakes or 
"playas," a tell-tale sign that water flows in and either seeps into the 
ground or evaporates. 

For most applications, it will be important to recognize that this area is 
"disconnected" and does not contribute to surface water flow at the basin outlet. 
On the other hand, if you are studying groundwater recharge, these areas may be important.

### Fill threshold

The script allows you to fill in small holes and keep big ones. 
If you have set `fill` to True, the configuration variable `fill_threshold` 
controls what size holes get filled in. 

The size threshold is roughly the number of pixels on a 3 arcsecond grid. 
In the source data, a pixels is a 0.000833° square. This is about 90m x 90m 
near the equator, or about 0.0081 km². The pixels get smaller in terms of 
surface area as you move north or south away from the equator.

For example, setting `fill_threshold=10` will fill in any donut holes with a 
width of 10 pixels or less and leave larger holes intact.

Setting `fill_threshold=0` will fill *all* donut holes, no matter their 
size. 

## Simplify Output

The output from the script may contain more detail than you need. 
To remove some vertices from the watershed boundary,
 set `simplify=True`. 
 
 You will also need to set `simplify_tolerance` to a value in decimal degrees. 
 This corresponds to the distance parameter in the 
 [Douglas-Peucker algorithm](https://cartography-playground.gitlab.io/playgrounds/douglas-peucker-algorithm/). 
 
Note that the vector polygons in the input data have been digitized from pixels 
with an edge length of 0.000833°. Setting `simplify_tolerance` to about half of this size or
higher will remove the jagged "staircase" appearance of the watershed boundary.

## The Importance of Reviewing All Output

Automated watershed delineation often makes mistakes. The good news is that 
these errors can often be fixed by slightly 
moving the location of your watershed outlet. But not always! 
I recommend carefully reviewing every watershed you delineate to ensure  
that it appears correct.

## Algorithm Details

This software has three aspects that make it fast and efficient in terms of 
memory use and processing time:

1. It uses a combination of raster and vector datasets. 
2. It accesses a set of pre-computed catchment geometries of various sizes, applying the concept of Hierarchical Spatial Aggregation. 
3. Geodata is stored in a relational sqlite database and accessed via queries.

For a detailed
description of the algorithm, see my manuscript, "[Fast, accurate watershed delineation 
with a hybrid of raster and vector methods](https://mghydro.com/pages/Heberger_delineation_2025.pdf)".

![](docs/method.jpg)

New in version 2.0.0 of `delineator` is the use of hierarchical aggregatrf catchments.
This new feature makes it much faster than previous versions for large watersheds. 
A large watershed like the Amazon or Mississippi contains tens of thousands of small
unit catchments. This is a well-known bottleneck in geospatial processing. 
In GIS, spatial operations--especially topological dissolves--are computationally expensive 
($O(n \log n)$ or worse, depending on vertex count. We solve this problem by 
pre-computing a set of nested catchments at 5 different sizes, which we 
refer to as Level 0 (for the smallest unit catchments) to Level 4. The goal is to fill the
"interior" of your delineated watershed with the largest possible pre-computed polygons, 
leaving only a thin outer shell of smaller unit catchments to be processed.

This map shows the nested catchments in the southern end of Madagascar. 

![](docs/nested_basins.jpg)

Finally, rather than using an old-fashioned format like shapefiles, vector geodata is
stored in sqlite, a modern relational database format. This allows us to use
the full power of SQL queries to efficiently access and manipulate the data.
This is *much* faster than using a traditional method of loading all the data into
memory and then doing processing on it. The datasets also have a spatial index, 
making it lightning fast to look up features by their location.

## Citation

I'm working on a paper on this project, but in the meantime, if you use this 
package in your research, please cite this GitHub repository:

```bibtex
@software{delineator,
  author = {Matthew Heberger},
  title = {delineator: Global Watershed Delineation with Python},
  year = {2026},
  publisher = {GitHub},
  version = {2.0.0},
  url = {https://github.com/mheberger/delineator}
}  
```

## Contributing

This project is open source and welcomes contributions. 
If you have comments or suggestions, please open an issue or pull request, 
or simply drop me an email.

