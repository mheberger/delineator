# delineator: Global Watershed Delineation with Python

Fast, accurate watershed delineation for any point on Earth's land surface, 
using a hybrid of vector- and raster-based methods with data from 
[MERIT-Hydro](https://global-hydrodynamics.github.io/MERIT_Hydro/) 
and [MERIT-Basins](https://www.reachhydro.org/home/params/merit-basins).

- Near-global coverage (excludes Greenland, Antarctica, and some small islands)
- Bundled sample data for Iceland; other regions download automatically on first use
- Returns watershed polygon, river network, and outlet points as GeoPandas GeoDataFrames

## Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Command line reference](#command-line-reference)
- [Configuration reference](#configuration-reference)
- [Data files](#data-files)
- [⚠️ Always review your results](#review-results)
- [Usage examples](#usage-examples)
- [Algorithm](#algorithm)
- [Citation](#citation)
- [Contributing](#contributing)


<a id="installation"></a>
## Installation

Requires Python ≥ 3.10. Python 3.11+ is recommended for speed. Also
recommended to install in a fresh virtual environment to avoid 
dependency conflicts.

**macOS/Linux:**
```bash
python3 -m venv venv
source venv/bin/activate
pip install delineator
```

**Windows:**
```bash
python -m venv venv
venv\Scripts\activate
pip install delineator
```

<a id="quick-start"></a>
## Quick start

The bundled Iceland data lets you run immediately after installation; no separate 
download required.

**Command line usage**

```bash
delineate --point 63.938 -21.004
```

This creates the watershed for the Ölfusá River at Route 1 in Iceland.
Output is written to `./output/watershed.gpkg` in your current directory. 
To create geodata for the river network and outlet points, run:

```bash
delineate --point 63.938 -21.004 --rivers --outlets
```

**Python script usage**

Alternatively, you can use the `delineate()` function in your 
own Python scripts or notebooks.

```python
from delineator import delineate, write_outputs

# The delineate function returns three GeoDataFrames
# Note the order of latitude, longitude!
watershed_gdf, rivers_gdf, outlets_gdf = delineate(63.938, -21.004)

# Do whatever you wish with the resulting GeoDataFrames.
# This utility function will write them to disk in one line. 
write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, id="olfusa")
```

Here is an example of the output displayed in QGIS:

![Example output](docs/example_output.png)

<a id="command-line-reference"></a>
## Command line reference

```bash
# Single point
delineate --point 63.938 -21.004

# Include rivers and outlet points
delineate --point 63.938 -21.004 --rivers --outlets

# Output different file formats
delineate --point 63.938 -21.004 --output-format geojson
delineate --point 63.938 -21.004 --output-format shp
delineate --point 63.938 -21.004 --output-format kml
delineate --point 63.938 -21.004 --output-format parquet

# Batch delineation of multiple outlet points in a CSV file
delineate --csv outlets.csv

# Custom output directory
delineate --csv outlets.csv --output-dir /path/to/output/

# List all the command line options
delineate --help
```

For batch delineation, the CSV file must contain at minimum `id`, `lat`, and `lon` columns.
Other columns are OK but will be ignored by the script. Example CSV file:

```
id,lat,lon,name
6401070,64.71072,-21.60337,Nordhura River at Stekkur
6401080,64.69229,-21.41046,Hvita River at Kljafoss
6401090,63.93796,-21.00666,Olfusa River at Selfoss
```

### Output files

When `--output-format gpkg` (the default), all layers are written to a single file 
(`watershed_<id>.gpkg`) with three layers: `watershed`, `rivers`, and `outlets`. 

For other formats like `shp`, each layer is written to a separate file, for example 
`rivers.shp`, `outlets.shp`, and `watershed.shp`.

### Environment variables

Instead of passing options to the command line, you can set environment variables
for the default data directory and the output director. There are three 
environment variables:

- `DELINEATOR_DATA_DIR`: directory where input data files are saved
- `DELINEATOR_OUTPUT_DIR`: directory where output files will be saved
- `DELINEATOR_AUTO_DOWNLOAD`: whether to automatically download data files as 
  they are needed

Environment variables add are useful when you want configuration 
that is global, repeatable, automatable, or sensitive, 
without forcing every CLI call or Python function call to spell 
everything out.

Environment variables work with the command-line interface or with 
the Python functions (`delineate()`, `downloader()`). Note that 
command line arguments will override environment variables, as will 
 the `DelineatorConfig` object passed to `delineate()`.

Set the three available environment variables as follows:

Mac/Linux:

```bash
export DELINEATOR_DATA_DIR=/mnt/data/delineator
export DELINEATOR_OUTPUT_DIR =/home/user/documents/watersheds
export DELINEATOR_AUTO_DOWNLOAD=false
delineator --csv outlets.csv
```

Windows CMD:

```cmd
set DELINEATOR_DATA_DIR=D:\Data\delineator
set DELINEATOR_OUTPUT_DIR=C:\Users\user\Documents\watersheds
set DELINEATOR_AUTO_DOWNLOAD=false
delineator --csv outlets.csv
```

Windows Powershell:
```powershell
$env:DELINEATOR_DATA_DIR = "D:\Data\delineator"
$env:DELINEATOR_OUTPUT_DIR = "C:\Users\user\Documents\watersheds"
$env:DELINEATOR_AUTO_DOWNLOAD = "false"
delineator --csv outlets.csv
```

<a id="configuration-reference"></a>
## Configuration reference

When using the Python function `delineate()`, options are passed via a 
`DelineatorConfig` object:

```python
from delineator import delineate, DelineatorConfig

config = DelineatorConfig(
    high_res=True,
    rivers=True,
    fill=True,
    output_format="gpkg",
    output_dir="/path/to/output",
)

watershed_gdf, rivers_gdf, outlets_gdf = delineate(63.938, -21.004, config)

# Config objects are mutable - update and reuse
config.rivers = False
config.outlets = False
config.output_format = "geojson"
watershed_gdf, _, _ = delineate(63.938, -21.59, config)
```

All options with their defaults:

| Option               | Default        | Description                                                                                                                                                                  |
|----------------------|----------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `auto_download`      | `True`         | Automatically download missing data files on first use.                                                                                                                      |
 | `calc_area`         | `True`         | Calculate the watershed area in km² and add to output geodata attribute table |
| `clean`              | `False`        | Apply a small buffer/unbuffer to repair seam artifacts in the watershed polygon.                                                                                             |
| `data_dir`           | system default | Override the data cache location.                                                                                                                                            |
| `fill`               | `True`         | Fill small interior holes caused by topological gaps in MERIT-Hydro data.                                                                                                    |
| `fill_threshold`     | `100`          | Maximum hole size to fill, in pixels on the 3″ grid (~90 m/pixel near the equator). Set `0` to fill all holes.                                                               |
| `high_res`           | `True`         | Refine the watershed boundary at the outlet using raster methods. More accurate but slower. Set `False` to skip (watershed will include some area downstream of the outlet). |
| `low_res_threshold`  | `6e6`          | Area in km² above which the script automatically falls back to low-res mode. The Amazon is ~5.9×10⁶ km².                                                                     |
| `rivers`             | `True`         | Include the upstream river network in output.                                                                                                                                |
| `num_stream_orders`  | `4`            | The number of Strahler stream orders to include in river network output. Set ≥ 9 for all available reaches.                                                                  |
| `outlets`            | `True`         | Include requested and snapped outlet points in output.                                                                                                                       |
| `output_format`      | `gpkg`         | Output format: `gpkg`, `geojson`, `shp`, `kml`, `parquet`, or any GeoPandas-supported driver.                                                                                |
| `output_dir`         | `./output/`    | Directory for output files.                                                                                                                                                  |
| `search_dist`        | `0.1`          | Search radius in decimal degrees when the outlet falls outside all unit catchments (~10 km at the equator). Set `0` to require an exact hit.                                 |
| `simplify`           | `False`        | Simplify output geometry using Douglas-Peucker. Reduces file size and removes staircase artifacts from raster-origin boundaries.                                             |
| `threshold_single`   | `3000`         | Number of upstream pixels that defines a stream for snapping the outlet, when the outlet is in a unit catchment with no upstream contributing catchments.                    |
| `threshold_multiple` | `5000`         | Number of upstream pixels that defines a stream for snapping the outlet, when the outlet is in a unit catchment wih upstream contributing catchments.                        |
| `verbose`            | `False`        | Print messages to the console to monitor the script's progress. |


### Notes on select options


#### Filling holes

Setting `fill=True` removes small interior gaps or "donut holes" in the watershed polygon. These 
arise from slivers between unit catchments in the source data and are usually unwanted. 
The `fill_threshold` parameter (in pixels) controls which holes are filled — 
larger holes representing genuine endorheic (internally draining) basins can 
be preserved by setting a threshold.

For example, the Rio Grande watershed contains a large endorheic basin between 
the main stem and the Pecos River that should probably *not* be filled, at least
for studies of surface drainage:

![Rio Grande Watershed](docs/rio_grande.jpg)


#### Search distance

If the outlet point falls just offshore, in an estuary, or in a gap between unit
 catchments, `search_dist` controls how far (in decimal degrees) the script 
 searches for the nearest catchment. A value of at least `0.005` is recommended 
 for coastal outlets.


#### Simplify

The watershed boundary inherits the staircase pattern of the underlying raster 
grid (pixel edge length ≈ 0.000833°). Setting `simplify=True` with 
`simplify_tolerance ≈ 0.0004` or higher removes this artifact and reduces file size.
The `simplify_tolerance` parameter is equivalent to the threshold for 
[Douglas-Peucker simplification](https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm).


#### Thresholds for snapping

The process of "snapping" the outlet point to a river centerline is where 
watershed delineation becomes both an art and a science. The `threshold_single` 
and `threshold_multiple` parameters control how many upstream pixels are 
required to define a stream for snapping the outlet point. The values for these 
parameters define how many upstream pixels are required to define a stream. 

![Accumulation raster](docs/accum11_screenshot2.jpg)


<a id="data-files"></a>
## Data files

The `delineator` package comes bundled with data for Iceland. Beyond this, 
you will need data files for other regions. 
The globe is divided into 59 **megabasins** (integer IDs 11–86, data for 
Greenland, megabasin 91, has been omitted):

![Megabasins map](docs/megabasins.jpg)

Each megabasin requires four data files (vector catchments, vector rivers, 
flow-direction raster, accumulation raster). These download automatically 
on first use and are saved in your system's default data directory:

- **Windows:** `C:\Users\<username>\AppData\Local\delineator`
- **Linux:** `~/.local/share/delineator`
- **macOS:** `~/Library/Application Support/delineator`

To pre-download data for a region:
```bash
delineator_download --basin 62    # e.g. basin 62 = Amazon
delineator_dir                    # show the cache location
```

You can also download these datasets manually by visiting:
[https://mghydro.com/watersheds/delineator-data.html](https://mghydro.com/watersheds/delineator-data.html).

Some regional datasets are up to 3 GB, so pre-downloading is recommended for 
large basins.

Override the default data directory with an environment variable:
```bash
# macOS/Linux
export DELINEATOR_DATADIR=~/gis/delineator_data

# Windows
set DELINEATOR_DATADIR=D:\GIS\delineator_data
```


<a id="review-results"></a>
## ⚠️ Always review your results

**No automated watershed delineation software can replace human judgment. Always visually inspect every watershed you create with this package — there is no guarantee the output is correct.**

Errors are common and often easy to miss without inspection. The good news is 
that many mistakes can be fixed by slightly adjusting the outlet coordinates 
and re-running. An experienced analyst can usually identify and resolve problems 
quickly, especially with an interactive map display.

### Where delineation is most likely to fail

Certain landscapes are inherently challenging for any automated tool:

- **Flat terrain** — where flow direction is ambiguous. Examples: Florida, 
  the Netherlands, the Ganges-Brahmaputra Delta.
- **Arid and semi-arid areas** — where channels are sparse or ephemeral. 
  Examples: North Africa, Central China, the American Southwest.
- **Frozen environments** — glaciers, tundra, and permafrost. Examples: Iceland,
   Greenland, northern Canada, northern Russia.
- **Karst and highly permeable terrain** — where surface drainage boundaries are
   poorly defined because water moves through the subsurface. Examples: the 
   Yucatán Peninsula, parts of the Deschutes basin in Oregon, the Karst Plateau 
   along the Italy–Slovenia border.
- **Urban areas** — where impervious surfaces, curbs, storm sewers, and drains 
  alter or override natural flow paths.
- **Heavily engineered basins** — irrigation canals, inter-basin transfers, and 
pipelines can reroute water in ways that no terrain-based algorithm can detect.

### The most common error: incorrect pour point snapping

Even in well-behaved terrain, the most frequent source of error is pour point 
snapping — the outlet being snapped to the wrong river reach, often a nearby 
tributary. This produces a watershed on a completely different branch of the 
river network. Such errors are not correlated with watershed size or geography 
and can be subtle if you are not looking carefully.

If the result looks wrong, try nudging the outlet coordinates toward the river 
centerline and re-running. Overlaying the MERIT-Basins river network on your 
map makes this much easier. The [`examples/demo_webapp.py`](examples/demo_webapp.py) 
interactive map is useful for this kind of iterative review.

### Areas with no data

MERIT-Hydro does not cover Greenland, Antarctica, or some small islands 
(e.g., Hawaii, the Azores). Delineation will fail silently for outlet points in 
these areas.


<a id="usage-examples"></a>
## Usage examples

The `examples/` directory on the project's GitHub page contains ready-to-run scripts.
The example scripts show how to use `delineator`and even how to set up a local, 
web-based point-and-click 
watershed delineation service similar to [Global Watersheds](https://mghydro.com/watersheds).


<a id="algorithm"></a>
## Algorithm

The `delineator` combines three techniques to achieve speed and low memory use
compared to traditional raster watershed delineation methods:

1. **Hybrid raster/vector approach**: vector unit catchments handle the bulk of 
  the upstream area; raster methods refine only the home catchment 
  around the outlet.
2. **Hierarchical Spatial Aggregation**: pre-computed nested catchments at five 
  size levels (L0–L4) minimize the number of polygons that must be dissolved at 
  runtime.
3. **SQLite-backed geodata**: vector data is stored in relational SQLite 
  databases with spatial indexes, enabling fast SQL lookups rather than loading 
  entire datasets into memory.

![Method diagram](docs/method.jpg)

The nested catchments at the southern end of Madagascar illustrate the 
aggregation levels:

![Nested basins](docs/nested_basins.jpg)

For a more detailed description, see the manuscript: [Fast, accurate watershed delineation 
with a hybrid of raster and vector methods]
(https://mghydro.com/pages/Heberger_delineation_2025.pdf).


<a id="citation"></a>
## Citation

If you use `delineator` in your research, please cite the project 
homepage, this GitHub repository. Here's a BibTeX entry:

```bibtex
@software{delineator,
  author    = {Matthew Heberger},
  title     = {delineator: Global Watershed Delineation with Python},
  year      = {2026},
  publisher = {GitHub},
  version   = {2.0},
  url       = {https://github.com/mheberger/delineator}
}
```

<a id="contributing"></a>
## Contributing

This project is open source and welcomes contributions. If you have comments 
or suggestions, please open an issue or drop the author an email.
