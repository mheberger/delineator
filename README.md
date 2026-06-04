# delineator: Global Watershed Delineation with Python

Fast, accurate watershed delineation for any point on Earth's land surface, using a hybrid of vector- and raster-based methods with data from [MERIT-Hydro](https://global-hydrodynamics.github.io/MERIT_Hydro/) and [MERIT-Basins](https://www.reachhydro.org/home/params/merit-basins).

- Near-global coverage (excludes Greenland, Antarctica, and some small islands)
- Bundled sample data for Iceland; other regions download automatically on first use
- Returns watershed polygon, river network, and outlet points as GeoPandas GeoDataFrames


## Installation

Requires Python ≥ 3.12. A virtual environment is strongly recommended.

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

## Quick start

The bundled Iceland data lets you run immediately after install; no download required.

**Command line:**
```bash
delineate --point 63.938 -21.004          # Ölfusá River at Route 1, Iceland
delineate --point 63.938 -21.004 --rivers --outlets
```

Output is written to `./output/watershed.gpkg` in your current directory.

**Python API:**

```python
from delineator import delineate, write_outputs

# Returns up to three GeoDataFrames: watershed, rivers, outlets
# Note: latitude comes before longitude
watershed_gdf, rivers_gdf, outlets_gdf = delineate(63.938, -21.004)

# Write all outputs to disk in one call
write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, id="olfusa")
```

Here is an example of the output displayed in QGIS:

![Example output](docs/example_output.png)


## Command line reference

```bash
# Single point
delineate --point 63.938 -21.004

# Include rivers and outlet points
delineate --point 63.938 -21.004 --rivers --outlets

# Batch delineation of multiple outlet points in a CSV file
delineate --csv outlets.csv

# Output format
delineate --point 63.938 -21.004 --output-format geojson
delineate --point 63.938 -21.004 --output-format shp
delineate --point 63.938 -21.004 --output-format kml
delineate --point 63.938 -21.004 --output-format parquet

# Custom output directory
delineate --csv outlets.csv --output-dir /path/to/output/

# List all the command line options
delineate --help
```

The CSV file must contain at minimum `id`, `lat`, and `lon` columns:

```
id,lat,lon,name
6401070,64.71072,-21.60337,Nordhura River at Stekkur
6401080,64.69229,-21.41046,Hvita River at Kljafoss
6401090,63.93796,-21.00666,Olfusa River at Selfoss
```

## Environment variables

Instead of passing options to the command line, you can set environment variables
for the default data directory and the output director. There are three environment variables:

- `DELINEATOR_DATA_DIR`: directory where input data files are saved
- `DELINEATOR_OUTPUT_DIR`: directory where output files will be saved
- `DELINEATOR_AUTO_DOWNLOAD`: whether to automatically download data files as they are needed

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

## Data

The globe is divided into 59 **megabasins** (integer IDs 11–86, data for Greenland,
#91, has been omitted):

![Megabasins map](docs/megabasins.jpg)

Each megabasin requires four data files (vector catchments, vector rivers, flow-direction raster, accumulation raster). These download automatically on first use and are cached in your system's default data directory:

- **Windows:** `C:\Users\<username>\AppData\Local\delineator`
- **Linux:** `~/.local/share/delineator`
- **macOS:** `~/Library/Application Support/delineator`

To pre-download data for a region:
```bash
delineator_download --basin 62    # e.g. basin 62 = Amazon
delineator_dir                    # show the cache location
```

Some regional datasets are up to 3 GB, so pre-downloading is recommended for large basins.

Override the cache location with an environment variable:
```bash
# macOS/Linux
export DELINEATOR_DATADIR=~/gis/delineator_data

# Windows
set DELINEATOR_DATADIR=D:\GIS\delineator_data
```

## Configuration reference

When using the Python API, options are passed via a `DelineatorConfig` object:

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
watershed_gdf, _, _ = delineate(63.938, -21.59, config)
```

All options with their defaults:

| Option | Default | Description |
|---|---|---|
| `high_res` | `True` | Refine the watershed boundary at the outlet using raster methods. More accurate but slower. Set `False` to skip (watershed will include some area downstream of the outlet). |
| `low_res_threshold` | `6e6` | Area in km² above which the script automatically falls back to low-res mode. The Amazon is ~5.9×10⁶ km². |
| `fill` | `True` | Fill small interior holes caused by topological gaps in MERIT-Hydro data. |
| `fill_threshold` | `100` | Maximum hole size to fill, in pixels on the 3″ grid (~90 m/pixel near the equator). Set `0` to fill all holes. |
| `rivers` | `True` | Include the upstream river network in output. |
| `num_stream_orders` | `4` | Strahler stream orders to include in river output. Set ≥ 9 for all available reaches. |
| `outlets` | `True` | Include requested and snapped outlet points in output. |
| `output_format` | `"gpkg"` | Output format: `gpkg`, `geojson`, `shp`, `kml`, `parquet`, or any GeoPandas-supported driver. |
| `output_dir` | `./output/` | Directory for output files. |
| `data_dir` | system default | Override the data cache location. |
| `search_dist` | `0.1` | Search radius in decimal degrees when the outlet falls outside all unit catchments (~10 km at the equator). Set `0` to require an exact hit. |
| `simplify` | `False` | Simplify output geometry using Douglas-Peucker. Reduces file size and removes staircase artifacts from raster-origin boundaries. |
| `simplify_tolerance` | `0.0008` | Tolerance in decimal degrees for simplification (~half a pixel edge length). |
| `clean` | `False` | Apply a small buffer/unbuffer to repair seam artifacts in the watershed polygon. |
| `auto_download` | `True` | Automatically download missing data files on first use. |


## Notes on select options

### Filling holes

Setting `fill=True` removes small interior gaps in the watershed polygon. These arise from slivers between unit catchments in the source data and are almost always unwanted. The `fill_threshold` parameter (in pixels) controls which holes are filled — larger holes representing genuine endorheic (internally draining) basins can be preserved by setting a threshold.

For example, the Rio Grande watershed contains a large endorheic basin between the main stem and the Pecos River that should *not* be filled:

![Rio Grande Watershed](docs/rio_grande.jpg)

### Search distance

If the outlet point falls just offshore, in an estuary, or in a gap between unit catchments, `search_dist` controls how far (in decimal degrees) the script searches for the nearest catchment. A value of at least `0.005` is recommended for coastal outlets.

### Simplify

The watershed boundary inherits the staircase pattern of the underlying raster grid (pixel edge length ≈ 0.000833°). Setting `simplify=True` with `simplify_tolerance ≈ 0.0004` or higher removes this artifact and reduces file size.


## Examples

The `examples/` directory contains ready-to-run scripts.

### Interactive web map (`examples/webapp.py`)

A self-contained Flask application that spins up a local delineation service and serves an interactive Leaflet map. Click anywhere on the map and the watershed, river network, and outlet points appear within seconds.

![Webapp screenshot](docs/webapp_screenshot.png)

**Install Flask and launch:**
```bash
pip install flask
python examples/webapp.py
```

Then open **http://localhost:5000** in your browser. The app is a single file — the Flask backend handles delineation via a `/delineate` POST endpoint, and the HTML/JS frontend is served inline.

The Iceland data is bundled, so clicking anywhere in Iceland works immediately. Clicking elsewhere in the world triggers an automatic data download for that megabasin on first use.

### Batch delineation from a CSV file with live results table and a Leaflet map viewer 
 (`examples/batch_app.py`)

Usage:

```python
pip install flask
python batch_app.py iceland_outlets.csv

# Or to use a different CSV:
python batch_app.py path/to/my_outlets.csv
```

Then open http://localhost:5001 in your browser.
 
Expected CSV columns: id, lat, lng, name, area

  - id   : unique identifier (will become a clickable link)
  - lat  : outlet latitude  (decimal degrees)
  - lng  : outlet longitude (decimal degrees)
  - name : human-readable label
  - area : a priori watershed area (km²), used to compute % difference
 
Results are saved to disk as each watershed completes:

  - `output/<id>/watershed.geojson`
  - `output/<id>/rivers.geojson`      (if rivers were returned)
  - `output/<id>/outlets.geojson`     (if outlets were returned)
  - `output/results.csv`              (summary table, appended row by row)
 
On restart, any outlet whose output directory already exists is skipped and its
previously computed results are loaded from results.csv instead of re-running
delineation. Delete the output directory (or a specific subdirectory) to force
a re-run.


### Python API demos (`examples/demo_core.py`)

Shows how to call `delineate()` directly for a single point or a batch of named outlets across different continents:

```python
from delineator.core import delineate
from delineator.settings import DelineatorConfig
from delineator.util import write_outputs

# Single point — Seine River at Paris
lat, lng = 48.834, 2.263
config = DelineatorConfig(high_res=True, output_format="geojson")
watershed_gdf, rivers_gdf, outlets_gdf = delineate(lat, lng, config)
write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, config)
```

A set of named test points (Amazon, Chattahoochee, Madagascar, Iceland, and others) is also defined in the script, useful for verifying that downloads and delineation work correctly across different megabasins.

### Sample outlet CSV (`examples/sample_outlets.csv`)

A ready-made CSV with ten outlet points spanning six continents, useful for testing batch delineation:

```bash
delineate --csv examples/sample_outlets.csv --rivers --output-dir output/
```

## Output files

In GeoPackage mode (default), all layers are written to a single file (`watershed_<id>.gpkg`) with three layers: `watershed`, `rivers`, and `outlets`. For other formats, each layer is written to a separate file.


## ⚠️ Always review your results

**No automated watershed delineation software can replace human judgment. Always visually inspect every watershed you create with this package — there is no guarantee the output is correct.**

Errors are common and often easy to miss without inspection. The good news is that many mistakes can be fixed by slightly adjusting the outlet coordinates and re-running. An experienced analyst can usually identify and resolve problems quickly, especially with an interactive map display.

### Where delineation is most likely to fail

Certain landscapes are inherently difficult for any automated tool:

- **Flat terrain** — where flow direction is ambiguous. Examples: Florida, the Netherlands, the Ganges-Brahmaputra Delta.
- **Arid and semi-arid areas** — where channels are sparse or ephemeral. Examples: North Africa, Central China, the American Southwest.
- **Frozen environments** — glaciers, tundra, and permafrost. Examples: Iceland, Greenland, northern Canada, northern Russia.
- **Karst and highly permeable terrain** — where surface drainage boundaries are poorly defined because water moves through the subsurface. Examples: the Yucatán Peninsula, parts of the Deschutes basin in Oregon, the Karst Plateau along the Italy–Slovenia border.
- **Urban areas** — where impervious surfaces, curbs, storm sewers, and drains alter or override natural flow paths.
- **Heavily engineered basins** — irrigation canals, inter-basin transfers, and pipelines can reroute water in ways that no terrain-based algorithm can detect.

### The most common error: incorrect pour point snapping

Even in well-behaved terrain, the most frequent source of error is pour point snapping — the outlet being snapped to the wrong river reach, often a nearby tributary. This produces a watershed on a completely different branch of the river network. Such errors are not correlated with watershed size or geography and can be subtle if you are not looking carefully.

If a result looks wrong, try nudging the outlet coordinates toward the river centerline and re-running. Overlaying the MERIT-Basins river network on your map makes this much easier. The [`examples/webapp.py`](examples/demo_webapp.py) interactive map is useful for this kind of iterative review.

### Areas with no data

MERIT-Hydro does not cover Greenland, Antarctica, or some small islands (e.g., Hawaii, the Azores). Delineation will fail silently for outlet points in these areas.

`delineator` combines three techniques to achieve speed and low memory use:

1. **Hybrid raster/vector approach**: vector unit catchments handle the bulk of the upstream area; raster flow-direction grids refine only the home catchment around the outlet.
2. **Hierarchical Spatial Aggregation**: pre-computed nested catchments at five size levels (L0–L4) minimize the number of polygons that must be dissolved at runtime.
3. **SQLite-backed geodata**: vector data is stored in relational SQLite databases with spatial indexes, enabling fast SQL lookups rather than loading entire datasets into memory.

![Method diagram](docs/method.jpg)

The nested catchments at the southern end of Madagascar illustrate the aggregation levels:

![Nested basins](docs/nested_basins.jpg)

For a full description, see the manuscript: [Fast, accurate watershed delineation with a hybrid of raster and vector methods](https://mghydro.com/pages/Heberger_delineation_2025.pdf).


## Citation

If you use `delineator` in your research, please cite:

```bibtex
@software{delineator,
  author    = {Matthew Heberger},
  title     = {delineator: Global Watershed Delineation with Python},
  year      = {2026},
  publisher = {GitHub},
  version   = {2.0.0},
  url       = {https://github.com/mheberger/delineator}
}
```


## Contributing

This project is open source and welcomes contributions. If you have comments or suggestions, please open an issue or pull request, or drop the author an email.


## License

[MIT LICENSE](LICENSE.txt).
