
### Interactive web map (`examples/demo_webapp.py`)

A self-contained Flask application that spins up a local delineation service 
and serves an interactive Leaflet map. Click anywhere on the map and the 
watershed, river network, and outlet points appear within seconds.

This is basically a simplified version of the (Global Watersheds web app)
[https://mghydro.com/watersheds]  that you can run locally.

![Webapp screenshot](docs/webapp_example.png)

**Install Flask and launch:**
```bash
pip install flask
python examples/webapp.py
```

Then open **http://localhost:5000** in your browser. The app is a single file 
— the Flask backend handles delineation via a `/delineate` POST endpoint, and 
the HTML/JS frontend is served inline.

The Iceland data is bundled, so clicking anywhere in Iceland works immediately. 
Clicking elsewhere in the world triggers an automatic data download for that 
megabasin on first use.

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

Shows how to call `delineate()` directly for a single point or a batch of named 
outlets across different continents:

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

A set of named test points is also defined in the script, useful for verifying 
that downloads and delineation work correctly across different megabasins.

### Sample outlet CSV (`examples/sample_outlets.csv`)

A ready-made CSV with ten outlet points spanning six continents, useful for 
testing batch delineation:

```bash
delineate --csv examples/sample_outlets.csv --rivers --output-dir output/
```
