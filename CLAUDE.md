# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`delineator` is a pip-installable Python package for fast, near-global watershed
delineation from any lat/lon point. It combines vector unit catchments (MERIT-Basins)
with raster refinement (MERIT-Hydro) so large watersheds can be delineated on a laptop.
The shippable package lives entirely under `src/delineator/`; everything else in the
repo is supporting material (data-prep scripts, tests, docs, benchmarks).

## Commands

The virtualenv is at `.venv/`. Use its tools directly (Windows layout):

```powershell
# Run the test suite
.venv\Scripts\pytest.exe

# Run a single test
.venv\Scripts\pytest.exe tests/test_core.py::test_delineate_iceland_returns_watershed

# Lint
.venv\Scripts\ruff.exe check src/

# Dependency hygiene (config in [tool.deptry] in pyproject.toml)
.venv\Scripts\python.exe -m deptry .

# Editable install after dependency/build changes
.venv\Scripts\pip.exe install -e .
```

`[tool.pytest.ini_options]` in `pyproject.toml` skips tests marked `network` by default.
Those download several GB into the default data directory; run them with
`.venv\Scripts\pytest.exe -m network`.
`[tool.ruff.lint]` pins ruff's classic default rules (`E4, E7, E9, F`), because
ruff 0.17 widened its defaults.

### CLI entry points (defined in `pyproject.toml` [project.scripts], all map to `cli.py`)

```bash
delineate --point 64.71072 -21.60337 --rivers --outlets   # single point
delineate --csv outlets.csv --output-format geojson        # batch (needs id,lat,lon cols)
delineator_download --basin 62                             # pre-download a megabasin
delineator_dir                                            # print the data directory
delineator_serve                                          # local Flask app at :5000
```

The command names don't all match their function targets: `delineate → cli.main` and
`delineator_serve → cli.serve_cli` (the other two match their command names).

## Data model — critical to understand before changing anything

The globe is split into **60 megabasins** (integer IDs 11–86; Greenland/91 excluded),
listed in `constants.VALID_MEGABASINS`. Each megabasin needs four data files:
`basins{NN}.db`, `rivers{NN}.db` (SQLite/SpatiaLite vector), and `flowdir{NN}.tif`,
`accum{NN}.tif` (rasters).

- **Only megabasin 27 (Iceland) is bundled** in `src/delineator/data/`. All other
  regions download on first use via `pooch`, verified against SHA-256 hashes in
  `constants.HASHES`. Tests rely on the bundled Iceland data (coordinates ~64.7, -21.6).
  The `.db`/`.tif` data files are gitignored, and hatchling drops VCS-ignored files from
  the build; `[tool.hatch.build] artifacts = ["src/delineator/data/**"]` in `pyproject.toml`
  is what forces the bundled Iceland data back into the wheel. Don't remove it.
- Inside each `basins{NN}.db`: `l0_basins` is the raw MERIT unit catchments (a `nextdown`
  column encodes the drainage graph); `l1`–`l4` are pre-aggregated nested catchments at
  coarser sizes ("Hierarchical Spatial Aggregation"). Upstream traversal is a recursive
  SQL CTE over `l0_basins.nextdown` (see `queries.get_upstream_comids`). The aggregation
  levels let the code dissolve a handful of large polygons instead of thousands of small
  ones at runtime.
- When `constants.HASHES` changes because data files were regenerated, bump
  `constants.DATA_DIR_NAME` (e.g. `"v2.2"`). Default-managed data lives in a
  version-named subfolder so corrected files land in a fresh directory; user-supplied
  data dirs (CLI `--data-dir`, `DELINEATOR_DATA_DIR`, or `DelineatorConfig(data_dir=...)`)
  use a flat layout. See `settings.DelineatorConfig._resolve_data_dir`.
- The runtime never loads the SpatiaLite extension (removed July 2026 after benchmarks
  showed no advantage once the DBs had attribute indexes). The `.db` files are still
  SpatiaLite-format — `data_mgmt/` uses SpatiaLite to build them, and `spatial.py`
  reads their R*Tree index tables (`idx_{table}_{geom}`) as plain SQLite virtual tables.

## Runtime architecture (`src/delineator/`)

`delineate(lat, lng, config)` in `core.py` is the pipeline entry point. It returns a
3-tuple of GeoDataFrames `(watershed, rivers, outlets)`, any of which may be `None`.
It **logs warnings and returns `None` on failure rather than raising** — preserve this.

Pipeline steps in `core.delineate`:
1. Validate coords (`util._validate_outlet_coordinates`).
2. Point-in-polygon to find the megabasin, then the home unit catchment
   (`spatial._point_in_polygon_analysis`, with `search_dist` fallback for near-misses).
3. Recursive-CTE query for all upstream comids (`queries.get_upstream_comids`), and —
   for multi-catchment watersheds — the hierarchical grouping query
   (`queries.get_hiearchical_basins`, note the typo), both while the scoped basins
   DB connection is open.
4. **High-res only** (`config.high_res`, on by default): split the home catchment at
   the outlet using raster methods (`merit_detailed._split_catchment`, backed by
   `pysheds` + the flowdir/accum rasters). This is also where pour-point **snapping**
   happens, gated by `config.stream_threshold_km2` (used as given; if no cell in the
   home catchment exceeds it, snapping fails with a warning). When the snap threshold is
   below `constants.MERIT_REACH_THRESHOLD_KM2` (25 km², the dataset's reach threshold),
   or snapping is off, `merit_detailed._is_on_mainstem` checks the outlet cell's
   accumulation against that reach threshold; a point snapped at ≥ 25 km² is on the
   reach by construction, so the check is skipped. Below the reach threshold,
   the outlet is off the reach, so no upstream catchment drains to it and the watershed
   is just the split.
5. Fetch upstream geometries by aggregation level and `unary_union` them
   (`queries.get_upstream_geometries`), splicing in the split home catchment.
6. Optional hole-filling (`util._close_holes`), simplify, `clean` (buffer/unbuffer),
   and smooth.
7. Optional rivers (`queries.get_rivers`) and outlet points.

Module map:
- `core.py` — orchestration (`delineate`, `downloader`).
- `settings.py` — `DelineatorConfig` dataclass: the single source of truth for all
  options + validation + env-var resolution (`DELINEATOR_DATA_DIR`,
  `DELINEATOR_OUTPUT_DIR`, `DELINEATOR_AUTO_DOWNLOAD`). Config objects are mutable and
  reused across calls; `core.delineate` does not mutate them (the web app in `serve.py`
  still copies per request before applying client option overrides).
- `queries.py` — all SQL against the basin/river DBs.
- `spatial.py` — point-in-polygon lookups.
- `merit_detailed.py` — the only raster/pysheds code; runs on a single catchment for speed.
- `smoothing.py` — optional output beautification (`config.smooth`): Catmull-Rom splines for
  rivers, Chaikin corner cutting for the watershed boundary. Runs last, after simplify/clean.
- `data.py` / `download.py` — locating and fetching data files (pooch, hashes, stale-data warnings).
- `constants.py` — megabasin IDs, pixel geometry, `HASHES`, `DATA_DIR_NAME`, `USE_SPATIALITE`.
- `util.py` — `write_outputs` and geometry helpers. `validation.py` — CSV normalization.
- `cli.py` — Click commands. `serve.py` — Flask web app.
- `__init__.py`'s `__version__` must be kept in sync with the version in `pyproject.toml`.

## `data_mgmt/` — offline data-prep pipeline (NOT shipped)

Numbered folders (`0_raster` → `6_illustrate`) reproduce the source `.db`/`.tif` files
from raw MERIT-Basins/MERIT-Hydro downloads: fix shapefile geometries, build the
SpatiaLite basin DBs with nested aggregation (`aggregate_basins.py` + `consolidate.py`),
build river DBs, compute SHA-256 hashes for `constants.HASHES`. Several scripts run
inside QGIS's Python (filenames prefixed `qgis_`). This code is excluded from deptry
(`extend_exclude`) and never imports as part of the package. Regenerating data here is
what forces a `DATA_DIR_NAME` bump and a new `HASHES` block.

## Conventions

- Coordinate order is **(lat, lng)** in the public `delineate()` API, but shapely
  `Point`s are constructed `(lng, lat)`. Watch this at every boundary.
- All geodata is EPSG:4326; areas are computed on the ellipsoid via `pyproj.Geod`.
- Private helpers are underscore-prefixed and imported across modules by name.
