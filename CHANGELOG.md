# Changelog

All notable changes to this project are documented here.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).


## [Unreleased]

### Changed

- The `delineate` command now writes rivers and outlet points by default, so a
  GeoPackage holds the watershed, rivers and outlets layers. Use
  `--rivers false` and `--outlets false` to leave them out. A bare `--rivers`
  or `--outlets` still works and means `true`.
- The `delineate` command now fills small holes in the watershed by default,
  matching `DelineatorConfig(fill=True)`, so the command line and Python use
  the same defaults. Use `--fill false` to keep the holes. A bare `--fill`
  still works and means `true`.
- `stream_threshold_km2` is now used exactly as given. Before, it was silently
  capped at 25% of the outlet's unit catchment area
  (`ADAPTIVE_THRESHOLD_FRACTION`), so in small catchments the outlet could
  snap to a channel much smaller than the threshold asked for. Removed the
  `ADAPTIVE_THRESHOLD_FRACTION` constant. If no cell in the unit catchment
  drains more than the threshold, snapping now fails with a warning that
  gives the catchment's largest drainage area. At the default of 25 km² (the
  MERIT-Basins reach threshold) this should not happen; above it, small
  headwater outlets may need a lower threshold.

### Fixed

- `write_outputs` now respects `round_coordinates=False`. Before, it always
  rounded the watershed and river coordinates to 6 decimal places, and it did
  so in place, changing the GeoDataFrames passed to it.
- `write_outputs` no longer raises when the watershed is `None`, as it is when
  `delineate` fails. It writes whichever layers are present, and if there are
  none, logs a warning and writes nothing.
- `delineator.__version__` reported 2.2.5 in the 2.2.6 release.
- Parquet output (`output_format="parquet"`) failed with an `ImportError` on a
  plain install, because nothing installed `pyarrow`. It is now the optional
  extra `pip install "delineator[parquet]"`. Without it, choosing parquet
  fails at once with that install hint, instead of after the delineation.
- `delineate` now reports a bad option value, such as an unknown
  `--output-format`, as a one-line error instead of a traceback.
- `delineate` no longer crashes with a traceback when it can't write an output
  file, for example because the file is open in QGIS or ArcGIS. It names the
  file, skips that outlet, and exits with status 1 if nothing was written.
- Delineating very large watersheds used several hundred megabytes more memory
  than it needed. The query that groups the upstream unit catchments wrote
  every upstream ID into the SQL text, and SQLite copied that list many times:
  about 700 MB extra for the Amazon at its mouth. The IDs now go through a
  temporary table. Peak memory for that watershed falls from about 820 MB to
  510 MB, the query runs in about half the time, and long-running processes
  such as the web app no longer grow by several megabytes after each large
  delineation.
- With `cache=True`, a cache full of continent-sized watersheds took about
  600 MB. It now takes about 240 MB, because the cached upstream IDs are kept
  as a compact array instead of a list of Python integers.


## [2.2.6] - 2026-11-22

### Fixed

- Watersheds could leave out an upstream branch that does drain to the outlet
  ([#13](https://github.com/mheberger/delineator/issues/13)), a regression in
  2.2.5. The new check for which upstream unit catchments reach the outlet
  missed a branch whose inlet pixel lay on the edge of the raster window,
  because pysheds clears the window's outer ring of cells before tracing. The
  check now also counts edge cells that flow into the outlet's catchment.
- `delineator_download --data-dir DIR` (and `downloader(basin, data_dir=DIR)`)
  now saves the files directly in `DIR`, which is where `delineate --data-dir DIR`
  looks for them. Before, the files went into a version subfolder (`DIR/v3/`),
  so `delineate` couldn't find them and downloaded them again, or failed with
  `--no-downloads`.

### Changed

- `delineator_dir` now prints two lines: the data directory, and the folder
  the data files are actually in. For the default location, the files are in
  a version-named subfolder (e.g. `v3/`), so the old one-line output led users
  to the wrong folder. Scripts that read its output need updating.
- The `--data-dir` help for `delineate` and `delineator_download` now says to
  run `delineator_dir` to see the default.
- Simpler test for whether upstream unit catchments drain to the outlet. A
  MERIT-Basins unit catchment begins at the confluence where its upstream
  catchments flow in, and every river reach drains more than 25 km². So the
  outlet receives the upstream flow exactly when its snapped cell drains at
  least 25 km²; otherwise the watershed is only the local catchment above it.
  This replaces a second pysheds trace over the neighbouring catchments
  (including the edge-cell workaround for #13) with one lookup in the
  accumulation raster. The partial case (an outlet above a confluence within
  its unit catchment) cannot occur in MERIT-Basins and is no longer handled.
  Results on a sample of Iceland outlets are unchanged. With snapping off,
  the accumulation raster is now also required.


## [2.2.5] - 2026-11-21

### Fixed

- Watersheds no longer include upstream unit catchments that do not drain to
  the outlet. This happened when the outlet was on a small tributary within
  its unit catchment (most often when snapping to a low `stream_threshold_km2`,
  or with snapping disabled), and the result included all the drainage area
  entering that unit catchment from upstream. The code now checks, using the
  flow-direction raster, which upstream unit catchments actually drain to the
  outlet, and builds the watershed, area, and rivers from those alone. It also
  handles an outlet above a confluence within its unit catchment, keeping only
  the branch it sits on. Previously this check ran only with snapping
  disabled, and it was all-or-nothing.
- Fixed a crash (`UnboundLocalError`) in high-res mode with snapping disabled,
  when the outlet was on the river network of a multi-catchment watershed.
- In low-res mode, `area_km2` now includes the outlet's own unit catchment.
  It was previously left out, although the polygon included it.

### Changed

- The web map no longer uses the Carto basemap, which now requires an API
  key. The default basemap is now OpenStreetMap, which does not require a key.
- Removed the `DISCONNECT_TOL_DEG` constant, which the old check used.


## [2.2.4] - 2026-07-11

### Added

- Added indexes to the SQLite databases for faster lookups. Vector data
  files (.db) have been bumped to version 3.

### Removed

- The runtime no longer loads the SpatiaLite extension for point-in-polygon
  lookups. Benchmarking showed no performance advantage over the GeoPandas
  path (which reads the same R*Tree spatial index directly).
  This removes the optional `mod_spatialite`
  system dependency, the `USE_SPATIALITE` constant, and the `use_spatialite`
  parameter of `spatial._point_in_polygon_analysis`. The data files are
  unchanged (still SpatiaLite-format), and the offline data-prep pipeline in
  `data_mgmt/` still uses SpatiaLite to build them.


## [2.2.3] - 2026-07-10

### Fixed

- Fixed a crash when rivers were requested but none exist for the watershed
  (for example, coastal catchments have no river geometries) and coordinate
  rounding was enabled — which is the default configuration.
- `delineate()` no longer raises an error when a rivers database is missing
  or contains no record for the outlet's catchment. In keeping with its
  contract (log a warning and return rather than raise), it now returns the
  watershed without the `area_km2` attribute and logs a warning.
- Fixed a latent `KeyError` in low-res mode that could occur when every
  upstream unit catchment was aggregated into higher-level groups.
- A whitespace-only `DELINEATOR_DATA_DIR` environment variable is now
  ignored (falling back to the platform default data directory) instead of
  creating a directory with a whitespace name.
- All SQLite database connections are now closed promptly after use.
  Previously, `delineate()` left connections open, which could keep the
  `.db` data files locked for the life of the process, especially on Windows.

### Changed

- Input coordinates are now validated before any database connection or
  other resource is opened.
- `rasterio` is now declared as a direct dependency. It was previously
  relied upon transitively through `pysheds`, even though the package
  imports it directly.
- Internal code cleanup: removed unused helper functions, consolidated
  duplicate definitions of `KM_PER_DEGREE` and of the data-directory
  resolution logic into single sources of truth, and updated stale
  docstrings, comments, and developer documentation.

### Added

- New configuration option `cache` (default `False`). When enabled, the most
  expensive steps of delineation, querying the upstream unit catchments and
  dissolving their polygons into a single boundary, are cached in memory,
  keyed on the outlet's unit catchment. Repeated delineations whose outlets
  fall in the same unit catchment reuse the cached result, which is much
  faster for large watersheds (for outlets near the mouth of the Rio Negro,
  about 8 seconds for the first call and about 2.5 seconds for subsequent
  ones). The cache holds up to 64 watersheds and can be emptied with the new
  `delineator.clear_cache()` function. The local web app now enables caching
  by default, since map users often click several points along the same river.
- New regression tests covering the fixes above, including a check that
  all database connections are closed after delineation.

## [2.2.1] - 2026-07-06

### Changed

- All the configuration options that were previously in units of decimal degrees
  are now in units of kilometers, and those with units of number of pixels
  are now in units of square kilometers.
- Several configuration options have been renamed or changed:

#### Removed / replaced fields

| Old field | Old default | Replaced by | New default | Notes |
|---|---|---|---|---|
| `fill_threshold: int` | `100` (pixels on 3-arcsec grid) | `fill_area_max: float` | `1.0` (km²) | Hole-fill cutoff is now a real area, not a pixel count |
| `threshold_multi: int` | `5000` (upstream pixels) | `stream_threshold_km2: float` | `25.0` (km²) | Stream-definition threshold for outlet snapping, now in km² of drainage area |
| `threshold_single: int` | `300` (upstream pixels) | *(folded into adaptive logic)* | — | No longer a config field; small-catchment threshold is derived adaptively (see `ADAPTIVE_THRESHOLD_FRACTION` in `constants.py`) |
| `low_res_threshold: float` | `6e6` (km²) | *(removed)* | — | Auto-switch-to-low-res threshold is gone; low-res is controlled by `high_res` |

#### Units / defaults changed (same field name)

| Field | Old | New              | Change                                                                |
|---|---|------------------|-----------------------------------------------------------------------|
| `search_dist: float` | `0.1` decimal degrees | `5.0` kilometers | Unit changed to km; now validated to the range 0–10 km                |
| `simplify_tolerance: float` | `0.0008` decimal degrees | `0.1` kilometers | Unit changed to km  |
| `clean: bool` | `False` | `True`           | Default flipped — seam-cleaning is now on by default                  |

  
### Fixed

- The watershed area is now correctly reported when the configuration option
  `fill = True`. The attribute `filled_area` is now calculated and added to 
  the geodata attributes. This reports the area in km² for any of the "donut
  holes" in the watershed that were filled in. 
- A variety of minor code optimizations and cleanups. 

## Added

- Added documentation for how to output to parquet and flatgeobuf formats.
- The watershed geodata now contains the two fields: 
  - `area_km2`: The actual area of the final watershed, as delivered
  - `fill_area`: The area of any internal holes that have been filled, in km²
    This will be `0` if you chose `fill=False`, or if the script did not find any 
    donut holes to fill in.
- Added smoothing option for watersheds and rivers, similar to the 
  feature on the Global Watersheds web app.
- Added the configuration options to the demo web app.  

## [2.2.0] - 2026-07-05

### Fixed

- Improved the handling of multipart polygon features throughout the code. 
  Previously, scripts were not correctly handling "bowtie" polygons that 
  are attached at a single point. Small gaps in the split unit catchment 
  were also being filled in, even when configuration option `fill=False`.
- Created new unit catchment spatialite databases, with correct handling 
  of multipart features. Dangling features attached at a single point are 
  no longer dropped, and holes are no longer filled in. Users of old 
  versions will be prompted to re-download them and delete old copies.

### Added

- Added a new configuration option `round_coordinates`. The new default is to  
  round the coordinates in geodata files to 6 decimal places, saving disk space
  with minimal loss of precision (around 10 cm near the equator).
- Added a new configuration option `snapping`, defaulting to `True`. In certain
  (rare) cases, the user may want to disable to pour point snapping algorithm,
  for example, if they have implemented their own snapping algorithm. 

### Changed

- Rounded the fields `uparea` and `lengthkm` in the rivers databases to 2 digits. 
  Former entries like 8.252159925700308 are now 8.25.
- Converted the docstrings for all functions to NumPy format, including 
  `Parameters` and `Returns` sections. It is my hope that consistent, thorough
  documentation will help those who wish to modify the code or contribute. 

## [2.1.2] - 2026-06-29

### Fixed

- Fixed a bug where the `delineate()` function would fail if the script 
  cannot find the default data directory, which is always the case after
  a new installation.

### Changed

- Minor edits to the README.

## [2.1.1] - 2026-06-26

### Fixed

- Fixed a mistake that caused the watershe area to be overestimated by double-
  counting the area of the home unit catchment.

### Added

- Made the demo web app a part of the package.

### Changed

- Removed some print statements in spatial.py that were purely for testing.


## [2.1.0] - 2026-06-25

### Fixed

- Updated the vector data files to use the MERIT-Basins 'bugfix' release 
- Fixed a bug where unit catchment lookups were slow because the field `comid` 
  was not indexed in the vector unit catchment sqlite databases.

### Changed

- Updated requirements for `geopandas` to version 1.0 and up
- Added requirement for `pyogrio` version 0.7 and up

### Added

- New versions of the river and basins vector data files are now available.
- Added a routine to bust the cache for old versions of vector data files,
  fires on initialization of the DelineatorConfig class. This will delete old
  versions of the vector data files and re-download them as needed.


## [2.0.6] - 2026-06-10

### Fixed

- Made logging behavior consistent across all modules, and 
  more consistent with Python standards. Now shows warnings
  and critical errors in the console.

### Added

- Added new options to the configuration class: `verbose` 
  will show more detailed information in the console as
  the script runs and `calc_area` gives the option to 
  turn off the area calculation if it is not needed
  (default is True, and area calculation is still done).

## [2.0.5] - 2026-06-08

### Fixed

- Fixed bug where `split_catchment()` fails when raster data not found. 
  Now fails more gracefully -- logs a warning and returns `None`.

## [2.0.4] - 2026-06-08

### Changed

- Made it the default to download necessary data files when using the command-line interface.
- Detailed raster delineation in `merit_detailed.py` now fails gracefully when the raster is not available
  and emits a warning instead of raising an error.

## Added

- Added tests for the command-line interface.

## [2.0.3] – First release on PyPI, 2026-06-06

### Changed

- Major updates to `spatial.py` to ensure the use of the R*tree index for
  spatial lookups.
- Added the parameters for pour point snapping to the configuration class.
- Updated the documentation and added a README to the `examples` folder.


## [2.0.2] – Released on Test PyPI, 2026-06-07

### Fixed

- Fixed a bug where spatial lookups failed when `spatialite` was not installed
  on the user's machine by updating the `_find_with_geopandas` function.

## [2.0.1] – Released on Test PyPI on 

### Changed

- Eased requirement to Python 3.10+. Python 3.11+ is still recommended for speed. 

### Fixed

- Made updates to the basins geodata; coastal catchment with `comid = 0`
  caused the script to hang.

## [2.0.0] — Unreleased

Version 2.0 is a major rewrite. The project has been restructured from a
collection of scripts you edit and run into a proper Python package you install,
import, and call, while preserving all the core algorithms and data from v1.

### Added

- **pip-installable package** — `pip install delineator`. No more manually
  downloading scripts or editing paths.
- **Python API** — `delineate(lat, lng, config)` returns three
  GeoDataFrames (watershed, rivers, outlets) that callers can inspect,
  transform, or persist however they like.
- **CLI entry point** — `delineate` command with `--point`, `--csv`,
  `--output-dir`, `--output-format`, `--rivers`, `--outlets`, `--fill`,
  `--simplify`, `--high-res`/`--low-res`, and `--verbose` flags.
- **`DelineatorConfig` dataclass** — all settings in one validated object,
  replacing the manually-edited `config.py` file. Invalid values raise a
  descriptive `TypeError` or `ValueError` at construction time.
- **Automatic data download** — data files are fetched on first use via
  `pooch` and cached in the platform's standard data directory
  (`platformdirs`). No more manual downloading and path configuration.
- **`delineator_download` CLI command** — pre-download data for a given
  megabasin before running delineation.
- **`delineator_dir` CLI command** — print the data cache location.
- **`DELINEATOR_DATADIR` environment variable** — override the cache
  location without touching code.
- **Hierarchical Spatial Aggregation** — pre-computed nested catchment
  polygons at five size levels (L0–L4) are stored in SQLite. Large
  watersheds are now assembled by dissolving a small number of large
  pre-aggregated polygons rather than thousands of unit catchments,
  dramatically reducing processing time for basins like the Amazon or
  Mississippi.
- **SQLite-backed geodata** — all vector data (catchments, rivers,
  megabasin boundaries) is stored in SQLite/SpatiaLite databases with
  spatial indexes, replacing the ESRI shapefiles used in v1. SQL queries
  replace in-memory GeoDataFrame operations for the most expensive lookups.
- **Bundled Iceland data** — megabasin 27 data ships with the package so
  first-time users can run immediately with no download step.
- **`write_outputs()` utility** — a single function to write watershed,
  rivers, and outlets to disk in any supported format.
- **GeoParquet output** — `output_format="parquet"` is now supported
  alongside GeoPackage, GeoJSON, Shapefile, and KML.
- **`clean` option** — small buffer/unbuffer pass to repair seam artifacts
  in the watershed polygon.
- **`outlets` output layer** — the snapped outlet point (and the original
  requested point) are returned and optionally written to disk.
- **Interactive web map example** (`examples/webapp.py`) — a single-file
  Flask + Leaflet application that serves a click-to-delineate map at
  `http://localhost:5000`.
- **`examples/` directory** — `demo_core.py` (API usage across multiple
  continents), `demo_cli.py` (CLI testing harness), and
  `sample_outlets.csv` (ten outlet points spanning six continents).
- Python 3.12+ is now required.

### Changed

- **Script → package structure** — source code moved to `src/delineator/`
  with a `pyproject.toml` build definition. The old `delineate.py`,
  `config.py`, and `py/` module directory are gone.
- **Configuration** — settings previously scattered across `config.py`
  uppercase constants (`HIGH_RES`, `FILL`, `OUTPUT_EXT`, etc.) are now
  fields on `DelineatorConfig` (`high_res`, `fill`, `output_format`, etc.).
- **Data format** — input data has moved from MERIT-Basins shapefiles to
  purpose-built SQLite databases. The new format is smaller, faster to
  query, and requires no separate spatial index creation step.
- **River network traversal** — the recursive `addnode()` function that
  walked the river network in v1 has been replaced with SQL-based upstream
  queries against the SQLite databases.
- **Pickle caching removed** — v1 cached GeoDataFrames as `.pkl` files to
  amortize the cost of reading shapefiles. The SQLite backend is fast
  enough that this workaround is no longer needed.
- **Map output removed** — the v1 `MAKE_MAP` option that produced a
  standalone HTML viewer (`.js` + `_viewer.html`) has been removed.
  The new `examples/webapp.py` provides a more capable replacement.
- **`OUTPUT_CSV` removed** — the per-run `OUTPUT.csv` summary file is no
  longer produced by default. Callers can inspect the returned GeoDataFrames
  directly.
- **Area matching removed** — the `MATCH_AREAS` / `AREA_MATCHING_THRESHOLD`
  pour-point relocation feature from v1 is not present in v2.
- **`low_res_threshold` default raised** — v1 defaulted to 50,000 km²;
  v2 defaults to 6,000,000 km² (effectively disabled for most use cases,
  since the SQLite/aggregation backend is fast enough for large basins).
- Minimum supported Python version raised from 3.9 to 3.12.

### Fixed

- Outlet snapping is now more robust: `search_dist` defaults to 0.1°
  (vs. 0 in v1), reducing silent failures for coastal or near-coastal
  outlets.
- Coordinate validation now raises clear errors for out-of-range or
  non-finite lat/lon values rather than failing silently downstream.



## [1.3] — 2023-11-16

### Changed

- Upgraded Python library dependencies.



## [1.2] — 2023-09-28

### Added

- Pickle file caching for GeoDataFrames — avoids re-reading shapefiles and
  rebuilding spatial indexes on repeated runs, improving performance.
- Automatic creation of output directories if they do not exist.

### Fixed

- Input shapefiles with missing `.prj` files (introduced in a MERIT-Basins
  upstream bugfix release) are now handled gracefully.



## [1.1] — 2023-09-21

### Fixed

- Rare bug where the dissolve step produced a `MultiPolygon` instead of
  the expected `Polygon`, causing downstream errors.



## [1.0] — 2022-11-11

Initial release.

- Watershed delineation from a CSV of outlet points using MERIT-Hydro and
  MERIT-Basins data.
- Hybrid raster/vector approach: vector unit catchments for the bulk of the
  upstream area; `pysheds` raster flow-direction for the terminal catchment.
- High-resolution and low-resolution modes.
- Output to GeoPackage, GeoJSON, or Shapefile via GeoPandas.
- Optional interactive HTML map viewer (`MAKE_MAP` in `config.py`).
- Sample data for Iceland bundled in the repository.
- Configuration via `config.py`.

[2.0.0]: https://github.com/mheberger/delineator/compare/v1.3...HEAD
[1.3]: https://github.com/mheberger/delineator/compare/v1.2...v1.3
[1.2]: https://github.com/mheberger/delineator/compare/v1.1...v1.2
[1.1]: https://github.com/mheberger/delineator/compare/v1.0...v1.1
[1.0]: https://github.com/mheberger/delineator/releases/tag/v1.0
