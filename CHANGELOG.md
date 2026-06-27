# Changelog

All notable changes to this project are documented here.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

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
