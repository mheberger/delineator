"""
Created with Claude.ai on 2026-06-19

Thoroughly inspect an ESRI Shapefile for common problems *before* you rely on it
in an analysis. The checker reads the file (and its sidecar files) and reports
issues grouped by severity, without ever modifying the data.

Severities
----------
ERROR    Something that will break or silently corrupt downstream work
         (unreadable file, missing required sidecar, non-finite coordinates,
         field-name collisions, format-limit violations).
WARNING  Likely to cause trouble or wrong results (missing .prj, invalid
         geometries, mixed types, suspicious coordinate ranges, duplicates).
INFO     Worth knowing but not a problem on its own (geographic vs projected
         CRS, 3D coordinates, per-column null counts).

Usage
-----
    from shapefile_check import check_shapefile
    report = check_shapefile("watershed.shp")
    print(report)                 # human-readable summary
    if not report.ok:             # True only when there are no ERRORs
        ...
    report.errors, report.warnings, report.infos   # structured access

Or from the command line:
    python shapefile_check.py watershed.shp [--expected-crs EPSG:4326]

Dependencies: geopandas (>=0.12), shapely (>=2.0), pyproj.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

import geopandas as gpd
import numpy as np
import pandas as pd
from pyproj import CRS
from shapely.validation import explain_validity

try:
    # PyShp reads polygon rings exactly as stored on disk, without the
    # ring-reorganization that GDAL applies, so it is the only reliable way to
    # audit on-disk ring orientation. Optional: the checker degrades gracefully.
    import shapefile as _pyshp
except ImportError:
    _pyshp = None

ERROR = "ERROR"
WARNING = "WARNING"
INFO = "INFO"

# Shapefile format hard limits
MAX_FILE_BYTES = 2 * 1024 ** 3   # .shp and .dbf are each capped at 2 GB
MAX_FIELD_NAME_LEN = 10          # DBF truncates field names to 10 characters
MAX_FIELDS = 255                 # DBF maximum number of fields


@dataclass
class ShapefileReport:
    """Container for the results of checking one shapefile."""

    path: str
    issues: list[tuple[str, str, str]] = field(default_factory=list)
    stats: dict = field(default_factory=dict)

    def add(self, severity: str, category: str, message: str) -> None:
        self.issues.append((severity, category, message))

    def _by(self, severity: str) -> list[tuple[str, str, str]]:
        return [i for i in self.issues if i[0] == severity]

    @property
    def errors(self) -> list[tuple[str, str, str]]:
        return self._by(ERROR)

    @property
    def warnings(self) -> list[tuple[str, str, str]]:
        return self._by(WARNING)

    @property
    def infos(self) -> list[tuple[str, str, str]]:
        return self._by(INFO)

    @property
    def ok(self) -> bool:
        """True when no ERROR-level issues were found."""
        return len(self.errors) == 0

    def __str__(self) -> str:
        lines = [f"Shapefile check: {self.path}"]
        if self.stats:
            for key, val in self.stats.items():
                lines.append(f"    {key}: {val}")
        lines.append(
            f"    -> {len(self.errors)} error(s), "
            f"{len(self.warnings)} warning(s), {len(self.infos)} info"
        )
        if not self.issues:
            lines.append("    No issues found.")
        for severity in (ERROR, WARNING, INFO):
            for _, category, message in self._by(severity):
                lines.append(f"  [{severity:<7}] ({category}) {message}")
        return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Individual checks. Each is wrapped by check_shapefile so a single failing
# check never aborts the whole report.
# --------------------------------------------------------------------------- #

def _sidecar_present(shp_path: Path, ext: str) -> bool:
    """Case-insensitive test for a sidecar file next to the .shp."""
    target = (shp_path.stem + ext).lower()
    try:
        return any(f.name.lower() == target for f in shp_path.parent.iterdir())
    except OSError:
        return False


def _check_sidecar_files(shp_path: Path, report: ShapefileReport) -> None:
    required = {".shp": "geometry", ".shx": "shape index", ".dbf": "attributes"}
    for ext, desc in required.items():
        if not _sidecar_present(shp_path, ext):
            report.add(ERROR, "files", f"Missing required {ext} file ({desc})")
    if not _sidecar_present(shp_path, ".prj"):
        report.add(WARNING, "files",
                   "Missing .prj file: the CRS is undefined and coordinates are ambiguous")
    if not _sidecar_present(shp_path, ".cpg"):
        report.add(INFO, "files",
                   "Missing .cpg file: text-field encoding is not declared (often fine, "
                   "but a source of mojibake for non-ASCII data)")


def _check_crs(gdf: gpd.GeoDataFrame, report: ShapefileReport,
               expected_crs: Optional[Union[str, int, CRS]]) -> None:
    crs = gdf.crs
    if crs is None:
        report.add(WARNING, "crs", "No CRS is defined on the layer")
        return
    report.stats["crs"] = crs.to_string()
    if crs.is_geographic:
        report.add(INFO, "crs",
                   f"Geographic CRS '{crs.name}' (units = degrees; area/length "
                   "calculations require a projected CRS)")
    else:
        report.add(INFO, "crs", f"Projected CRS '{crs.name}'")
    if crs.to_epsg() is None:
        report.add(WARNING, "crs",
                   "CRS has no EPSG code (custom or non-standard); it may not "
                   "round-trip cleanly through other tools")
    if expected_crs is not None:
        try:
            if not crs.equals(CRS.from_user_input(expected_crs)):
                report.add(WARNING, "crs",
                           f"CRS '{crs.to_string()}' does not match expected "
                           f"'{CRS.from_user_input(expected_crs).to_string()}'")
        except Exception as exc:  # pragma: no cover - bad user input
            report.add(WARNING, "crs", f"Could not interpret expected_crs: {exc}")


def _check_geometry(gdf: gpd.GeoDataFrame, report: ShapefileReport,
                    max_examples: int) -> None:
    geom = gdf.geometry
    null_mask = geom.isna()

    if null_mask.any():
        rows = list(geom.index[null_mask][:max_examples])
        report.add(ERROR, "geometry",
                   f"{int(null_mask.sum())} null/missing geometries (e.g. rows {rows})")

    present = geom[~null_mask]
    if present.empty:
        return

    empty_mask = present.is_empty
    if empty_mask.any():
        report.add(WARNING, "geometry", f"{int(empty_mask.sum())} empty geometries")

    # Geometry-type consistency
    types = sorted(present.geom_type.dropna().unique())
    report.stats["geom_types"] = types
    base_types = {t.replace("Multi", "") for t in types}
    if len(base_types) > 1:
        report.add(WARNING, "geometry", f"Mixed geometry types present: {types}")
    elif any(t.startswith("Multi") for t in types) and \
            any(not t.startswith("Multi") for t in types):
        report.add(INFO, "geometry",
                   "Mix of single-part and multi-part geometries of the same base type")

    # Validity (ignore empties, which report as valid anyway)
    solid = present[~empty_mask]
    invalid_mask = ~solid.is_valid
    if invalid_mask.any():
        examples = []
        for idx in solid.index[invalid_mask][:max_examples]:
            examples.append(f"row {idx}: {explain_validity(solid.loc[idx])}")
        report.add(WARNING, "geometry",
                   f"{int(invalid_mask.sum())} invalid geometries "
                   f"(fixable with make_valid / buffer(0)). Examples: "
                   + "; ".join(examples))

    # 3D coordinates
    has_z = solid.has_z
    if has_z.any():
        report.add(INFO, "geometry",
                   f"{int(has_z.sum())} geometries carry Z (3D) coordinates")

    # Zero / negative area: a classic symptom of reversed ring orientation or
    # mis-nested holes (works even without pyshp installed).
    poly_mask = solid.geom_type.isin(["Polygon", "MultiPolygon"])
    if poly_mask.any():
        bad_area = poly_mask & (solid.area <= 0)
        if bad_area.any():
            report.add(WARNING, "geometry",
                       f"{int(bad_area.sum())} polygon(s) with zero or negative area "
                       "(often reversed ring orientation or mis-nested holes)")


def _check_coordinates(gdf: gpd.GeoDataFrame, report: ShapefileReport) -> None:
    geom = gdf.geometry
    usable = ~(geom.isna() | geom.is_empty)
    bounds = geom.bounds[["minx", "miny", "maxx", "maxy"]].to_numpy()

    finite = np.isfinite(bounds).all(axis=1)
    bad = (~finite) & usable.to_numpy()
    if bad.any():
        report.add(ERROR, "coordinates",
                   f"{int(bad.sum())} geometries contain non-finite (NaN/Inf) coordinates")

    total = gdf.total_bounds
    if np.isfinite(total).all():
        report.stats["bounds"] = [round(float(v), 6) for v in total]

    crs = gdf.crs
    if crs is not None and crs.is_geographic and np.isfinite(total).all():
        minx, miny, maxx, maxy = total
        if minx < -180 or maxx > 180 or miny < -90 or maxy > 90:
            report.add(WARNING, "coordinates",
                       f"CRS is geographic but bounds {report.stats.get('bounds')} fall "
                       "outside lon[-180,180]/lat[-90,90] - the data may actually be "
                       "projected or the .prj may be wrong")


def _check_attributes(gdf: gpd.GeoDataFrame, report: ShapefileReport) -> None:
    cols = [c for c in gdf.columns if c != gdf.geometry.name]
    report.stats["n_fields"] = len(cols)

    long_names = [c for c in cols if len(c) > MAX_FIELD_NAME_LEN]
    if long_names:
        report.add(WARNING, "attributes",
                   f"Field name(s) over {MAX_FIELD_NAME_LEN} chars will be truncated "
                   f"by the DBF format: {long_names}")

    truncated: dict[str, list[str]] = {}
    for c in cols:
        truncated.setdefault(c[:MAX_FIELD_NAME_LEN], []).append(c)
    collisions = {k: v for k, v in truncated.items() if len(v) > 1}
    if collisions:
        report.add(ERROR, "attributes",
                   f"Field names collide after {MAX_FIELD_NAME_LEN}-char truncation: "
                   f"{collisions}")

    empty_cols, partial_nulls = [], {}
    for c in cols:
        n_null = int(gdf[c].isna().sum())
        if n_null == len(gdf) and len(gdf) > 0:
            empty_cols.append(c)
        elif n_null > 0:
            partial_nulls[c] = n_null
    if empty_cols:
        report.add(WARNING, "attributes", f"Entirely null/empty column(s): {empty_cols}")
    if partial_nulls:
        report.add(INFO, "attributes", f"Columns with some null values: {partial_nulls}")

    # Replacement character usually signals a text-encoding problem
    for c in cols:
        col = gdf[c]
        if col.dtype == object or pd.api.types.is_string_dtype(col):
            try:
                if col.astype("string").str.contains("\ufffd", na=False).any():
                    report.add(WARNING, "attributes",
                               f"Column '{c}' contains the Unicode replacement character "
                               "(likely a text-encoding problem)")
            except Exception:
                pass


def _check_duplicates(gdf: gpd.GeoDataFrame, report: ShapefileReport) -> None:
    geom = gdf.geometry
    present = geom[~geom.isna()]
    if present.empty:
        return
    wkb = present.to_wkb()
    dup_geom = wkb.duplicated(keep=False)
    if dup_geom.any():
        report.add(WARNING, "duplicates",
                   f"{int(dup_geom.sum())} features share an identical geometry with "
                   "another feature")

    # Fully duplicate rows: compare attributes alongside the geometry's WKB.
    # WKB bytes (and None) are hashable, so duplicated() handles them directly.
    attr_cols = [c for c in gdf.columns if c != gdf.geometry.name]
    tmp = gdf[attr_cols].copy()
    tmp["_wkb"] = [g.wkb if g is not None else None for g in geom]
    dup_full = tmp.duplicated(keep=False)
    if dup_full.any():
        report.add(WARNING, "duplicates",
                   f"{int(dup_full.sum())} fully duplicate rows (identical geometry "
                   "and attributes)")


def _check_format_limits(shp_path: Path, gdf: gpd.GeoDataFrame,
                         report: ShapefileReport) -> None:
    for ext in (".shp", ".dbf"):
        f = shp_path.with_suffix(ext)
        if f.exists():
            size = f.stat().st_size
            if size > MAX_FILE_BYTES:
                report.add(ERROR, "format",
                           f"{ext} is {size / 1024**3:.2f} GB, over the 2 GB shapefile limit")
    n_fields = len([c for c in gdf.columns if c != gdf.geometry.name])
    if n_fields > MAX_FIELDS:
        report.add(ERROR, "format",
                   f"{n_fields} fields exceed the DBF {MAX_FIELDS}-field limit")


def _iter_polygon_rings(geom):
    """Yield every LinearRing (exterior + interiors) of a (Multi)Polygon."""
    if geom.geom_type == "Polygon":
        polys = [geom]
    elif geom.geom_type == "MultiPolygon":
        polys = list(geom.geoms)
    else:
        return
    for p in polys:
        if p.exterior is not None:
            yield p.exterior
        yield from p.interiors


def _check_bowties(gdf: gpd.GeoDataFrame, report: ShapefileReport,
                   max_examples: int) -> None:
    """Flag 'bow-tie' features: a ring that crosses itself (is not simple)."""
    geom = gdf.geometry
    bad = []
    for idx, g in geom[~geom.isna()].items():
        if g is None or g.is_empty:
            continue
        try:
            if any(r is not None and not r.is_simple for r in _iter_polygon_rings(g)):
                bad.append(idx)
        except Exception:
            continue
    if bad:
        report.add(WARNING, "geometry",
                   f"{len(bad)} feature(s) have self-intersecting 'bow-tie' rings "
                   f"(a ring crosses itself; e.g. rows {bad[:max_examples]})")


def _signed_area(ring) -> float:
    """Shoelace signed area of a coordinate ring.

    Positive => counter-clockwise, negative => clockwise (standard x-east,
    y-north orientation). Handles rings that are not explicitly closed.
    """
    n = len(ring)
    if n < 3:
        return 0.0
    s = 0.0
    for i in range(n - 1):
        x1, y1 = ring[i][0], ring[i][1]
        x2, y2 = ring[i + 1][0], ring[i + 1][1]
        s += x1 * y2 - x2 * y1
    if ring[0][0] != ring[-1][0] or ring[0][1] != ring[-1][1]:
        s += ring[-1][0] * ring[0][1] - ring[0][0] * ring[-1][1]
    return 0.5 * s


def _ring_roles(rings) -> list[str]:
    """Classify each raw ring as 'outer' or 'hole' using the even-odd rule.

    A ring is a hole when it is nested inside an odd number of the other rings.
    Nesting is decided by testing a *vertex* of each ring against the filled
    area of the others: a vertex of an inner ring sits inside its outer ring,
    while a vertex of an outer ring sits outside an inner one. (Testing an
    interior point instead would misfire for concentric rings, where the outer
    ring's interior point lands inside the hole.) Orientation itself is read
    from the raw signed area, so a self-intersecting ring is still judged from
    its real coordinates.
    """
    from shapely.geometry import Point as _Point, Polygon as _Polygon

    polys, test_points = [], []
    for coords in rings:
        try:
            p = _Polygon(coords)
            if not p.is_valid:
                p = p.buffer(0)
            polys.append(None if p.is_empty else p)
        except Exception:
            polys.append(None)
        try:
            test_points.append(_Point(coords[0][0], coords[0][1]) if coords else None)
        except Exception:
            test_points.append(None)

    roles = []
    for i in range(len(rings)):
        depth = 0
        if test_points[i] is not None:
            for j, pj in enumerate(polys):
                if j != i and pj is not None:
                    try:
                        if pj.contains(test_points[i]):
                            depth += 1
                    except Exception:
                        pass
        roles.append("outer" if depth % 2 == 0 else "hole")
    return roles


def _check_ring_orientation(shp_path: Path, report: ShapefileReport,
                            max_examples: int) -> None:
    """Audit on-disk ring orientation against the ESRI shapefile convention:
    outer rings clockwise, holes counter-clockwise. Reads raw rings with pyshp
    because GDAL silently reorganizes (and can corrupt) orientation on read.
    """
    if _pyshp is None:
        report.add(INFO, "orientation",
                   "On-disk ring-orientation check skipped: install 'pyshp' (pip "
                   "install pyshp) for a definitive audit. Reversed orientation also "
                   "tends to surface as invalid geometry or zero/negative area above.")
        return

    polygon_types = {5, 15, 25}  # Polygon, PolygonZ, PolygonM
    try:
        reader = _pyshp.Reader(str(shp_path))
    except Exception as exc:
        report.add(INFO, "orientation",
                   f"Could not raw-read the shapefile for orientation check: {exc}")
        return

    try:
        if reader.shapeType not in polygon_types:
            return  # ring orientation is meaningful only for polygons

        outer_bad, hole_bad, degenerate = [], [], 0
        for fid, shp in enumerate(reader.iterShapes()):
            if not shp.points or not shp.parts:
                continue
            starts = list(shp.parts) + [len(shp.points)]
            rings = [shp.points[starts[k]:starts[k + 1]]
                     for k in range(len(shp.parts))]
            areas = [_signed_area(r) for r in rings]
            roles = _ring_roles(rings)
            for area, role in zip(areas, roles):
                if abs(area) < 1e-12:
                    degenerate += 1
                    continue
                is_cw = area < 0
                if role == "outer" and not is_cw:
                    outer_bad.append(fid)
                elif role == "hole" and is_cw:
                    hole_bad.append(fid)
    finally:
        try:
            reader.close()
        except Exception:
            pass

    if outer_bad:
        report.add(WARNING, "orientation",
                   f"{len(set(outer_bad))} feature(s) have an OUTER ring wound "
                   "counter-clockwise; the ESRI spec requires outer rings clockwise "
                   f"(e.g. rows {sorted(set(outer_bad))[:max_examples]})")
    if hole_bad:
        report.add(WARNING, "orientation",
                   f"{len(set(hole_bad))} feature(s) have a HOLE ring wound clockwise; "
                   "the ESRI spec requires holes counter-clockwise "
                   f"(e.g. rows {sorted(set(hole_bad))[:max_examples]})")
    if degenerate:
        report.add(WARNING, "orientation",
                   f"{degenerate} ring(s) have ~zero area (degenerate; orientation "
                   "undefined)")
    if not (outer_bad or hole_bad or degenerate):
        report.add(INFO, "orientation",
                   "Polygon rings follow the shapefile convention (outer clockwise, "
                   "holes counter-clockwise)")


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #

def check_shapefile(
    path: Union[str, Path],
    *,
    expected_crs: Optional[Union[str, int, CRS]] = None,
    max_examples: int = 5,
    check_duplicates: bool = True,
    check_bowties: bool = True,
    check_ring_orientation: bool = True,
) -> ShapefileReport:
    """Run a battery of validity checks on a shapefile and return a report.

    Parameters
    ----------
    path
        Path to the ``.shp`` file.
    expected_crs
        If given (e.g. ``"EPSG:4326"`` or ``4326``), warn when the file's CRS
        does not match it.
    max_examples
        How many offending row indices / validity reasons to list per issue.
    check_duplicates
        Whether to check for duplicate geometries and rows (set False for very
        large layers where the WKB comparison is costly).
    check_bowties
        Whether to flag self-intersecting ('bow-tie') rings.
    check_ring_orientation
        Whether to audit on-disk ring orientation against the ESRI spec
        (outer clockwise, holes counter-clockwise). Requires the optional
        'pyshp' package for a definitive result.

    Returns
    -------
    ShapefileReport
        Use ``report.ok`` for a pass/fail, ``print(report)`` for a summary, or
        ``report.errors`` / ``.warnings`` / ``.infos`` for structured access.
    """
    shp_path = Path(path)
    report = ShapefileReport(path=str(shp_path))

    _check_sidecar_files(shp_path, report)

    try:
        gdf = gpd.read_file(shp_path)
    except Exception as exc:
        report.add(ERROR, "read", f"Could not open the shapefile: {exc}")
        return report

    report.stats["n_features"] = len(gdf)
    if len(gdf) == 0:
        report.add(WARNING, "read", "The shapefile contains no features")

    for check in (
        lambda: _check_crs(gdf, report, expected_crs),
        lambda: _check_geometry(gdf, report, max_examples),
        lambda: _check_bowties(gdf, report, max_examples) if check_bowties else None,
        lambda: _check_coordinates(gdf, report),
        lambda: _check_attributes(gdf, report),
        lambda: _check_duplicates(gdf, report) if check_duplicates else None,
        lambda: _check_format_limits(shp_path, gdf, report),
        lambda: (_check_ring_orientation(shp_path, report, max_examples)
                 if check_ring_orientation else None),
    ):
        try:
            check()
        except Exception as exc:  # a broken check should not sink the report
            report.add(WARNING, "internal", f"A check raised an exception: {exc!r}")

    return report


def _main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Check a shapefile for problems.")
    parser.add_argument("shapefile", help="path to the .shp file")
    parser.add_argument("--expected-crs", default=None,
                        help="warn if the CRS differs from this (e.g. EPSG:4326)")
    parser.add_argument("--no-duplicates", action="store_true",
                        help="skip the (slower) duplicate-feature checks")
    parser.add_argument("--no-orientation", action="store_true",
                        help="skip the on-disk ring-orientation check")
    args = parser.parse_args()

    report = check_shapefile(
        args.shapefile,
        expected_crs=args.expected_crs,
        check_duplicates=not args.no_duplicates,
        check_ring_orientation=not args.no_orientation,
    )
    print(report)
    return 0 if report.ok else 1


if __name__ == "__main__":
    raise SystemExit(_main())
