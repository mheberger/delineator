"""
Optional smoothing ("beautification") of output geometries.

Two algorithms, chosen per geometry type (both tuned by eye in the
Global Watersheds web app, https://mghydro.com/watersheds):

- Rivers (polylines): centripetal Catmull-Rom splines. The curve passes
  through every input vertex, so the river stays glued to its data while the
  angular joints between segments are rounded off. Core math adapted from
  https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline
- Watershed boundary (polygons): Chaikin corner cutting. Each iteration
  replaces every vertex with two points at 1/4 and 3/4 of the adjoining
  edges, slicing corners off closed rings. (Equivalent to PostGIS
  ST_ChaikinSmoothing, implemented here in numpy.)

Both are geometric refinement, not generalization: apply them AFTER
Douglas-Peucker simplification (config.simplify). Smoothing an unsimplified
geometry faithfully traces the raster staircase instead of removing it.
"""
import logging

import numpy as np
from shapely.geometry import (
    GeometryCollection,
    LineString,
    MultiLineString,
    MultiPolygon,
    Polygon,
)

from delineator.constants import (
    CATMULL_ROM_ALPHA,
    CATMULL_ROM_POINTS_PER_SEGMENT,
    CHAIKIN_ITERATIONS,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Catmull-Rom (rivers)
# ---------------------------------------------------------------------------

def _catmull_rom_segment(P0, P1, P2, P3, num_points: int, alpha: float) -> np.ndarray:
    """
    Interpolate one spline segment between control points P1 and P2.

    P0 and P3 shape the curve's tangents at the segment ends. Returns
    `num_points` points from P1 to P2 inclusive, as an (N, 2) array.
    alpha: 0.0 uniform, 0.5 centripetal, 1.0 chordal parameterization.
    """
    def tj(ti: float, pi: np.ndarray, pj: np.ndarray) -> float:
        return ti + float(np.hypot(*(pj - pi))) ** alpha

    t0 = 0.0
    t1 = tj(t0, P0, P1)
    t2 = tj(t1, P1, P2)
    t3 = tj(t2, P2, P3)
    t = np.linspace(t1, t2, num_points).reshape(num_points, 1)

    A1 = (t1 - t) / (t1 - t0) * P0 + (t - t0) / (t1 - t0) * P1
    A2 = (t2 - t) / (t2 - t1) * P1 + (t - t1) / (t2 - t1) * P2
    A3 = (t3 - t) / (t3 - t2) * P2 + (t - t2) / (t3 - t2) * P3
    B1 = (t2 - t) / (t2 - t0) * A1 + (t - t0) / (t2 - t0) * A2
    B2 = (t3 - t) / (t3 - t1) * A2 + (t - t1) / (t3 - t1) * A3
    return (t2 - t) / (t2 - t1) * B1 + (t - t1) / (t2 - t1) * B2


def _catmull_rom_coords(
    coords,
    points_per_segment: int = CATMULL_ROM_POINTS_PER_SEGMENT,
    alpha: float = CATMULL_ROM_ALPHA,
) -> np.ndarray:
    """
    Smooth a polyline's coordinate sequence with a Catmull-Rom spline chain.

    The first and last input vertices are preserved exactly (a phantom
    control point is extrapolated past each end so the spline spans the full
    line). Consecutive duplicate vertices are dropped first: a zero-length
    chord makes the knot sequence non-increasing and produces NaNs.

    Returns an (N, 2) array; if fewer than 2 distinct points remain, the
    deduplicated input is returned unchanged.
    """
    pts = np.asarray(coords, dtype=float)
    if pts.ndim != 2:
        pts = pts.reshape(-1, 2)

    # Drop consecutive duplicates (see docstring).
    keep = np.ones(len(pts), dtype=bool)
    keep[1:] = np.any(pts[1:] != pts[:-1], axis=1)
    pts = pts[keep]

    if len(pts) < 2:
        return pts

    # Phantom control points mirrored past each end.
    first = 2.0 * pts[0] - pts[1]
    last = 2.0 * pts[-1] - pts[-2]
    pts = np.vstack([first, pts, last])

    pieces = []
    n_segments = len(pts) - 3
    for i in range(n_segments):
        seg = _catmull_rom_segment(pts[i], pts[i + 1], pts[i + 2], pts[i + 3],
                                   points_per_segment, alpha)
        # Segments share endpoints; keep the joint vertex only once.
        pieces.append(seg if i == n_segments - 1 else seg[:-1])
    return np.vstack(pieces)


def _smooth_river_geometry(geom):
    """
    Smooth a river reach geometry with Catmull-Rom splines.

    Accepts LineString, MultiLineString, or a GeometryCollection (the split
    of the home reach by the catchment polygon can produce one with stray
    Points, which are dropped). Anything else, or a degenerate result, is
    returned unchanged.
    """
    if geom is None or geom.is_empty:
        return geom

    if isinstance(geom, LineString):
        smoothed = _catmull_rom_coords(geom.coords)
        if len(smoothed) < 2:
            return geom
        return LineString(smoothed)

    if isinstance(geom, (MultiLineString, GeometryCollection)):
        parts = [_smooth_river_geometry(g) for g in geom.geoms
                 if isinstance(g, (LineString, MultiLineString))]
        lines = []
        for p in parts:
            if isinstance(p, LineString):
                lines.append(p)
            elif isinstance(p, MultiLineString):
                lines.extend(p.geoms)
        if not lines:
            return geom
        return lines[0] if len(lines) == 1 else MultiLineString(lines)

    return geom


# ---------------------------------------------------------------------------
# Chaikin (watershed boundary)
# ---------------------------------------------------------------------------

def _chaikin_closed(coords, iterations: int = CHAIKIN_ITERATIONS) -> np.ndarray:
    """
    Smooth a closed ring with Chaikin's corner-cutting algorithm.

    Each iteration replaces every vertex with two points at 1/4 and 3/4
    along the edges that meet there. The ring is treated as cyclic (no
    pinned endpoints). Returns a closed (first == last) (N, 2) array.
    """
    pts = np.asarray(coords, dtype=float)
    if pts.ndim != 2:
        pts = pts.reshape(-1, 2)
    # Work on the open ring.
    if len(pts) > 1 and np.array_equal(pts[0], pts[-1]):
        pts = pts[:-1]

    for _ in range(iterations):
        nxt = np.roll(pts, -1, axis=0)
        cut = np.empty((2 * len(pts), pts.shape[1]), dtype=float)
        cut[0::2] = 0.75 * pts + 0.25 * nxt
        cut[1::2] = 0.25 * pts + 0.75 * nxt
        pts = cut

    return np.vstack([pts, pts[:1]])


def _smooth_watershed_geometry(geom, iterations: int = CHAIKIN_ITERATIONS):
    """
    Smooth a watershed boundary (Polygon or MultiPolygon) with Chaikin
    corner cutting. Interior rings (unfilled holes) are smoothed too.

    Chaikin on a simple ring stays simple, but as a belt-and-suspenders
    measure an invalid result falls back to the input geometry.
    """
    if geom is None or geom.is_empty:
        return geom

    if isinstance(geom, MultiPolygon):
        return MultiPolygon([_smooth_watershed_geometry(p, iterations) for p in geom.geoms])

    if not isinstance(geom, Polygon):
        return geom

    exterior = _chaikin_closed(geom.exterior.coords, iterations)
    interiors = [_chaikin_closed(ring.coords, iterations) for ring in geom.interiors]
    smoothed = Polygon(exterior, interiors)
    if not smoothed.is_valid:
        logger.warning("Chaikin smoothing produced an invalid polygon; leaving it unsmoothed.")
        return geom
    return smoothed
