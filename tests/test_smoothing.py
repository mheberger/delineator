import numpy as np
import pytest
from shapely.geometry import (
    GeometryCollection,
    LineString,
    MultiLineString,
    MultiPolygon,
    Point,
    Polygon,
)

from delineator.core import delineate
from delineator.settings import DelineatorConfig
from delineator.smoothing import (
    _catmull_rom_coords,
    _chaikin_closed,
    _smooth_river_geometry,
    _smooth_watershed_geometry,
)


# ---------------------------------------------------------------------------
# Catmull-Rom (rivers)
# ---------------------------------------------------------------------------

def test_catmull_rom_preserves_endpoints_and_densifies():
    """
    The Catmull-Rom spline adds points, keeps the river's two ends fixed and produces no NaNs.

    The ends must not move because a river reach has to stay joined to the
    reaches above and below it, and to the outlet.
    """
    coords = [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0)]
    out = _catmull_rom_coords(coords)
    assert len(out) > len(coords)
    np.testing.assert_allclose(out[0], coords[0], atol=1e-9)
    np.testing.assert_allclose(out[-1], coords[-1], atol=1e-9)
    assert np.all(np.isfinite(out))


def test_catmull_rom_passes_through_input_vertices():
    """
    The smoothed line passes through every original vertex.

    This is the property that makes Catmull-Rom suitable for rivers: the curve
    rounds off the corners between segments but stays on the data, unlike
    corner cutting, which moves vertices.
    """
    coords = [(0.0, 0.0), (1.0, 2.0), (3.0, 1.0), (4.0, 3.0)]
    out = _catmull_rom_coords(coords)
    for vertex in coords:
        dists = np.hypot(out[:, 0] - vertex[0], out[:, 1] - vertex[1])
        assert dists.min() < 1e-9


def test_catmull_rom_handles_consecutive_duplicate_points():
    """
    Repeated consecutive vertices don't produce NaNs.

    A zero-length segment would make the spline's knot sequence
    non-increasing, which produces NaNs. _catmull_rom_coords drops consecutive
    duplicates first, so the output must be finite and still start and end at
    the original endpoints.
    """
    coords = [(0.0, 0.0), (0.0, 0.0), (1.0, 1.0), (1.0, 1.0), (2.0, 0.0)]
    out = _catmull_rom_coords(coords)
    assert np.all(np.isfinite(out))
    np.testing.assert_allclose(out[0], (0.0, 0.0), atol=1e-9)
    np.testing.assert_allclose(out[-1], (2.0, 0.0), atol=1e-9)


def test_catmull_rom_degenerate_input_returned_unchanged():
    """
    A single point, or two identical points, is returned as one point instead of raising.

    After duplicates are removed, both inputs have fewer than two distinct
    points, so there is no line to smooth.
    """
    assert len(_catmull_rom_coords([(1.0, 1.0)])) == 1
    assert len(_catmull_rom_coords([(1.0, 1.0), (1.0, 1.0)])) == 1


def test_smooth_river_geometry_linestring():
    """
    A LineString river comes back as a LineString with more vertices.
    """
    line = LineString([(0, 0), (1, 1), (2, 0)])
    out = _smooth_river_geometry(line)
    assert isinstance(out, LineString)
    assert len(out.coords) > len(line.coords)


def test_smooth_river_geometry_multilinestring():
    """
    A MultiLineString river comes back as a MultiLineString with the same number of parts.
    """
    mls = MultiLineString([[(0, 0), (1, 1)], [(2, 2), (3, 1), (4, 2)]])
    out = _smooth_river_geometry(mls)
    assert isinstance(out, MultiLineString)
    assert len(out.geoms) == 2


def test_smooth_river_geometry_collection_drops_points():
    """
    Stray Points in a GeometryCollection are dropped, leaving a single LineString.

    Splitting the home reach with the catchment polygon can produce a
    GeometryCollection with stray Points in it. Smoothing keeps only the line
    parts, and returns a plain LineString when only one is left.
    """
    gc = GeometryCollection([LineString([(0, 0), (1, 1), (2, 0)]), Point(5, 5)])
    out = _smooth_river_geometry(gc)
    assert isinstance(out, LineString)


def test_smooth_river_geometry_passthrough():
    """
    Geometry that isn't a line, and None, are returned unchanged.

    The same object must come back (not a copy), so callers can apply the
    function to every row without checking the type first.
    """
    pt = Point(1, 2)
    assert _smooth_river_geometry(pt) is pt
    assert _smooth_river_geometry(None) is None


# ---------------------------------------------------------------------------
# Chaikin (watershed)
# ---------------------------------------------------------------------------

def test_chaikin_ring_stays_closed_and_densifies():
    """
    Chaikin corner cutting keeps a ring closed and adds vertices.

    The first and last points must still be equal, or the ring can't be used
    as a polygon boundary.
    """
    square = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
    out = _chaikin_closed(square, iterations=2)
    np.testing.assert_allclose(out[0], out[-1])
    assert len(out) > len(square)


def test_smooth_watershed_polygon_valid_and_close_to_original():
    """
    A smoothed polygon is valid, a bit smaller than the original, and inside it.

    Corner cutting only ever removes area, so the result must be within the
    original shape and somewhat smaller. A square is the worst case, so the
    lower bound on the area ratio is set for it; see the comments below.
    """
    square = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    out = _smooth_watershed_geometry(square)
    assert isinstance(out, Polygon)
    assert out.is_valid
    # Corner cutting shrinks the shape. A 4-vertex square is the worst case
    # (its corners ARE the shape; the ratio converges to 5/6). Real watershed
    # rings have hundreds of shallow-angle vertices and lose far less.
    assert 0.83 < out.area / square.area < 1.0
    # The smoothed shape stays inside the original (corner cutting is convex-hull-contained).
    assert square.buffer(1e-9).contains(out)


def test_smooth_watershed_preserves_holes():
    """
    Holes in a watershed polygon survive smoothing.

    Holes left unfilled (fill=False, or larger than fill_area_max) are smoothed
    too, and must stay as holes in a valid polygon.
    """
    outer = [(0, 0), (10, 0), (10, 10), (0, 10)]
    hole = [(4, 4), (6, 4), (6, 6), (4, 6)]
    poly = Polygon(outer, [hole])
    out = _smooth_watershed_geometry(poly)
    assert isinstance(out, Polygon)
    assert out.is_valid
    assert len(out.interiors) == 1


def test_smooth_watershed_multipolygon():
    """
    Each part of a MultiPolygon is smoothed, and the result is still a valid MultiPolygon with the same number of parts.
    """
    mp = MultiPolygon([
        Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
        Polygon([(5, 5), (6, 5), (6, 6), (5, 6)]),
    ])
    out = _smooth_watershed_geometry(mp)
    assert isinstance(out, MultiPolygon)
    assert len(out.geoms) == 2
    assert out.is_valid


# ---------------------------------------------------------------------------
# Integration: delineate() with smooth=True (bundled Iceland data)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("simplify", [True, False])
def test_delineate_smooth_iceland(tmp_path, simplify):
    """
    delineate() with smooth=True returns a valid watershed and valid rivers.

    Runs on the Nordhura River outlet with simplify both on and off. Smoothing
    runs after simplification, so this checks it on both a simplified boundary
    and the full-detail boundary that follows the raster's pixel edges.
    """
    config = DelineatorConfig(
        output_dir=tmp_path,
        simplify=simplify,
        smooth=True,
    )
    watershed_gdf, rivers_gdf, outlets_gdf = delineate(64.71072, -21.60337, config)

    assert watershed_gdf is not None
    assert watershed_gdf.geometry.iloc[0].is_valid
    assert rivers_gdf is not None
    assert rivers_gdf.geometry.is_valid.all()
