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
    coords = [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0)]
    out = _catmull_rom_coords(coords)
    assert len(out) > len(coords)
    np.testing.assert_allclose(out[0], coords[0], atol=1e-9)
    np.testing.assert_allclose(out[-1], coords[-1], atol=1e-9)
    assert np.all(np.isfinite(out))


def test_catmull_rom_passes_through_input_vertices():
    coords = [(0.0, 0.0), (1.0, 2.0), (3.0, 1.0), (4.0, 3.0)]
    out = _catmull_rom_coords(coords)
    for vertex in coords:
        dists = np.hypot(out[:, 0] - vertex[0], out[:, 1] - vertex[1])
        assert dists.min() < 1e-9


def test_catmull_rom_handles_consecutive_duplicate_points():
    # Zero-length chords would make the knot sequence non-increasing -> NaNs.
    coords = [(0.0, 0.0), (0.0, 0.0), (1.0, 1.0), (1.0, 1.0), (2.0, 0.0)]
    out = _catmull_rom_coords(coords)
    assert np.all(np.isfinite(out))
    np.testing.assert_allclose(out[0], (0.0, 0.0), atol=1e-9)
    np.testing.assert_allclose(out[-1], (2.0, 0.0), atol=1e-9)


def test_catmull_rom_degenerate_input_returned_unchanged():
    assert len(_catmull_rom_coords([(1.0, 1.0)])) == 1
    assert len(_catmull_rom_coords([(1.0, 1.0), (1.0, 1.0)])) == 1


def test_smooth_river_geometry_linestring():
    line = LineString([(0, 0), (1, 1), (2, 0)])
    out = _smooth_river_geometry(line)
    assert isinstance(out, LineString)
    assert len(out.coords) > len(line.coords)


def test_smooth_river_geometry_multilinestring():
    mls = MultiLineString([[(0, 0), (1, 1)], [(2, 2), (3, 1), (4, 2)]])
    out = _smooth_river_geometry(mls)
    assert isinstance(out, MultiLineString)
    assert len(out.geoms) == 2


def test_smooth_river_geometry_collection_drops_points():
    gc = GeometryCollection([LineString([(0, 0), (1, 1), (2, 0)]), Point(5, 5)])
    out = _smooth_river_geometry(gc)
    assert isinstance(out, LineString)


def test_smooth_river_geometry_passthrough():
    pt = Point(1, 2)
    assert _smooth_river_geometry(pt) is pt
    assert _smooth_river_geometry(None) is None


# ---------------------------------------------------------------------------
# Chaikin (watershed)
# ---------------------------------------------------------------------------

def test_chaikin_ring_stays_closed_and_densifies():
    square = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
    out = _chaikin_closed(square, iterations=2)
    np.testing.assert_allclose(out[0], out[-1])
    assert len(out) > len(square)


def test_smooth_watershed_polygon_valid_and_close_to_original():
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
    outer = [(0, 0), (10, 0), (10, 10), (0, 10)]
    hole = [(4, 4), (6, 4), (6, 6), (4, 6)]
    poly = Polygon(outer, [hole])
    out = _smooth_watershed_geometry(poly)
    assert isinstance(out, Polygon)
    assert out.is_valid
    assert len(out.interiors) == 1


def test_smooth_watershed_multipolygon():
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
