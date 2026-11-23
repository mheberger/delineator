"""
Test suite for the util module and the DelineatorConfig dataclass.
Written by the PyCharm AI Assistant.
"""

import math
from pathlib import Path
import geopandas as gpd
import pytest
from shapely.geometry import LineString, MultiPolygon, Polygon

from delineator.constants import KM_PER_DEGREE
from delineator.settings import DelineatorConfig
from delineator.util import (
    _close_holes,
    _ring_area_km2,
    _validate_outlet_coordinates,
    write_outputs,
)


def test_delineator_config_coerces_string_paths_and_normalizes_output_format(tmp_path):
    """
    DelineatorConfig converts its inputs to standard types.

    Checks three conversions:
    - string paths become pathlib.Path
    - ".GeoJSON" becomes "geojson" (leading dot dropped, lowercased)
    - integer values for float options become floats
    """
    data_dir = tmp_path / "data"
    output_dir = tmp_path / "outputs"

    config = DelineatorConfig(
        data_dir=str(data_dir),
        output_dir=str(output_dir),
        output_format=".GeoJSON",
        fill_area_max=2,
        search_dist=1,
        simplify_tolerance=0,
        stream_threshold_km2=30,
    )

    assert config.data_dir == data_dir
    assert config.output_dir == output_dir
    assert isinstance(config.data_dir, Path)
    assert isinstance(config.output_dir, Path)
    assert config.output_format == "geojson"
    assert config.fill_area_max == 2.0
    assert config.search_dist == 1.0
    assert config.simplify_tolerance == 0.0
    assert config.stream_threshold_km2 == 30.0


@pytest.mark.parametrize(
    "field_name",
    [
        "auto_download",
        "clean",
        "fill",
        "high_res",
        "outlets",
        "rivers",
        "round_coordinates",
        "simplify",
        "snapping",
    ],
)
def test_delineator_config_rejects_non_bool_boolean_fields(field_name):
    """
    Boolean options reject non-bool values such as "yes" with a TypeError.

    A truthy string would otherwise quietly switch an option on. Runs once for
    each boolean option listed.
    """
    kwargs = {field_name: "yes"}

    with pytest.raises(TypeError, match=f"{field_name} must be a bool"):
        DelineatorConfig(**kwargs)


@pytest.mark.parametrize(
    ("kwargs", "error_type", "message"),
    [
        ({"data_dir": 123}, TypeError, "data_dir must be a pathlib.Path or string"),
        ({"output_dir": 123}, TypeError, "output_dir must be a pathlib.Path or string"),
        ({"fill_area_max": "big"}, TypeError, "fill_area_max must be a number"),
        ({"fill_area_max": True}, TypeError, "fill_area_max must be a number, not bool"),
        ({"fill_area_max": -1}, ValueError, "fill_area_max must be greater than or equal to 0"),
        ({"num_stream_orders": 1.5}, TypeError, "num_stream_orders must be an int"),
        ({"num_stream_orders": 0}, ValueError, "num_stream_orders must be greater than or equal to 1"),
        ({"search_dist": True}, TypeError, "search_dist must be a number, not bool"),
        ({"search_dist": 0}, ValueError, r"search_dist \(in km\) must be greater than 0"),
        ({"search_dist": 100}, ValueError, r"search_dist \(in km\) must be less than or equal to 50"),
        ({"simplify_tolerance": True}, TypeError, "simplify_tolerance must be a number, not bool"),
        ({"simplify_tolerance": -0.1}, ValueError, "simplify_tolerance must be greater than or equal to 0"),
        ({"stream_threshold_km2": 0}, ValueError, "stream_threshold_km2 must be greater than 0"),
        ({"stream_threshold_km2": 60000}, ValueError, "must be less than or equal to 50000"),
        ({"output_format": "definitely_not_supported"}, ValueError, "Unsupported output_format"),
    ],
)
def test_delineator_config_validation_errors(kwargs, error_type, message):
    """
    Invalid option values raise the right error type with a helpful message.

    Each case sets one bad value. Wrong types raise TypeError, and values out
    of range raise ValueError. bool is rejected for numeric options even though
    it is a subclass of int in Python.
    """
    with pytest.raises(error_type, match=message):
        DelineatorConfig(**kwargs)


def test_delineator_config_treats_data_dir_env_var_as_custom(monkeypatch, tmp_path):
    """
    A data_dir taken from DELINEATOR_DATA_DIR counts as custom, so it uses the flat layout.

    The spaces around the path are stripped, as in
    test_get_data_dir_strips_whitespace_around_env_var in test_queries.py.
    """
    monkeypatch.setenv("DELINEATOR_DATA_DIR", f"  {tmp_path}  ")

    config = DelineatorConfig()

    assert config.custom_data_dir is True
    assert config.data_dir == tmp_path


def test_delineator_config_ignores_whitespace_only_data_dir_env_var(monkeypatch):
    """
    A DELINEATOR_DATA_DIR containing only spaces doesn't count as a custom data directory.
    """
    monkeypatch.setenv("DELINEATOR_DATA_DIR", "   ")

    config = DelineatorConfig()

    # A blank env var falls back to the managed platform default
    assert config.custom_data_dir is False


def test_delineator_config_uses_output_dir_environment_variable(monkeypatch, tmp_path):
    """
    DelineatorConfig takes its default output_dir from DELINEATOR_OUTPUT_DIR.
    """
    output_dir = tmp_path / "env-output"
    monkeypatch.setenv("DELINEATOR_OUTPUT_DIR", str(output_dir))

    config = DelineatorConfig()

    assert config.output_dir == output_dir


def test_validate_outlet_coordinates_accepts_valid_float_coordinates():
    """
    Valid coordinates are returned unchanged, in (lat, lng) order.
    """
    assert _validate_outlet_coordinates(45.0, -120.0) == (45.0, -120.0)


@pytest.mark.parametrize(
    ("lat", "lng", "error_type", "message"),
    [
        (math.inf, -120.0, ValueError, "Latitude must be finite"),
        (45.0, math.nan, ValueError, "Longitude must be finite"),
        (-60.0, -120.0, ValueError, "Latitude must be greater than -60"),
        (85.0, -120.0, ValueError, "Latitude must be less than 85"),
        (45.0, -180.0, ValueError, "Longitude must be greater than -180"),
        (45.0, 180.0, ValueError, "Longitude must be less than 180"),
    ],
)
def test_validate_outlet_coordinates_rejects_invalid_coordinates(lat, lng, error_type, message):
    """
    Infinite, NaN and out-of-range coordinates raise ValueError with a message naming the problem.

    delineate() accepts latitudes between -60 and 85 and longitudes between
    -180 and 180, all exclusive, so each boundary value must be rejected.
    """
    with pytest.raises(error_type, match=message):
        _validate_outlet_coordinates(lat, lng)


def test_close_holes_fills_all_polygon_holes_when_area_max_is_zero():
    """
    With area_max_km2=0, _close_holes fills every hole in a Polygon.

    0 is the special value meaning "fill all holes, whatever their size". The
    result must be the outer ring alone, with the full area of the square.
    """
    polygon_with_hole = Polygon(
        shell=[(0, 0), (5, 0), (5, 5), (0, 5)],
        holes=[[(1, 1), (2, 1), (2, 2), (1, 2)]],
    )

    result = _close_holes(polygon_with_hole, area_max_km2=0)

    assert isinstance(result, Polygon)
    assert len(result.interiors) == 0
    assert result.area == 25



def test_close_holes_applies_to_each_polygon_in_multipolygon():
    """
    _close_holes fills holes in every part of a MultiPolygon and keeps all the parts.
    """
    polygon_with_hole = Polygon(
        shell=[(0, 0), (5, 0), (5, 5), (0, 5)],
        holes=[[(1, 1), (2, 1), (2, 2), (1, 2)]],
    )
    polygon_without_hole = Polygon([(10, 10), (12, 10), (12, 12), (10, 12)])
    multipolygon = MultiPolygon([polygon_with_hole, polygon_without_hole])

    result = _close_holes(multipolygon, area_max_km2=0)

    assert isinstance(result, MultiPolygon)
    assert len(result.geoms) == 2
    assert all(len(poly.interiors) == 0 for poly in result.geoms)


def test_ring_area_km2_approximates_equatorial_square():
    """
    _ring_area_km2 gives roughly the right area for a 1° square at the equator.

    At the equator, a degree is about the same distance in either direction,
    so the area should be close to KM_PER_DEGREE squared. This checks the
    degrees-to-km² conversion that _close_holes uses to compare holes with
    area_max_km2.
    """
    ring = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]).exterior

    # A 1° x 1° square near the equator is about KM_PER_DEGREE² in km²
    assert _ring_area_km2(ring) == pytest.approx(KM_PER_DEGREE ** 2, rel=0.01)


def test_close_holes_keeps_holes_larger_than_area_max_km2():
    """
    With a positive area_max_km2, only holes at or below that size are filled.

    The polygon has one hole of about 0.3 km² and one of about 12,000 km². With
    a 1 km² limit, the small one must be filled and the large one kept exactly
    as it was.
    """
    # Near the equator, a 0.005° square hole is ~0.3 km²; a 1° square is ~12,000 km²
    small_hole = [(0.2, 0.2), (0.205, 0.2), (0.205, 0.205), (0.2, 0.205)]
    large_hole = [(1, 1), (2, 1), (2, 2), (1, 2)]
    polygon = Polygon(
        shell=[(0, 0), (3, 0), (3, 3), (0, 3)],
        holes=[small_hole, large_hole],
    )

    result = _close_holes(polygon, area_max_km2=1.0)

    assert len(result.interiors) == 1
    assert Polygon(result.interiors[0]).equals(Polygon(large_hole))


def test_close_holes_rejects_unsupported_geometry_type():
    """
    _close_holes raises ValueError for geometry that isn't a Polygon or MultiPolygon.
    """
    line = LineString([(0, 0), (1, 1)])

    with pytest.raises(ValueError, match="Unsupported geometry type"):
        _close_holes(line, area_max_km2=0)


def test_write_outputs_parquet_passes_path_per_layer(tmp_path, monkeypatch):
    """
    For parquet output, write_outputs passes a single file path for each layer.

    A regression once passed a tuple instead of a path. to_parquet needs
    pyarrow, which may not be installed, so the test replaces it with a stub
    that records the path it was given. The only layer here is the watershed,
    which must go to watershed_<id>.parquet in output_dir.
    """
    paths = []

    def fake_to_parquet(self, path, **kwargs):
        paths.append(path)

    monkeypatch.setattr(gpd.GeoDataFrame, "to_parquet", fake_to_parquet)

    watershed_gdf = gpd.GeoDataFrame(
        geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
        crs="EPSG:4326",
    )
    config = DelineatorConfig(output_dir=tmp_path, output_format="parquet")

    write_outputs(watershed_gdf, None, None, config, id="test")

    assert paths == [tmp_path / "watershed_test.parquet"]
    assert all(isinstance(p, Path) for p in paths)
