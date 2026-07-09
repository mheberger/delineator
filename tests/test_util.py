"""
Test suite for the util module and the DelineatorConfig dataclass.
Written by the PyCharm AI Assistant.
"""

import math
from pathlib import Path
import pytest
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon

from delineator.settings import DelineatorConfig
from delineator.util import (
    _close_holes,
    _get_overlap_lines,
    _validate_outlet_coordinates,
)


def test_delineator_config_coerces_string_paths_and_normalizes_output_format(tmp_path):
    data_dir = tmp_path / "data"
    output_dir = tmp_path / "outputs"

    config = DelineatorConfig(
        data_dir=str(data_dir),
        output_dir=str(output_dir),
        output_format=".GeoJSON",
        low_res_threshold=100,
        search_dist=1,
        simplify_tolerance=0,
    )

    assert config.data_dir == data_dir
    assert config.output_dir == output_dir
    assert isinstance(config.data_dir, Path)
    assert isinstance(config.output_dir, Path)
    assert config.output_format == "geojson"
    assert config.low_res_threshold == 100.0
    assert config.search_dist == 1.0
    assert config.simplify_tolerance == 0.0


@pytest.mark.parametrize(
    "field_name",
    [
        "auto_download",
        "clean",
        "fill",
        "high_res",
        "outlets",
        "rivers",
        "simplify",
    ],
)
def test_delineator_config_rejects_non_bool_boolean_fields(field_name):
    kwargs = {field_name: "yes"}

    with pytest.raises(TypeError, match=f"{field_name} must be a bool"):
        DelineatorConfig(**kwargs)


@pytest.mark.parametrize(
    ("kwargs", "error_type", "message"),
    [
        ({"data_dir": 123}, TypeError, "data_dir must be a pathlib.Path or string"),
        ({"output_dir": 123}, TypeError, "output_dir must be a pathlib.Path or string"),
        ({"fill_threshold": 1.5}, TypeError, "fill_threshold must be an int"),
        ({"fill_threshold": -1}, ValueError, "fill_threshold must be greater than or equal to 0"),
        ({"num_stream_orders": 1.5}, TypeError, "num_stream_orders must be an int"),
        ({"num_stream_orders": 0}, ValueError, "num_stream_orders must be greater than or equal to 1"),
        ({"low_res_threshold": True}, TypeError, "low_res_threshold must be a number, not bool"),
        ({"low_res_threshold": 0}, ValueError, "low_res_threshold must be greater than 0"),
        ({"search_dist": True}, TypeError, "search_dist must be a number, not bool"),
        ({"search_dist": 0}, ValueError, "search_dist must be greater than 0"),
        ({"simplify_tolerance": True}, TypeError, "simplify_tolerance must be a number, not bool"),
        ({"simplify_tolerance": -0.1}, ValueError, "simplify_tolerance must be greater than or equal to 0"),
        ({"output_format": "definitely_not_supported"}, ValueError, "Unsupported output_format"),
    ],
)
def test_delineator_config_validation_errors(kwargs, error_type, message):
    with pytest.raises(error_type, match=message):
        DelineatorConfig(**kwargs)


def test_validate_outlet_coordinates_accepts_valid_float_coordinates():
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
    with pytest.raises(error_type, match=message):
        _validate_outlet_coordinates(lat, lng)


def test_get_overlap_lines_returns_line_segments_inside_polygon():
    polygon = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
    line = LineString([(-1, 1), (3, 1)])

    result = _get_overlap_lines(line, polygon)

    assert len(result) == 1
    assert result[0].equals(LineString([(0, 1), (2, 1)]))


def test_get_overlap_lines_returns_empty_list_when_no_overlap():
    polygon = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
    line = LineString([(3, 3), (4, 4)])

    result = _get_overlap_lines(line, polygon)

    assert result == []


def test_get_overlap_lines_handles_multiline_overlap():
    polygon = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
    multiline = MultiLineString(
        [
            [(-1, 1), (1, 1)],
            [(1, 1.5), (3, 1.5)],
        ]
    )

    result = _get_overlap_lines(multiline, polygon)

    assert len(result) == 2
    assert all(segment.geom_type == "LineString" for segment in result)


def test_get_overlap_lines_ignores_point_only_touches():
    polygon = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
    line = LineString([(-1, -1), (0, 0)])

    result = _get_overlap_lines(line, polygon)

    assert result == []


def test_close_holes_fills_all_polygon_holes_when_area_max_is_zero():
    polygon_with_hole = Polygon(
        shell=[(0, 0), (5, 0), (5, 5), (0, 5)],
        holes=[[(1, 1), (2, 1), (2, 2), (1, 2)]],
    )

    result = _close_holes(polygon_with_hole, area_max_km2=0)

    assert isinstance(result, Polygon)
    assert len(result.interiors) == 0
    assert result.area == 25



def test_close_holes_applies_to_each_polygon_in_multipolygon():
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


def test_close_holes_rejects_unsupported_geometry_type():
    line = LineString([(0, 0), (1, 1)])

    with pytest.raises(ValueError, match="Unsupported geometry type"):
        _close_holes(line, area_max_km2=0)
