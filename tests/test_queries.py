"""
Tests for the queries module: graceful degradation when data files are
missing or a catchment has no record, instead of raising.
"""

from pathlib import Path

from platformdirs import user_data_dir

import delineator.queries
from delineator.data import _get_data_dir
from delineator.queries import get_rivers, get_upstream_area
from delineator.settings import DelineatorConfig


def test_get_upstream_area_returns_none_when_rivers_db_unavailable(monkeypatch):
    # _find_data_file returns None when the file is missing and cannot be
    # downloaded; get_upstream_area must not crash on it
    monkeypatch.setattr(delineator.queries, "_find_data_file", lambda *args: None)

    result = get_upstream_area(27000119, DelineatorConfig())

    assert result is None


def test_get_upstream_area_returns_none_for_unknown_comid():
    # 27999999 is a valid-looking Iceland comid with no record in the
    # bundled rivers27.db
    result = get_upstream_area(27999999, DelineatorConfig())

    assert result is None


def test_get_upstream_area_returns_value_from_bundled_data():
    result = get_upstream_area(27000119, DelineatorConfig())

    assert result is not None
    assert result > 0


def test_get_rivers_returns_none_when_rivers_db_unavailable(monkeypatch):
    monkeypatch.setattr(delineator.queries, "_find_data_file", lambda *args: None)

    result = get_rivers([27000119], None, DelineatorConfig())

    assert result is None


def test_get_data_dir_ignores_whitespace_only_env_var(monkeypatch):
    monkeypatch.setenv("DELINEATOR_DATA_DIR", "   ")

    result = _get_data_dir()

    assert result == Path(user_data_dir("delineator", appauthor=False))


def test_get_data_dir_strips_whitespace_around_env_var(monkeypatch, tmp_path):
    monkeypatch.setenv("DELINEATOR_DATA_DIR", f"  {tmp_path}  ")

    assert _get_data_dir() == tmp_path
