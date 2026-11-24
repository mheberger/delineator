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
    """
    get_upstream_area() returns None when the rivers database can't be found.

    _find_data_file returns None when a file is missing and can't be
    downloaded. get_upstream_area() must return None instead of crashing, so
    that delineate() can still return a watershed, just without area columns
    (see test_delineate_omits_area_when_upstream_area_unknown in test_core.py).
    """
    monkeypatch.setattr(delineator.queries, "_find_data_file", lambda *args: None)

    result = get_upstream_area(27000119, DelineatorConfig())

    assert result is None


def test_get_upstream_area_returns_none_for_unknown_comid():
    """
    get_upstream_area() returns None for a comid that isn't in the rivers database.

    27999999 looks like a valid Iceland comid but has no record in the bundled
    rivers27.db.
    """
    result = get_upstream_area(27999999, DelineatorConfig())

    assert result is None


def test_get_upstream_area_returns_value_from_bundled_data():
    """
    get_upstream_area() returns a positive area for a real Iceland comid.

    The positive case for the two tests above, using the bundled rivers27.db.
    """
    result = get_upstream_area(27000119, DelineatorConfig())

    assert result is not None
    assert result > 0


def test_get_rivers_returns_none_when_rivers_db_unavailable(monkeypatch):
    """
    get_rivers() returns None when the rivers database can't be found.

    Same situation as test_get_upstream_area_returns_none_when_rivers_db_unavailable.
    delineate() relies on getting None here so it can return the watershed
    without rivers.
    """
    monkeypatch.setattr(delineator.queries, "_find_data_file", lambda *args: None)

    result = get_rivers([27000119], None, DelineatorConfig())

    assert result is None


def test_get_data_dir_ignores_whitespace_only_env_var(monkeypatch):
    """
    A DELINEATOR_DATA_DIR containing only spaces is treated as unset.

    _get_data_dir() must fall back to the platform default rather than use a
    blank path.
    """
    monkeypatch.setenv("DELINEATOR_DATA_DIR", "   ")

    result = _get_data_dir()

    assert result == Path(user_data_dir("delineator", appauthor=False))


def test_get_data_dir_strips_whitespace_around_env_var(monkeypatch, tmp_path):
    """
    Spaces around the DELINEATOR_DATA_DIR path are stripped.

    Stray spaces are easy to pick up when setting an environment variable, for
    example when copying a path.
    """
    monkeypatch.setenv("DELINEATOR_DATA_DIR", f"  {tmp_path}  ")

    assert _get_data_dir() == tmp_path
