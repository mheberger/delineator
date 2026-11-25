"""
Tests for the queries module: graceful degradation when data files are
missing or a catchment has no record, instead of raising.
"""

import sqlite3
from contextlib import closing
from pathlib import Path

from platformdirs import user_data_dir

import delineator.queries
from delineator.data import _find_data_file, _get_data_dir
from delineator.queries import get_hiearchical_basins, get_rivers, get_upstream_area, get_upstream_comids
from delineator.settings import DelineatorConfig

# The outlet catchment of Iceland's largest watershed in the bundled data. Its
# 184 upstream unit catchments group into a mix of L1 and L0 catchments.
ICELAND_LARGEST_OUTLET = 27001383


def _iceland_basins_db():
    return _find_data_file("basins27.db", DelineatorConfig())


def test_get_hiearchical_basins_partitions_the_upstream_catchments():
    """
    The grouped catchments cover every upstream unit catchment exactly once.

    Each L1-L4 catchment is expanded back to its unit catchments through the
    catchment_hierarchy table. Together with the L0 catchments, they must be
    exactly the upstream list, with nothing missing and nothing counted twice.
    At least one aggregated (L1-L4) catchment must be used, or the grouping
    would not be saving any work.
    """
    with closing(sqlite3.connect(_iceland_basins_db())) as conn:
        upstream = get_upstream_comids(conn, ICELAND_LARGEST_OUTLET)[1:]
        groups = get_hiearchical_basins(conn, upstream)

        covered = list(groups.get("L0", []))
        for level, ids in groups.items():
            if level == "L0":
                continue
            column = f"{level.lower()}_id"
            for group_id in ids:
                covered += [row[0] for row in conn.execute(
                    f"SELECT l0_id FROM catchment_hierarchy WHERE {column} = ?", (group_id,))]

    assert set(groups) - {"L0"}
    assert sorted(covered) == sorted(upstream)


def test_get_hiearchical_basins_leaves_its_connection_reusable():
    """
    get_hiearchical_basins() can run twice on one connection and leaves no temp tables.

    It loads the upstream IDs into a temporary table for the query. That table
    must be dropped afterwards, so a second call on the same connection works
    and gives the same answer.
    """
    with closing(sqlite3.connect(_iceland_basins_db())) as conn:
        upstream = get_upstream_comids(conn, ICELAND_LARGEST_OUTLET)[1:]
        first = get_hiearchical_basins(conn, upstream)
        second = get_hiearchical_basins(conn, upstream)
        leftover_tables = conn.execute("SELECT name FROM sqlite_temp_master").fetchall()

    assert first == second
    assert leftover_tables == []


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
