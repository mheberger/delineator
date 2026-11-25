import logging
import sqlite3
from array import array
from contextlib import closing
from pathlib import Path
import pandas as pd
import pytest
import shapely

import delineator.core
from delineator.data import _find_data_file
from delineator.queries import get_upstream_comids
from delineator.util import write_outputs
from delineator.core import delineate, downloader
from delineator.settings import DelineatorConfig

TESTS_DIR = Path(__file__).parent


def test_delineate_iceland_returns_watershed(tmp_path):
    """
    Basic smoke test: a known Iceland outlet gives a non-empty watershed with a CRS.

    The outlet is the Nordhura River at Stekkur, run in high-res mode with the
    bundled megabasin 27 data, so it needs no network access.
    """
    config = DelineatorConfig(
        high_res=True,
        output_dir=tmp_path,
        output_format="geojson",
    )

    watershed_gdf, rivers_gdf, outlets_gdf = delineate(
        64.71072,
        -21.60337,
        config,
    )
    assert watershed_gdf is not None
    assert not watershed_gdf.empty
    assert watershed_gdf.crs is not None


def test_delineate_invalid_coordinates_returns_none(tmp_path):
    """
    Out-of-range coordinates return (None, None, None) instead of raising an error.

    999 is not a valid latitude or longitude, so the first validation step of
    delineate() rejects it. delineate() must log a warning and return None for
    all three outputs rather than raise.
    """
    config = DelineatorConfig(
        high_res=True,
        output_dir=tmp_path,
        output_format="geojson",
    )

    watershed_gdf, rivers_gdf, outlets_gdf = delineate(
        999.0,
        999.0,
        config,
    )

    assert watershed_gdf is None
    assert rivers_gdf is None
    assert outlets_gdf is None


def test_write_outputs_creates_file_for_iceland(tmp_path):
    """
    write_outputs() writes a GeoJSON file named with the given id to output_dir.

    Delineates the Nordhura River outlet, then checks that a file whose name
    contains the id ("iceland") and ends in .geojson appears in tmp_path.
    """
    config = DelineatorConfig(
        high_res=True,
        output_dir=tmp_path,
        output_format="geojson",
    )

    watershed_gdf, rivers_gdf, outlets_gdf = delineate(
        64.71072,
        -21.60337,
        config,
    )

    assert watershed_gdf is not None

    write_outputs(
        watershed_gdf,
        rivers_gdf,
        outlets_gdf,
        config,
        id="iceland",
    )

    output_files = list(tmp_path.glob("*iceland*.geojson"))
    assert output_files


def test_iceland_csv_outlets_can_be_delineated(tmp_path):
    """
    Every outlet in tests/iceland_outlets.csv gives a non-empty watershed.

    The CSV holds 23 gauged outlets across Iceland, so this covers a range of
    watershed sizes and locations. It passes no config, so it also checks that
    delineate() works with the defaults.
    """
    df = pd.read_csv(TESTS_DIR / "iceland_outlets.csv")

    for row in df.itertuples(index=False):
        watershed_gdf, rivers_gdf, outlets_gdf = delineate(
            float(row.lat),
            float(row.lng),
        )

        assert watershed_gdf is not None
        assert not watershed_gdf.empty


def test_delineate_reports_area_and_snapped_outlet(tmp_path):
    """
    With default options, the watershed has the expected area, a snapped outlet, rivers and outlet points.

    Checks several defaults at once:
    - the snapped outlet is reported and stays near the requested point
    - the area is within 20% of the published 507 km² for the Nordhura River
      at Stekkur
    - the filled_area column is present
    - rivers are returned
    - the outlet points include both the requested and the snapped location
    """
    config = DelineatorConfig(
        high_res=True,
        output_dir=tmp_path,
        output_format="geojson",
    )

    # Nordhura River at Stekkur; published drainage area is 507 km²
    watershed_gdf, rivers_gdf, outlets_gdf = delineate(
        64.71072,
        -21.60337,
        config,
    )

    assert watershed_gdf is not None
    watershed = watershed_gdf.iloc[0]

    # In high-res mode the outlet is snapped to a nearby stream cell
    assert watershed["outlet_lat"] == pytest.approx(64.71072, abs=0.05)
    assert watershed["outlet_lng"] == pytest.approx(-21.60337, abs=0.05)

    # calc_area and fill are on by default
    assert watershed["area_km2"] == pytest.approx(507, rel=0.2)
    assert "filled_area" in watershed_gdf.columns

    # rivers and outlets are returned by default
    assert rivers_gdf is not None
    assert not rivers_gdf.empty
    assert set(outlets_gdf["type"]) == {"requested", "snapped"}


# A point on a small tributary inside unit catchment 27000373 (93 km²). Its two
# upstream neighbours join the mainstem, not this tributary; the tributary's
# local catchment is about 21.6 km², versus 80.8 km² if the 59.2 km² of
# upstream catchments were wrongly included.
TRIBUTARY_OUTLET = (65.75417, -22.29167)


@pytest.mark.parametrize("snapping, threshold", [(True, 1.0), (False, 25.0)])
def test_delineate_tributary_outlet_excludes_upstream_catchments(tmp_path, snapping, threshold):
    """
    An outlet on a tributary leaves out upstream catchments that only drain into the mainstem.

    Regression test for the bug fixed in 2.2.5. The outlet's home unit catchment
    has two upstream neighbours, but they join the mainstem, not this tributary,
    so the watershed must be about 21.6 km² rather than 80.8 km². The rivers must
    be just the home catchment's reach. Runs with snapping on (a 1 km² threshold,
    small enough to keep the snapped point on the tributary) and with snapping
    off (the default 25 km² threshold, which then has no effect).
    """
    config = DelineatorConfig(snapping=snapping, stream_threshold_km2=threshold,
                              output_dir=tmp_path)
    watershed_gdf, rivers_gdf, _ = delineate(*TRIBUTARY_OUTLET, config)

    assert watershed_gdf is not None
    assert watershed_gdf.iloc[0]["area_km2"] == pytest.approx(21.6, abs=2)
    assert list(rivers_gdf["comid"]) == [27000373]


def test_delineate_threshold_above_catchment_drainage_fails_to_snap(tmp_path, caplog):
    """
    A stream threshold larger than anything in the home catchment is not lowered.

    The home unit catchment 27000373 drains about 152 km² at its outlet, so a
    1000 km² threshold leaves no snap targets. Snapping must fail with a warning
    rather than quietly using a smaller threshold.
    """
    config = DelineatorConfig(stream_threshold_km2=1000, output_dir=tmp_path)
    with caplog.at_level(logging.WARNING):
        watershed_gdf, rivers_gdf, outlets_gdf = delineate(*TRIBUTARY_OUTLET, config)

    assert (watershed_gdf, rivers_gdf, outlets_gdf) == (None, None, None)
    assert "No stream cells in the unit catchment" in caplog.text


def test_delineate_mainstem_outlet_keeps_all_upstream_branches(tmp_path):
    """
    An outlet on the mainstem below a confluence includes every upstream branch.

    The outlet is on the reach of 27000015, which both upstream neighbours
    (27000035, 174.8 km², and 27000058, 31.3 km²) drain into, so the watershed
    is about 211.5 km² and the rivers include all three reaches.
    """
    config = DelineatorConfig(output_dir=tmp_path)
    watershed_gdf, rivers_gdf, _ = delineate(66.06333, -15.75833, config)

    assert watershed_gdf.iloc[0]["area_km2"] == pytest.approx(211.5, abs=2)
    assert {27000015, 27000035, 27000058} <= set(rivers_gdf["comid"])


@pytest.mark.parametrize("snapping", [True, False])
def test_delineate_keeps_neighbour_whose_inlet_is_on_the_window_edge(tmp_path, snapping):
    """
    An upstream neighbour whose only inlet pixel is on the raster window's edge is still included.

    Regression test for GitHub issue #13. The only pixel through which 27001625
    drains into 27001458 is in the one-pixel pad around the raster window, where
    the upstream-neighbour test of version 2.2.5 could not see it, so that
    version dropped the branch and reported 746.0 km²; the correct area is
    809.7 km². Runs with snapping both on and off.
    """
    config = DelineatorConfig(snapping=snapping, output_dir=tmp_path)
    watershed_gdf, rivers_gdf, _ = delineate(64.1425, -21.06833333333333, config)

    assert watershed_gdf.iloc[0]["area_km2"] == pytest.approx(809.7, abs=1)
    assert 27001625 in set(rivers_gdf["comid"])


def test_delineate_low_res_returns_watershed_without_snapping(tmp_path):
    """
    Low-res mode builds the watershed from whole unit catchments, with no raster split or snapped outlet.

    Uses the Olfusa River at Selfoss, a large watershed spanning many unit
    catchments. Without the raster step there is no snapped outlet, so the
    outlet_lat/outlet_lng columns are absent and the only outlet point is the
    requested one.
    """
    config = DelineatorConfig(
        high_res=False,
        output_dir=tmp_path,
        output_format="geojson",
    )

    # Olfusa River at Selfoss: a large watershed spanning many unit catchments
    watershed_gdf, rivers_gdf, outlets_gdf = delineate(
        63.93796,
        -21.00666,
        config,
    )

    assert watershed_gdf is not None
    assert not watershed_gdf.empty
    # No raster split in low-res mode, so no snapped outlet is reported
    assert "outlet_lat" not in watershed_gdf.columns
    assert list(outlets_gdf["type"]) == ["requested"]


def test_delineate_simplify_reduces_vertex_count(tmp_path):
    """
    simplify=True gives a watershed boundary with fewer vertices than the default.

    Delineates the same outlet twice, once with simplify off and once with it
    on, and compares the number of vertices in the two polygons.
    """
    detailed_config = DelineatorConfig(
        high_res=True,
        output_dir=tmp_path,
        output_format="geojson",
    )
    simplified_config = DelineatorConfig(
        high_res=True,
        simplify=True,
        output_dir=tmp_path,
        output_format="geojson",
    )

    detailed_gdf, _, _ = delineate(64.71072, -21.60337, detailed_config)
    simplified_gdf, _, _ = delineate(64.71072, -21.60337, simplified_config)

    assert detailed_gdf is not None
    assert simplified_gdf is not None

    detailed_vertices = len(shapely.get_coordinates(detailed_gdf.geometry.iloc[0]))
    simplified_vertices = len(shapely.get_coordinates(simplified_gdf.geometry.iloc[0]))
    assert simplified_vertices < detailed_vertices


def test_delineate_outlet_beyond_snap_distance_returns_none(tmp_path):
    """
    An outlet farther from a stream than search_dist returns (None, None, None).

    search_dist is also the maximum distance the pour point may be snapped. The
    nearest MERIT-Hydro stream cell to the Blanda River gauge at Langamyri is
    about 3 km away, so with search_dist set to 1 km the snap fails and
    delineate() must return None for all three outputs instead of snapping to
    the wrong place.
    """
    config = DelineatorConfig(
        high_res=True,
        search_dist=1.0,  # km
        output_dir=tmp_path,
        output_format="geojson",
    )

    # Blanda River at Langamyri: nearest MERIT-Hydro stream cell is ~3 km away
    watershed_gdf, rivers_gdf, outlets_gdf = delineate(
        65.49006,
        -20.11284,
        config,
    )

    assert watershed_gdf is None
    assert rivers_gdf is None
    assert outlets_gdf is None


def test_delineate_survives_rivers_query_returning_none(tmp_path, monkeypatch):
    """
    delineate() still returns a watershed when get_rivers() finds no rivers.

    get_rivers() returns None when the catchments have no river geometries (for
    example, some coastal catchments). With rivers=True and
    round_coordinates=True (both defaults), delineate() used to crash while
    rounding the river coordinates. The rivers output must be None, and the
    watershed and outlets must still be returned.
    """
    monkeypatch.setattr(delineator.core, "get_rivers", lambda *args: None)

    config = DelineatorConfig(
        rivers=True,
        round_coordinates=True,
        output_dir=tmp_path,
    )

    watershed_gdf, rivers_gdf, outlets_gdf = delineate(
        64.71072,
        -21.60337,
        config,
    )

    assert watershed_gdf is not None
    assert rivers_gdf is None
    assert outlets_gdf is not None


def test_delineate_omits_area_when_upstream_area_unknown(tmp_path, monkeypatch):
    """
    When the upstream area is unknown, the watershed is returned without area columns.

    get_upstream_area() returns None when the rivers database is unavailable.
    The watershed should still be delineated, just without the area_km2 and
    filled_area columns, even though calc_area=True.
    """
    monkeypatch.setattr(delineator.core, "get_upstream_area", lambda *args: None)

    config = DelineatorConfig(
        calc_area=True,
        output_dir=tmp_path,
    )

    watershed_gdf, rivers_gdf, outlets_gdf = delineate(
        64.71072,
        -21.60337,
        config,
    )

    assert watershed_gdf is not None
    assert "area_km2" not in watershed_gdf.columns
    assert "filled_area" not in watershed_gdf.columns


def test_delineate_closes_all_database_connections(tmp_path, monkeypatch):
    """
    Every sqlite connection that delineate() opens is closed by the time it returns.

    delineate() opens connections to the megabasins, basins and rivers databases.
    If any is left open, the .db file stays locked, which is a real problem on
    Windows. The test wraps sqlite3.connect to record each connection, then
    checks that running a query on each one raises ProgrammingError, which only
    happens on a closed connection.
    """
    opened = []
    real_connect = sqlite3.connect

    def tracking_connect(*args, **kwargs):
        conn = real_connect(*args, **kwargs)
        opened.append(conn)
        return conn

    monkeypatch.setattr(sqlite3, "connect", tracking_connect)

    config = DelineatorConfig(output_dir=tmp_path)
    watershed_gdf, _, _ = delineate(64.71072, -21.60337, config)

    assert watershed_gdf is not None
    assert opened  # the tracking wrapper actually saw connections
    for conn in opened:
        # executing on a closed connection raises ProgrammingError
        with pytest.raises(sqlite3.ProgrammingError):
            conn.execute("SELECT 1")


def test_cache_skips_repeat_upstream_queries(tmp_path, monkeypatch):
    """
    With cache=True, a second call for the same outlet reuses cached results and gives the same watershed.

    Counts calls to get_upstream_comids() and get_upstream_geometries(). The
    first delineation must query once; the second, for an outlet in the same
    home unit catchment, must not query again. Both watersheds must have equal
    geometry. The cache is cleared before and after so the test neither depends
    on nor affects other tests.
    """
    calls = {"comids": 0, "geoms": 0}
    real_comids = delineator.core.get_upstream_comids
    real_geoms = delineator.core.get_upstream_geometries

    def counting_comids(*args, **kwargs):
        calls["comids"] += 1
        return real_comids(*args, **kwargs)

    def counting_geoms(*args, **kwargs):
        calls["geoms"] += 1
        return real_geoms(*args, **kwargs)

    monkeypatch.setattr(delineator.core, "get_upstream_comids", counting_comids)
    monkeypatch.setattr(delineator.core, "get_upstream_geometries", counting_geoms)

    delineator.core.clear_cache()
    try:
        config = DelineatorConfig(
            cache=True,
            rivers=False,
            outlets=False,
            output_dir=tmp_path,
        )
        first_gdf, _, _ = delineate(64.71072, -21.60337, config)
        counts_after_first = dict(calls)
        second_gdf, _, _ = delineate(64.71072, -21.60337, config)
    finally:
        delineator.core.clear_cache()

    assert first_gdf is not None
    assert second_gdf is not None
    assert counts_after_first["comids"] == 1
    assert calls == counts_after_first
    assert first_gdf.geometry.iloc[0].equals(second_gdf.geometry.iloc[0])


def test_cache_stores_upstream_comids_as_compact_array(tmp_path):
    """
    The cache keeps upstream comids as an int64 array, not a list of Python ints.

    A large watershed has over 100,000 upstream comids. As a list, each one is a
    separate Python int object, and a full cache of them took more than twice
    the memory of the same IDs in an array. The cached array must hold the same
    comids that get_upstream_comids() returns.
    """
    delineator.core.clear_cache()
    try:
        config = DelineatorConfig(cache=True, rivers=False, outlets=False, output_dir=tmp_path)
        delineate(64.71072, -21.60337, config)
        (key, (cached_ids, _)), = [(k, v) for k, v in delineator.core._CACHE.items()
                                   if k[0] == "upstream"]
    finally:
        delineator.core.clear_cache()

    assert isinstance(cached_ids, array)
    assert cached_ids.typecode == "q"
    db_path = _find_data_file("basins27.db", DelineatorConfig())
    with closing(sqlite3.connect(db_path)) as conn:
        assert list(cached_ids) == get_upstream_comids(conn, key[1])


def test_downloader_accepts_bundled_iceland_data(tmp_path, capsys):
    """
    downloader(27) finds the bundled Iceland files and reports where they are.

    Iceland's data ships with the package, so even with an empty data_dir,
    downloader() must find all four files (basins, rivers, flow direction,
    accumulation) without downloading anything, and print each one's path.
    """
    downloader(27, data_dir=str(tmp_path))

    captured = capsys.readouterr()

    assert "Unit catchments file is at:" in captured.out
    assert "Rivers file is at:" in captured.out
    assert "Flow Direction file is at:" in captured.out
    assert "Flow Accumulation file is at:" in captured.out
    assert "basins27.db" in captured.out
    assert "rivers27.db" in captured.out
    assert "flowdir27.tif" in captured.out
    assert "accum27.tif" in captured.out


def _capture_downloader_paths(monkeypatch, **kwargs):
    """
    Run downloader() for megabasin 11 without network access and return the local paths it would use.

    Replaces _find_data_file with a stub that records the path for each file
    instead of checking for it or downloading it. Extra keyword arguments,
    such as data_dir, are passed through to downloader().
    """
    from delineator.download import _local_path
    paths = []
    monkeypatch.setattr(delineator.core, "_find_data_file",
                        lambda rel, config: paths.append(_local_path(rel, config)))
    downloader(11, **kwargs)
    return paths


def test_downloader_custom_data_dir_is_flat(monkeypatch, tmp_path):
    """
    With a custom data_dir, downloader() puts files directly in that folder.

    User-supplied data directories use a flat layout, with no version subfolder.
    This must match where delineate() looks when given the same data_dir,
    otherwise downloaded files would never be found.
    """
    paths = _capture_downloader_paths(monkeypatch, data_dir=str(tmp_path))

    assert paths
    assert all(p.parent == tmp_path for p in paths)


def test_downloader_default_data_dir_is_versioned(monkeypatch, tmp_path):
    """
    With the default data directory, downloader() puts files in a DATA_DIR_NAME subfolder.

    Default-managed data lives in a version-named subfolder (for example "v3"),
    so regenerated data files land in a fresh folder instead of mixing with old
    ones. The test clears DELINEATOR_DATA_DIR and points the default user data
    directory at tmp_path.
    """
    from delineator.constants import DATA_DIR_NAME
    monkeypatch.delenv("DELINEATOR_DATA_DIR", raising=False)
    monkeypatch.setattr("delineator.settings.user_data_dir", lambda *a, **k: str(tmp_path))

    paths = _capture_downloader_paths(monkeypatch)

    assert paths
    assert all(p.parent == tmp_path / DATA_DIR_NAME for p in paths)
