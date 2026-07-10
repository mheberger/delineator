import sqlite3
from pathlib import Path
import pandas as pd
import pytest
import shapely

import delineator.core
from delineator.util import write_outputs
from delineator.core import delineate, downloader
from delineator.settings import DelineatorConfig

TESTS_DIR = Path(__file__).parent


def test_delineate_iceland_returns_watershed(tmp_path):
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

    df = pd.read_csv(TESTS_DIR / "iceland_outlets.csv")

    for row in df.itertuples(index=False):
        watershed_gdf, rivers_gdf, outlets_gdf = delineate(
            float(row.lat),
            float(row.lng),
        )

        assert watershed_gdf is not None
        assert not watershed_gdf.empty


def test_delineate_reports_area_and_snapped_outlet(tmp_path):
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


def test_delineate_low_res_returns_watershed_without_snapping(tmp_path):
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
    # get_rivers returns None when no river geometries exist (e.g. coastal
    # catchments). With rivers=True and round_coordinates=True (both defaults)
    # this used to crash in the coordinate-rounding step.
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
    # get_upstream_area returns None when the rivers database is unavailable;
    # the watershed should still be delineated, just without area attributes
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
    # delineate opens sqlite connections to the megabasins, basins, and
    # rivers databases; all of them must be closed by the time it returns,
    # or the .db files stay locked (a real problem on Windows)
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
    # With cache=True, a second delineation with an outlet in the same home
    # unit catchment must not re-run the upstream comid query or re-read the
    # upstream geometries, and must produce the same watershed.
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


def test_downloader_accepts_bundled_iceland_data(tmp_path, capsys):
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
