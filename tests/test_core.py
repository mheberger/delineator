from pathlib import Path
import pandas as pd

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
    config = DelineatorConfig(
        high_res=True,
        output_dir=tmp_path,
        output_format="geojson",
    )

    df = pd.read_csv(TESTS_DIR / "iceland_outlets.csv")

    for row in df.itertuples(index=False):
        watershed_gdf, rivers_gdf, outlets_gdf = delineate(
            float(row.lat),
            float(row.lng),
            config,
        )

        assert watershed_gdf is not None
        assert not watershed_gdf.empty


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
