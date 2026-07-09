from click.testing import CliRunner
from pathlib import Path

from delineator import DelineatorConfig
from delineator.cli import delineator_download, main, delineator_dir


TESTS_DIR = Path(__file__).parent


def test_show_data_dir_exits_successfully():
    runner = CliRunner()

    result = runner.invoke(delineator_dir)

    assert result.exit_code == 0
    assert result.output.strip()


def test_download_data_for_bundled_iceland_data_succeeds():
    runner = CliRunner()

    with runner.isolated_filesystem():
        result = runner.invoke(delineator_download, ["--basin", "27", "--data-dir", "data"])

    assert result.exit_code == 0
    assert "Unit catchments file is at:" in result.output
    assert "Rivers file is at:" in result.output
    assert "Flow Direction file is at:" in result.output
    assert "Flow Accumulation file is at:" in result.output
    assert "basins27.db" in result.output
    assert "rivers27.db" in result.output
    assert "flowdir27.tif" in result.output
    assert "accum27.tif" in result.output


def test_download_data_rejects_invalid_basin():
    runner = CliRunner()

    result = runner.invoke(delineator_download, ["--basin", "999"])

    assert result.exit_code != 0
    assert isinstance(result.exception, ValueError)
    assert "Invalid megabasin ID" in str(result.exception)


def test_main_requires_input_source():
    runner = CliRunner()

    result = runner.invoke(main)

    assert result.exit_code != 0
    assert "Specify either --csv or --point" in result.output


def test_main_rejects_csv_and_point_together():
    runner = CliRunner()
    csv_path = TESTS_DIR / "iceland_outlets.csv"

    result = runner.invoke(
        main,
        [
            "--csv",
            str(csv_path),
            "--point",
            "64.71072",
            "-21.60337",
        ],
    )

    assert result.exit_code != 0
    assert "Specify either --csv or --point, not both" in result.output


def test_main_rejects_missing_csv_file():
    runner = CliRunner()

    result = runner.invoke(main, ["--csv", "does-not-exist.csv"])

    assert result.exit_code != 0
    assert "Invalid value for '--csv'" in result.output
    assert "does-not-exist.csv" in result.output


def test_main_delineates_iceland_csv_to_output_directory(tmp_path):
    runner = CliRunner()
    csv_path = TESTS_DIR / "iceland_outlets.csv"
    output_dir = tmp_path / "output"

    result = runner.invoke(
        main,
        [
            "--csv",
            str(csv_path),
            "--output-dir",
            str(output_dir),
            "--output-format",
            "geojson",
        ],
    )

    assert result.exit_code == 0
    assert "Done. Output written to" in result.output
    assert output_dir.exists()
    assert list(output_dir.glob("*.geojson"))


def test_main_delineates_single_iceland_point_to_output_directory(tmp_path):
    runner = CliRunner()
    output_dir = tmp_path / "output"

    result = runner.invoke(
        main,
        [
            "--point",
            "64.71072",
            "-21.60337",
            "--id",
            "iceland-point",
            "--output-dir",
            str(output_dir),
            "--output-format",
            "geojson",
        ],
    )

    assert result.exit_code == 0
    assert "Done. Output written to" in result.output
    assert output_dir.exists()
    assert list(output_dir.glob("*.geojson"))


def test_delineator_config_uses_output_dir_environment_variable(monkeypatch, tmp_path):
    output_dir = tmp_path / "env-output"

    monkeypatch.setenv("DELINEATOR_OUTPUT_DIR", str(output_dir))

    config = DelineatorConfig()

    assert config.output_dir == output_dir

def test_delineator_config_uses_data_dir_environment_variable(monkeypatch, tmp_path):
    data_dir = tmp_path / "env-data"

    monkeypatch.setenv("DELINEATOR_DATA_DIR", str(data_dir))

    config = DelineatorConfig()

    assert config.data_dir == data_dir

def test_cli_csv(tmp_path):
    runner = CliRunner()
    output_dir = tmp_path / "env-output"

    result = runner.invoke(
        main,
        [
            "--csv",
            TESTS_DIR  / "sample_outlets.csv",
            "--output-dir",
            str(output_dir),
        ],
    )

    assert result.exit_code == 0
    assert output_dir.exists()
    assert list(output_dir.glob("*.gpkg"))
