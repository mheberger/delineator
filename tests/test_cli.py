import geopandas as gpd
import pytest
from click.testing import CliRunner
from pathlib import Path
from shapely.geometry import Polygon

import delineator.cli
from delineator.cli import delineator_download, main, delineator_dir, serve_cli
from delineator.settings import DelineatorConfig


TESTS_DIR = Path(__file__).parent

# Nordhura River at Stekkur, in the bundled Iceland data
ICELAND_POINT = ["64.71072", "-21.60337"]


@pytest.fixture
def fake_delineation(monkeypatch):
    """
    Replace delineate() and write_outputs() in the CLI with fast stubs, and return the configs they were called with.

    Lets the option tests check how each command-line flag is passed to
    DelineatorConfig without running a real delineation or writing files.
    The stub returns a small square watershed, so the CLI sees a success
    and exits 0.
    """
    configs = []
    square = gpd.GeoDataFrame(geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])], crs="EPSG:4326")

    def stub_delineate(lat, lng, config):
        configs.append(config)
        return square, None, None

    monkeypatch.setattr(delineator.cli, "delineate", stub_delineate)
    monkeypatch.setattr(delineator.cli, "write_outputs", lambda *args, **kwargs: None)
    return configs


# ── delineator_dir ───────────────────────────────

def test_show_data_dir_default_shows_versioned_files_folder(monkeypatch, tmp_path):
    """
    With the default data directory, delineator_dir shows the DATA_DIR_NAME subfolder as where files live.

    Clears DELINEATOR_DATA_DIR and points the platform default at tmp_path. The
    "Data directory" line must be the base folder, and the "Data files" line
    must be its version subfolder (for example "v3"), matching where
    downloader() and delineate() put and look for files.
    """
    from delineator.constants import DATA_DIR_NAME
    monkeypatch.delenv("DELINEATOR_DATA_DIR", raising=False)
    monkeypatch.setattr("delineator.settings.user_data_dir", lambda *a, **k: str(tmp_path))

    result = CliRunner().invoke(delineator_dir)

    assert result.exit_code == 0
    assert f"Data directory: {tmp_path}\n" in result.output
    assert f"Data files:     {tmp_path / DATA_DIR_NAME}\n" in result.output


def test_show_data_dir_env_var_is_flat(monkeypatch, tmp_path):
    """
    With DELINEATOR_DATA_DIR set, delineator_dir says so and shows a flat layout.

    A user-supplied data directory has no version subfolder, so the "Data
    files" line must be the same folder, and the "Data directory" line must
    note that it came from the environment variable.
    """
    monkeypatch.setenv("DELINEATOR_DATA_DIR", str(tmp_path))

    result = CliRunner().invoke(delineator_dir)

    assert result.exit_code == 0
    assert f"Data directory: {tmp_path} (set by DELINEATOR_DATA_DIR)\n" in result.output
    assert f"Data files:     {tmp_path}\n" in result.output


# ── delineator_download ──────────────────────────

def test_download_data_for_bundled_iceland_data_succeeds(tmp_path):
    """
    delineator_download --basin 27 finds the bundled Iceland files and prints their paths.

    Iceland's data ships with the package, so this must succeed without
    downloading anything, even though the --data-dir folder doesn't exist.
    This is the command-line counterpart of
    test_downloader_accepts_bundled_iceland_data in test_core.py.
    """
    result = CliRunner().invoke(
        delineator_download, ["--basin", "27", "--data-dir", str(tmp_path / "data")]
    )

    assert result.exit_code == 0
    assert "Unit catchments file is at:" in result.output
    assert "Rivers file is at:" in result.output
    assert "Flow Direction file is at:" in result.output
    assert "Flow Accumulation file is at:" in result.output
    assert "basins27.db" in result.output
    assert "rivers27.db" in result.output
    assert "flowdir27.tif" in result.output
    assert "accum27.tif" in result.output


def test_download_data_rejects_invalid_basin(tmp_path):
    """
    delineator_download with an ID that isn't a megabasin fails with a usage error.

    The ID is checked before anything is downloaded. Click's usage errors exit
    with status 2, and the message must name the bad value.
    """
    result = CliRunner().invoke(
        delineator_download, ["--basin", "999", "--data-dir", str(tmp_path)]
    )

    assert result.exit_code == 2
    assert "999 is not a valid megabasin ID" in result.output
    assert not list(tmp_path.iterdir())


def test_download_data_exits_with_error_when_files_unavailable(monkeypatch, tmp_path):
    """
    delineator_download exits with status 1 when the data files can't be found or downloaded.

    Downloads are switched off with DELINEATOR_AUTO_DOWNLOAD=0 and the data
    directory is empty, so none of megabasin 11's four files can be found.
    The command must report which files are missing and exit 1, so scripts
    that pre-download data can tell that it failed.
    """
    monkeypatch.setenv("DELINEATOR_AUTO_DOWNLOAD", "0")

    result = CliRunner().invoke(delineator_download, ["--basin", "11", "--data-dir", str(tmp_path)])

    assert result.exit_code == 1
    assert "Could not find or download: Unit catchments, Rivers, Flow Direction, Flow Accumulation" in result.output


# ── delineate: input checks and failures ────────

def test_main_requires_input_source():
    """
    delineate with neither --csv nor --point fails with a usage error.
    """
    result = CliRunner().invoke(main)

    assert result.exit_code == 2
    assert "Specify either --csv or --point" in result.output


def test_main_rejects_csv_and_point_together():
    """
    delineate with both --csv and --point fails with a usage error.

    The two input modes are mutually exclusive, so the error message must say
    to use one or the other.
    """
    csv_path = TESTS_DIR / "iceland_outlets.csv"

    result = CliRunner().invoke(main, ["--csv", str(csv_path), "--point", *ICELAND_POINT])

    assert result.exit_code == 2
    assert "Specify either --csv or --point, not both" in result.output


def test_main_rejects_missing_csv_file():
    """
    delineate --csv with a file that doesn't exist fails before any delineation starts.

    Click's own path validation rejects it, so the error names the --csv option
    and the missing file.
    """
    result = CliRunner().invoke(main, ["--csv", "does-not-exist.csv"])

    assert result.exit_code == 2
    assert "Invalid value for '--csv'" in result.output
    assert "does-not-exist.csv" in result.output


def test_main_rejects_csv_without_required_columns(tmp_path):
    """
    delineate --csv with a file that lacks the id, lat and lon columns exits with status 1.

    The CSV exists, so click accepts it, but the CLI's own validation rejects
    it. This used to print "Invalid input" and still exit 0.
    """
    csv_path = tmp_path / "bad.csv"
    csv_path.write_text("foo,bar\n1,2\n")

    result = CliRunner().invoke(main, ["--csv", str(csv_path), "--output-dir", str(tmp_path / "output")])

    assert result.exit_code == 1
    assert "Invalid input" in result.output


def test_main_exits_with_error_when_no_watershed_is_delineated(tmp_path):
    """
    delineate exits with status 1 when none of the outlets can be delineated.

    Null Island (0, 0) is in the ocean, in no megabasin. With every outlet
    failing, the CLI must say that nothing was written and exit 1. This used
    to print "No outputs written" and still exit 0.
    """
    output_dir = tmp_path / "output"

    result = CliRunner().invoke(main, ["--point", "0", "0", "--output-dir", str(output_dir)])

    assert result.exit_code == 1
    assert "No outputs written" in result.output
    assert not list(output_dir.glob("*"))


def test_main_parquet_without_pyarrow_fails_before_delineating(fake_delineation, no_pyarrow, tmp_path):
    """
    --output-format parquet without pyarrow exits 1 with the install hint, before delineating anything.

    pyarrow is the optional extra delineator[parquet]. Without the early
    check, a long delineation would finish and then fail at the write.
    """
    result = CliRunner().invoke(main, ["--point", *ICELAND_POINT, "--output-dir", str(tmp_path),
                                       "--output-format", "parquet"])

    assert result.exit_code == 1
    assert isinstance(result.exception, SystemExit)
    assert 'pip install "delineator[parquet]"' in result.output
    assert fake_delineation == []


def test_main_rejects_unknown_output_format_without_traceback(fake_delineation, tmp_path):
    """An unknown --output-format is a one-line error and exit 1, not a ValueError traceback."""
    result = CliRunner().invoke(main, ["--point", *ICELAND_POINT, "--output-dir", str(tmp_path),
                                       "--output-format", "foo"])

    assert result.exit_code == 1
    assert isinstance(result.exception, SystemExit)
    assert "Unsupported output_format: 'foo'" in result.output
    assert fake_delineation == []


def _locked_file_error(filename):
    """The PermissionError pyogrio raises on Windows when an output file is open in another program."""
    return PermissionError(13, "The process cannot access the file because it is being used by another process",
                           filename)


def test_main_reports_locked_output_file_without_traceback(fake_delineation, monkeypatch, tmp_path, caplog):
    """
    A locked output file gives a short warning naming the file and exits 1, not a traceback.

    On Windows, a shapefile or GeoJSON open in QGIS or ArcGIS can't be
    overwritten, and pyogrio raises PermissionError. This used to crash the
    CLI with a long traceback.
    """
    locked = tmp_path / "outlets.shp"

    def locked_write(*args, **kwargs):
        raise _locked_file_error(str(locked))

    monkeypatch.setattr(delineator.cli, "write_outputs", locked_write)

    result = CliRunner().invoke(main, ["--point", *ICELAND_POINT, "--output-dir", str(tmp_path)])

    assert result.exit_code == 1
    assert isinstance(result.exception, SystemExit)
    assert "No outputs written" in result.output
    assert str(locked) in caplog.text
    assert "close it and try again" in caplog.text


def test_main_csv_continues_after_locked_output_file(fake_delineation, monkeypatch, tmp_path, caplog):
    """In a batch, a locked file for one outlet is skipped and the other outlets are still written."""
    csv_path = tmp_path / "outlets.csv"
    csv_path.write_text("id,lat,lon\na,64.7,-21.6\nb,64.8,-21.7\n")
    written = []

    def write_unless_a(*args, id=None, **kwargs):
        if id == "a":
            raise _locked_file_error(str(tmp_path / "watershed_a.gpkg"))
        written.append(id)

    monkeypatch.setattr(delineator.cli, "write_outputs", write_unless_a)

    result = CliRunner().invoke(main, ["--csv", str(csv_path), "--output-dir", str(tmp_path)])

    assert result.exit_code == 0, result.output
    assert written == ["b"]
    assert "watershed_a.gpkg" in caplog.text


def test_main_no_downloads_fails_instead_of_downloading(tmp_path):
    """
    With --no-downloads and an empty --data-dir, an outlet outside Iceland fails without downloading anything.

    The outlet is the Fulda River in Germany (megabasin 23), whose data isn't
    in the empty data directory. With downloads off, delineation must fail
    with exit status 1 and nothing may be written to the data directory.
    """
    data_dir = tmp_path / "data"
    data_dir.mkdir()

    result = CliRunner().invoke(main, [
        "--point", "51.004", "9.72",
        "--data-dir", str(data_dir),
        "--no-downloads",
        "--output-dir", str(tmp_path / "output"),
    ])

    assert result.exit_code == 1
    assert "No outputs written" in result.output
    assert not list(data_dir.iterdir())


# ── delineate: options ───────────────────────────

def test_main_default_options(fake_delineation, tmp_path):
    """
    With no options, delineate uses the same defaults as DelineatorConfig.

    This test pins the command-line defaults so that any change to them is
    deliberate, and checks each against DelineatorConfig's own default.
    """
    result = CliRunner().invoke(main, ["--point", *ICELAND_POINT, "--output-dir", str(tmp_path)])

    assert result.exit_code == 0
    config = fake_delineation[0]
    assert config.high_res is True
    assert config.rivers is True
    assert config.outlets is True
    assert config.fill is True
    assert config.simplify is False
    assert config.verbose is False
    assert config.output_format == "gpkg"

    python_defaults = DelineatorConfig()
    for name in ("high_res", "rivers", "outlets", "fill", "simplify", "verbose", "output_format"):
        assert getattr(config, name) == getattr(python_defaults, name), name


@pytest.mark.parametrize(
    ("args", "attribute", "expected"),
    [
        (["--low-res"], "high_res", False),
        (["--fill"], "fill", True),
        (["--fill", "false"], "fill", False),
        (["--fill", "true"], "fill", True),
        (["--rivers"], "rivers", True),
        (["--rivers", "false"], "rivers", False),
        (["--rivers", "true"], "rivers", True),
        (["--outlets"], "outlets", True),
        (["--outlets", "false"], "outlets", False),
        (["--outlets", "true"], "outlets", True),
        (["--simplify"], "simplify", True),
        (["--verbose"], "verbose", True),
        (["--no-downloads"], "auto_download", False),
        (["--auto-download"], "auto_download", True),
        (["--output-format", "geojson"], "output_format", "geojson"),
    ],
)
def test_main_passes_option_to_config(fake_delineation, tmp_path, args, attribute, expected):
    """
    Each command-line option sets the matching DelineatorConfig field.

    Runs once per option. --verbose is included because it used to have no
    effect: it was never passed to DelineatorConfig.
    """
    result = CliRunner().invoke(main, ["--point", *ICELAND_POINT, "--output-dir", str(tmp_path), *args])

    assert result.exit_code == 0
    assert getattr(fake_delineation[0], attribute) == expected


def test_main_passes_data_dir_and_output_dir_to_config(fake_delineation, tmp_path):
    """
    --data-dir and --output-dir set the config's data_dir and output_dir.

    A data directory given on the command line counts as custom, so it uses
    the flat layout with no version subfolder.
    """
    data_dir = tmp_path / "data"
    output_dir = tmp_path / "output"

    result = CliRunner().invoke(main, [
        "--point", *ICELAND_POINT, "--data-dir", str(data_dir), "--output-dir", str(output_dir),
    ])

    assert result.exit_code == 0
    config = fake_delineation[0]
    assert config.data_dir == data_dir
    assert config.custom_data_dir is True
    assert config.output_dir == output_dir


# ── delineate: real runs on the bundled Iceland data ──

def test_main_delineates_iceland_csv_to_output_directory(tmp_path):
    """
    Batch mode writes one GeoJSON watershed file for each outlet in the Iceland CSV.

    Uses the 23 outlets in tests/iceland_outlets.csv, which only need the
    bundled data. With --rivers false and --outlets false, there must be
    exactly 23 files, all watershed_<id>.geojson. The output
    directory doesn't exist beforehand, so this also checks that the CLI
    creates it.
    """
    csv_path = TESTS_DIR / "iceland_outlets.csv"
    output_dir = tmp_path / "output"
    ids = [line.split(",")[0] for line in csv_path.read_text().splitlines()[1:] if line]

    result = CliRunner().invoke(main, [
        "--csv", str(csv_path), "--output-dir", str(output_dir), "--output-format", "geojson",
        "--rivers", "false", "--outlets", "false",
    ])

    assert result.exit_code == 0
    assert "Done. Output written to" in result.output
    assert sorted(p.name for p in output_dir.iterdir()) == sorted(f"watershed_{i}.geojson" for i in ids)


def test_main_delineates_single_iceland_point_to_output_directory(tmp_path):
    """
    Single-point mode writes one watershed file named with the --id value.

    Uses the Nordhura River outlet, which only needs the bundled Iceland data.
    Rivers and outlets are turned off, so only the watershed is written.
    """
    output_dir = tmp_path / "output"

    result = CliRunner().invoke(main, [
        "--point", *ICELAND_POINT, "--id", "iceland-point",
        "--output-dir", str(output_dir), "--output-format", "geojson",
        "--rivers", "false", "--outlets", "false",
    ])

    assert result.exit_code == 0
    assert "Done. Output written to" in result.output
    assert [p.name for p in output_dir.iterdir()] == ["watershed_iceland-point.geojson"]


def test_main_writes_rivers_and_outlets_by_default(tmp_path):
    """
    With no --rivers or --outlets option, rivers and outlets each get their own output file.

    With GeoJSON output, each layer is written to its own file. The outlets
    file must hold both the requested and the snapped outlet, since high-res
    mode snaps the outlet to a stream.
    """
    output_dir = tmp_path / "output"

    result = CliRunner().invoke(main, [
        "--point", *ICELAND_POINT, "--id", "x",
        "--output-dir", str(output_dir), "--output-format", "geojson",
    ])

    assert result.exit_code == 0
    assert sorted(p.name for p in output_dir.iterdir()) == [
        "outlets_x.geojson", "rivers_x.geojson", "watershed_x.geojson",
    ]
    assert not gpd.read_file(output_dir / "rivers_x.geojson").empty
    assert set(gpd.read_file(output_dir / "outlets_x.geojson")["type"]) == {"requested", "snapped"}


# ── delineator_serve ─────────────────────────────

def test_serve_help():
    """
    delineator_serve --help runs without starting the server.

    Catches import errors in the web app and checks that the command is
    wired up, without opening a port.
    """
    result = CliRunner().invoke(serve_cli, ["--help"])

    assert result.exit_code == 0
    assert "Launches a local web server" in result.output


# ── outlets around the world (downloads data) ───

@pytest.mark.network
def test_main_delineates_global_sample_csv(tmp_path):
    """
    Batch mode on outlets around the world writes a watershed for each one it can delineate.

    tests/sample_outlets.csv has outlets in eight megabasins (12, 23, 24, 56,
    63, 71, 73 and 82), so this downloads several GB into the default data
    directory on its first run. It is skipped unless you run
    `pytest -m network`. The CSV also has one outlet on Null Island (0, 0),
    which is in no megabasin: the run must carry on past it, exit 0, and
    write one GeoPackage for each of the other nine outlets.
    """
    csv_path = TESTS_DIR / "sample_outlets.csv"
    output_dir = tmp_path / "output"
    ids = [line.split(",")[0] for line in csv_path.read_text().splitlines()[1:] if line]
    expected = sorted(f"watershed_{i}.gpkg" for i in ids if i != "9999999")

    result = CliRunner().invoke(main, ["--csv", str(csv_path), "--output-dir", str(output_dir)])

    assert result.exit_code == 0
    assert sorted(p.name for p in output_dir.iterdir()) == expected
