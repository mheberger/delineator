"""
Command Line Interface for delineator
"""
from pathlib import Path
import click
import pandas as pd
import logging

from delineator.core import delineate, downloader
from delineator.settings import DelineatorConfig, _resolve_default_data_dir
from delineator.util import write_outputs
from delineator.constants import DATA_DIR_NAME
from delineator.validation import _validate_and_normalize_df
from delineator.serve import serve

logger = logging.getLogger(__name__)


@click.command()
def serve_cli():
    """
    Launches a local web server to serve the delineator app.
    with a simple point-and-click map interface.
    """
    serve()


@click.command()
@click.option("--point", nargs=2, type=float, metavar="LAT LON", help="Creates watershed for a single outlet point")
@click.option("--id", default=None, help="ID for the watershed (used with --point)")
@click.option("--csv", "outlets_csv", type=click.Path(exists=True), help="CSV file of outlet points")
@click.option("--data-dir", default=None, help="Directory for delineator's input data files. "
                                               "Run `delineator_dir` to see the default.")
@click.option("--auto-download/--no-downloads", default=None, help="Automatically download data files if needed")
@click.option("--high-res/--low-res", default=True, help="In higher-resolution mode, uses raster methods"
                                                         " to find the watershed boundary around the outlet point."
                                                         " In lower-resolution mode, the script runs faster but"
                                                         " produces boundary with less accuracy near the outlet point."
                                                         " Default is True (high-res mode).")
@click.option("--output-dir", default=None, help="Directory for output files. Defaults to ./output/")
@click.option("--output-format", default="gpkg",
              help="Output format: geojson, shp, kml, gpkg. Default is gpkg (Geopackage).")
@click.option("--fill", is_flag=True, help="Fill interior gaps or \"donut holes\" in the watershed")
@click.option("--rivers", is_flag=True, help="Output the river network")
@click.option("--outlets", is_flag=True, help="Output the outlet points")
@click.option("--simplify", is_flag=True, help="Simplify the geometry of the watershed and rivers")
@click.option("--verbose", is_flag=True, default=False, help="Show more information on the script's progress")
def main(point, id, outlets_csv, data_dir, auto_download, high_res, output_dir,
         output_format, fill, rivers, outlets, simplify, verbose):
    """Delineate watersheds from a CSV file or a single lat/lon point.

    Usage:
        delineate --help  # List options

        delineate --point 48.834 2.263  # Create watershed for a single outlet

        # Single outlet with a variety of options
        delineate --point 48.863, 2.314 --id alexandreIII --rivers --outlets --fill --simplify --output-dir D:/Projects/watersheds --output-format geojson

        # Batch delienation of a list of outlet points in a CSV file; see docs for correct format
        delineate --csv outlets.csv
    """

    # Pipe the logs to stdout
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(levelname)s %(name)s: %(message)s"))
    logger.addHandler(handler)
    logger.propagate = False
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    # Check that user specified either --csv or --point, not both
    if outlets_csv and point:
        raise click.UsageError("Specify either --csv or --point, not both.")
    if not outlets_csv and not point:
        raise click.UsageError("Specify either --csv or --point.")

    # Normalize input
    if point:
        lat, lon = point
        outlets_df = pd.DataFrame([{"id": id, "lat": lat, "lon": lon}])
    else:
        # Verify that the CSV file exists
        if not Path(outlets_csv).exists():
            click.echo(f"Error: CSV file not found: {outlets_csv}")
            return
        outlets_df = pd.read_csv(outlets_csv)

    # Validate input
    outlets_df = _validate_and_normalize_df(outlets_df)
    if outlets_df is None:
        click.echo("Invalid input. See above for details.")
        return

    # Get the options for the DelineatorConfig object
    config_kwargs = {
        "high_res": high_res,
        "output_format": output_format,
        "fill": fill,
        "rivers": rivers,
        "outlets": outlets,
        "simplify": simplify,
    }

    if data_dir is not None:
        config_kwargs["data_dir"] = Path(data_dir)
    if output_dir is not None:
        config_kwargs["output_dir"] = Path(output_dir)
    if auto_download is not None:
        config_kwargs["auto_download"] = auto_download

    config = DelineatorConfig(**config_kwargs)

    # Set the output directory (may not have been set by user)
    output_path = config.output_dir

    if not output_path.exists():
        click.echo(f"Creating output directory: {output_path}")
        try:
            output_path.mkdir(parents=True, exist_ok=True)
        except FileNotFoundError:
            click.echo(f"ERROR: Could not find or create the directory: {output_path}")
            click.echo("***No outputs written.***")
            return

    # Run watershed delineation
    success = _delineate_dataframe(outlets_df, config)
    if success:
        click.echo(f"Done. Output written to {output_path}/")
    else:
        click.echo("***No outputs written.***")


@click.command()
def delineator_dir():
    """Show where the delineator Python package keeps its data files.

    Prints two paths. "Data directory" is the default location, which depends
    on the operating system, or the location in the DELINEATOR_DATA_DIR
    environment variable, if it is set. "Data files" is the folder the data
    files are actually in. For the default location, this is a version-named
    subfolder of the data directory (e.g. v3/). For DELINEATOR_DATA_DIR, it is
    the same folder. The data for Iceland (megabasin 27) is bundled with the
    package and is not in either folder.

    Usage:
        > delineator_dir
    """
    data_dir, is_custom = _resolve_default_data_dir()
    data_dir.mkdir(parents=True, exist_ok=True)
    # Same rule as download._local_path: custom dirs are flat, the default is versioned.
    files_dir = data_dir if is_custom else data_dir / DATA_DIR_NAME
    source = " (set by DELINEATOR_DATA_DIR)" if is_custom else ""
    click.echo(f"Data directory: {data_dir}{source}")
    click.echo(f"Data files:     {files_dir}")


@click.command()
@click.option("--basin", type=int, required=True, help="The megabasin ID, an integer from 11 to 86")
@click.option("--data-dir", default=None, help="Directory for delineator's data files. "
                                               "Defaults to the system default data directory; "
                                               "run `delineator_dir` to see it.")
def delineator_download(basin, data_dir):
    """
    Download the data files needed for delineation in a megabasin.
    Use the command `delineator_dir` to get the location of the data directory.

    Usage:
        # Download the data files for megabasin 11 (E. Africa) to the default data directory
        delineator_download --basin 11

        # Download data files to your preferred location
        delineator_download --basin 11  --data-dir D:/Projects/watersheds
    """
    # Leave data_dir as None for the default, so the files go in its version
    # subfolder. DelineatorConfig creates the directory if it doesn't exist.
    downloader(int(basin), data_dir)


def _delineate_dataframe(outlets_df, config: DelineatorConfig) -> bool:
    """
    For the command line interface, delineates watersheds based on a
    dataframe of outlet points, derived from the user's CSV file or lat/lon point.

    Parameters
    ----------
    outlets_df : pandas.DataFrame
        DataFrame containing the outlet points. Must have columns "lat", "lon", and "id".

    config : DelineatorConfig

    Returns
    -------
    boolean
        indicating whether any outputs were written to disk.

    """
    num_outputs_written = 0
    for _, row in outlets_df.iterrows():
        lat, lon = row["lat"], row["lon"]
        id = row["id"]
        watershed_gdf, rivers_gdf, outlets_gdf = delineate(lat, lon, config)
        if watershed_gdf is None:
            logging.warning(f"No watershed found for point {id} ({lat}, {lon}).")
            continue
        write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, config, id=id)
        num_outputs_written += 1

    if num_outputs_written > 0:
        return True
    else:
        return False
