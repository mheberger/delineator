from pathlib import Path
import click
import pandas as pd
import logging

from delineator.core import delineate, downloader
from delineator.settings import DelineatorConfig
from delineator.util import write_outputs
from delineator.data import _get_data_dir
from delineator.validation import _validate_and_normalize_df


@click.command()
@click.option("--point", nargs=2, type=float, metavar="LAT LON", help="Single outlet point")
@click.option("--id", default=None, help="ID for the watershed (used with --point)")
@click.option("--csv", "outlets_csv", type=click.Path(exists=True), help="CSV file of outlet points")
@click.option("--data-dir", default=None, help="Directory for delineator's input data files. ")
@click.option("--auto-download/--no-downloads", is_flag=True, default=True, help="Automatically download data files if needed")
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
        delineate --help
        delineate --point 48.834 2.263
        delineate --csv outlets.csv
    """

    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(levelname)s %(name)s: %(message)s"))
    logger = logging.getLogger("delineator")
    logger.addHandler(handler)
    logger.propagate = False
    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

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

    config_kwargs = {
        "high_res": high_res,
        "output_format": output_format,
        "fill": fill,
        "rivers": rivers,
        "outlets": outlets,
        "auto_download": True,
        "simplify": simplify,
    }
    if data_dir is not None:
        config_kwargs["data_dir"] = Path(data_dir)

    if output_dir is not None:
        config_kwargs["output_dir"] = Path(output_dir)
    if auto_download is not None:
        config_kwargs["auto_download"] = auto_download

    config = DelineatorConfig(**config_kwargs)

    output_path = config.output_dir
    if not output_path.exists():
        click.echo(f"Creating output directory: {output_path}")
        output_path.mkdir(parents=True, exist_ok=True)

    # Run delineation
    _delineate_dataframe(outlets_df, config)
    click.echo(f"Done. Output written to {output_path}/")


@click.command()
def delineator_dir():
    """Show the location of the data directory used by the delineator Python package.

    Will display the path to the data directory, which is either the default
    and depends on the operating system, or the custom location specified by the user,
    which can be set by setting the DELINEATOR_DATA environment variable.
    """
    click.echo(_get_data_dir())


@click.command()
@click.option("--basin", type=int, required=True, help="The megabasin ID, an integer from 11 to 86")
@click.option("--data-dir", default=None, help="Directory for delineator's data files. "
                                               "Defaults to the system default data directory.")
def delineator_download(basin, data_dir):
    """
    Download the data files needed for delineation in a megabasin.
    Use the command `delineator_dir` to get the location of the data directory.
    """
    # Resolve output directory
    if data_dir is None:
        data_dir = _get_data_dir()

    output_path = Path(data_dir) if data_dir else Path.cwd() / "output"
    if not output_path.exists():
        click.echo(f"Creating output directory: {output_path}")
        output_path.mkdir(parents=True, exist_ok=True)

    downloader(int(basin), str(output_path))


def _delineate_dataframe(outlets_df, config: DelineatorConfig):
    """
    For the command line interface, delineates watersheds based on a
    dataframe of outlet points, derived from the user's CSV file or lat/lon point.

    Parameters
    ----------
    outlets_df : DataFrame
        DataFrame containing the outlet points. Must have columns "lat", "lon", and "id".

    config : DelineatorConfig

    Returns
    -------
    Noting is returned, but the function writes the output files to disk.

    """

    for _, row in outlets_df.iterrows():
        lat, lon = row["lat"], row["lon"]
        id = row["id"]
        watershed_gdf, rivers_gdf, outlets_gdf = delineate(lat, lon, config)
        if watershed_gdf is None:
            logging.warning(f"No watershed found for point {id} ({lat}, {lon}).")
            continue
        write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, config, id=id)
