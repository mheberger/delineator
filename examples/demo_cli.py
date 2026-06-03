"""
Debugging scripts for the command line interface.

Named the scripts `try_###()` instead of `test_###()` so they don't invoke pytest.

Just for development. Add to .gitignore
 
"""

from click.testing import CliRunner
from delineator.cli import main, delineator_download
import shutil


def try_downloader():
    runner = CliRunner()
    result = runner.invoke(delineator_download, ["--basin", "54"])
    print(result.output)
    if result.exception:
        import traceback
        traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)


def try_csv():
    """
    Try running the command line function with the CSV option
    """
    runner = CliRunner()
    result = runner.invoke(main, ["--csv", "sample_outlets.csv", "--output-dir",  r"C:\Users\mheberger\Desktop"])
    print(result.output)
    if result.exception:
        import traceback
        traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)


def try_point():
    """
    Run the command line function with the point option
    """
    runner = CliRunner()
    result = runner.invoke(main, ["--point", "63.938", "-21.004", "--output-format", "gpkg",
                                  "--rivers", "--outlets", "--output-dir", r"C:\Users\mheberger\Desktop"])
    print(result.output)
    if result.exception:
        import traceback
        traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)


def clear_outputs():
    # Delete the folder output in the current directory
    # (This gets created by default when there is no
    shutil.rmtree('../tests/output', ignore_errors=True)


if __name__ == "__main__":
    #test_downloader()
    #test_delineator()
    clear_outputs()
    try_point()
