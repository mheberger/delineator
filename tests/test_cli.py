"""
Debug script for the command line interface.
"""

from click.testing import CliRunner
from delineator.cli import main

runner = CliRunner()
result = runner.invoke(main, ["--csv", "sample_outlets.csv"])
print(result.output)
if result.exception:
    import traceback
    traceback.print_exception(type(result.exception), result.exception, result.exception.__traceback__)