from click.testing import CliRunner
from delineator.cli import main

runner = CliRunner()

result = runner.invoke(
    main,
    [
        "--csv",
        "sample_outlets.csv"
    ],
)

print(result.exit_code)
print(result.output)
print(result.exception)