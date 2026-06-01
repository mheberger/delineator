
## Global Datasets

To delineate a watershed outside Iceland, you will need datasets covering 
your region. The `delineator` package can download these datasets for you as they are needed,
if your computer is connected to the internet. Alternatively, you can download 
the datasets in advance so that they are staged and ready. The `delineator` package has a
utility function to download the datasets: `download_data()`. For example, to download 
the datasets for the Amazon Basin:

```bash
download_data(62)
```

By default, the datasets will be downloaded to your system's default data directory.
The usual locations are: 

  - Windows: `C:\Users\<username>\AppData\Local\delineator`.
  - Linux: `~/.local/share/delineator`.
  - Mac: `~/Library/Application Support/delineator`.

You can change this by setting the `data_dir` option. See below for more options.

## Batch Processing

To delineate multiple watersheds, you can use the `delineate` command from the command line, 
and pass the path to a CSV file containing the coordinates of your desired watershed outlet.

For example:
```bash
delineate --csv outlets.csv
```

The CSV file should contain at least three columns: `id`, `lat`, and `lon`. 
In this example, the fields `name` and `area` are not required and are ignored.

```
id,lat,lng,name,area
6401070,64.71072,-21.60337,"Nordhura River at Stekkur",507
6401080,64.69229,-21.41046,"Hvita River at Kljafoss",1574
6401090,63.93796,-21.00666,"Olfusa River at Selfoss",5662
```

## Configuration Options

The `delineate` function accepts a number of options to customize the output.

To see a list of all options for the command line interface, run:

```bash
delineate --help
```

For a list of all options for the Python API, see the `DelineatorConfig` dataclass
in `delineator.settings`. Any of the options can be passed to the `delineate` function. 
For example:

```python
from delineator import DelineatorConfig, delineate

config = DelineatorConfig(
    data_dir="/path/to/data",
    rivers=False,
    fill=False
)

watershed_gdf, rivers_gdf, outlets_gdf = delineate(63.938, -21.004, config)

```
