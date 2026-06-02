

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

The `delineate` function accepts a number of options to customize the output. A
limited set of options are available for the command line interface, and a more
complete set of options are available for the Python API.

To see a list of all options for the command line interface, run:

```bash
delineate --help
```

When using the Python API, options are set via the `DelineatorConfig` object, 
which is passed to the `delineate` function, as shown below.

```python
from delineator import DelineatorConfig, delineate

config = DelineatorConfig(
    data_dir="/path/to/data",
    rivers=False,  # the delineate function will not return rivers
    fill=False,    # do not fill in the gaps in the watersheds
    output_dir=r"C:\Users\user\Desktop\output"  # Set a custom output directory
    # Default output directory is the current working directory + "output"
)

watershed_gdf, rivers_gdf, outlets_gdf = delineate(63.938, -21.004, config)

# Update the config object to change options
config.rivers = True
watershed_gdf, rivers_gdf, outlets_gdf = delineate(63.938, -21.59, config)
```

## The Importance of Reviewing the Output

Automated watershed delineation often makes mistakes. The good news is that 
these errors can often be fixed by slightly 
moving the location of your watershed outlet. But not always! 
I recommend carefully reviewing every watershed you delineate to ensure  
that it appears correct.