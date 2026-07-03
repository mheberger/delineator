"""
Demonstration methods for delineating single coordinates,
multiple coordinates, processing CSV data, and using the
utility functions to download data files and write outputs.

Matthew Heberger, June 2026
"""

import pandas as pd
from delineator.core import delineate, downloader
from delineator.settings import DelineatorConfig
from delineator.util import write_outputs


def tryme():
    config = DelineatorConfig(verbose=True, fill=False, rivers=False, outlets=False, snapping=False)
    #w, r, o = delineate(44.705, 4.206, config)  # France (23,
    #w, r, o = delineate(63.938, -21.004, config)  # Iceland (27, built in)
    #w, r, o = delineate(-0.854, 41.350, config)  # East Africa (11)
    #w, r, o = delineate(48.982, 7.011, config)
    w, r, o = delineate(65.85166667, -17.8975, config)

    write_outputs(w, r, o, config)


def try_one():
    # Delineate a single watershed
    lat, lng = 48.863, 2.314  # Seine River at the Pont Alexandre III, Paris

    config = DelineatorConfig(high_res=True,
                              output_format="geojson",
                              outlets=False, # Skip outlets; will create watershed and rivers only
                              simplify=True,
                              verbose=True
                              )

    watershed_gdf, rivers_gdf, outlets_gdf = delineate(lat, lng, config)

    if watershed_gdf is None:
        print("No watershed found.")
    else:
        write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, config)


def try_several():
    test_coords = {
        "Amazon": (-0.63717, -51.01179),
        "Chatahoochee": (31.123, -85.055),
        "uMngeni": (-29.763, 30.934),
        "Ikopa": (-18.829, 47.333),
        "Iceland": (64.71072, -21.60337),
    }

    for name, (lat, lng) in test_coords.items():
        print(f"Delineating {name}")
        watershed_gdf, rivers_gdf, outlets_gdf = delineate(lat, lng)
        if watershed_gdf is None:
            print("No watershed found.")
        else:
            write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, id=name)


def try_csv():
    # Outlets all in Iceland -- should work out of the box
    fname = "iceland_outlets.csv"

    # Contains outlets around the world -- may trigger downloads if data files not already present
    #fname = "sample_outlets.csv"

    config = DelineatorConfig(high_res=True)

    df = pd.read_csv(fname)

    for i in range(len(df)):
        lat = df.iloc[i]['lat']
        lng = df.iloc[i]['lng']
        id = df.iloc[i]['id']
        name = df.iloc[i]['name']
        print(f"Delineating {name} ({id})")

        watershed_gdf, rivers_gdf, outlets_gdf = delineate(lat, lng, config)
        if watershed_gdf is None:
            print("No watershed found.")
        else:
            write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, config, id=id)


def try_downloader():
    """
    Download the data files for megabasin 11 (East Africa)
    and save them to the default data directory.
    """
    downloader(11)

    # To download to a specific directory:
    downloader(11, r"D:\data\delineator")


if __name__ == "__main__":

    tryme()
    #try_one()
    #try_downloader()
    #try_several()
    #try_csv()
