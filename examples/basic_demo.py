"""
Demonstration methods for delineating single coordinates,
multiple coordinates, processing CSV data, and using the
utility functions to download data files and write outputs.

Matthew Heberger, June 2026
"""

import pandas as pd
import logging
from delineator.core import delineate, downloader
from delineator.settings import DelineatorConfig
from delineator.util import write_outputs


def tryme():
    config = DelineatorConfig(
        calc_area=False,
        data_dir=r"C:\Users\mheberger\Documents\watershed_app\static"
    )
    w, r, o = delineate(17.03333333, 18.67333334, config)
    write_outputs(w, r, o)


def try_one():
    "Delineate a single watershed"
    lat, lng = 48.863, 2.314  # Seine River at the Pont Alexandre III, Paris

    config = DelineatorConfig(high_res=True,
                              output_format="geojson",
                              outlets=False, # Skip outlets; will create watershed and rivers only
                              simplify=True
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
