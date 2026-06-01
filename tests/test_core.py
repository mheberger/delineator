"""
Simple test script to run the delineator.
Plan to build this out into a proper test suite later.
"""

import pandas as pd
import logging
from delineator.core import delineate, downloader
from delineator.settings import DelineatorConfig
from delineator.util import write_outputs

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def test_one():
    test_coords = {
        "Iceland": (64.71072, -21.60337),
    }

    config = DelineatorConfig(high_res=True, output_format="geojson")
    for name, (lat, lng) in test_coords.items():
        print(f"Delineating {name}")
        watershed_gdf, rivers_gdf, outlets_gdf = delineate(lat, lng, config)
        if watershed_gdf is None:
            print("No watershed found.")
        else:
            write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, config, id=name)


def test_several():
    test_coords = {
        "Amazon": (-0.63717, -51.01179),
        "Chatahoochee": (31.123, -85.055),
        "uMngeni": (-29.763, 30.934),
        "Ikopa": (-18.829, 47.333),
        "Iceland": (64.71072, -21.60337),
    }

    config = DelineatorConfig(high_res=True, output_format="geojson")
    for name, (lat, lng) in test_coords.items():
        print(f"Delineating {name}")
        watershed_gdf, rivers_gdf, outlets_gdf = delineate(lat, lng, config)
        if watershed_gdf is None:
            print("No watershed found.")
        else:
            write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, config, id=name)


def test_csv():
    # Outlets all in Iceland -- should work out of the box
    fname = "iceland_outlets.csv"

    # Contains outlets around the world -- may trigger downloads if data files not already present
    #fname = "sample_outlets.csv"
    config = DelineatorConfig(high_res=True)

    df = pd.read_csv(fname)

    for i in range(len(df)):
        lat = df.iloc[i]['lat']
        lng = df.iloc[i]['lon']
        id = df.iloc[i]['id']
        name = df.iloc[i]['name']
        print(f"Delineating {name} ({id})")

        watershed_gdf, rivers_gdf, outlets_gdf = delineate(lat, lng, config)
        if watershed_gdf is None:
            print("No watershed found.")
        else:
            write_outputs(watershed_gdf, rivers_gdf, outlets_gdf, config, id=id)


def test_downloader():
    downloader(11)


if __name__ == "__main__":
    test_one()
    test_several()
    test_csv()
    test_downloader()
