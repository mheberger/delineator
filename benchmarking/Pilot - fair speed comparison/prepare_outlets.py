"""
Choose the pilot outlets and pre-snap each one to the centre of a stream cell.

The timed runs give every tool the same pixel-centre coordinate with
snapping turned off, so all three delineate exactly the same watershed and
the comparison is of delineation speed alone. (Each tool snaps differently:
delineator uses a drainage-area threshold within the home unit catchment,
and pysheds and TauDEM use a stream mask over the whole raster.)

Snapping is done once here with delineator's default threshold (25 km²),
then moved to the exact centre of the snapped cell of the flow-direction
grid, which all three tools share.

Writes pilot_outlets.csv. Run once before run_pilot.py.
"""

import csv
import json
import math
from pathlib import Path

import rasterio
from rasterio.windows import Window

from delineator import DelineatorConfig, delineate

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "Experiment 3 - benchmark speed, vary basin area" / "expt3_outlets.csv"

# Ten Experiment 3 outlets: three megabasins, ~100 km² to the whole Amazon.
PILOT_IDS = [
    6401310,  # Iceland, ~100 km²
    6401120,  # Iceland, ~7,300 km²
    6604635,  # W. Europe, ~100 km²
    6118800,  # W. Europe, Somme, ~5,600 km²
    6335070,  # W. Europe, Rhine at Andernach, ~140,000 km²
    3649511,  # Amazon basin, ~230 km²
    3618730,  # Amazon basin, ~14,000 km²
    3618052,  # Amazon basin, Rio Negro, ~190,000 km²
    3627035,  # Amazon basin, Madeira at Humaita, ~1.1 million km²
    9,        # Amazon at its mouth, ~5.9 million km²
]


def main():
    cfg = json.loads((HERE / "bench_config.json").read_text())
    data_dir = Path(cfg["delineator_data_dir"])
    rows = {int(r["id"]): r for r in csv.DictReader(SOURCE.open(newline="", encoding="utf-8"))}
    config = DelineatorConfig(rivers=False, outlets=False, fill=False, clean=False,
                              auto_download=False, data_dir=str(data_dir))
    out = []
    for oid in PILOT_IDS:
        r = rows[oid]
        mb = int(r["megabasin"])
        watershed, _, _ = delineate(float(r["lat"]), float(r["lng"]), config)
        if watershed is None:
            raise SystemExit(f"Could not delineate outlet {oid}")
        snapped_lat = float(watershed["outlet_lat"].iloc[0])
        snapped_lng = float(watershed["outlet_lng"].iloc[0])
        with rasterio.open(data_dir / f"flowdir{mb}.tif") as src:
            col_f, row_f = ~src.transform * (snapped_lng, snapped_lat)
            col, row = math.floor(col_f), math.floor(row_f)
            lng, lat = src.transform * (col + 0.5, row + 0.5)
        with rasterio.open(data_dir / f"accum{mb}.tif") as src:
            accum = int(src.read(1, window=Window(col, row, 1, 1))[0, 0])
        out.append({
            "id": oid, "megabasin": mb, "name": r["name"],
            "area_grdc_km2": r["area"] or "",
            "lat_requested": r["lat"], "lng_requested": r["lng"],
            "lat": f"{lat:.10f}", "lng": f"{lng:.10f}", "row": row, "col": col,
            "accum_pixels": accum,
            "area_delineator_km2": round(float(watershed["area_km2"].iloc[0]), 1),
        })
        print(f"{oid}: cell ({row}, {col}), {accum:,} upstream cells, "
              f"{out[-1]['area_delineator_km2']:,} km²")

    with (HERE / "pilot_outlets.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(out[0]))
        writer.writeheader()
        writer.writerows(out)
    print(f"Wrote {HERE / 'pilot_outlets.csv'}")


if __name__ == "__main__":
    main()
