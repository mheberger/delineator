"""
When I created my sqlite catchments databases, I forgot to add the coastal catchments.
Go back in and append them!
"""
import os

import shapely.wkb as shapely_wkb
from shapely.geometry import Polygon, MultiPolygon
import sqlite3
import geopandas as gpd
from typing import Optional


def _to_multipolygon_wkb(geom) -> Optional[bytes]:
    """Convert a Shapely geometry to WKB bytes, promoting Polygon → MultiPolygon."""
    if geom is None or geom.is_empty:
        return None
    if isinstance(geom, Polygon):
        geom = MultiPolygon([geom])
    return shapely_wkb.dumps(geom)


BASINS = [11]  #,12,13,14,15,16,17,18,21,22,23,24,25,26,28,29,31,32,33,34,35,36,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,61,62,63,64,65,66,67,71,72,73,74,75,76,77,78,81,82,83,84,85,86]


for basin in BASINS:
    print(f"Processing basin {basin}...")
    fname = rf'C:\Data\GIS\MERITBasins\coastal_hillslopes\fixed\fixed_{basin}.shp'
    if not os.path.exists(fname):
        print(f"  {fname} does not exist. Skipping...")
        continue

    gdf = gpd.read_file(fname)

    conn = sqlite3.connect(rf"C:\Users\mheberger\AppData\Local\delineator\2026-06-24\vector\basins{basin}.db")
    conn.enable_load_extension(True)
    conn.load_extension("mod_spatialite")
    cur = conn.cursor()

    rows = []

    for _, row in gdf.iterrows():
        geom_bytes = _to_multipolygon_wkb(row.geometry)
        if geom_bytes is None:
            continue
        rows.append((
            int(row["comid"]),
            1,
            round(float(row["unitarea"]), 2),
            0,
            geom_bytes,
        ))

    cur.executemany(
        """
        INSERT OR REPLACE INTO l0_basins
            (comid, member_count, area_km2, nextdown, geometry)
        VALUES (?, ?, ?, ?, GeomFromWKB(?, 4326))
        """,
        rows,
    )

    # Note the spatial index is updated automatically because the db has "triggers" enabled

    conn.commit()
    conn.close()
    print(f"  Wrote {len(rows)} rows → table l0_basins")
