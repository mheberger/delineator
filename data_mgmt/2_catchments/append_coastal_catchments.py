"""
When I created my sqlite catchments databases, I forgot to add the coastal catchments.
Go back in and append them!

The coastal catchments had comids like 0, 1, 2, ...
and this occasionally caused trouble because in the MERIT-Basins
network structure, 0 is used as the terminal node (ocean or
inland sink). So when we try to find ALL the upstream nodes, the
script adds every node in the network, which takes a long time
and is NOT the desired result. So, I renumbered all the comids
so that they begin with the 2-digit megabasin ID (consistent with
the other MERIT-Basins catchments) and appended the existing
comid. For example, 0 in basin 23 becomes 230 000
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


BASINS = [51,52,53,54,
          55,56,57,61,62,63,64,65,66,67,71,72,73,74,75,76,77,78,81,82,83,84,85,86]


for basin in BASINS:
    print(f"Processing basin {basin}...")
    fname = rf'C:\Data\GIS\MERITBasins\coastal_hillslopes\fixed\fixed_{basin}.shp'
    if not os.path.exists(fname):
        print(f"  {fname} does not exist. Skipping...")
        continue

    gdf = gpd.read_file(fname)

    # Update the COMID to a more useful value (see note above)
    val = basin * 10000
    gdf['comid'] = gdf['comid'] + val

    conn = sqlite3.connect(rf"C:\Users\mheberger\AppData\Local\delineator\vector\basins{basin}.db")
    conn.enable_load_extension(True)
    conn.load_extension("mod_spatialite")
    cur = conn.cursor()

    rows = []

    for comid, row in gdf.iterrows():
        geom_bytes = _to_multipolygon_wkb(row.geometry)
        if geom_bytes is None:
            continue
        rows.append((
            int(comid),
            1,
            round(float(row.get("area_km2", 0.0)), 2),
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
