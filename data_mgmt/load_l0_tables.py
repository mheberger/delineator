"""
Going in and fixing a previous mistake where I accidentally cleared the l0_basins table
"""

import sqlite3
from typing import Optional, Tuple
from shapely.geometry import Polygon, MultiPolygon
from shapely import wkb as shapely_wkb
import geopandas as gpd


basins = [15,16,17,18,21,22,24,25,26,28,29,31,32,33,34,35,36,
          41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,61,62,63,64,65,66,67,
          71,72,73,74,75,76,77,78,81,82,83,84,85,86]

basins = [13]

def _find_col(df, needle: str) -> str:
    """Case/substring-tolerant column lookup (e.g. 'nextdown' -> 'NextDownID')."""
    for c in df.columns:
        if c.lower() == needle.lower():
            return c
    for c in df.columns:
        if needle.lower() in c.lower():
            return c
    raise KeyError(f"No column matching {needle!r} in {list(df.columns)}")


def get_geodata() -> gpd.GeoDataFrame:
    """
    Load catchments and rivers for the current BASIN from shapefiles.

    Topology (`nextdown`) is taken from the
    rivers layer's NextDownID and attached to the catchment polygons on COMID.

    Returns
    -------
    basins_gdf : GeoDataFrame
        Indexed by comid (int). Columns: 'unitarea', 'nextdown'; active
        geometry column named 'geom' (polygons), CRS as read from the .prj.
    rivers_gdf : GeoDataFrame | None
        Only if RIVERS is True. Indexed by comid; columns 'lengthkm';
        active geometry 'geom' (lines).
    """
    cat_path = CAT_SHP
    riv_path = RIV_SHP

    # --- catchments (polygons): comid, unitarea, geometry ---
    basins_gdf = gpd.read_file(cat_path).rename(columns={"COMID": "comid"})
    basins_gdf["comid"] = basins_gdf["comid"].astype(int)

    # --- rivers: carry the topology (and lengthkm for the rivers path) ---
    rivers_raw = gpd.read_file(riv_path).rename(columns={"COMID": "comid"})
    rivers_raw["comid"] = rivers_raw["comid"].astype(int)

    # New 2026-06-21 In the MERIT-Basins 'bugfix' version, endorheic basins have
    # the value of their own COMID in the NextDownID field, which creates a loop
    rivers_raw = rivers_raw[(rivers_raw["NextDownID"] != rivers_raw["comid"])].copy()
    # If rivers have zero length, drop (some non-existent rivers have lengthkm = 9e-5 (!)
    rivers_gdf = rivers_raw[rivers_raw["lengthkm"] > 1e-4].copy()
    print("null river geoms:", rivers_raw.geometry.isna().sum())
    print("empty river geoms:", rivers_gdf.geometry.is_empty.sum())

    nextdown_col = _find_col(rivers_gdf, "nextdown")          # e.g. 'NextDownID'
    topo = rivers_gdf[["comid", nextdown_col]].rename(
        columns={nextdown_col: "nextdown"}
    )

    # attach nextdown onto the catchments (the old r JOIN b ON r.comid = b.comid)
    basins_gdf = basins_gdf.merge(topo, on="comid", how="left")
    # catchments with no matching river row are terminal sinks or coastal catchments
    basins_gdf["nextdown"] = basins_gdf["nextdown"].fillna(0).astype(int)

    # match the layout the rest of the script expects: geom col 'geom', index 'comid'
    basins_gdf = basins_gdf.rename_geometry("geom").set_index("comid")

    return basins_gdf


def _to_multipolygon_wkb(geom) -> Optional[bytes]:
    """Convert a Shapely geometry to WKB bytes, promoting Polygon → MultiPolygon."""
    if geom is None or geom.is_empty:
        return None
    if isinstance(geom, Polygon): 
        geom = MultiPolygon([geom])
    return shapely_wkb.dumps(geom)


def write_basin_table(
    cur: sqlite3.Cursor,
    gdf: gpd.GeoDataFrame,
    table: str,
    member_count_col: Optional[str] = None,
) -> None:
    """
    Write a GeoDataFrame to a pre-created SpatiaLite basin table.

    Parameters
    ----------
    cur : sqlite3.Cursor

    gdf : GeoDataFrame
        Source data. Index must be the comid. Must contain
        'unitarea' and the geometry column; optionally 'nextdown' and
        the column named by member_count_col.
    table : str
        Target table name, e.g. 'l0_basins' or 'l2_basins'.
    member_count_col : str | None
        Column in gdf holding the merge count. If None, defaults to 1.
    """


    rows = []
    for comid, row in gdf.iterrows():
        geom_bytes = _to_multipolygon_wkb(row.geom)
        if geom_bytes is None:
            continue
        rows.append((
            int(comid),
            int(row[member_count_col]) if member_count_col else 1,
            round(float(row.get("unitarea", 0.0)), 2),
            int(row.get("nextdown", 0)),
            geom_bytes,
        ))

    cur.executemany(
        f"""
        INSERT OR REPLACE INTO {table}
            (comid, member_count, area_km2, nextdown, geometry)
        VALUES (?, ?, ?, ?, GeomFromWKB(?, 4326))
        """,
        rows,
    )
    print(f"  Wrote {len(rows)} rows → table '{table}'")


# MAIN

for basin in basins:

    DB_PATH = f"C:/Users/mheberger/AppData/Local/delineator/vector/basins{basin}.db"
    CAT_SHP = rf"C:\Data\GIS\MERITBasins\catchments\bugfix_repaired\cat_pfaf_{basin}_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"
    RIV_SHP = rf"C:\Data\GIS\MERITBasins\rivers\bugfix_repaired\riv_pfaf_{basin}_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"

    con = sqlite3.connect(DB_PATH)
    con.enable_load_extension(True)
    con.load_extension("mod_spatialite")

    cur = con.cursor()
    sql = "DROP TABLE IF EXISTS l0_basins;"
    cur.execute(sql)

    sql = "SELECT DisableSpatialIndex('l0_basins', 'geometry');"
    cur.execute(sql)
    sql = "DROP TABLE IF EXISTS idx_l0_basins_geometry;"
    cur.execute(sql)

    gdf = get_geodata()

    sql = """       
    CREATE TABLE IF NOT EXISTS l0_basins (
        comid        INTEGER, 
        member_count INTEGER NOT NULL DEFAULT 1,
        area_km2     REAL,
        nextdown     INTEGER NOT NULL DEFAULT 0
    );
    """
    cur.execute(sql)


    sql = "SELECT AddGeometryColumn('l0_basins', 'geometry', 4326, 'MULTIPOLYGON', 'XY')"
    cur.execute(sql)

    write_basin_table(cur, gdf, "l0_basins")

    sql = "SELECT CreateSpatialIndex('l0_basins', 'geometry')"
    cur.execute(sql)

    con.commit()
    con.close()




