"""
Convert MERIT-Basins rivers geodata from shapefiles to SpatiaLite databases.
Matthew Heberger, June 2026
"""
import os
import sqlite3
from pathlib import Path

import geopandas as gpd


SHAPEFILE_PATH = Path(r'C:\Data\GIS\MERITBasins\rivers\bugfix_repaired')

OUTPUT_DIR = Path(r'C:\Users\mheberger\AppData\Local\delineator\vector')

BASINS = [11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 41, 42, 43, 44,
          45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 56, 57, 61, 62, 63, 64, 65, 66, 67, 71, 72, 73, 74, 75, 76, 77, 78,
          81, 82, 83, 84, 85, 86, 91]

# Map source shapefile field (lower-cased) -> rivers table column.
# MERIT-Basins stores the Strahler stream order in a field named "order",
# which is a SQL reserved word, so we store it in the column "sorder".
ATTR_COLUMNS = {
    'comid':      'comid',
    'lengthkm':   'lengthkm',
    'uparea':     'uparea',
    'order':      'sorder',
    'nextdownid': 'nextdownid',
    'maxup':      'maxup',
    'up1':        'up1',
    'up2':        'up2',
    'up3':        'up3',
    'up4':        'up4',
}

INT_FIELDS = {'comid', 'order', 'nextdownid', 'maxup', 'up1', 'up2', 'up3', 'up4'}


def init_spatialite_db(path):
    """Create and initialize a fresh SpatiaLite database with the rivers table."""
    if os.path.exists(path):
        os.remove(path)  # start clean so re-runs don't fail on CREATE TABLE
    con = sqlite3.connect(path)
    con.enable_load_extension(True)
    con.load_extension("mod_spatialite")
    con.execute("SELECT InitSpatialMetaData(1)")
    con.execute("""
        CREATE TABLE rivers (
            basin INTEGER,
            comid INTEGER PRIMARY KEY,
            lengthkm FLOAT, 
            uparea FLOAT,
            sorder INTEGER,
            nextdownid INTEGER,
            maxup INTEGER,
            up1 INTEGER,
            up2 INTEGER,
            up3 INTEGER,
            up4 INTEGER
        )
    """)
    con.execute("SELECT AddGeometryColumn('rivers', 'geometry', 4326, 'LINESTRING', 'XY')")
    return con


def load_rivers(con, gdf, basin):
    """Insert every feature from a rivers GeoDataFrame into the rivers table."""
    # Grab the geometry as WKB *before* touching the columns, so renaming can't
    # disturb GeoPandas' active-geometry tracking. Vectorized; returns bytes.
    wkb = gdf.geometry.to_wkb().tolist()

    gdf.columns = [c.lower() for c in gdf.columns]

    missing = [c for c in ATTR_COLUMNS if c not in gdf.columns]
    if missing:
        raise KeyError(f"Basin {basin}: shapefile is missing expected field(s) {missing}. "
                       f"Found: {list(gdf.columns)}")

    geom_types = set(gdf.geom_type.unique())
    if geom_types - {'LineString'}:
        # The geometry column is declared LINESTRING, so a MultiLineString would
        # be rejected by GeomFromWKB. Flag it loudly rather than failing mid-insert.
        raise ValueError(f"Basin {basin}: expected only LineString geometries, got {geom_types}")

    n = len(gdf)

    # Build each column as a list of *native* Python objects. .tolist() converts
    # numpy int64/float64 to int/float; some sqlite3 builds reject numpy scalars.
    col = {'basin': [basin] * n}
    for src in ATTR_COLUMNS:
        dtype = 'int64' if src in INT_FIELDS else 'float64'
        col[src] = gdf[src].astype(dtype).tolist()

    records = list(zip(
        col['basin'],
        col['comid'],
        col['lengthkm'],
        col['uparea'],
        col['order'],       # -> sorder
        col['nextdownid'],
        col['maxup'],
        col['up1'],
        col['up2'],
        col['up3'],
        col['up4'],
        wkb,
    ))

    con.executemany("""
        INSERT INTO rivers
            (basin, comid, lengthkm, uparea, sorder,
             nextdownid, maxup, up1, up2, up3, up4, geometry)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, GeomFromWKB(?, 4326))
    """, records)
    con.commit()
    return n


def main():

    for basin in BASINS:
        print(f"Processing {basin}...")

        shapefile = SHAPEFILE_PATH / f'riv_pfaf_{basin}_MERIT_Hydro_v07_Basins_v01_bugfix1.shx'
        db_path = OUTPUT_DIR / f'rivers{basin}.db'

        gdf = gpd.read_file(shapefile)
        gdf = gdf[gdf.geometry.notna()].copy()

        con = init_spatialite_db(db_path)
        n = load_rivers(con, gdf, basin)

        # Build the spatial index *after* the bulk insert -- far faster than
        # maintaining the R*Tree row by row during the load.
        con.execute("SELECT CreateSpatialIndex('rivers', 'geometry')")
        con.commit()
        con.close()

        print(f"  inserted {n:,} reaches -> {db_path.name}")


if __name__ == "__main__":
    main()
