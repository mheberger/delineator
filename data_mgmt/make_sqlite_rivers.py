"""
Exports MERIT basin geometries from PostgreSQL/PostGIS to per-continent SpatiaLite databases.
Uses a server-side cursor to stream rows and avoid loading large result sets into memory.
"""
import psycopg2
import sqlite3
from shapely.geometry import MultiPolygon
from shapely import wkb
from credentials import *

CHUNK_SIZE = 500
OUTPUT_DIR = r'C:\Users\mheberger\AppData\Local\delineator\vector'

BASINS = [11,12,13,14,15,16,17,18,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,41,42,43,44,45,46,47,48,49,51,52,53,54,
          55,56,57,61,62,63,64,65,66,67,71,72,73,74,75,76,77,78,81,82,83,84,85,86,91]

def iter_basins(conn, basin: int, chunk_size=CHUNK_SIZE):
    """Server-side cursor that streams rows from PostGIS without loading all into memory."""
    sql = f"""
        SELECT basin, comid, lengthkm, uparea, sorder, nextdownid, maxup, up1, up2, up3, up4, ST_AsEWKB(geom_detail) AS geometry
        FROM merit_rivers
        WHERE basin = {basin}
        AND geom_detail IS NOT NULL
    """
    with conn.cursor(name='basin_cursor') as cur:
        cur.itersize = chunk_size
        cur.execute(sql)
        for basin, comid, lengthkm, uparea, sorder, nextdownid, maxup, up1, up2, up3, up4, geom_ewkb in cur:
            geom = wkb.loads(bytes(geom_ewkb))
            #geom = to_multi(geom)
            if geom is None:
                print(f"  WARNING: Could not coerce geometry for comid={comid}, skipping.")
                continue
            yield (basin, comid, float(lengthkm), float(uparea), sorder, nextdownid, maxup, up1, up2, up3, up4, wkb.dumps(geom))


def init_spatialite_db(path):
    """Create and initialize a fresh SpatiaLite database with the rivers table."""
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


def main():
    postgres_conn = psycopg2.connect(
        host="localhost",
        database=DBNAME,
        user=DBUSER,
        password=DBPW,
    )

    try:
        for basin in BASINS:
            print(f"Processing {basin}...")

            fpath = fr'{OUTPUT_DIR}\rivers_{basin}.db'
            sqlite_con = init_spatialite_db(fpath)

            try:
                rows = iter_basins(postgres_conn, basin)
                sqlite_con.executemany(
                    "INSERT INTO rivers VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, GeomFromWKB(?, 4326))",
                    rows,
                )
                sqlite_con.execute("SELECT CreateSpatialIndex('rivers', 'geometry')")
                sqlite_con.commit()
                print(f"  Done: {fpath}")
            except Exception as e:
                sqlite_con.rollback()
                print(f"  ERROR processing basin {basin}: {e}")
                raise
            finally:
                sqlite_con.close()

    finally:
        postgres_conn.close()


if __name__ == "__main__":
    main()