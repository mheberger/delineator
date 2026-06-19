"""
This version works well for putting shapefile data into a sqlite database.
"""

import sqlite3
import shapely
import geopandas as gpd
from shapely.geometry import MultiPolygon


def to_multi(geom):
    if geom.geom_type == 'Polygon':
        return MultiPolygon([geom])
    return geom


fname = r'C:\Users\mheberger\Documents\delineator\data\shp\basins_level2\cleaned\merit_hydro_vect_level2.shp'
gdf = gpd.read_file(fname)


fname = r'C:\Users\mheberger\AppData\Local\delineator\vector\megabasins.db'
con = sqlite3.connect(fname)
con.enable_load_extension(True)
con.load_extension("mod_spatialite")
con.execute("SELECT InitSpatialMetaData(1)")
con.execute("""
    CREATE TABLE megabasins (
        basin INTEGER
    )
""")
con.execute("SELECT AddGeometryColumn('megabasins', 'geometry', 4326, 'MULTIPOLYGON', 'XY')")

rows = [
    (row['BASIN'], shapely.to_wkb(to_multi(row.geometry)))
    for _, row in gdf.iterrows()
]
con.executemany("INSERT INTO megabasins VALUES (?, GeomFromWKB(?, 4326))", rows)
con.execute("SELECT CreateSpatialIndex('megabasins', 'geometry')")
con.commit()
