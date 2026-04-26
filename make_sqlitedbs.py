import sqlite3
import geopandas as gpd

con = sqlite3.connect(":memory:")
con.enable_load_extension(True)

try:
    con.load_extension("mod_spatialite")
    version = con.execute("SELECT spatialite_version()").fetchone()[0]
    print(f"SpatiaLite loaded successfully. Version: {version}")
except Exception as e:
    print(f"Failed to load SpatiaLite: {e}")


fname = r'C:\Users\mheberger\AppData\Local\delineator\vector\merit_hydro_vect_level2.shp'

gdf = gpd.read_file(fname)
fname = r'C:\Users\mheberger\AppData\Local\delineator\vector\megabasins.db'

gdf.to_file(fname, driver='SQLite', engine='fiona', spatialite=True, layer='megabasins')
