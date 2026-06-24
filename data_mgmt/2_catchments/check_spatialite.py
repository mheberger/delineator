# Verifying my spatialite installation

import sqlite3

# 1
DB_PATH = r'C:\Users\mheberger\Documents\delineator\output\basins23.db'
conn = sqlite3.connect(DB_PATH)
conn.enable_load_extension(True)
conn.load_extension("mod_spatialite")
row = conn.execute("SELECT spatialite_version()").fetchone()
print(row[0])

# 2
conn = sqlite3.connect(":memory:")
conn.enable_load_extension(True)
try:
    conn.load_extension("mod_spatialite")
    print("SpatiaLite loaded OK")
    print(conn.execute("SELECT spatialite_version()").fetchone()[0])
except Exception as e:
    print(f"Failed: {e}")

# 3
import subprocess
# Windows
subprocess.run(["where", "mod_spatialite.dll"], shell=True)
