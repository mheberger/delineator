"""
Updates the spatial index for the unit catchments or `l0_basins` 
in the vector geodata files for the `delineator` Python package.

Matthew Heberger
May 2026

"""

import sqlite3

basins = [11, 12,13,14,15,16,17,18,21,22,23,24,25,26,28,29,31,32,33,
34,35,36,41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,61,62,63,
64,65,66,67,71,72,73,74,75,76,77,78,81,82,83,84,85,86]



for basin in basins:
    print(basin)
    fname = rf"C:\Users\mheberger\AppData\Local\delineator\vector\basins{basin}.db"

    con = sqlite3.connect(fname)

    con.enable_load_extension(True)
    con.load_extension("mod_spatialite")

    con.execute("SELECT DisableSpatialIndex('l0_basins', 'geometry')")
    con.execute("DROP TABLE IF EXISTS idx_l0_basins_geometry")
    con.execute("SELECT CreateSpatialIndex('l0_basins', 'geometry')")
    con.commit()
