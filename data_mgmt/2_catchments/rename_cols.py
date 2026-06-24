"""
Temp. utility script to rename columns
"""
import sqlite3

basins = [11,12,13,14,15,16,17,18,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,41,42,43,44,45,46,47,48,49,51,52,53,54,
          55,56,57,61,62,63,64,65,66,67,71,72,73,74,75,76,77,78,81,82,83,84,85,86]

if False:
    tables = [
        'l0_basins',
        'l1_basins',
        'l2_basins',
        'l3_basins',
        'l4_basins',
    ]


    for basin in basins:
        con = sqlite3.connect(rf"C:\Users\mheberger\AppData\Local\delineator\vector\basins{basin}.db")
        cur = con.cursor()
        for table in tables:
            sql = f"ALTER TABLE {table} RENAME COLUMN poly_id TO comid;"
            cur.execute(sql)
        con.commit()
        con.close()

for basin in basins:
    con = sqlite3.connect(rf"C:\Users\mheberger\AppData\Local\delineator\vector\basins{basin}.db")
    cur = con.cursor()
    sql = "ALTER TABLE catchment_hierarchy RENAME unit_id to l0_id;"
    cur.execute(sql)
    con.commit()
    con.close()
