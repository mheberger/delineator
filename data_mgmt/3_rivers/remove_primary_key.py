import sqlite3


basins = [13,14,15,16,17,18,21,22,23,24,25,26,28,29,31,32,33,34,35,36,41,42,
          43,44,45,46,47,48,49,51,52,53,54,55,56,57,61,62,63,64,65,66,67,71,72,
          73,74,75,76,77,78,81,82,83,84,85,86]

for basin in basins:
    print(basin)
    db = rf"C:\Users\mheberger\AppData\Local\delineator\vector\basins{basin}.db"

    conn = sqlite3.connect(db)
    cur = conn.cursor()

    cur.executescript("""
        PRAGMA foreign_keys = OFF;
        BEGIN;
    
        CREATE TABLE l0_basins_new (
            comid   INTEGER,
            member_count    INTEGER,
            area_km2   REAL,
            nextdown  INTEGER,
            geometry  MULTIPOLYGON
        );
    
        INSERT INTO l0_basins_new (comid, member_count, area_km2, nextdown, geometry)
            SELECT comid, member_count, area_km2, nextdown, geometry FROM l0_basins_new;
    
        DROP TABLE l0_basins;
        ALTER TABLE l0_basins_new RENAME TO l0_basins;
    
        COMMIT;
        PRAGMA foreign_keys = ON;
    """)
    conn.commit()
    conn.execute("VACUUM")
    conn.commit()
    conn.close()