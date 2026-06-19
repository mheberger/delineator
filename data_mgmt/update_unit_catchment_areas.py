"""
I used this script to update the area of the unit catchments for
 the delienator Python package.

Matthew Heberger

"""

import sqlite3

basins = [11,12,13,14,15,16,17,18,21,22,23,24,25,26,28,29,31,32,33,34,35,36,41,42,43,44,45,46,47,48,49,51,52,53,54,
          55,56,57,61,62,63,64,65,66,67,71,72,73,74,75,76,77,78,81,82,83,84,85,86]

for basin in basins[1:]:
    print(basin)
    db_path = rf"C:\Users\mheberger\AppData\Local\delineator\vector\basins{basin}.db"

    con = sqlite3.connect(db_path)

    table = "l0_basins"

    # Load the spatialite extension
    con.enable_load_extension(True)
    con.load_extension('mod_spatialite')

    # Test that the extension is loaded
    con.execute("SELECT ST_Point(0, 0)").fetchone()

    sql = """
    UPDATE l0_basins
    SET area_km2 = (
        ROUND( 
            ST_Area( 
                ST_Transform(geometry, 6933) -- World Cylindrical Equal Area (meters) 
            ) / 1000000.0 -- convert m² → km²
        , 2)  -- round to 2 decimals
    )
    WHERE area_km2 = 0;
    """

    cursor = con.execute(sql)
    print(cursor.rowcount, "rows updated")

    con.commit()
    con.close()
