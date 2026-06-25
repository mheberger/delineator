import sqlite3

basins = [11,12,13,14,15,16,17,18,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,
          41,42,43,44,45,46,47,48,49,51,52,53,54,55,56,57,61,62,63,64,65,66,67,
          71,72,73,74,75,76,77,78,81,82,83,84,85,86]

for basin in basins:
    fname = rf"C:\Users\mheberger\AppData\Local\delineator\vector\basins{basin}.db"
    con = sqlite3.connect(fname)
    cursor = con.cursor()

    # Number of rows in l0_basins
    #cursor.execute("SELECT COUNT(*) FROM l0_basins")
    #print(f"Basin {basin}: {cursor.fetchone()[0]}")

    # Does l0_basins have a primary key?
    cursor.execute("SELECT name FROM pragma_table_info('l0_basins') WHERE pk >= 1;")
    print(f"Basin {basin}: {len(cursor.fetchall())} ")
