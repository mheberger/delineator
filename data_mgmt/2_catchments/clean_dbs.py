"""
Compress database files

Matthew Heberger, May 2026

I found that after working with the database files, they had grown
in size. A simple `VACUUM` query cleaned things up.
"""

import sqlite3

basins = [11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33,
          34, 35, 36, 41, 42, 43, 44, 45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 56, 57, 61, 62, 63,
          64, 65, 66, 67, 71, 72, 73, 74, 75, 76, 77, 78, 81, 82, 83, 84, 85, 86]


for basin in basins[6:]:
    print(basin)

    #db_path = rf"C:\Users\mheberger\AppData\Local\delineator\vector\rivers{basin}.db"
    db_path = rf"C:\Users\mheberger\AppData\Local\delineator\vector\basins{basin}.db"

    con = sqlite3.connect(db_path)
    con.execute("VACUUM")
    con.commit()
    con.close()
