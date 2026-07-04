"""
Round the fields `uparea` and `lengthkm` in the rivers databases.
The source has entries like 8.252159925700308. Better to round them to 8.25
"""
import sqlite3

SQL = "UPDATE rivers SET lengthkm = ROUND(lengthkm, 2), uparea = ROUND(uparea, 2);"

basins = [11, 12, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 41, 42, 43, 44,
          45, 46, 47, 48, 49, 51, 52, 53, 54, 55, 56, 57, 61, 62, 63, 64, 65, 66, 67, 71, 72, 73, 74, 75, 76, 77, 78,
          81, 82, 83, 84, 85, 86]

folder = r"C:\Users\mheberger\Documents\watershed_app\static\data"

for basin in basins:
    print(f"Processing basin {basin}")
    fname = f'{folder}/rivers{basin}.db'
    con = sqlite3.connect(fname)
    con.execute(SQL)
    con.commit()
    con.close()
