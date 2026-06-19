"""
Update the comids of coastal unit catchments

Matthew Heberger, May 2026

The coastal catchments had comids like 0, 1, 2, ...
and this occasionally caused trouble because in the MERIT-Basins
network structure, 0 is used as the terminal node (ocean or 
inland sink). So when we try to find ALL the upstream nodes, the 
script adds every node in the network, which takes a long time
and is NOT the desired result. So, I renumbered all the comids
so that they begin with the 2-digit megabasin ID (consistent with
the other MERIT-Basins catchments) and appended the existing 
comid. For example, 0 in basin 23 becomes 230 000
"""

import sqlite3

BASINS = [28,29,31,32,33,34,35,36,41,42,43,44,45,46,47,48,49,51,52,
53,54,55,56,57,61,62,63,64,65,66,67,71,72,73,74,75,76,77,78,81,82,
83,84,85,86]

dir = r"C:\Users\mheberger\AppData\Local\delineator\vector"

for basin in BASINS:
    DB = fr"{dir}\basins{basin}.db"
    print(DB)
    con = sqlite3.connect(DB)
    cur = con.cursor()

    val = basin * 10000
    sql = f"UPDATE l0_basins SET comid = comid + {val} WHERE comid < 10000"
    cur.execute(sql)
    con.commit()
    con.close()
