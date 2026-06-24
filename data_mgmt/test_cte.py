"""
Test a complicated query using CTE to find upstream basins
efficiently, using my new nested hierarhical geodata format

Matthew Heberger, June 2026
"""

from collections import defaultdict
import sqlite3
from shapely import unary_union
from shapely.ops import unary_union
from shapely.geometry import mapping
import geopandas as gpd
import folium
import json
from py.fast_dissolve import close_holes

# Home unit catchment (contains the watershed outlet)
comid = 23014799

# Megabasin is always the first two digits of the COMID
megabasin = str(comid)[0:2]

# Connect to the SQLite database
DB_PATH = fr'C:\Users\mheberger\AppData\Local\delineator\vector\basins{megabasin}.db'
con = sqlite3.connect(DB_PATH)
cur = con.cursor()

# Query #1, gets a complete list of upstream basins comids
sql = """
WITH RECURSIVE upstream(comid) AS (
    SELECT comid FROM l0_basins WHERE comid = ?
    UNION ALL
    SELECT uc.comid
    FROM l0_basins uc
    INNER JOIN upstream u ON uc.nextdown = u.comid
)
SELECT comid FROM upstream;
"""

cur.execute(sql, (comid,))
rows = cur.fetchall()

# Get the upstream comids as a Python list
upstream_comids = [row[0] for row in rows]

# Conver the list of integer comids to a comma-delimited string
placeholders = ", ".join([str(i) for i in upstream_comids])

# Query #2, gets the nested upstream basins
sql = f"""
WITH upstream(l0_id) AS (
    SELECT l0_id
    FROM catchment_hierarchy
    WHERE l0_id IN ({placeholders})
),
l4_matches AS (
    SELECT h.l4_id AS comid
    FROM upstream t
    JOIN catchment_hierarchy h ON t.l0_id = h.l0_id
    JOIN expected_counts ec ON ec.parent_id = h.l4_id AND ec.level_name = 'L4'
    WHERE h.l4_id IS NOT NULL
    GROUP BY h.l4_id, ec.expected_count
    HAVING COUNT(t.l0_id) = ec.expected_count
),
l3_matches AS (
    SELECT h.l3_id AS comid
    FROM upstream t
    JOIN catchment_hierarchy h ON t.l0_id = h.l0_id
    JOIN expected_counts ec ON ec.parent_id = h.l3_id AND ec.level_name = 'L3'
    WHERE h.l3_id IS NOT NULL
      AND NOT EXISTS (SELECT 1 FROM l4_matches lm WHERE lm.comid = h.l4_id)
    GROUP BY h.l3_id, ec.expected_count
    HAVING COUNT(t.l0_id) = ec.expected_count
),
l2_matches AS (
    SELECT h.l2_id AS comid
    FROM upstream t
    JOIN catchment_hierarchy h ON t.l0_id = h.l0_id
    JOIN expected_counts ec ON ec.parent_id = h.l2_id AND ec.level_name = 'L2'
    WHERE h.l2_id IS NOT NULL
      AND NOT EXISTS (SELECT 1 FROM l4_matches lm WHERE lm.comid = h.l4_id)
      AND NOT EXISTS (SELECT 1 FROM l3_matches lm WHERE lm.comid = h.l3_id)
    GROUP BY h.l2_id, ec.expected_count
    HAVING COUNT(t.l0_id) = ec.expected_count
),
l1_matches AS (
    SELECT h.l1_id AS comid
    FROM upstream t
    JOIN catchment_hierarchy h ON t.l0_id = h.l0_id
    JOIN expected_counts ec ON ec.parent_id = h.l1_id AND ec.level_name = 'L1'
    WHERE h.l1_id IS NOT NULL
      AND NOT EXISTS (SELECT 1 FROM l4_matches lm WHERE lm.comid = h.l4_id)
      AND NOT EXISTS (SELECT 1 FROM l3_matches lm WHERE lm.comid = h.l3_id)
      AND NOT EXISTS (SELECT 1 FROM l2_matches lm WHERE lm.comid = h.l2_id)
    GROUP BY h.l1_id, ec.expected_count
    HAVING COUNT(t.l0_id) = ec.expected_count
),
residuals AS (
    SELECT t.l0_id AS comid
    FROM upstream t
    JOIN catchment_hierarchy h ON t.l0_id = h.l0_id
    WHERE NOT EXISTS (SELECT 1 FROM l4_matches lm WHERE lm.comid = h.l4_id)
      AND NOT EXISTS (SELECT 1 FROM l3_matches lm WHERE lm.comid = h.l3_id)
      AND NOT EXISTS (SELECT 1 FROM l2_matches lm WHERE lm.comid = h.l2_id)
      AND NOT EXISTS (SELECT 1 FROM l1_matches lm WHERE lm.comid = h.l1_id)
)
    SELECT 'L4' AS match_level, comid
    FROM l4_matches

    UNION ALL

    SELECT 'L3' AS match_level, comid
    FROM l3_matches

    UNION ALL

    SELECT 'L2' AS match_level, comid
    FROM l2_matches

    UNION ALL

    SELECT 'L1' AS match_level, comid
    FROM l1_matches

    UNION ALL

    SELECT 'L0' AS match_level, comid
    FROM residuals;
"""

cur.execute(sql)
rows = cur.fetchall()

# Convert the results to a dictionary (mostly for my convenience -- easier to inspect and debug!)
basin_collection = defaultdict(list)
for match_level, comid in rows:
    basin_collection[match_level].append(comid)

basin_collection = dict(basin_collection)
print(basin_collection)

# Query #3, gets the geometries of the nested basins

TABLE_FOR_LEVEL = {
    'L0': 'l0_basins',
    'L1': 'l1_basins',
    'L2': 'l2_basins',
    'L3': 'l3_basins',
    'L4': 'l4_basins',
}

def dissolve_from_matches(db_path: str, matches: dict[str, list[int]]):
    geometries = []

    for level_name, comids in matches.items():
        if not comids:
            continue
        table = TABLE_FOR_LEVEL[level_name]
        id_list = ", ".join(str(c) for c in comids)
        gdf = gpd.read_file(db_path, layer=table, where=f"comid IN ({id_list})")
        geometries.extend(gdf.geometry.tolist())

    return unary_union(geometries)

# The result is a shapely geometry in WKB format (well-known binary)
watershed = dissolve_from_matches(DB_PATH, basin_collection)

watershed = close_holes(watershed, 100000)

# We now have the watershed as a shapely geometry. Let's visualize it

# Folium and web maps expect GeoJSON format. Shapely's `mapping()` does this conversion easily.
geojson_data = mapping(watershed)

# Folium maps require a starting location [Latitude, Longitude].
# We can extract the centroid of our geometry to center the map.
start_location = [watershed.centroid.y, watershed.centroid.x]
m = folium.Map(location=start_location, zoom_start=8, tiles='CartoDB positron')

# Add the GeoJSON representation of our WKB to the Folium map
folium.GeoJson(
    geojson_data,
    name="Watershed"
).add_to(m)

# Add layer control (optional, but helpful)
folium.LayerControl().add_to(m)

# Save the map to an HTML file
m.save("test_webmap.html")

# Export the GeoJSON for inspection
with open("test_watershed.geojson", "w") as f:
    json.dump(geojson_data, f, indent=2)