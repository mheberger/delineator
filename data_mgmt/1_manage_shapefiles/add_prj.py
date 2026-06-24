"""
Some of the MERIT-Basins shapefiles don't have a .prj file to define their coordinate system.
This script adds one if it is missing.
"""

from pathlib import Path

# ESRI-flavored WKT for EPSG:4326 — the form ArcGIS/QGIS/GDAL all read cleanly
WKT = ('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",'
       'SPHEROID["WGS_1984",6378137.0,298.257223563]],'
       'PRIMEM["Greenwich",0.0],'
       'UNIT["Degree",0.0174532925199433]]')

folder = Path(r"C:\Data\GIS\MERITBasins\rivers\bugfix")

for shp in folder.glob("*.shp"):          # use rglob("*.shp") to recurse subfolders
    prj = shp.with_suffix(".prj")
    if not prj.exists():
        prj.write_text(WKT)
        print(f"Wrote {prj.name}")