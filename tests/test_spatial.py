import sqlite3
from delineator.spatial import _find_with_geopandas, _point_in_polygon_analysis
from shapely.geometry import Point
from pathlib import Path


MYDIR = Path(__file__).parent.parent


def test_find_with_geopandas():
    dbfile = MYDIR / r"src\delineator\data\basins27.db"
    conn = sqlite3.connect(dbfile)
    mypoint = Point(-21.004, 63.938)
    result = _find_with_geopandas(conn, mypoint, "l0_basins", "geometry", "comid", 1000)
    assert result == 27001411


def test_find_with_spatialite():
    dbfile = MYDIR / r"src\delineator\data\basins27.db"
    conn = sqlite3.connect(dbfile)
    mypoint = Point(-21.004, 63.938)
    result = _point_in_polygon_analysis(conn, mypoint, "l0_basins", "geometry", "comid",
                                        1000, use_spatialite=True )
    print(result)
    assert result == 27001411


def test_find_without_spatialite_no_spatialite():
    dbfile = MYDIR / r"src\delineator\data\basins27.db"
    conn = sqlite3.connect(dbfile)
    mypoint = Point(-21.004, 63.938)
    result = _point_in_polygon_analysis(conn, mypoint, "l0_basins", "geometry", "comid",
                                        1000, use_spatialite=False )
    assert result == 27001411


def test_find_bad_point():
    dbfile = MYDIR / r"src\delineator\data\basins27.db"
    conn = sqlite3.connect(dbfile)
    mypoint = Point(-100.004, 73.938)
    result = _point_in_polygon_analysis(
        conn, mypoint, "l0_basins", "geometry", "comid", 1000, use_spatialite=True
    )
    print(result)
    assert result is None
