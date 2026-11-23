import sqlite3
from delineator.spatial import _find_with_geopandas, _point_in_polygon_analysis
from shapely.geometry import Point
from pathlib import Path


MYDIR = Path(__file__).parent.parent


def test_find_with_geopandas():
    """
    The lower-level lookup finds the unit catchment that contains a point.

    Uses the bundled Iceland data (megabasin 27). basins27.db has an R*Tree
    index, so this covers the index pre-filter followed by the exact
    containment test. No search distance is passed, so only exact containment
    is tried. Note that shapely points are (lng, lat).
    """
    dbfile = MYDIR / r"src\delineator\data\basins27.db"
    conn = sqlite3.connect(dbfile)
    mypoint = Point(-21.004, 63.938)
    result = _find_with_geopandas(conn, mypoint, "l0_basins", "geometry", "comid", 1000)
    assert result == 27001411


def test_point_in_polygon_analysis():
    """
    The public wrapper returns the same unit catchment as the lower-level lookup.

    Same Iceland point as test_find_with_geopandas, but called the way
    core.delineate calls it, with the default 10 km nearest-feature fallback
    switched on. A point inside a catchment must still give an exact match.
    """
    dbfile = MYDIR / r"src\delineator\data\basins27.db"
    conn = sqlite3.connect(dbfile)
    mypoint = Point(-21.004, 63.938)
    result = _point_in_polygon_analysis(conn, mypoint, "l0_basins", "geometry", "comid", 1000)
    assert result == 27001411


def test_find_bad_point():
    """
    A point outside every catchment returns None instead of raising.

    The point is in the Canadian Arctic (lng -100, lat 74), far outside the
    Iceland data and beyond the 10 km nearest-feature fallback. core.delineate
    relies on getting None here so it can log a warning and return no result.
    """
    dbfile = MYDIR / r"src\delineator\data\basins27.db"
    conn = sqlite3.connect(dbfile)
    mypoint = Point(-100.004, 73.938)
    result = _point_in_polygon_analysis(
        conn, mypoint, "l0_basins", "geometry", "comid", 1000
    )
    assert result is None


# ── megabasins.db regression tests ───────────────
# Commit 112018c dropped Greenland from megabasins.db and VACUUMed it, which
# renumbered megabasin 85's rowid without rebuilding the R*Tree. Every outlet
# in megabasin 85 then failed with "not in any megabasin".

def _megabasin_rows():
    """
    Return (rowid, basin ID, shapely geometry) for every row in the bundled megabasins.db.

    The rowid is included because the R*Tree index refers to rows by rowid
    (its pkid column), not by basin ID. Geometries are decoded from SpatiaLite
    BLOBs the same way the package does at runtime, without loading SpatiaLite.
    """
    import shapely
    from delineator.constants import MEGABASINS_DB_FILE
    from delineator.queries import _spatialite_to_wkb
    with sqlite3.connect(MEGABASINS_DB_FILE) as conn:
        rows = conn.execute("SELECT rowid, basin, geometry FROM megabasins").fetchall()
    return [(rowid, basin, shapely.from_wkb(_spatialite_to_wkb(blob))) for rowid, basin, blob in rows]


def test_megabasins_match_valid_list():
    """
    The basin IDs in megabasins.db are exactly constants.VALID_MEGABASINS (60 of them).

    Catches a megabasin added to or dropped from the database or the list
    without the other being updated. This test passed even on the broken
    file from 112018c, because all 60 IDs were still present. The two tests
    below are the ones that catch that bug.
    """
    from delineator.constants import VALID_MEGABASINS
    assert sorted(basin for _, basin, _ in _megabasin_rows()) == sorted(VALID_MEGABASINS)


def test_megabasins_rtree_matches_table():
    """
    Every megabasin row has one R*Tree entry, and that entry's box covers its polygon.

    Checks three things:
    - every rowid in the table has an index entry (pkid)
    - each entry's box contains the polygon's real bounds
    - no index entry is left over without a matching row

    This catches an index that was not rebuilt after rows were deleted or
    renumbered, as happened in 112018c. The 1e-4 degree tolerance is there
    because R*Tree stores boxes as float32, rounded outward, so they are
    slightly larger than the exact float64 bounds.
    """
    from delineator.constants import MEGABASINS_DB_FILE
    with sqlite3.connect(MEGABASINS_DB_FILE) as conn:
        rtree = {pkid: (xmin, xmax, ymin, ymax) for pkid, xmin, xmax, ymin, ymax in
                 conn.execute("SELECT pkid, xmin, xmax, ymin, ymax FROM idx_megabasins_geometry")}
    for rowid, basin, geom in _megabasin_rows():
        assert rowid in rtree, f"megabasin {basin} (rowid {rowid}) has no R*Tree entry"
        minx, miny, maxx, maxy = geom.bounds
        xmin, xmax, ymin, ymax = rtree.pop(rowid)
        # R*Tree boxes are float32, rounded outward
        assert xmin <= minx + 1e-4 and xmax >= maxx - 1e-4, f"megabasin {basin}: stale x bounds"
        assert ymin <= miny + 1e-4 and ymax >= maxy - 1e-4, f"megabasin {basin}: stale y bounds"
    assert not rtree, f"R*Tree entries with no matching row: {sorted(rtree)}"


def test_every_megabasin_is_found():
    """
    A point inside each megabasin resolves to that megabasin through the runtime lookup.

    This is the end-to-end version of the check above. It uses the same
    lookup call as step 2 of core.delineate. The test point is shapely's
    representative_point(), which is always inside the polygon, unlike the
    centroid. All failures are collected and reported together rather than
    stopping at the first one. With the broken 112018c file, this reports
    [(85, None)].
    """
    from delineator.constants import MEGABASINS_DB_FILE
    failures = []
    with sqlite3.connect(MEGABASINS_DB_FILE) as conn:
        for _, basin, geom in _megabasin_rows():
            found = _point_in_polygon_analysis(conn, geom.representative_point(), table="megabasins",
                                               geom_col="geometry", id_col="basin", search_dist_km=5)
            if found != basin:
                failures.append((basin, found))
    assert not failures, f"(expected, found): {failures}"
