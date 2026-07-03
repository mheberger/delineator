"""
Geospatial database utilities for working with SQLite3 and SpatiaLite.

Main function is to do "point in polygon" queries on the vector catchment tables
in sqlite databases.  If the mod_spatialite extension is available it is used for
a database-side, index-accelerated query.  Otherwise the lookup falls back to
GeoPandas, which uses the same R*Tree spatial index when one is present.
"""

import logging
import numbers
import sqlite3
import warnings
import geopandas as gpd
from shapely.geometry import Point

logger = logging.getLogger(__name__)

# A geographic CRS makes shapely/GeoPandas .distance() emit this warning; we use
# planar distance deliberately (consistent with the SpatiaLite path), so it is
# suppressed where we do that on purpose.
_GEOGRAPHIC_CRS_WARNING = "Geometry is in a geographic CRS"


# ── small helpers ─

def _to_py_scalar(value):
    """
    Normalise a possibly-numpy scalar (e.g. a GeoDataFrame index label) to a
    plain Python ``int``/``str`` so both code paths return identical types.
    """
    item = getattr(value, "item", None)
    return item() if callable(item) else value


def _sql_in_list(values) -> str:
    """
    Render *values* as a SQL ``IN`` list literal, quoting non-integers safely.

    pyogrio's ``where`` filter is raw SQL with no parameter binding, so the
    values must be formatted here.  Integers (including numpy integers) are
    emitted bare; everything else is treated as text and single-quoted with
    embedded quotes escaped, so string IDs neither break the query nor open an
    injection hole.
    """
    parts = []
    for v in values:
        if isinstance(v, numbers.Integral):
            parts.append(str(int(v)))
        else:
            parts.append("'" + str(v).replace("'", "''") + "'")
    return ",".join(parts)


# ── SpatiaLite extension loading ──

def _load_spatialite(conn: sqlite3.Connection) -> bool:
    """
    Attempt to load the SpatiaLite extension.  Returns True on success.

    Tries the two most common extension names; returns False (rather than
    raising) if neither loads, if extension loading is disabled, or if the
    sqlite3 module was compiled without extension support.
    """
    try:
        conn.enable_load_extension(True)
    except (AttributeError, sqlite3.Error):
        # AttributeError: sqlite3 built without load-extension support.
        # sqlite3.Error: loading extensions disabled at build/runtime.
        return False

    try:
        for name in ("mod_spatialite", "spatialite"):
            try:
                conn.load_extension(name)
                conn.execute("SELECT ST_Point(0, 0)").fetchone()
                return True
            except sqlite3.Error:
                continue
        return False
    finally:
        try:
            conn.enable_load_extension(False)
        except (AttributeError, sqlite3.Error):
            pass


# ── R*Tree helpers ────────────────

def _get_rtree_table(
    conn: sqlite3.Connection,
    table: str,
    geom_col: str,
) -> str | None:
    """
    Return the R*Tree virtual-table name for (table, geom_col), or None.

    SpatiaLite names it idx_{table}_{geom_col}.  Confirms existence in
    sqlite_master before returning.
    """
    candidate = f"idx_{table}_{geom_col}"
    row = conn.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name=?",
        (candidate,),
    ).fetchone()
    return candidate if row else None


def _query_rtree(
    conn: sqlite3.Connection,
    point: Point,
    rtree_table: str,
    main_table: str,
    index_col: str,
) -> list[tuple]:
    """
    Return ``(id,)`` rows whose bounding boxes contain *point*, using the
    R*Tree for the envelope pre-filter.

    SpatiaLite convention: rtree.pkid → main_table.rowid.
    """
    sql = f"""
        SELECT m.{index_col}
        FROM   {main_table} m
        JOIN   {rtree_table} r ON m.rowid = r.pkid
        WHERE  r.xmin <= ? AND r.xmax >= ?
          AND  r.ymin <= ? AND r.ymax >= ?
    """
    return conn.execute(sql, (point.x, point.x, point.y, point.y)).fetchall()


def _query_rtree_bbox(
    conn: sqlite3.Connection,
    rtree_table: str,
    main_table: str,
    index_col: str,
    minx: float,
    miny: float,
    maxx: float,
    maxy: float,
) -> list[tuple]:
    """
    Return ``(id,)`` rows whose bounding boxes intersect the box
    [minx, miny, maxx, maxy], using the R*Tree.

    Used for the nearest-feature search: any geometry within distance *d* of a
    point necessarily has an envelope intersecting the point's +/- d box, so
    this pre-filter is exact (no false negatives) for a subsequent
    distance <= d test.
    """
    sql = f"""
        SELECT m.{index_col}
        FROM   {main_table} m
        JOIN   {rtree_table} r ON m.rowid = r.pkid
        WHERE  r.xmin <= ? AND r.xmax >= ?
          AND  r.ymin <= ? AND r.ymax >= ?
    """
    return conn.execute(sql, (maxx, minx, maxy, miny)).fetchall()


# ── Core lookup ───────────────────

def _get_db_path(conn: sqlite3.Connection) -> str | None:
    """
    Return the file path of the main attached database, or None for an
    in-memory database.

    Uses the ``pragma_database_list`` table-valued function so the *main*
    schema is selected explicitly rather than by row position.  In-memory
    databases report an empty file path, which is normalised to None.
    """
    row = conn.execute(
        "SELECT file FROM pragma_database_list WHERE name = 'main'"
    ).fetchone()
    path = row[0] if row else None
    return path or None  # '' (in-memory) -> None


def _find_with_spatialite(
    conn: sqlite3.Connection,
    point: Point,
    table: str,
    geom_col: str,
    id_col: str,
    search_dist: float | None = None,
) -> int | str | None:
    """
    Use SpatiaLite for a database-side point-in-polygon query.

    If no polygon contains the point and search_dist is provided, return the
    nearest feature whose geometry is within search_dist of the point.

    Note
    ----
    SpatiaLite does *not* use the R*Tree automatically for ``ST_Within`` or
    ``ST_Distance``; a bare predicate scans the whole table.  When a spatial
    index exists we restrict candidates with the ``SpatialIndex`` meta-table
    first (a point frame for containment, a +/- search_dist box for the nearest
    search).  We only do this when an index is actually present: on an
    unindexed table the index sub-query returns no rows and would silently hide
    a genuine match.  search_dist is in the geometry CRS units (degrees for
    EPSG:4326 — not metres).
    """
    rtree = _get_rtree_table(conn, table, geom_col)

    # ── containment ──
    if rtree:
        containment_sql = f"""
            SELECT {id_col}
            FROM   {table}
            WHERE  ST_Within(ST_Point(?, ?), {geom_col}) = 1
              AND  ROWID IN (
                  SELECT ROWID FROM SpatialIndex
                  WHERE  f_table_name = ?
                    AND  f_geometry_column = ?
                    AND  search_frame = ST_Point(?, ?)
              )
            LIMIT  1
        """
        containment_params = (point.x, point.y, table, geom_col, point.x, point.y)
    else:
        containment_sql = f"""
            SELECT {id_col}
            FROM   {table}
            WHERE  ST_Within(ST_Point(?, ?), {geom_col}) = 1
            LIMIT  1
        """
        containment_params = (point.x, point.y)

    row = conn.execute(containment_sql, containment_params).fetchone()
    if row:
        return row[0]

    # ── nearest-feature fallback ──
    if search_dist is None or search_dist <= 0:
        return None

    minx, miny = point.x - search_dist, point.y - search_dist
    maxx, maxy = point.x + search_dist, point.y + search_dist

    # ST_Distance is computed once per candidate via the subquery alias.
    if rtree:
        nearest_sql = f"""
            SELECT {id_col} FROM (
                SELECT {id_col},
                       ST_Distance({geom_col}, ST_Point(?, ?)) AS _dist
                FROM   {table}
                WHERE  ROWID IN (
                    SELECT ROWID FROM SpatialIndex
                    WHERE  f_table_name = ?
                      AND  f_geometry_column = ?
                      AND  search_frame = BuildMbr(?, ?, ?, ?)
                )
            )
            WHERE _dist <= ?
            ORDER BY _dist
            LIMIT 1
        """
        nearest_params = (
            point.x, point.y,
            table, geom_col,
            minx, miny, maxx, maxy,
            search_dist,
        )
    else:
        nearest_sql = f"""
            SELECT {id_col} FROM (
                SELECT {id_col},
                       ST_Distance({geom_col}, ST_Point(?, ?)) AS _dist
                FROM   {table}
            )
            WHERE _dist <= ?
            ORDER BY _dist
            LIMIT 1
        """
        nearest_params = (point.x, point.y, search_dist)

    row = conn.execute(nearest_sql, nearest_params).fetchone()
    return row[0] if row else None


def _containing_id(gdf, point: Point, index_col: str):
    """Exact containment test on a loaded GeoDataFrame; returns a Python scalar."""
    if gdf.empty:
        return None

    gdf = gdf.set_index(index_col)

    matches = gdf[gdf.geometry.contains(point)]
    if matches.empty:
        return None
    # If several polygons contain the point (overlapping geometries) the first
    # is returned; for a clean tessellation this is unambiguous.
    return _to_py_scalar(matches.index[0])


def _nearest_id(gdf, point: Point, index_col: str, search_dist: float):
    """Nearest feature within search_dist on a loaded GeoDataFrame, or None."""
    if gdf.empty:
        return None
    gdf = gdf.set_index(index_col)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=_GEOGRAPHIC_CRS_WARNING)
        dist = gdf.geometry.distance(point)
    within = dist[dist <= search_dist]
    if within.empty:
        return None
    return _to_py_scalar(within.idxmin())


def _find_with_geopandas(
    conn: sqlite3.Connection,
    point: Point,
    table: str,
    geom_col: str,
    index_col: str,
    scan_threshold: int,
    search_dist: float | None = None,
) -> int | str | None:
    """
    GeoPandas-based fallback for the containment (and nearest) lookup.

    If a SpatiaLite R*Tree spatial index is present, the R*Tree is queried for
    candidate IDs and only those features are loaded via ``gpd.read_file()`` for
    the exact ``contains()`` test.  If no spatial index is present, the whole
    table is read once (with a warning for large tables) and reused for both the
    containment and the nearest tests.

    If no polygon contains the point and search_dist is provided, the nearest
    feature within search_dist is returned, mirroring the SpatiaLite path.
    search_dist is in the geometry CRS units (degrees for EPSG:4326).

    Parameters:
    -----------
    conn: sqlite3.Connection to a database file
    point: shapely.geometry.Point, the location to search for
    table: str, the name of the table to search
    geom_col: str, the name of the geometry column, usually "geometry" or "geom"
    index_col: str, the name of the index column, usually "comid"
    scan_threshold: int, the number of rows above which a warning is emitted
    search_dist: optional nearest-feature tolerance in CRS units

    Returns:
    --------
    index: int or string, the ID of the matching feature, or None if no match

    Notes:
    ------
    Requires a file-backed database.  In-memory databases (conn opened with
    ':memory:') are not supported and will raise RuntimeError.
    """
    db_path = _get_db_path(conn)
    if not db_path:
        raise RuntimeError(
            "_find_with_geopandas requires a file-backed database. "
            "In-memory SQLite connections are not supported."
        )

    rtree_table = _get_rtree_table(conn, table, geom_col)
    has_dist = search_dist is not None and search_dist > 0

    if rtree_table:
        # ── containment via point-in-envelope candidates ──
        candidates = _query_rtree(conn, point, rtree_table, table, index_col)
        if candidates:
            ids = [row[0] for row in candidates]
            gdf = gpd.read_file(
                db_path, layer=table,
                where=f"{index_col} IN ({_sql_in_list(ids)})",
            )
            hit = _containing_id(gdf, point, index_col)
            if hit is not None:
                return hit

        if not has_dist:
            return None

        # ── nearest via envelope-intersects-box candidates ──
        bbox_ids = [
            row[0] for row in _query_rtree_bbox(
                conn, rtree_table, table, index_col,
                point.x - search_dist, point.y - search_dist,
                point.x + search_dist, point.y + search_dist,
            )
        ]
        if not bbox_ids:
            return None
        ngdf = gpd.read_file(
            db_path, layer=table,
            where=f"{index_col} IN ({_sql_in_list(bbox_ids)})",
        )
        return _nearest_id(ngdf, point, index_col, search_dist)

    # ── no spatial index: one full read, reused for containment + nearest ──
    count = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
    if count > scan_threshold:
        warnings.warn(
            f"find_home_catchment: no spatial index found on {table!r} "
            f"({count:,} rows). GDAL will scan the entire table. "
            f"Add a SpatiaLite spatial index for faster lookups: "
            f"SELECT CreateSpatialIndex('{table}', '{geom_col}');",
            stacklevel=3,
        )
    gdf = gpd.read_file(db_path, layer=table)

    hit = _containing_id(gdf, point, index_col)
    if hit is not None:
        return hit
    if not has_dist:
        return None
    return _nearest_id(gdf, point, index_col, search_dist)


def _point_in_polygon_analysis(
    conn: sqlite3.Connection,
    point: Point,
    table: str = "l0_basins",
    geom_col: str = "geometry",
    id_col: str = "comid",
    scan_threshold: int = 50_000,
    use_spatialite: bool = True,
    search_dist: float = 0.1,
) -> int | str | None:
    """
    Return the ID of the unit catchment that contains *point*.

    Tries the fastest available method in this order:

    1. **SpatiaLite** — a database-side ``ST_Within`` query, requiring the
       mod_spatialite extension.  When the table has a spatial index, candidates
       are pre-filtered through the R*Tree (the ``SpatialIndex`` meta-table)
       before the exact test; without an index it scans the table.
    2. **GeoPandas fallback** — when SpatiaLite is unavailable or disabled.  If a
       spatial index is present, the R*Tree is queried directly for candidate
       IDs and only those features are loaded via ``gpd.read_file()`` for the
       exact ``contains()`` test; otherwise the whole table is loaded (and a
       ``UserWarning`` is emitted past ``scan_threshold`` rows).

    Both paths honour ``search_dist`` and return Python scalars (``int``/``str``).

    Parameters:
    ----------
    conn:
        Open ``sqlite3.Connection``.  The connection is not closed here.
    point:
        Outlet location as a Shapely ``Point(x, y)`` (lon, lat).
    table:
        Name of the unit-catchment table.
    geom_col:
        Column holding the geometry.  The SpatiaLite path needs a SpatiaLite
        geometry BLOB; the GeoPandas fallback reads whatever GDAL recognises.
    id_col:
        Column whose value is returned for the matching catchment.
    scan_threshold:
        Row count above which a missing spatial index emits a ``UserWarning``.
    use_spatialite:
        ``True`` (default) – use SpatiaLite if the mod_spatialite extension
        loads, and fall back to the GeoPandas path if it cannot be loaded.
        ``False`` – always use the GeoPandas path; never attempt SpatiaLite.
    search_dist:
        Optional nearest-feature tolerance.  If no polygon contains the point,
        the nearest feature within this distance is returned (supported on both
        paths).  The value is in the geometry CRS units — for EPSG:4326 that is
        **degrees, not metres** (≈ 0.01 deg ≈ 1.1 km near the equator).  Project
        to a metric CRS if you need metre-based tolerances.

    Returns:
    -------
    index: The value of *id_col* for the matching catchment, or ``None`` if no
    catchment contains the point (and none is within ``search_dist``).

    Notes:
    -----
    If more than one polygon contains the point (overlapping geometries) the
    first match is returned.
    """

    # ── SpatiaLite path ───────────
    if use_spatialite:
        try:
            spatialite_loaded = _load_spatialite(conn)
        except Exception as e:  # defensive: _load_spatialite shouldn't raise
            spatialite_loaded = False
            logger.warning("SpatiaLite not available; falling back to GeoPandas: %s", e)

        if spatialite_loaded:
            return _find_with_spatialite(conn, point, table, geom_col, id_col, search_dist)

    # ── GeoPandas / GDAL path ─────
    return _find_with_geopandas(
        conn, point, table, geom_col, id_col, scan_threshold, search_dist
    )
