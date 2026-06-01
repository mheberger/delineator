"""
Point-in-polygon analysis
Geospatial database utilities for working with SQLite3 and SpatiaLite.

Provides helper functions to load SpatiaLite extensions, query R*Tree spatial
indices, and perform geospatial containment lookups using either SpatiaLite
or GeoPandas as fallback. Includes safeguards for handling both memory-based
and file-backed databases, along with performance warning mechanisms for
large tables without spatial indices.
"""


import sqlite3
import warnings
from typing import Optional
import geopandas as gpd
from shapely.geometry import Point


# ── SpatiaLite extension loading ──────────────────────────────────────────────

def _load_spatialite(conn: sqlite3.Connection) -> bool:
    """
    Attempt to load the SpatiaLite extension.  Returns True on success.
    Tries the two most common extension names; silently returns False if
    neither loads or if the sqlite3 module was compiled without extension
    support.
    """
    try:
        conn.enable_load_extension(True)
        for name in ("mod_spatialite", "spatialite"):
            try:
                conn.load_extension(name)
                conn.execute("SELECT ST_Point(0, 0)").fetchone()
                return True
            except sqlite3.OperationalError:
                continue
    except AttributeError:
        pass  # sqlite3 built without extension support
    finally:
        try:
            conn.enable_load_extension(False)
        except Exception:
            pass
    return False


# ── R*Tree helpers ────────────────────────────────────────────────────────────

def get_rtree_table(
    conn: sqlite3.Connection,
    table: str,
    geom_col: str,
) -> Optional[str]:
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


def query_rtree(
    conn: sqlite3.Connection,
    point: Point,
    rtree_table: str,
    main_table: str,
    geom_col: str,
    id_col: str,
) -> list[tuple]:
    """
    Return (id, geom_blob) rows whose bounding boxes contain *point*,
    using the R*Tree for the envelope pre-filter.

    SpatiaLite convention: rtree.pkid → main_table.rowid.
    """
    sql = f"""
        SELECT m.{id_col}, m.{geom_col}
        FROM   {main_table} m
        JOIN   {rtree_table} r ON m.rowid = r.pkid
        WHERE  r.xmin <= ? AND r.xmax >= ?
          AND  r.ymin <= ? AND r.ymax >= ?
    """
    return conn.execute(sql, (point.x, point.x, point.y, point.y)).fetchall()


# ── Core lookup ───────────────────────────────────────────────────────────────

def _get_db_path(conn: sqlite3.Connection) -> Optional[str]:
    """
    Return the file path of the main attached database, or None for
    in-memory databases.

    SQLite exposes this via PRAGMA database_list, which returns rows of
    (seq, name, file).  The main database is always seq=0.
    """
    row  = conn.execute("PRAGMA database_list").fetchone()
    path = row[2] if row else None
    return path if path else None  # empty string for :memory:


def _find_with_spatialite(
    conn: sqlite3.Connection,
    point: Point,
    table: str,
    geom_col: str,
    id_col: str,
    search_dist: float | None = None,
) -> Optional[int | str]:
    """
    Use SpatiaLite for a database-side point-in-polygon query.

    If no polygon contains the point and search_dist is provided, return the
    nearest feature whose geometry is within search_dist of the point.
    """
    containment_sql = f"""
        SELECT {id_col}
        FROM   {table}
        WHERE  ST_Within(ST_Point(?, ?), {geom_col}) = 1
        LIMIT  1
    """
    row = conn.execute(containment_sql, (point.x, point.y)).fetchone()
    if row:
        return row[0]

    if search_dist is None or search_dist <= 0:
        return None

    nearest_sql = f"""
        SELECT {id_col}
        FROM   {table}
        WHERE  ST_Distance({geom_col}, ST_Point(?, ?)) <= ?
        ORDER BY ST_Distance({geom_col}, ST_Point(?, ?))
        LIMIT  1
    """
    row = conn.execute(
        nearest_sql,
        (
            point.x,
            point.y,
            search_dist,
            point.x,
            point.y,
        ),
    ).fetchone()
    return row[0] if row else None


def _find_with_geopandas(
    conn: sqlite3.Connection,
    point: Point,
    table: str,
    geom_col: str,
    scan_threshold: int,
) -> Optional[int | str]:
    """
    GeoPandas-based fallback for containment lookup.

    Passes a bounding-box spatial filter to ``gpd.read_file()``, which
    translates it to GDAL's ``SetSpatialFilterRect()``.  The SQLite/SpatiaLite
    GDAL driver uses the R*Tree spatial index to resolve that filter, so no
    manual index querying is needed here.  All geometry parsing is handled by
    GDAL, which reads SpatiaLite blobs natively.

    The existing ``conn`` is used only to check whether a spatial index is
    present (for the large-table warning).  GeoPandas opens the database file
    independently via GDAL for the actual read.

    Notes
    -----
    Requires a file-backed database.  In-memory databases (conn opened with
    ':memory:') are not supported and will raise RuntimeError.
    """
    db_path = _get_db_path(conn)
    if not db_path:
        raise RuntimeError(
            "_find_with_geopandas requires a file-backed database. "
            "In-memory SQLite connections are not supported."
        )

    # Warn if there is no spatial index and the table is large, since GDAL
    # will fall back to a full scan.
    if not get_rtree_table(conn, table, geom_col):
        count = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        if count > scan_threshold:
            warnings.warn(
                f"find_home_catchment: no spatial index found on {table!r} "
                f"({count:,} rows). GDAL will scan the entire table. "
                f"Add a SpatiaLite spatial index for faster lookups: "
                f"SELECT CreateSpatialIndex('{table}', '{geom_col}');",
                stacklevel=3,
            )

    # GDAL resolves the bbox against the R*Tree (if present), returning only
    # candidates whose envelopes intersect the point.  All columns
    # are present in the resulting GeoDataFrame.
    px, py = point.x, point.y
    gdf = gpd.read_file(db_path, layer=table, bbox=(px, py, px, py), fid_as_index=True)

    matches = gdf[gdf.geometry.contains(point)]
    if matches.empty:
        return None
    return matches.index.values[0]


def point_in_polygon_analysis(
    conn: sqlite3.Connection,
    point: Point,
    table: str = "l0_basins",
    geom_col: str = "geometry",
    id_col: str = "comid",
    scan_threshold: int = 50_000,
    use_spatialite: Optional[bool] = None,
    search_dist: float | None = None,
) -> Optional[int | str]:
    """
    Return the ID of the unit catchment that contains *point*.

    Tries the fastest available method in this order:

    1. **SpatiaLite** — a single ``ST_Within`` query; requires the
       mod_spatialite extension.
    2. **GeoPandas + GDAL bbox filter** — passes the point coordinates as a
       bounding-box spatial filter to ``gpd.read_file()``.  GDAL's SQLite
       driver uses the R*Tree spatial index automatically when one is present,
       then GeoPandas tests the small candidate set with ``contains()``.

    Parameters
    ----------
    conn:
        Open ``sqlite3.Connection``.  The connection is not closed here.
    point:
        Outlet location as a Shapely ``Point(x, y)`` (lon, lat).
    table:
        Name of the unit-catchment table.
    geom_col:
        Column holding the geometry (SpatiaLite blob or plain WKB).
    id_col:
        Column whose value is returned for the matching catchment.
    scan_threshold:
        Row count above which a missing spatial index emits a
        ``UserWarning``.
    use_spatialite:
        ``True``  – require SpatiaLite (raise if unavailable).
        ``False`` – skip SpatiaLite even if installed.
        ``None``  – auto-detect (default).
    search_dist:
        Optional distance tolerance for the SpatiaLite path. If no polygon
        contains the point, the nearest feature within this distance is
        returned. Units are the units of the geometry CRS.

    Returns
    -------
    The value of *id_col* for the matching catchment, or ``None`` if no
    catchment contains the point.

    Raises
    ------
    RuntimeError
        If ``use_spatialite=True`` and SpatiaLite cannot be loaded, or if
        the database is in-memory and SpatiaLite is unavailable.
    """
    # ── SpatiaLite path ───────────────────────────────────────────────────────
    if use_spatialite is not False:
        spatialite_ok = _load_spatialite(conn)
        if use_spatialite is True and not spatialite_ok:
            raise RuntimeError(
                "use_spatialite=True but mod_spatialite could not be loaded. "
                "Install SpatiaLite or set use_spatialite=False."
            )
        if spatialite_ok:
            return _find_with_spatialite(
                conn,
                point,
                table,
                geom_col,
                id_col,
                search_dist,
            )

    # ── GeoPandas / GDAL path ─────────────────────────────────────────────────
    return _find_with_geopandas(
        conn, point, table, geom_col, id_col, scan_threshold
    )