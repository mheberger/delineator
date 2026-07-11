"""
Geospatial database utilities for "point in polygon" queries on the vector
catchment tables in the SQLite databases.

Candidates are pre-filtered through the R*Tree spatial index (created by
SpatiaLite during data prep, but readable as a plain SQLite virtual table),
then the exact containment test runs in GeoPandas on just those features.
"""

import math
import numbers
import sqlite3
import warnings
import geopandas as gpd
from shapely.geometry import Point

from delineator.constants import KM_PER_DEGREE


def _search_box_degrees(point: Point, search_dist_km: float) -> tuple[float, float]:
    """
    Convert a search distance in km to degree half-widths (dx, dy) at the
    latitude of *point*.

    A symmetric radius in km is an *anisotropic* box in degrees: away from
    the equator, one km spans more degrees of longitude than of latitude, by
    a factor of 1/cos(latitude). Using a single degree value for both axes
    would either miss features east-west or over-include them north-south.

    The longitude half-width is clamped to 180 degrees so the box remains
    valid arbitrarily close to the poles.
    """
    dy = search_dist_km / KM_PER_DEGREE
    cos_lat = abs(math.cos(math.radians(point.y)))
    dx = min(dy / cos_lat, 180.0) if cos_lat > 1e-12 else 180.0
    return dx, dy

# A geographic CRS makes shapely/GeoPandas .distance() emit this warning; we use
# planar distance deliberately (see _nearest_id), so it is suppressed where we
# do that on purpose.
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


def _nearest_id(gdf, point: Point, index_col: str, search_dist_km: float):
    """
    Nearest feature within search_dist_km on a loaded GeoDataFrame, or None.

    Distances are measured in kilometers using a local equirectangular
    approximation: longitudes are scaled by cos(latitude) so that planar
    distance in scaled degree space, multiplied by KM_PER_DEGREE, gives km.
    A raw planar degree distance would be anisotropic (biased toward
    north-south neighbors away from the equator); at search-tolerance scales
    (a few km) the equirectangular error is negligible.
    """
    if gdf.empty:
        return None
    gdf = gdf.set_index(index_col)
    cos_lat = abs(math.cos(math.radians(point.y)))
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=_GEOGRAPHIC_CRS_WARNING)
        scaled = gdf.geometry.scale(xfact=cos_lat, yfact=1.0, origin=(0, 0))
        scaled_point = Point(point.x * cos_lat, point.y)
        dist_km = scaled.distance(scaled_point) * KM_PER_DEGREE
    within = dist_km[dist_km <= search_dist_km]
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
    search_dist_km: float | None = None,
) -> int | str | None:
    """
    GeoPandas-based containment (and nearest) lookup.

    If an R*Tree spatial index is present, the R*Tree is queried for candidate
    IDs and only those features are loaded via ``gpd.read_file()`` for the
    exact ``contains()`` test. If no spatial index is present, the whole table
    is read once (with a warning for large tables) and reused for both the
    containment and the nearest tests.

    If no polygon contains the point and ``search_dist_km`` is provided, the
    nearest feature within ``search_dist_km`` kilometers is returned.

    Parameters
    ----------
    conn : sqlite3.Connection
        Connection to a database file.
    point : shapely.geometry.Point
        The location to search for.
    table : str
        The name of the table to search.
    geom_col : str
        The name of the geometry column, usually "geometry" or "geom".
    index_col : str
        The name of the index column, usually "comid".
    scan_threshold : int
        The number of rows above which a warning is emitted.
    search_dist_km : float, optional
        Nearest-feature tolerance in kilometers.

    Returns
    -------
    index : int or str or None
        The ID of the matching feature, or None if no match.

    Notes
    -----
    Requires a file-backed database. In-memory databases (``conn`` opened with
    ':memory:') are not supported and will raise RuntimeError.
    """
    db_path = _get_db_path(conn)
    if not db_path:
        raise RuntimeError(
            "_find_with_geopandas requires a file-backed database. "
            "In-memory SQLite connections are not supported."
        )

    rtree_table = _get_rtree_table(conn, table, geom_col)
    has_dist = search_dist_km is not None and search_dist_km > 0

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
        # km converted to an anisotropic degree box at this latitude
        dx, dy = _search_box_degrees(point, search_dist_km)
        bbox_ids = [
            row[0] for row in _query_rtree_bbox(
                conn, rtree_table, table, index_col,
                point.x - dx, point.y - dy,
                point.x + dx, point.y + dy,
            )
        ]
        if not bbox_ids:
            return None
        ngdf = gpd.read_file(
            db_path, layer=table,
            where=f"{index_col} IN ({_sql_in_list(bbox_ids)})",
        )
        return _nearest_id(ngdf, point, index_col, search_dist_km)

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
    return _nearest_id(gdf, point, index_col, search_dist_km)


def _point_in_polygon_analysis(
    conn: sqlite3.Connection,
    point: Point,
    table: str = "l0_basins",
    geom_col: str = "geometry",
    id_col: str = "comid",
    scan_threshold: int = 50_000,
    search_dist_km: float = 10.0,
) -> int | str | None:
    """
    Return the ID of the unit catchment that contains *point*.

    If a spatial index is present, the R*Tree is queried directly for candidate
    IDs and only those features are loaded via ``gpd.read_file()`` for the
    exact ``contains()`` test; otherwise the whole table is loaded (and a
    ``UserWarning`` is emitted past ``scan_threshold`` rows).

    Parameters:
    ----------
    conn:
        Open ``sqlite3.Connection``.  The connection is not closed here.
    point:
        Outlet location as a Shapely ``Point(x, y)`` (lon, lat).
    table:
        Name of the unit-catchment table.
    geom_col:
        Column holding the geometry, in any format GDAL recognises.
    id_col:
        Column whose value is returned for the matching catchment.
    scan_threshold:
        Row count above which a missing spatial index emits a ``UserWarning``.
    search_dist_km:
        Optional nearest-feature tolerance, in **kilometers**.  If no polygon
        contains the point, the nearest feature within this distance is
        returned.  Distances use a local equirectangular approximation,
        accurate to well under 1% at tolerance scales.  Set to 0 or None to
        disable the nearest-feature fallback.

    Returns
    -------
    index : int or str
        The value of *id_col* for the matching catchment, or ``None`` if no
        catchment contains the point (and none is within ``search_dist``).

    Notes
    -----
    If more than one polygon contains the point (overlapping geometries) the
    first match is returned.
    """
    return _find_with_geopandas(
        conn, point, table, geom_col, id_col, scan_threshold, search_dist_km
    )
