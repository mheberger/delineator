import logging
import sqlite3
import struct
from collections import defaultdict
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon
from pyproj import Geod

from delineator.data import _find_data_file
from delineator.settings import DelineatorConfig

logger = logging.getLogger(__name__)


def get_upstream_comids(conn: sqlite3.Connection, home_unit_catchment: int) -> list:
    """
    Get the COMIDs of every unit catchment upstream of the home unit catchment.

    Walks the drainage graph (the ``nextdown`` column of ``l0_basins``) with a
    recursive SQL query.

    Parameters
    ----------
    conn : sqlite3.Connection
        SQLite connection to the basin database.
    home_unit_catchment : int
        COMID of the home unit catchment, e.g. 23016159.

    Returns
    -------
    list[int]
        Upstream basin COMIDs. The first element is the home unit catchment.
        Example: [23016159, 23016160, 23016603, ..., 23020142]
    """
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
    cur = conn.cursor()
    cur.execute(sql, (home_unit_catchment,))
    rows = cur.fetchall()

    # Get the upstream comids as a Python list
    upstream_comids = [row[0] for row in rows]

    return upstream_comids


def get_hiearchical_basins(conn: sqlite3.Connection, upstream_catchment_list: list[int]) -> dict:
    """
    Group upstream unit catchments into the largest pre-aggregated catchments.

    Replaces each complete set of unit catchments with the coarsest nested
    catchment (levels L1-L4) that covers exactly that set, so the caller can
    dissolve a handful of large polygons instead of thousands of small ones.
    Unit catchments that no aggregate covers are returned at level L0.

    Parameters
    ----------
    conn : sqlite3.Connection
        SQLite connection to the basin database.
    upstream_catchment_list : list[int]
        COMIDs of the unit catchments upstream of the home unit catchment,
        not including the home unit catchment itself.

    Returns
    -------
    dict[str, list[int]]
        Catchment IDs keyed by level name (``'L0'`` to ``'L4'``). Only levels
        with at least one catchment appear. For example::

            {
                'L2': [23016162, 23016170, 23016215, 23016306, 23017387],
                'L0': [23016160, 23016161, 23016603, 23017101, 23017136],
            }
    """
    placeholders = ", ".join([str(i) for i in upstream_catchment_list])

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

    cur = conn.cursor()
    cur.execute(sql)
    rows = cur.fetchall()

    # Convert the results to a dictionary (mostly for my convenience -- easier to inspect and debug)
    basin_collection = defaultdict(list)
    for match_level, comid in rows:
        basin_collection[match_level].append(comid)

    basin_collection = dict(basin_collection)

    return basin_collection


def get_home_unit_catchment_geom_and_area(db_path: str, unit_catchment: int) -> tuple[Polygon | MultiPolygon, float]:
    """
    Get the geometry and area of one unit catchment.

    Reads from the ``l0_basins`` table of the basin database ``basins{NN}.db``.

    Parameters
    ----------
    db_path : str or Path
        Path to the basin database.
    unit_catchment : int
        COMID of the unit catchment.

    Returns
    -------
    geom : shapely Polygon or MultiPolygon
        The geometry of the unit catchment.
    area : float
        Area of the unit catchment, in km².
    """
    table = "l0_basins"
    catchment_gdf = gpd.read_file(db_path, layer=table, where=f"comid = {unit_catchment}")
    multi_poly = catchment_gdf.geometry.iloc[0]
    area = catchment_gdf.area_km2.iloc[0]
    return multi_poly, area


def _spatialite_to_wkb(blob: bytes) -> bytes:
    """
    Convert an uncompressed XY SpatiaLite Polygon or MultiPolygon BLOB to WKB.

    A SpatiaLite BLOB is a 39-byte header (start byte 0x00, byte order, SRID,
    MBR, marker 0x7C), then the geometry type and body, then an end byte 0xFE.
    The body is laid out as in WKB, except that each part of a collection
    starts with a 0x69 marker where WKB has a byte-order byte. Raises
    ValueError for anything else (compressed, Z/M, or other geometry types).
    """
    if len(blob) < 44 or blob[0] != 0x00 or blob[38] != 0x7C or blob[-1] != 0xFE:
        raise ValueError("not a SpatiaLite geometry BLOB")
    byte_order = blob[1]
    fmt = "<i" if byte_order == 1 else ">i"
    body = blob[39:-1]
    geom_type = struct.unpack_from(fmt, body, 0)[0]
    if geom_type == 3:
        return bytes([byte_order]) + body
    if geom_type != 6:
        raise ValueError(f"unsupported SpatiaLite geometry type {geom_type}")

    n_parts = struct.unpack_from(fmt, body, 4)[0]
    out = bytearray([byte_order])
    out += body[:8]  # type, number of polygons
    pos = 8
    for _ in range(n_parts):
        if body[pos] != 0x69 or struct.unpack_from(fmt, body, pos + 1)[0] != 3:
            raise ValueError("malformed SpatiaLite MultiPolygon")
        start = pos + 5  # past the marker and the part's type
        n_rings = struct.unpack_from(fmt, body, start)[0]
        end = start + 4
        for _ in range(n_rings):
            n_points = struct.unpack_from(fmt, body, end)[0]
            end += 4 + 16 * n_points  # XY doubles
        out.append(byte_order)
        out += body[pos + 1:end]
        pos = end
    if pos != len(body):
        raise ValueError("malformed SpatiaLite MultiPolygon")
    return bytes(out)


def get_upstream_area(home_unit_catchment: int, config: DelineatorConfig) -> float | None:
    """
    Get the total drainage area at the outlet of a unit catchment.

    Reads the pre-computed ``uparea`` column of the ``rivers`` table in the
    river database ``rivers{NN}.db``. The area includes the unit catchment
    itself.

    Parameters
    ----------
    home_unit_catchment : int
        COMID of the unit catchment. The first two digits identify the
        megabasin, which determines the river database to read.
    config : DelineatorConfig
        Configuration used to locate (and, if allowed, download) the river
        database.

    Returns
    -------
    float or None
        The upstream area, in km². Returns 0 for a COMID shorter than seven
        digits, and None if the river database is unavailable or has no
        record for this catchment.
    """
    if len(str(home_unit_catchment)) < 7:
        return 0

    megabasin = str(home_unit_catchment)[0:2]
    rivers_db_file = f"rivers{megabasin}.db"
    rivers_db_path = _find_data_file(rivers_db_file, config)
    if rivers_db_path is None:
        # _find_data_file has already logged a warning
        return None

    conn = sqlite3.connect(rivers_db_path)
    try:
        row = conn.execute(
            "SELECT uparea FROM rivers WHERE comid = ?", (home_unit_catchment,)
        ).fetchone()
    finally:
        conn.close()

    if row is None:
        logger.warning(
            f"No record for catchment {home_unit_catchment} in {rivers_db_file}; "
            f"upstream area is unknown."
        )
        return None
    return row[0]


def get_upstream_geometries(db_path: str, basins_dict: dict) -> list:
    """
    Get the geometries of the upstream catchments, by aggregation level.

    Reads each level's catchments from the matching table (``l0_basins`` to
    ``l4_basins``), using the pre-computed nested catchments so that far fewer
    polygons need to be fetched and dissolved.

    Parameters
    ----------
    db_path : str or Path
        Path to the basin database.
    basins_dict : dict[str, list[int]]
        Catchment IDs keyed by level name (``'L0'`` to ``'L4'``), as returned
        by ``get_hiearchical_basins``. For example::

            {
                'L2': [23014809, 23014837, 23014870, 23014919],
                'L1': [23014801, 23014805],
                'L0': [23014799, 23014800, 23016113],
            }

    Returns
    -------
    list of shapely Polygon or MultiPolygon
        One geometry per catchment, not yet dissolved.
    """
    TABLE_FOR_LEVEL = {
        'L0': 'l0_basins',
        'L1': 'l1_basins',
        'L2': 'l2_basins',
        'L3': 'l3_basins',
        'L4': 'l4_basins',
    }

    geometries = []

    for level_name, comids in basins_dict.items():
        if not comids:
            continue
        table = TABLE_FOR_LEVEL[level_name]
        id_list = ", ".join(str(c) for c in comids)
        gdf = gpd.read_file(db_path, layer=table, where=f"comid IN ({id_list})")
        geometries.extend(gdf.geometry.tolist())

    return geometries


def get_rivers(upstream_catchment_list: list,
               split_catchment_polygon: Polygon | MultiPolygon | None,
               config: DelineatorConfig) -> gpd.GeoDataFrame | None:
    """
    Get the river reaches in the watershed.

    Keeps only the larger rivers: ``config.num_stream_orders`` stream orders,
    counting down from the order of the outlet reach (with the default of 4
    and an outlet reach of order 6, reaches of order 3 to 6). If a split
    catchment polygon is given, the outlet reach is clipped to it and its
    ``lengthkm`` recalculated, so the river ends at the outlet.

    Parameters
    ----------
    upstream_catchment_list : list[int]
        COMIDs of the unit catchments in the watershed. The first element
        must be the home unit catchment, which contains the outlet.
    split_catchment_polygon : shapely Polygon or MultiPolygon, or None
        The part of the home unit catchment upstream of the outlet, used to
        clip the outlet reach. Pass None to keep the whole reach (low-res mode).
    config : DelineatorConfig
        Configuration used to locate the river database and for
        ``num_stream_orders``.

    Returns
    -------
    GeoDataFrame or None
        The river reaches, with the attributes from the river database.
        None if the home unit catchment is coastal (no rivers), the river
        database is unavailable, or no reaches match.
    """
    home_unit_catchment = upstream_catchment_list[0]

    # Coastal catchments have the digit 9 in the 3rd place; they do not have river
    # geometries associated with them.
    if str(home_unit_catchment)[2] == '9':
        return None

    megabasin = str(upstream_catchment_list[0])[0:2]
    rivers_db_file = f"rivers{megabasin}.db"
    rivers_db_path = _find_data_file(rivers_db_file, config)
    if rivers_db_path is None:
        # _find_data_file has already logged a warning
        return None

    if len(upstream_catchment_list) == 1:
        sql = f"SELECT * FROM rivers WHERE comid = {upstream_catchment_list[0]}"
    else:
        placeholders = ", ".join([str(i) for i in upstream_catchment_list])
        sql = f"""
            SELECT *
            FROM rivers
            WHERE comid IN ({placeholders})  -- your list of comids
              AND sorder > (
                  SELECT sorder FROM rivers WHERE comid = {upstream_catchment_list[0]}
              ) - {config.num_stream_orders}
        """

    rivers_gdf = gpd.read_file(rivers_db_path, sql=sql)

    if len(rivers_gdf) == 0:
        return None

    # The result looks nicer if we split the final river reach using the split catchment polygon
    if split_catchment_polygon is not None:
        # Get the geometry of the most downstream river reach
        segment_geom = rivers_gdf[rivers_gdf['comid'] == upstream_catchment_list[0]].geometry.iloc[0]

        # Use shapely to split the river reach
        split_reach = segment_geom.intersection(split_catchment_polygon)

        rivers_gdf.loc[rivers_gdf['comid'] == upstream_catchment_list[0], 'geometry'] = split_reach

        # Update the length attribute for the split river reach
        # Initialize a Geod object using the WGS84 ellipsoid (CRS 4326)
        geod = Geod(ellps="WGS84")

        # Calculate the true geodetic length in meters
        length_meters = geod.geometry_length(split_reach)

        # Convert to kilometers
        length_km = round(length_meters / 1000, 2)

        rivers_gdf.loc[rivers_gdf["comid"] == upstream_catchment_list[0], "lengthkm"] = length_km

    return rivers_gdf
