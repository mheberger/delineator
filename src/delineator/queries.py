import sqlite3
from collections import defaultdict
import geopandas as gpd
import shapely
from shapely.geometry import Polygon

from delineator.data import _find_data_file
from delineator.settings import DelineatorConfig


def get_upstream_comids(conn: sqlite3.Connection, home_unit_catchment: int) -> list:
    """
    Query #1, gets a complete list of upstream basins comids

    Parameters
    ----------
    conn: a sqlite3 connection to the database
    home_unit_catchment: the comid of the home unit catchment, e.g. 23016159

    Returns
    -------
    the list of upstream basins comids. The first element is the home unit catchment.

    Example:

        [23016159, 23016160, 23016603, .. , 23020142]
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
    Converts the list of unit catchment comids to a dictionary of nested hiearchical catchments

    Parameters
    ---------
    conn: a sqlite3 connection to the database
    upstream_catchment_list: integer list of unit catchment comids that are upstream of the home unit catchment

    Returns
    -------
    a dictionary of nested hiearchical catchments, where the keys are the level names ('L0', 'L1', ... 'L4')
    and the values are lists of basin comids.

    Example:

        {
            'L2': [23016162, 23016170, 23016215, 23016306, 23017387, 23017390, 23018599],
            'L0': [23016160, 23016161, 23016603, 23017101, 23017136, 23017189]
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


def get_home_unit_catchment_geom_and_area(db_path: str, unit_catchment: int) -> tuple[shapely.geometry.Polygon, float]:
    """
    Retrieves the geometry of a specified unit catchment based on its `comid`.
    Gets the geometry from the 'l0_basins' table in the sqlite database `basins##.db`
    params:
        db_path: str, path to the sqlite database
        unit_catchment: int, the comid of the unit catchment to retrieve the geometry for
    returns:
        geom: a shapely Polygon object representing the geometry of the unit catchment
        area: the area of the unit catchment, in km²
    """
    table = "l0_basins"
    catchment_gdf = gpd.read_file(db_path, layer=table, where=f"comid = {unit_catchment}")
    multi_poly = catchment_gdf.geometry.iloc[0]
    area = catchment_gdf.area_km2.iloc[0]
    return multi_poly, area


def get_upstream_area(home_unit_catchment: int, config: DelineatorConfig) -> float:
    """
    Retrieves the upstream area of a specified home unit catchment.
    Connects to the _rivers_ sqlite database and queries the table 'rivers'
    and the 'uparea' column.

    Parameters:
        home_unit_catchment (int): The catchment ID for which upstream area needs
                                   to be determined. The first two digits of this
                                   value identify the megabasin.
        config (DelineatorConfig): Configuration object containing resource paths
                                   and relevant settings.

    Returns:
        float: The upstream area of the specified catchment, in km²
    """
    if len(str(home_unit_catchment)) < 7:
        return 0

    megabasin = str(home_unit_catchment)[0:2]
    rivers_db_file = f"vector/rivers{megabasin}.db"
    rivers_db_path = _find_data_file(rivers_db_file, config)
    sql = f"SELECT uparea FROM rivers WHERE comid = {home_unit_catchment}"
    conn = sqlite3.connect(rivers_db_path)
    cur = conn.cursor()
    cur.execute(sql)
    rows = cur.fetchall()
    return rows[0][0]


def get_upstream_geometries(db_path: str, basins_dict: dict) -> list:
    """
    Get the geometries of the upstream basins. Uses an efficient approach that
    is based on nested catchment hierarchies that have been pre-computed.

    Parameters
    ----------
    basins_dict : dict
        A dictionary of basins, where the keys are the level names ('L0', 'L1', ... 'L4')
        and the values are lists of basin comids. For example:

        basins_dict = {
            'L2': [23014809, 23014837, 23014870, 23014919, 23016153, 23016154, 23016162, 23016170, 23016215, 23016306, 23017387, 23017390, 23018599],
            'L1': [23014801, 23014805],
            'L0': [23014799, 23014800, 23016113]
        }

    Returns
    -------
    A list of Shapely MultiPolygons (which we will dissolve later to create a single Polygon).

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
               split_catchment_polygon: Polygon,
               config: DelineatorConfig) -> gpd.GeoDataFrame | None:
    """
    Retrieves river data for the specified catchments, optionally splitting the geometry of the most downstream
    river reach using a provided polygon.

    Parameters:
        upstream_catchment_list (list): A list of catchment identifiers for which river data is queried.
        split_catchment_polygon (Polygon): A polygon used to split the geometry of the final river reach. Can be None if no splitting is required.
        config (DelineatorConfig): Configuration object containing settings for data processing, including the
            number of stream orders and file paths.

    Returns:
        gpd.GeoDataFrame: A GeoDataFrame containing the queried river records, optionally with updated geometry for
        the downstream river reach. Returns None if no river data is found.

    """
    home_unit_catchment = upstream_catchment_list[0]

    # Coastal catchments have a short comid; they have no river polylines associated with them
    if len(str(home_unit_catchment)) < 6:
        return None

    megabasin = str(upstream_catchment_list[0])[0:2]
    rivers_db_file = f"vector/rivers{megabasin}.db"
    rivers_db_path = _find_data_file(rivers_db_file, config)

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

    return rivers_gdf
