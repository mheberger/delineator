"""
Aggregate MERIT-Basins unit catchments to larger sizes and write to SpatiaLite.

Geoprocessing (dissolve) is performed with PostGIS via a SQLAlchemy engine.
The previous GeoPandas dissolve of rivers produced MultiLineStrings that
were not merged at shared endpoints; ST_LineMerge(ST_Union(...)) in PostGIS
collapses connecting LineStrings into single LineStrings.

Output is a single SpatiaLite .db file containing one table per level:

  Spatial tables (each with a geometry column and spatial index):
    l0_basins       — original base unit polygons
    l1_basins       — aggregated at ~1 000 km² threshold
    l2_basins       — aggregated at ~10 000 km² threshold
    l3_basins       — aggregated at ~100 000 km² threshold
    l4_basins       — aggregated at ~1 000 000 km² threshold

  Non-spatial tables:
    catchment_hierarchy — maps each base unit to its L1–L4 parent comid
    expected_counts     — number of base units inside each pre-staged polygon
"""

import sqlite3
from pathlib import Path
from typing import Dict, Optional, Tuple

import geopandas as gpd
from geopandas import read_postgis
from matplotlib import pyplot as plt
from shapely.geometry import MultiPolygon, Polygon
from shapely import wkb as shapely_wkb
from sqlalchemy import create_engine, text
from sqlalchemy.engine import Engine

from consolidate import consolidate_network, show_area_stats
from delineator.util import _close_holes
from graph_tools import make_river_network
from credentials import *


# ---------------------------------------------------------------------------
# SpatiaLite helpers
# ---------------------------------------------------------------------------

def get_db_conn(db_path: str) -> sqlite3.Connection:
    """
    Establishes a connection to an SQLite database and loads the SpatiaLite extension.
    """
    conn = sqlite3.connect(db_path)
    conn.enable_load_extension(True)
    conn.load_extension("mod_spatialite")
    return conn


def _to_multipolygon_wkb(geom) -> Optional[bytes]:
    """Convert a Shapely geometry to WKB bytes, promoting Polygon → MultiPolygon."""
    if geom is None or geom.is_empty:
        return None
    if isinstance(geom, Polygon):
        geom = MultiPolygon([geom])
    return shapely_wkb.dumps(geom)


def _create_basin_table(cur: sqlite3.Cursor, table: str) -> None:
    """
    Create one spatial basin table and register it with SpatiaLite.

    Each level gets its own table: l0_basins, l1_basins, l2_basins, …
    All share the same column layout:
        poly_id      — comid from MERIT-Basins (PK)
        member_count — number of unit catchments dissolved into this polygon
        area_km2     — total area
        geometry     — MULTIPOLYGON, EPSG:4326
    """
    cur.execute(f"""
        CREATE TABLE IF NOT EXISTS {table} (
            poly_id      INTEGER PRIMARY KEY,
            member_count INTEGER NOT NULL DEFAULT 1,
            area_km2     REAL,
            nextdown     INTEGER NOT NULL DEFAULT 0
        )
    """)
    cur.execute(
        f"SELECT AddGeometryColumn('{table}', 'geometry', 4326, 'MULTIPOLYGON', 'XY')"
    )
    cur.execute(f"SELECT CreateSpatialIndex('{table}', 'geometry')")


def init_db(db_path: str, n_levels: int = 4) -> None:
    """
    Create all tables and spatial indexes in a fresh SpatiaLite database.
    Deletes any existing file at db_path first.

    Tables created:
        l0_basins     — base unit polygons
        l1_basins … lN_basins — one per aggregation level
        catchment_hierarchy — relational hierarchy
        expected_counts     — base-unit counts per pre-staged polygon
    """
    Path(db_path).unlink(missing_ok=True)
    conn = get_db_conn(db_path)
    cur = conn.cursor()

    # Required SpatiaLite metadata (1 = suppress PROJ warnings)
    cur.execute("SELECT InitSpatialMetaData(1)")

    # -- lN_basins (one per level) ----------------------------------------
    for lv in range(0, n_levels + 1):
        _create_basin_table(cur, f"l{lv}_basins")

    # -- catchment_hierarchy ---------------------------------------------
    # One row per base unit; lN_id is the comid of its parent at level N.
    # NULL means no fully-contained parent exists at that level.
    level_cols = "\n".join(
        f"    l{lv}_id  INTEGER," for lv in range(1, n_levels + 1)
    ).rstrip(",")
    cur.execute(f"""
        CREATE TABLE IF NOT EXISTS catchment_hierarchy (
            unit_id  INTEGER PRIMARY KEY,
{level_cols}
        )
    """)
    for lv in range(1, n_levels + 1):
        cur.execute(
            f"CREATE INDEX IF NOT EXISTS idx_hierarchy_l{lv} "
            f"ON catchment_hierarchy(l{lv}_id)"
        )

    # -- expected_counts -------------------------------------------------
    # Number of base unit catchments inside each pre-staged parent polygon.
    # Used at query time (HAVING clause) to confirm full containment.
    cur.execute("""
        CREATE TABLE IF NOT EXISTS expected_counts (
            parent_id      INTEGER NOT NULL,
            level_name     TEXT    NOT NULL,
            expected_count INTEGER NOT NULL CHECK (expected_count > 0),
            PRIMARY KEY (parent_id, level_name)
        )
    """)
    cur.execute("""
        CREATE INDEX IF NOT EXISTS idx_expected_parent
            ON expected_counts(parent_id, level_name)
    """)

    conn.commit()
    conn.close()
    print(f"Initialised SpatiaLite database: {db_path}")


def write_basin_table(
    db_path: str,
    gdf: gpd.GeoDataFrame,
    table: str,
    member_count_col: Optional[str] = None,
) -> None:
    """
    Write a GeoDataFrame to a pre-created SpatiaLite basin table.

    Parameters
    ----------
    db_path : str
        Path to the SpatiaLite .db file.
    gdf : GeoDataFrame
        Source data. Index must be the comid (poly_id). Must contain
        'unitarea' and the geometry column; optionally 'nextdown' and
        the column named by member_count_col.
    table : str
        Target table name, e.g. 'l0_basins' or 'l2_basins'.
    member_count_col : str | None
        Column in gdf holding the merge count. If None, defaults to 1.
    """
    conn = get_db_conn(db_path)
    cur = conn.cursor()


    rows = []
    for poly_id, row in gdf.iterrows():
        geom_bytes = _to_multipolygon_wkb(row.geom)
        if geom_bytes is None:
            continue
        rows.append((
            int(poly_id),
            int(row[member_count_col]) if member_count_col else 1,
            round(float(row.get("unitarea", 0.0)), 2),
            int(row.get("nextdown", 0)),
            geom_bytes,
        ))

    cur.executemany(
        f"""
        INSERT OR REPLACE INTO {table}
            (poly_id, member_count, area_km2, nextdown, geometry)
        VALUES (?, ?, ?, ?, GeomFromWKB(?, 4326))
        """,
        rows,
    )
    conn.commit()
    conn.close()
    print(f"  Wrote {len(rows)} rows → table '{table}'")


def write_hierarchy_and_counts(
    db_path: str,
    hierarchy: Dict[int, Dict[str, Optional[int]]],
    n_levels: int,
) -> None:
    """
    Write catchment_hierarchy and expected_counts to the database.

    hierarchy maps: base_unit_id → {level_name → parent_id | None}
    expected_counts is derived by counting base units per parent at each level.
    Also inserts L0 expected_counts (always 1) for the l0_basins table (smallest unit catchments).
    """
    level_names = [f"L{i}" for i in range(1, n_levels + 1)]
    col_names = ", ".join(f"l{i}_id" for i in range(1, n_levels + 1))
    placeholders = ", ".join(["?"] * (n_levels + 1))   # unit_id + N level cols

    conn = get_db_conn(db_path)

    # -- catchment_hierarchy ---------------------------------------------
    hier_rows = []
    for unit_id, lmap in hierarchy.items():
        vals = [lmap.get(ln) for ln in level_names]
        while len(vals) < n_levels:     # pad any missing levels to NULL
            vals.append(None)
        hier_rows.append((int(unit_id), *vals))

    conn.executemany(
        f"INSERT OR REPLACE INTO catchment_hierarchy (unit_id, {col_names}) "
        f"VALUES ({placeholders})",
        hier_rows,
    )

    # -- expected_counts — L0 (unit catchments, always 1) ----------------
    conn.executemany(
        "INSERT OR REPLACE INTO expected_counts (parent_id, level_name, expected_count) "
        "VALUES (?, ?, ?)",
        [(int(uid), "L0", 1) for uid in hierarchy],
    )

    # -- expected_counts — L1 through LN ---------------------------------
    for level_name in level_names:
        parent_counts: Dict[int, int] = {}
        for lmap in hierarchy.values():
            pid = lmap.get(level_name)
            if pid is not None:
                parent_counts[pid] = parent_counts.get(pid, 0) + 1

        conn.executemany(
            "INSERT OR REPLACE INTO expected_counts (parent_id, level_name, expected_count) "
            "VALUES (?, ?, ?)",
            [(int(pid), level_name, int(cnt)) for pid, cnt in parent_counts.items()],
        )

    conn.commit()
    conn.close()
    print(
        f"  Wrote catchment_hierarchy ({len(hierarchy)} base units) "
        f"and expected_counts (L0–L{n_levels})."
    )


# ---------------------------------------------------------------------------
# PostGIS / SQLAlchemy helpers (unchanged)
# ---------------------------------------------------------------------------

def get_engine() -> Engine:
    """Create the SQLAlchemy engine for the project's PostGIS database."""
    return create_engine(
        f"postgresql+psycopg2://{DBUSER}:{DBPW}@localhost:5432/{DBNAME}"
    )


def postgis_dissolve(
    gdf: gpd.GeoDataFrame,
    engine: Engine,
    by: str,
    agg_funcs: Dict[str, str],
    geom_type: str = "polygon",
    temp_table: str = "tmp_dissolve",
) -> gpd.GeoDataFrame:
    """
    Dissolve a GeoDataFrame using PostGIS.

    The GeoDataFrame is staged in a temporary PostGIS table, aggregated by
    `by`, and read back. ST_Union (polygons) and ST_LineMerge(ST_Union(...))
    (lines) handle invalid input and shared-endpoint line merging better
    than GeoPandas' dissolve.

    Parameters
    ----------
    gdf : GeoDataFrame
        Input data. Must contain the grouping column and every column listed
        in `agg_funcs`. CRS is preserved on output.
    engine : SQLAlchemy Engine
        Engine for a PostGIS-enabled PostgreSQL database.
    by : str
        Grouping column. If the GeoDataFrame's index has this name, it is
        promoted to a column first. On output it becomes the index, renamed
        to 'comid'.
    agg_funcs : dict[str, str]
        Map of non-geometry column names to SQL aggregate functions, e.g.
        {'unitarea': 'SUM', 'merge_count': 'SUM', 'nextdown': 'MAX'}.
    geom_type : {'polygon', 'line'}
        - 'polygon': ST_Union(ST_MakeValid(geom))
        - 'line':    ST_LineMerge(ST_Union(geom))
    temp_table : str
        Name of the staging table (replaced if it exists). Use distinct
        names when dissolving basins and rivers in the same run.
    """
    geom_col = gdf.geometry.name

    if gdf.index.name == by:
        gdf = gdf.reset_index()

    gdf.to_postgis(temp_table, engine, if_exists="replace", index=False)

    if geom_type == "line":
        geom_expr = f'ST_LineMerge(ST_Union("{geom_col}"))'
    elif geom_type == "polygon":
        geom_expr = f'ST_Union(ST_MakeValid("{geom_col}"))'
    else:
        raise ValueError(f"geom_type must be 'polygon' or 'line', got {geom_type!r}")

    select_clauses = [f'"{by}" AS comid']
    for col, func in agg_funcs.items():
        select_clauses.append(f'{func}("{col}") AS {col}')
    select_clauses.append(f"{geom_expr} AS geom")
    select_sql = ",\n        ".join(select_clauses)

    sql = f"""
    SELECT
        {select_sql}
    FROM "{temp_table}"
    GROUP BY "{by}";
    """

    out = read_postgis(sql, engine, geom_col="geom", crs=gdf.crs, index_col="comid")

    with engine.begin() as conn:
        conn.execute(text(f'DROP TABLE IF EXISTS "{temp_table}"'))

    return out


# ---------------------------------------------------------------------------
# Aggregation logic
# ---------------------------------------------------------------------------

def aggregate(
    basins_gdf: gpd.GeoDataFrame,
    rivers_gdf: Optional[gpd.GeoDataFrame],
    threshold_area: float,
    engine: Engine,
) -> Tuple[
    Optional[gpd.GeoDataFrame],
    Optional[gpd.GeoDataFrame],
    Dict[int, int],
]:
    """
    Perform one round of aggregation on the basins/rivers network.

    Returns
    -------
    agg_gdf : GeoDataFrame | None
        Dissolved basin polygons. None if no change occurred.
    rivers_agg_gdf : GeoDataFrame | None
        Dissolved river lines. None if no change occurred.
    unit_to_target : dict[int, int]
        Maps every input node ID to its parent (target) ID in agg_gdf.
        Used by main() to compose the per-level hierarchy.
    """
    basins_df = basins_gdf.drop(columns=basins_gdf.geometry.name)
    total_area = basins_df["unitarea"].sum()
    print(f"BEFORE SUM of unitarea: {total_area} km²")

    G = make_river_network(basins_df)

    print("  Consolidating the network...")
    n_before = G.number_of_nodes()
    S, merges, rivers2merge, rivers2delete = consolidate_network(
        G, threshold_area=threshold_area
    )
    n_after = G.number_of_nodes()

    if n_after == n_before:
        print("  No change in the number of nodes. Exiting.")
        return None, None, {}

    # Assign targets — every node gets one (identity if not in merges)
    basins_gdf["target"] = basins_gdf.index
    basins_gdf["merge_count"] = 1
    for node, target in merges.items():
        basins_gdf.at[node, "target"] = target

    # Capture the full input→target mapping before the dissolve consumes it
    unit_to_target: Dict[int, int] = {
        int(uid): int(tgt)
        for uid, tgt in zip(basins_gdf.index, basins_gdf["target"])
    }

    print("  Dissolving catchment polygons (PostGIS)...")
    agg_gdf = postgis_dissolve(
        basins_gdf,
        engine,
        by="target",
        agg_funcs={
            "unitarea": "SUM",
            "merge_count": "SUM",
            "nextdown": "MAX",
        },
        geom_type="polygon",
        temp_table="tmp_basins_dissolve",
    )

    total_area = agg_gdf["unitarea"].sum()
    print(f"AFTER SUM unitarea = {total_area}")

    # Drop singleton catchments that drain to nothing
    catchments_to_drop = (agg_gdf["merge_count"] == 1) & (agg_gdf["nextdown"] == 0)
    removed_catchments = set(catchments_to_drop[catchments_to_drop].index.tolist())
    agg_gdf = agg_gdf[~catchments_to_drop].copy()

    agg_gdf.geometry = agg_gdf.geometry.apply(lambda p: _close_holes(p, 0))

    agg_gdf["nextdown"] = 0
    for idx in agg_gdf.index:
        try:
            nextdown = list(G.successors(idx))[0]
            if nextdown in removed_catchments:
                nextdown = 0
        except Exception:
            nextdown = 0
        agg_gdf.at[idx, "nextdown"] = nextdown

    if MAP:
        fig, axes = plt.subplots(1, 2)
        basins_gdf.plot(ax=axes[0])
        agg_gdf.plot(ax=axes[1])

        if COUNTRY is not None:
            shp = r"C:\Data\GIS\Natural_Earth\ne_010m\countries\ne_10m_admin_0_countries.shp"
            countries = gpd.read_file(shp)
            country = countries[countries["NAME"] == COUNTRY]
            country.boundary.plot(ax=axes[0], color="black", linewidth=1)
            country.boundary.plot(ax=axes[1], color="black", linewidth=1)

        plt.show()

    show_area_stats(agg_gdf)

    # Process rivers
    if RIVERS:
        print("  Processing rivers...")
        rivers_gdf = rivers_gdf.drop(rivers2delete, errors="ignore")
        rivers_agg_gdf = rivers_gdf

        if len(rivers2merge) > 0:
            rivers_gdf["target"] = rivers_gdf.index
            for target, node_list in rivers2merge.items():
                for node in node_list:
                    if node in rivers_gdf.index:
                        rivers_gdf.at[node, "target"] = target

            print("  Dissolving rivers (PostGIS, ST_LineMerge)...")
            rivers_agg_gdf = postgis_dissolve(
                rivers_gdf,
                engine,
                by="target",
                agg_funcs={"lengthkm": "SUM"},
                geom_type="line",
                temp_table="tmp_rivers_dissolve",
            )

            for idx in rivers_agg_gdf.index:
                try:
                    nextdown = list(G.successors(idx))[0]
                except Exception:
                    nextdown = 0
                rivers_agg_gdf.at[idx, "nextdown"] = nextdown

                if idx in G.nodes:
                    rivers_agg_gdf.at[idx, "strahler_order"] = G.nodes[idx].get("strahler_order", 0)
                    rivers_agg_gdf.at[idx, "shreve_order"] = G.nodes[idx].get("shreve_order", 0)
                else:
                    rivers_agg_gdf.at[idx, "strahler_order"] = 0
                    rivers_agg_gdf.at[idx, "shreve_order"] = 0

            rivers_agg_gdf = rivers_agg_gdf[
                rivers_agg_gdf.index.isin(agg_gdf.index)
            ].copy()

            if MAP:
                fig, axes = plt.subplots(1, 2)
                rivers_gdf.plot(ax=axes[0])
                rivers_agg_gdf.plot(ax=axes[1])
                plt.show()
    else:
        rivers_agg_gdf = None

    return agg_gdf, rivers_agg_gdf, unit_to_target


def get_geodata(engine: Engine) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Load basins and rivers for BASIN from PostGIS."""
    sql_basins = f"""
    SELECT r.comid, r.nextdownid AS nextdown, b.unitarea, b.geom_detail AS geom
    FROM merit_rivers r
    JOIN merit_basins b ON r.comid = b.comid
    WHERE r.basin = {BASIN};
    """
    basins_gdf = read_postgis(
        sql_basins, engine, geom_col="geom", crs="EPSG:4326", index_col="comid"
    )
    if RIVERS:
        sql_rivers = f"""
        SELECT comid, lengthkm, geom_simple AS geom
        FROM merit_rivers
        WHERE basin = {BASIN};
        """
        rivers_gdf = read_postgis(
            sql_rivers, engine, geom_col="geom", crs="EPSG:4326", index_col="comid"
        )
    else:
        rivers_gdf = None

    return basins_gdf, rivers_gdf


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Area thresholds (km²) for each aggregation level.
    # Add or remove entries to change the number of levels (up to 4).
    THRESHOLDS = [1e3, 1e4, 1e5, 1e6]   # -> L1 through L4

    DB_PATH = f"output/basins{BASIN}.db"
    Path("output").mkdir(exist_ok=True)

    engine = get_engine()
    basins_gdf, rivers_gdf = get_geodata(engine)

    n_levels = len(THRESHOLDS)

    # Create all tables up front
    init_db(DB_PATH, n_levels=n_levels)

    # Write base unit catchments to unit catchments table `l0_basins`
    write_basin_table(DB_PATH, basins_gdf, "l0_basins")

    # hierarchy[unit_id][level_name] = parent_id (or None if no valid parent)
    # Keyed on base unit comids; a new entry is added each iteration.
    hierarchy: Dict[int, Dict[str, Optional[int]]] = {
        int(uid): {} for uid in basins_gdf.index
    }

    level_num = 0
    for threshold in THRESHOLDS:
        level_num += 1
        level_name = f"L{level_num}"
        table_name = f"l{level_num}_basins"
        print(f"\n*** {level_name}: area threshold = {threshold:,.0f} km² ***")

        agg_gdf, rivers_gdf, unit_to_target = aggregate(
            basins_gdf, rivers_gdf, threshold, engine
        )

        if agg_gdf is None:
            print("  No change. Stopping early.")
            level_num -= 1   # don't count a level that produced no output
            break

        valid_ids = set(agg_gdf.index.astype(int))

        # Compose the mapping: base_unit -> parent at this level.
        # At L1 the input IDs are base unit comids.
        # At L2+ the input IDs are the previous level's poly comids,
        # so we follow the chain through hierarchy[uid][prev_level].
        for uid in hierarchy:
            if level_num == 1:
                prev_id = uid
            else:
                prev_level = f"L{level_num - 1}"
                prev_id = hierarchy[uid].get(prev_level)

            if prev_id is None:
                # Chain already broken at a prior level; propagate None
                hierarchy[uid][level_name] = None
                continue

            target = unit_to_target.get(int(prev_id))

            # A target that ended up in removed_catchments is not in agg_gdf
            if target is not None and int(target) in valid_ids:
                hierarchy[uid][level_name] = int(target)
            else:
                hierarchy[uid][level_name] = None

        # Write this level's pre-dissolved polygons to its own table
        write_basin_table(
            DB_PATH, agg_gdf, table_name, member_count_col="merge_count"
        )

        # Advance to the next level
        basins_gdf = agg_gdf

    # Write the full hierarchy and expected_counts derived from it
    write_hierarchy_and_counts(DB_PATH, hierarchy, level_num)

    print(f"\nDone. SpatiaLite database written to: {DB_PATH}")


# *** GLOBAL PARAMETERS ***
RIVERS = False
COUNTRY = None
MAP = False


if __name__ == "__main__":
    basins = [55,56,57,61,62,63,64,65,66,67,71,72,73,74,
              75,76,78,81,82,83,84,85,86]

    for basin in basins:
        BASIN = basin
        main()
