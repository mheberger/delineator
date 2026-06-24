"""
Aggregate MERIT-Basins unit catchments to larger sizes and write to SpatiaLite.

Geoprocessing (dissolve) is performed via

Note:
The previous GeoPandas dissolve of rivers produced MultiLineStrings that
were not merged at shared endpoints; ST_LineMerge(ST_Union(...)) in PostGIS
collapses connecting LineStrings into single LineStrings.

Output is a single SpatiaLite .db file containing one table per level:

  Spatial tables (each with a geometry column and spatial index):
    l0_basins  — original base unit polygons
    l1_basins  — aggregated at ~1 000 km² threshold
    l2_basins  — aggregated at ~10 000 km² threshold
    l3_basins  — aggregated at ~100 000 km² threshold
    l4_basins  — aggregated at ~1 000 000 km² threshold

  Non-spatial tables:
    catchment_hierarchy — maps each base unit to its L1–L4 parent comid
    expected_counts     — number of base units inside each pre-staged polygon
"""

import sqlite3
from pathlib import Path
from typing import Dict, Optional, Tuple

import geopandas as gpd
from matplotlib import pyplot as plt
import shapely
from shapely.geometry import MultiPolygon, Polygon
from shapely import wkb as shapely_wkb

from consolidate import consolidate_network, show_area_stats
from delineator.util import _close_holes
from graph_tools import make_river_network


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
        comid        — comid from MERIT-Basins (Primary Key)
        member_count — number of unit catchments dissolved into this polygon
        area_km2     — total area
        geometry     — MULTIPOLYGON, EPSG:4326
    """
    cur.execute(f"""
        CREATE TABLE IF NOT EXISTS {table} (
            comid      INTEGER, 
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
        l1_basins ... lN_basins — one per aggregation level
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
            l0_id  INTEGER,
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
        Source data. Index must be the comid. Must contain
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
    for comid, row in gdf.iterrows():
        geom_bytes = _to_multipolygon_wkb(row.geom)
        if geom_bytes is None:
            continue
        rows.append((
            int(comid),
            int(row[member_count_col]) if member_count_col else 1,
            round(float(row.get("unitarea", 0.0)), 2),
            int(row.get("nextdown", 0)),
            geom_bytes,
        ))

    cur.executemany(
        f"""
        INSERT OR REPLACE INTO {table}
            (comid, member_count, area_km2, nextdown, geometry)
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

    hierarchy maps: l0_id → {level_name → parent_id | None}
    expected_counts is derived by counting base units per parent at each level.
    Also inserts L0 expected_counts (always 1) for the l0_basins table (smallest unit catchments).
    """
    level_names = [f"L{i}" for i in range(1, n_levels + 1)]
    col_names = ", ".join(f"l{i}_id" for i in range(1, n_levels + 1))
    placeholders = ", ".join(["?"] * (n_levels + 1))   # l0_id + N level cols

    conn = get_db_conn(db_path)

    # -- catchment_hierarchy ---------------------------------------------
    hier_rows = []
    for unit_id, lmap in hierarchy.items():
        vals = [lmap.get(ln) for ln in level_names]
        while len(vals) < n_levels:     # pad any missing levels to NULL
            vals.append(None)
        hier_rows.append((int(unit_id), *vals))

    conn.executemany(
        f"INSERT OR REPLACE INTO catchment_hierarchy (l0_id, {col_names}) "
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
# Aggregation logic
# ---------------------------------------------------------------------------

def aggregate(
    basins_gdf: gpd.GeoDataFrame,
    rivers_gdf: Optional[gpd.GeoDataFrame],
    threshold_area: float,
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

    # DEBUG code: 2026-06-21
    import networkx as nx

    if not nx.is_directed_acyclic_graph(G):
        cycle = nx.find_cycle(G)
        print("CYCLE:", cycle)
        for u, v in cycle:
            row = basins_df.loc[u] if u in basins_df.index else None
            print(
                f"  {u} -> nextdown={row['nextdown'] if row is not None else '?'}, "
                f"area={row.get('unitarea') if row is not None else '?'}"
            )
        raise SystemExit("stopping to inspect")

    # END DEBUG code

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

    print("  Dissolving catchment polygons...")
    agg_gdf = dissolve(
        basins_gdf,
        by="target",
        agg_funcs={
            "unitarea": "SUM",
            "merge_count": "SUM",
            "nextdown": "MAX",
        },
        geom_type="polygon",
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

            print("  Dissolving rivers...")
            rivers_agg_gdf = dissolve(
                rivers_gdf,
                by="target",
                agg_funcs={"lengthkm": "SUM"},
                geom_type="line",
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


# ---------------------------------------------------------------------------
# Data source: read basins + rivers from shapefiles
# ---------------------------------------------------------------------------

def _find_col(df, needle: str) -> str:
    """Case/substring-tolerant column lookup (e.g. 'nextdown' -> 'NextDownID')."""
    for c in df.columns:
        if c.lower() == needle.lower():
            return c
    for c in df.columns:
        if needle.lower() in c.lower():
            return c
    raise KeyError(f"No column matching {needle!r} in {list(df.columns)}")


def get_geodata() -> Tuple[gpd.GeoDataFrame, Optional[gpd.GeoDataFrame]]:
    """
    Load catchments and rivers for the current BASIN from shapefiles.

    Topology (`nextdown`) is taken from the
    rivers layer's NextDownID and attached to the catchment polygons on COMID.

    Returns
    -------
    basins_gdf : GeoDataFrame
        Indexed by comid (int). Columns: 'unitarea', 'nextdown'; active
        geometry column named 'geom' (polygons), CRS as read from the .prj.
    rivers_gdf : GeoDataFrame | None
        Only if RIVERS is True. Indexed by comid; columns 'lengthkm';
        active geometry 'geom' (lines).
    """
    cat_path = CAT_SHP.format(basin=BASIN)
    riv_path = RIV_SHP.format(basin=BASIN)

    # --- catchments (polygons): comid, unitarea, geometry ---
    basins_gdf = gpd.read_file(cat_path).rename(columns={"COMID": "comid"})
    basins_gdf["comid"] = basins_gdf["comid"].astype(int)

    # --- rivers: carry the topology (and lengthkm for the rivers path) ---
    rivers_raw = gpd.read_file(riv_path).rename(columns={"COMID": "comid"})
    rivers_raw["comid"] = rivers_raw["comid"].astype(int)

    # New 2026-06-21 In the MERIT-Basins 'bugfix' version, endorheic basins have
    # the value of their own COMID in the NextDownID field, which creates a loop
    rivers_raw = rivers_raw[(rivers_raw["NextDownID"] != rivers_raw["comid"])].copy()
    # If rivers have zero length, drop (some non-existent rivers have lengthkm = 9e-5 (!)
    rivers_gdf = rivers_raw[rivers_raw["lengthkm"] > 1e-4].copy()
    print("null river geoms:", rivers_raw.geometry.isna().sum())
    print("empty river geoms:", rivers_gdf.geometry.is_empty.sum())

    nextdown_col = _find_col(rivers_gdf, "nextdown")          # e.g. 'NextDownID'
    topo = rivers_gdf[["comid", nextdown_col]].rename(
        columns={nextdown_col: "nextdown"}
    )

    # attach nextdown onto the catchments (the old r JOIN b ON r.comid = b.comid)
    basins_gdf = basins_gdf.merge(topo, on="comid", how="left")
    # catchments with no matching river row are terminal sinks or coastal catchments
    basins_gdf["nextdown"] = basins_gdf["nextdown"].fillna(0).astype(int)

    # match the layout the rest of the script expects: geom col 'geom', index 'comid'
    basins_gdf = basins_gdf.rename_geometry("geom").set_index("comid")

    if RIVERS:
        geom_name = rivers_gdf.geometry.name
        rivers_gdf = (
            rivers_gdf[["comid", "lengthkm", geom_name]]
            .rename_geometry("geom")
            .set_index("comid")
        )
    else:
        rivers_gdf = None

    return basins_gdf, rivers_gdf


# ---------------------------------------------------------------------------
# Geoprocessing: in-process dissolve (drop-in for postgis_dissolve)
# ---------------------------------------------------------------------------

_SQL_TO_PANDAS = {
    "SUM": "sum", "MAX": "max", "MIN": "min",
    "MEAN": "mean", "AVG": "mean", "COUNT": "count",
}


def dissolve(
    gdf: gpd.GeoDataFrame,
    by: str,
    agg_funcs: Dict[str, str],
    geom_type: str = "polygon",
    **_ignored,            # tolerate a stray temp_table=... from old call sites
) -> gpd.GeoDataFrame:
    """
    Dissolve a GeoDataFrame in-process with shapely 2.0 (GEOS).

    Same contract as the former postgis_dissolve(): group by `by`, aggregate
    the attribute columns per `agg_funcs`, combine geometry, and return a
    GeoDataFrame indexed by 'comid' with the active geometry column 'geom'.

    Geometry combine:
        polygon -> shapely.union_all(shapely.make_valid(geoms))
        line    -> shapely.line_merge(shapely.union_all(geoms))
    """
    geom_col = gdf.geometry.name

    if gdf.index.name == by:
        gdf = gdf.reset_index()

    # --- attribute aggregation ---
    agg_map = {col: _SQL_TO_PANDAS[func.upper()] for col, func in agg_funcs.items()}
    attrs = gdf.groupby(by, sort=False).agg(agg_map)

    # --- geometry aggregation ---
    if geom_type == "polygon":
        def _combine(s):
            return shapely.union_all(shapely.make_valid(s.to_numpy()))
    elif geom_type == "line":
        def _combine(s):
            return shapely.line_merge(shapely.union_all(s.to_numpy()))
    else:
        raise ValueError(f"geom_type must be 'polygon' or 'line', got {geom_type!r}")

    geoms = gdf.groupby(by, sort=False)[geom_col].apply(_combine)

    merged = attrs.join(geoms.rename("geom"))
    merged["geom"] = gpd.GeoSeries(
        merged["geom"].values, index=merged.index, crs=gdf.crs
    )
    out = gpd.GeoDataFrame(merged, geometry="geom", crs=gdf.crs)
    out.index.name = "comid"
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Area thresholds (km²) for each aggregation level.
    # Add or remove entries to change the number of levels (up to 4).
    THRESHOLDS = [1e3, 1e4, 1e5, 1e6]   # -> L1 through L4

    basins_gdf, rivers_gdf = get_geodata()

    n_levels = len(THRESHOLDS)

    # Create all tables up front
    init_db(DB_PATH, n_levels=n_levels)

    # Write base unit catchments to unit catchments table `l0_basins`
    write_basin_table(DB_PATH, basins_gdf, "l0_basins")

    # hierarchy[l0_id][level_name] = parent_id (or None if no valid parent)
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
            basins_gdf, rivers_gdf, threshold)

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
COUNTRY = 'France'
MAP = False


if __name__ == "__main__":
    basins = [11]

    for BASIN in basins:
        DB_PATH = f"C:/Users/mheberger/AppData/Local/delineator/vector/basins{BASIN}.db"
        CAT_SHP = rf"C:\Data\GIS\MERITBasins\catchments\bugfix_repaired\cat_pfaf_{BASIN}_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"
        RIV_SHP = rf"C:\Data\GIS\MERITBasins\rivers\bugfix_repaired\riv_pfaf_{BASIN}_MERIT_Hydro_v07_Basins_v01_bugfix1.shp"

        main()
