# =============================================================================
# create_merit_sqlite_indexes.py
#
# Purpose:
#   Create normal SQLite attribute indexes once for all delineator MERIT
#   vector databases.
#
# Why:
#   The RTree spatial indexes already exist and work. The slow part was normal
#   attribute lookup, especially l0_basins.nextdown and l0_basins.comid.
# =============================================================================

from pathlib import Path
import sqlite3
import time


DATA_DIR = Path(
    r"C:\Users\mheberger\AppData\Local\delineator"
)

VECTOR_DIR = DATA_DIR / "vector"


def table_exists(conn, table):
    row = conn.execute(
        """
        SELECT name
        FROM sqlite_master
        WHERE type = 'table'
          AND name = ?
        """,
        (table,),
    ).fetchone()

    return row is not None


def create_index(conn, table, column, index_name):
    if not table_exists(conn, table):
        print(f"    skip: table missing {table}")
        return

    t0 = time.perf_counter()

    conn.execute(
        f"CREATE INDEX IF NOT EXISTS {index_name} ON {table}({column})"
    )
    conn.commit()

    print(f"    {index_name:<45} {time.perf_counter() - t0:>8.2f} s")


def index_basins_db(db_file):
    print("")
    print("=" * 90)
    print(f"Basins DB: {db_file}")
    print("=" * 90)

    conn = sqlite3.connect(db_file)

    try:
        conn.execute("PRAGMA temp_store = MEMORY")
        conn.execute("PRAGMA cache_size = -200000")
        conn.execute("PRAGMA automatic_index = ON")

        create_index(
            conn,
            "l0_basins",
            "comid",
            "idx_l0_basins_comid_attr",
        )

        create_index(
            conn,
            "l0_basins",
            "nextdown",
            "idx_l0_basins_nextdown_attr",
        )

        create_index(
            conn,
            "catchment_hierarchy",
            "l0_id",
            "idx_catchment_hierarchy_l0_id_attr",
        )

        create_index(
            conn,
            "expected_counts",
            "parent_id",
            "idx_expected_counts_parent_id_attr",
        )

    finally:
        conn.close()


def index_rivers_db(db_file):
    print("")
    print("=" * 90)
    print(f"Rivers DB: {db_file}")
    print("=" * 90)

    conn = sqlite3.connect(db_file)

    try:
        conn.execute("PRAGMA temp_store = MEMORY")
        conn.execute("PRAGMA cache_size = -200000")
        conn.execute("PRAGMA automatic_index = ON")

        create_index(
            conn,
            "rivers",
            "comid",
            "idx_rivers_comid_attr",
        )

    finally:
        conn.close()


def main():
    if not VECTOR_DIR.exists():
        raise FileNotFoundError(VECTOR_DIR)

    basins_files = sorted(VECTOR_DIR.glob("basins*.db"))
    rivers_files = sorted(VECTOR_DIR.glob("rivers*.db"))

    print(f"Vector dir: {VECTOR_DIR}")
    print(f"Basins DB files: {len(basins_files)}")
    print(f"Rivers DB files: {len(rivers_files)}")

    t_all = time.perf_counter()

    for db_file in basins_files:
        index_basins_db(db_file)

    for db_file in rivers_files:
        index_rivers_db(db_file)

    print("")
    print("=" * 90)
    print(f"Indexing complete in {time.perf_counter() - t_all:.2f} s")
    print("=" * 90)


if __name__ == "__main__":
    main()