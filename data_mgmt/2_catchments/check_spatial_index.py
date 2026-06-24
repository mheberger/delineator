"""
Verify that every spatial layer in a SpatiaLite .db file has a spatial index
that covers ALL of its features.

For each layer registered in geometry_columns it reports:
  features   - row count in the table
  indexed    - entries in the R*Tree (idx_<table>_<geom>)
  null geom  - rows with NULL geometry (legitimately NOT indexable)
  check      - SpatiaLite's own CheckSpatialIndex() verdict

A layer is "covered" when CheckSpatialIndex() == 1, which guarantees one
correct R*Tree entry per non-NULL geometry. The count columns are shown so
any shortfall is easy to interpret.

Usage:
    python check_spatial_index.py output/basins27.db
    python check_spatial_index.py output/basins*.db      # check many at once
"""
import sqlite3


def get_db_conn(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.enable_load_extension(True)
    conn.load_extension("mod_spatialite")
    return conn


def check_db(db_path: str) -> bool:
    conn = get_db_conn(db_path)
    cur = conn.cursor()

    layers = cur.execute(
        """SELECT f_table_name, f_geometry_column, spatial_index_enabled
           FROM geometry_columns
           ORDER BY f_table_name"""
    ).fetchall()

    print(f"\n{db_path}")
    if not layers:
        print("  (no spatial layers registered in geometry_columns)")
        conn.close()
        return False

    print(f"  {'layer':<16}{'features':>10}{'indexed':>10}{'null geom':>11}{'check':>8}")
    all_ok = True
    for tbl, geom, sie in layers:
        n_feat = cur.execute(f'SELECT count(*) FROM "{tbl}"').fetchone()[0]
        n_null = cur.execute(
            f'SELECT count(*) FROM "{tbl}" WHERE "{geom}" IS NULL'
        ).fetchone()[0]

        if sie != 1:
            print(f"  {tbl:<16}{n_feat:>10}{'--':>10}{n_null:>11}{'NO IDX':>8}")
            all_ok = False
            continue

        n_idx = cur.execute(f'SELECT count(*) FROM "idx_{tbl}_{geom}"').fetchone()[0]
        chk = cur.execute(
            "SELECT CheckSpatialIndex(?, ?)", (tbl, geom)
        ).fetchone()[0]

        if chk == 1:
            status = "OK"
        elif chk == 0:
            status = "FAIL"      # index exists but is out of sync / incomplete
            all_ok = False
        else:
            status = "ERR"       # NULL -> could not be evaluated
            all_ok = False

        # cross-check the counts (redundant when chk==1, but informative)
        if chk == 1 and n_idx != n_feat - n_null:
            status = "MISS"
            all_ok = False

        print(f"  {tbl:<16}{n_feat:>10}{n_idx:>10}{n_null:>11}{status:>8}")

    print("  => " + ("ALL LAYERS COVERED" if all_ok else "PROBLEMS FOUND (see above)"))
    if not all_ok:
        raise Exception("Some layers are not covered by an index")

    return all_ok


if __name__ == "__main__":
    basins = [11,12,13,14,15,16,17,18,21,22,23,24,25,26,27,28,29,31,32,33,34,35,36,41,42,43,44,45,46,47,48,49,51,52,53,54,
          55,56,57,61,62,63,64,65,66,67,71,72,73,74,75,76,77,78,81,82,83,84,85,86]
    for basin in basins:
        db = rf"C:\Users\mheberger\AppData\Local\delineator\vector\basins{basin}.db"
        check_db(db)
