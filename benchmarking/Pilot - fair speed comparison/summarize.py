"""
Summarize a pilot run: medians and spread, setup costs, speedups, area
agreement, resource use, and how many repetitions the main run needs.

    python summarize.py results/<timestamp>

Writes summary.md and CSV tables next to results.jsonl.

Metrics
- warm, per outlet: total_s, the delineation alone (setup excluded; the
  warm-up run is dropped).
- cold, per outlet: wall_since_launch_s, from starting a fresh process to
  having the polygon, after purging that tool's data from the file cache.
  Also shown without Python start-up and imports, which matter for the two
  Python tools but are only harness glue for TauDEM.
- Spread: IQR and the coefficient of variation (CV = sd / mean) over
  repetitions.
- Repetitions needed: N such that the standard error of the median,
  approximated as 1.2533 * sd / sqrt(N), is at most TARGET_REL_SE of the
  median. Conservative, because the worst outlet decides.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TARGET_REL_SE = 0.02
TOOL_ORDER = ["delineator", "pysheds", "taudem"]


def md_table(df: pd.DataFrame, floatfmt: str = "{:,.3g}") -> str:
    def fmt(v):
        if isinstance(v, (float, np.floating)):
            if math.isnan(v):
                return ""
            return f"{v:,.0f}" if abs(v) >= 1000 else floatfmt.format(v)
        return str(v)
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines)


def load(results_dir: Path):
    records = [json.loads(line) for line in (results_dir / "results.jsonl").open(encoding="utf-8")]
    meta = json.loads((results_dir / "metadata.json").read_text())
    df = pd.DataFrame(records)
    return df, meta


def spread(values: pd.Series) -> pd.Series:
    v = values.dropna()
    mean = v.mean()
    return pd.Series({"n": len(v), "median": v.median(), "q25": v.quantile(0.25),
                      "q75": v.quantile(0.75), "min": v.min(), "max": v.max(),
                      "cv": v.std(ddof=1) / mean if len(v) > 1 and mean else np.nan})


def main(results_dir: Path) -> None:
    df, meta = load(results_dir)
    out = []

    outlets = df[(df.kind == "outlet")].copy()
    outlets["id"] = outlets["id"].astype("int64")
    ok = outlets[(outlets.status == "ok") & (~outlets.warmup.astype(bool))].copy()
    setups = df[df.kind == "setup"].copy()

    # Cold runs: subtract each worker's start-up and imports to get the data-dependent part
    if "scenario" in ok and (ok.scenario == "cold").any():
        s = (setups[["block", "stdlib_import_s", "import_s"]]
             .rename(columns={"stdlib_import_s": "worker_stdlib_s", "import_s": "worker_import_s"}))
        ok = ok.merge(s, on="block", how="left")
        ok["cold_excl_imports_s"] = np.where(
            ok.scenario == "cold",
            ok.wall_since_launch_s - ok.worker_stdlib_s.fillna(0) - ok.worker_import_s.fillna(0), np.nan)

    # Reference area per outlet: delineator's, else the first tool's
    area = ok.groupby(["id", "tool"]).area_km2.median().unstack()
    ref_area = area["delineator"] if "delineator" in area else area.iloc[:, 0]

    def metric(scenario):
        return "total_s" if scenario == "warm" else "wall_since_launch_s"

    # 1. Per-outlet timing
    rows = []
    for (scenario, tool, oid), g in ok.groupby(["scenario", "tool", "id"]):
        st = spread(g[metric(scenario)])
        row = {"scenario": scenario, "tool": tool, "id": oid, "megabasin": g.megabasin.iloc[0],
               "area_km2": ref_area.get(oid, np.nan), "metric": metric(scenario), **st.to_dict()}
        if scenario == "cold":
            row["median_excl_imports"] = g.cold_excl_imports_s.median()
            row["purge_ok"] = bool(g.purge_ok.all())
        rows.append(row)
    timings = pd.DataFrame(rows)
    timings["tool"] = pd.Categorical(timings.tool, TOOL_ORDER)
    timings = timings.sort_values(["scenario", "tool", "area_km2"])
    timings.to_csv(results_dir / "timings.csv", index=False)

    out.append(f"# Pilot summary: {results_dir.name}\n")
    cpu, storage = meta.get("cpu", {}), meta.get("storage", {})
    out.append(f"- Machine: {cpu.get('name')}, {cpu.get('physical_cores')} cores, "
               f"{meta.get('memory_mb', {}).get('total_mb', 0) / 1024:.0f} GB RAM; "
               + "; ".join(f"{d}: {s.get('Disk')} ({s.get('BusType')} {s.get('MediaType')})"
                           for d, s in storage.items() if s))
    out.append(f"- Power: AC={meta.get('on_ac_power')}, {meta.get('power_scheme')}")
    av = meta.get("antivirus") or {}
    if av:
        folders = av.get("data_folders_excluded") or {}
        if av.get("exclusions_readable"):
            not_excluded = [f for f, excluded in folders.items() if not excluded]
            exclusions = (f"{len(folders) - len(not_excluded)} of {len(folders)} data folders excluded"
                          + (f"; not excluded: {', '.join(not_excluded)}" if not_excluded else ""))
        else:
            exclusions = "exclusions unknown (not run as administrator)"
        out.append(f"- Antivirus: {', '.join(av.get('antivirus_products') or ['none found'])}; "
                   f"real-time protection {'on' if av.get('RealTimeProtectionEnabled') else 'off'}; "
                   f"{exclusions}")
    out.append(f"- Packages: {meta.get('packages')}")
    end = df[df.kind == "end"]
    if "tool_version" in end and end.tool_version.notna().any():
        out.append(f"- TauDEM: {end.tool_version.dropna().iloc[0]}")
    out.append(f"- delineator repo: {meta.get('delineator_repo', {}).get('describe')}")
    out.append(f"- Elapsed: {meta.get('elapsed_s', 0) / 60:.0f} min\n")

    for scenario in ("warm", "cold"):
        t = timings[timings.scenario == scenario]
        if t.empty:
            continue
        title = ("Warm batch: seconds per outlet, setup excluded" if scenario == "warm" else
                 "Cold one-off: seconds from process launch to polygon, data purged from cache")
        cols = ["tool", "id", "megabasin", "area_km2", "n", "median", "q25", "q75", "cv"]
        if scenario == "cold":
            cols += ["median_excl_imports", "purge_ok"]
        out.append(f"## {title}\n\n{md_table(t[cols])}\n")

    # 2. Setup costs
    if not setups.empty:
        good = setups[setups.status == "ok"]
        setup_tab = (good.groupby(["scenario", "tool", "megabasin"])
                     .agg(n=("import_s", "size"), import_s=("import_s", "median"),
                          setup_s=("setup_s", "median"), read_mb=("read_bytes", lambda x: x.median() / 1e6))
                     .reset_index())
        setup_tab.to_csv(results_dir / "setup.csv", index=False)
        out.append("## One-off costs per worker (medians)\n\nimport_s: importing the tool's "
                   "Python libraries. setup_s: pysheds reading the whole flow-direction raster; "
                   "nothing for the others. read_mb includes the imports.\n\n"
                   f"{md_table(setup_tab)}\n")

    # 3. Speedups, paired by outlet
    rows = []
    for scenario, col in (("warm", "total_s"), ("cold", "wall_since_launch_s"), ("cold", "cold_excl_imports_s")):
        sub = ok[ok.scenario == scenario] if "scenario" in ok else ok.iloc[0:0]
        if sub.empty or col not in sub:
            continue
        med = sub.groupby(["id", "tool"])[col].median().unstack()
        if "delineator" not in med:
            continue
        for tool in [t for t in ("pysheds", "taudem") if t in med]:
            ratio = (med[tool] / med["delineator"]).dropna()
            for oid, r in ratio.items():
                rows.append({"scenario": scenario, "metric": col, "tool": tool, "id": oid,
                             "area_km2": ref_area.get(oid, np.nan), "times_slower_than_delineator": r})
    if rows:
        speed = pd.DataFrame(rows).sort_values(["scenario", "metric", "tool", "area_km2"])
        speed.to_csv(results_dir / "speedups.csv", index=False)
        agg = (speed.groupby(["scenario", "metric", "tool"]).times_slower_than_delineator
               .agg(n="size", median="median", q25=lambda x: x.quantile(0.25),
                    q75=lambda x: x.quantile(0.75), min="min", max="max").reset_index())
        out.append("## How many times slower than delineator (other tool's median / delineator's "
                   "median, paired by outlet)\n\n"
                   f"{md_table(agg)}\n\nPer outlet: speedups.csv\n")

    # 4. Area agreement
    if not area.empty:
        areas = area.copy()
        areas["max_rel_diff_vs_delineator"] = (
            area.sub(ref_area, axis=0).abs().max(axis=1) / ref_area if "delineator" in area else np.nan)
        areas = areas.reset_index().sort_values("max_rel_diff_vs_delineator", na_position="last")
        areas.to_csv(results_dir / "areas.csv", index=False)
        out.append("## Area agreement (km², ellipsoidal, same method for every tool)\n\n"
                   f"{md_table(areas, '{:,.6g}')}\n")

    # 5. Phases (warm medians)
    warm = ok[ok.scenario == "warm"] if "scenario" in ok else ok.iloc[0:0]
    if not warm.empty:
        ph = pd.concat([warm[["tool", "id"]].reset_index(drop=True),
                        pd.json_normalize(warm.phases.tolist())], axis=1)
        ph = ph.groupby(["tool", "id"]).median().reset_index()
        ph.insert(2, "area_km2", ph.id.map(ref_area))
        ph = ph.sort_values(["tool", "area_km2"])
        ph.to_csv(results_dir / "phases.csv", index=False)
        for tool, g in ph.groupby("tool"):
            g = g.dropna(axis=1, how="all")
            out.append(f"## Phases, warm medians (s): {tool}\n\n{md_table(g.drop(columns='tool'))}\n")

    # 6. Resource use
    res = ok.assign(cpu_s=ok.cpu_user_s + ok.cpu_kernel_s, read_mb=ok.read_bytes / 1e6,
                    write_mb=ok.write_bytes / 1e6)
    res["cpu_per_wall"] = res.cpu_s / res.total_s
    res_tab = (res.groupby(["scenario", "tool", "id"])[["read_mb", "write_mb", "cpu_per_wall"]]
               .median().reset_index())
    res_tab.insert(3, "area_km2", res_tab.id.map(ref_area))
    peak = end.groupby(["scenario", "tool", "megabasin"]).peak_memory_mb.max().reset_index()
    res_tab.to_csv(results_dir / "resources.csv", index=False)
    peak.to_csv(results_dir / "peak_memory.csv", index=False)
    out.append("## Resource use per delineation (medians)\n\nread_mb: logical bytes read by the "
               "process tree, cache hits included. cpu_per_wall above 1 means several cores were "
               "busy (TauDEM's MPI ranks also spin while waiting).\n\n"
               f"{md_table(res_tab.sort_values(['scenario', 'tool', 'area_km2']))}\n")
    out.append("Peak committed memory per worker (MB; warm workers include setup):\n\n"
               f"{md_table(peak)}\n")

    # 7. Repetitions needed
    rows = []
    for (scenario, tool), g in timings.groupby(["scenario", "tool"], observed=True):
        cv = g.cv.dropna()
        if cv.empty:
            continue
        need = math.ceil((1.2533 * cv.max() / TARGET_REL_SE) ** 2)
        rows.append({"scenario": scenario, "tool": tool, "median_cv": cv.median(), "max_cv": cv.max(),
                     "worst_outlet": g.loc[g.cv.idxmax(), "id"], "n_for_2pct_se": need})
    if rows:
        out.append("## Repetitions needed for the main run\n\nN so that the standard error of each "
                   f"outlet's median is within {TARGET_REL_SE:.0%} of it, set by the noisiest outlet.\n\n"
                   f"{md_table(pd.DataFrame(rows))}\n")

    # 8. Problems
    bad = outlets[outlets.status != "ok"]
    failed_blocks = df[df.kind == "block_failed"] if "kind" in df else df.iloc[0:0]
    if not bad.empty or not failed_blocks.empty:
        out.append("## Failures\n")
        for r in bad.itertuples():
            out.append(f"- {r.scenario} {r.tool} outlet {r.id}: {getattr(r, 'error', '')}")
        for r in failed_blocks.itertuples():
            out.append(f"- block {r.block} failed: return code {r.returncode}, timed out {r.timed_out}")
        out.append("")
    if "purge_ok" in ok and (ok.scenario == "cold").any():
        cold = ok[ok.scenario == "cold"]
        out.append(f"Cache purge verified for {int(cold.purge_ok.sum())} of {len(cold)} cold runs.\n")

    text = "\n".join(out)
    (results_dir / "summary.md").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    main(Path(sys.argv[1]))
