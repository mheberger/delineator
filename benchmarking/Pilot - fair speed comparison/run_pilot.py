"""
Pilot benchmark: delineator vs. pysheds vs. TauDEM, warm and cold, repeated.

Two scenarios:

warm   Batch use. A fresh worker process per (tool, megabasin) block
       imports the tool, does its one-off setup (pysheds loads the whole
       flow-direction raster), delineates one discarded warm-up outlet, then
       every outlet of that megabasin in shuffled order. Just before each
       block, its data files are read end to end (untimed), so they are in
       the OS file cache. Setup and per-outlet times are reported separately.

cold   One-off use. Before each run, the files that tool reads are purged
       from the Windows file cache (and the purge is checked); then a fresh
       worker delineates a single outlet. wall_since_launch_s is the time
       from starting the process to having the polygon in memory, including
       Python start-up and imports. Python and its libraries stay in the file
       cache, as they would on a machine that has used them recently.

Each round/repetition visits every (tool, outlet) once in a freshly shuffled,
seeded order, so drift (heat, background tasks) spreads over all tools.

Results go to results/<timestamp>/: results.jsonl (one record per setup or
delineation), metadata.json (machine, software, data files, settings), and
logs/ (each worker's stderr, including TauDEM's own output).

Examples:
    python run_pilot.py                              # full pilot: 5 warm rounds, 3 cold reps
    python run_pilot.py --quick                      # smoke test: 2 Iceland outlets, 1 round, 1 rep
    python run_pilot.py --scenarios cold --cold-reps 5 --tools delineator pysheds
Then: python summarize.py results/<timestamp>
"""

import argparse
import csv
import json
import random
import subprocess
import sys
import threading
from datetime import datetime
from pathlib import Path
from time import perf_counter

HERE = Path(__file__).resolve().parent
HARNESS = HERE.parent / "harness"
REPO = HERE.parent.parent
sys.path.insert(0, str(HARNESS))

import machine_info  # noqa: E402
import winperf  # noqa: E402
from tools import TOOLS  # noqa: E402

QUICK_IDS = [6401310, 6401120]


def load_config(path: Path) -> dict:
    cfg = json.loads(path.read_text())
    scratch = Path(cfg["scratch_dir"])
    cfg["scratch_dir"] = str(scratch if scratch.is_absolute() else (path.parent / scratch).resolve())
    return cfg


def load_outlets(path: Path, ids) -> list[dict]:
    rows = list(csv.DictReader(path.open(newline="", encoding="utf-8")))
    outlets = [{"id": int(r["id"]), "megabasin": int(r["megabasin"]),
                "lat": float(r["lat"]), "lng": float(r["lng"])} for r in rows]
    if ids:
        outlets = [o for o in outlets if o["id"] in set(ids)]
    return outlets


class Runner:
    def __init__(self, args, cfg: dict, out: Path):
        self.args, self.cfg, self.out = args, cfg, out
        self.config_path = out / "bench_config.resolved.json"
        self.config_path.write_text(json.dumps(cfg, indent=2))
        (out / "logs").mkdir()
        (out / "jobs").mkdir()
        Path(cfg["scratch_dir"]).mkdir(parents=True, exist_ok=True)
        self.results = (out / "results.jsonl").open("a", encoding="utf-8")
        self.tools = {name: TOOLS[name](cfg) for name in args.tools}
        self.block_count = 0

    def files(self, tool: str, megabasin: int) -> list[Path]:
        return self.tools[tool].data_files(megabasin)

    def write(self, record: dict) -> None:
        self.results.write(json.dumps(record) + "\n")
        self.results.flush()

    def run_worker(self, tool: str, megabasin: int, outlets: list[dict], warmup: bool,
                   tags: dict, timeout_s: float) -> None:
        self.block_count += 1
        block = f"{self.block_count:04d}_{tags['scenario']}_{tool}_{megabasin}"
        job = self.out / "jobs" / f"{block}.json"
        job.write_text(json.dumps(outlets))
        cmd = [sys.executable, str(HARNESS / "worker.py"), "--tool", tool,
               "--config", str(self.config_path), "--megabasin", str(megabasin), "--outlets", str(job)]
        if warmup:
            cmd.append("--warmup")
        base = {**tags, "tool": tool, "megabasin": megabasin, "block": block}
        done_at = {}
        with (self.out / "logs" / f"{block}.log").open("w", encoding="utf-8") as log:
            t0 = perf_counter()
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log, text=True,
                                    bufsize=1, cwd=self.cfg["scratch_dir"])
            timed_out = threading.Event()

            def kill():
                timed_out.set()
                proc.kill()
            timer = threading.Timer(timeout_s, kill)
            timer.start()
            try:
                for line in proc.stdout:
                    if line.startswith("@@DONE "):
                        done_at[int(line.split()[1])] = perf_counter() - t0
                    elif line.startswith("@@RESULT "):
                        record = {**base, **json.loads(line[len("@@RESULT "):])}
                        if record["kind"] == "outlet":
                            record["wall_since_launch_s"] = done_at.get(record["order"])
                            self.report(record)
                        self.write(record)
                    else:
                        log.write(line)
                proc.wait()
            finally:
                timer.cancel()
        if proc.returncode != 0 or timed_out.is_set():
            self.write({**base, "kind": "block_failed", "returncode": proc.returncode,
                        "timed_out": timed_out.is_set(), "elapsed_s": perf_counter() - t0})
            print(f"    !! {block} failed (return code {proc.returncode}"
                  f"{', timed out' if timed_out.is_set() else ''}); see logs/{block}.log")

    @staticmethod
    def report(r: dict) -> None:
        flag = " (warm-up)" if r.get("warmup") else ""
        status = "" if r["status"] == "ok" else f"  {r['status'].upper()}: {r.get('error', '')[:120]}"
        wall = f", {r['wall_since_launch_s']:.2f} s from launch" if r["scenario"] == "cold" else ""
        print(f"    {r['tool']:<10} {r['id']:>8}  {r['total_s']:8.2f} s{wall}{flag}{status}")

    @staticmethod
    def prewarm(paths) -> None:
        """Read files end to end so they are in the file cache (not timed)."""
        for p in paths:
            with open(p, "rb", buffering=0) as f:
                while f.read(16 << 20):
                    pass

    def warm(self, outlets, rng: random.Random) -> None:
        # Each block re-reads its own data files first: a previous block's
        # large allocations (pysheds peaks at ~24 GB on the Amazon) can push
        # them out of the cache, and the warm scenario promises a warm cache.
        megabasins = sorted({o["megabasin"] for o in outlets})
        for rnd in range(self.args.warm_rounds):
            blocks = [(t, mb) for t in self.tools for mb in megabasins]
            rng.shuffle(blocks)
            for tool, mb in blocks:
                batch = [o for o in outlets if o["megabasin"] == mb]
                rng.shuffle(batch)
                print(f"[warm round {rnd + 1}/{self.args.warm_rounds}] {tool}, megabasin {mb}")
                self.prewarm(self.files(tool, mb))
                self.run_worker(tool, mb, batch, warmup=True,
                                tags={"scenario": "warm", "round": rnd},
                                timeout_s=self.args.block_timeout)

    def cold(self, outlets, rng: random.Random) -> dict:
        megabasins = sorted({o["megabasin"] for o in outlets})
        candidates = sorted({p for t in self.tools for mb in megabasins for p in self.files(t, mb)},
                            key=lambda p: p.stat().st_size, reverse=True)
        calibration = winperf.calibrate(candidates[0])
        print(f"Cache purge calibration on {candidates[0].name}: warm {calibration['warm_us']:.0f} µs, "
              f"cold {calibration['cold_us']:.0f} µs per random 4 KiB read")
        if calibration["ratio"] < 3:
            print("  WARNING: cold and warm reads are hard to tell apart on this storage; "
                  "purge checks are unreliable")
        for rep in range(self.args.cold_reps):
            items = [(t, o) for t in self.tools for o in outlets]
            rng.shuffle(items)
            for i, (tool, outlet) in enumerate(items, 1):
                checks = [winperf.purge_and_verify(p, calibration["threshold_us"])
                          for p in self.files(tool, outlet["megabasin"])]
                purge_ok = all(c["ok"] for c in checks)
                if not purge_ok:
                    print(f"    WARNING: purge did not take for "
                          f"{[Path(c['file']).name for c in checks if not c['ok']]} "
                          f"(is another program holding the file open?)")
                print(f"[cold rep {rep + 1}/{self.args.cold_reps}, {i}/{len(items)}] "
                      f"{tool}, outlet {outlet['id']}")
                self.run_worker(tool, outlet["megabasin"], [outlet], warmup=False,
                                tags={"scenario": "cold", "rep": rep, "purge_ok": purge_ok,
                                      "purge_checks": checks},
                                timeout_s=self.args.cold_timeout)
        return calibration


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--outlets", type=Path, default=HERE / "pilot_outlets.csv")
    p.add_argument("--config", type=Path, default=HERE / "bench_config.json")
    p.add_argument("--tools", nargs="+", default=list(TOOLS), choices=list(TOOLS))
    p.add_argument("--scenarios", nargs="+", default=["warm", "cold"], choices=["warm", "cold"])
    p.add_argument("--ids", nargs="+", type=int, help="only these outlet ids")
    p.add_argument("--warm-rounds", type=int, default=5)
    p.add_argument("--cold-reps", type=int, default=3)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--block-timeout", type=float, default=7200, help="seconds per warm block")
    p.add_argument("--cold-timeout", type=float, default=3600, help="seconds per cold run")
    p.add_argument("--quick", action="store_true",
                   help=f"smoke test: outlets {QUICK_IDS}, 1 warm round, 1 cold rep")
    p.add_argument("--label", default="", help="appended to the results folder name")
    args = p.parse_args()
    if args.quick:
        args.ids = args.ids or QUICK_IDS
        args.warm_rounds, args.cold_reps = 1, 1
        args.label = args.label or "quick"

    cfg = load_config(args.config)
    outlets = load_outlets(args.outlets, args.ids)
    if not outlets:
        raise SystemExit("No outlets selected")
    stamp = datetime.now().strftime("%Y-%m-%d_%H%M")
    out = HERE / "results" / (stamp + (f"_{args.label}" if args.label else ""))
    out.mkdir(parents=True)
    runner = Runner(args, cfg, out)
    rng = random.Random(args.seed)

    megabasins = sorted({o["megabasin"] for o in outlets})
    all_files = sorted({p for t in runner.tools for mb in megabasins for p in runner.files(t, mb)})
    missing = [str(p) for p in all_files if not p.exists()]
    if missing:
        raise SystemExit(f"Missing data files: {missing}")

    print("Collecting machine information")
    meta = machine_info.collect(all_files + [Path(cfg["scratch_dir"])], REPO,
                                cfg.get("mpiexec", "mpiexec"), cfg.get("gagewatershed", "GageWatershed"))
    meta.update(arguments={k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
                config=cfg, outlets=outlets, data_files=machine_info.file_inventory(all_files))
    if not meta["on_ac_power"]:
        print("WARNING: not on mains power; laptops throttle on battery")
    available_gb = meta["memory_mb"]["available_mb"] / 1024
    if "pysheds" in runner.tools and 62 in megabasins and available_gb < 26:
        print(f"WARNING: {available_gb:.0f} GB of RAM free. pysheds needs about 24 GB for the "
              f"Amazon at its mouth (TauDEM about 16 GB); close other programs, or the runs will "
              f"page to disk and the timings will be slow and noisy")
    (out / "metadata.json").write_text(json.dumps(meta, indent=2, default=str))

    started = perf_counter()
    if "warm" in args.scenarios:
        runner.warm(outlets, rng)
    if "cold" in args.scenarios:
        meta["purge_calibration"] = runner.cold(outlets, rng)
    meta["elapsed_s"] = perf_counter() - started
    (out / "metadata.json").write_text(json.dumps(meta, indent=2, default=str))
    print(f"Done in {meta['elapsed_s'] / 60:.1f} min. Results in {out}")
    print(f"Summarize with:  python summarize.py \"{out}\"")
    return 0


if __name__ == "__main__":
    sys.exit(main())
