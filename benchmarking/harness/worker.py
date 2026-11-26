"""
Benchmark worker: one process that imports one tool, sets it up for one
megabasin, and delineates a list of outlets. The controller starts a fresh
worker for every block (warm batch) or every single outlet (cold run).

Protocol on stdout (everything else the worker or the tools print goes to
stderr, which the controller saves as the block's log):
    @@DONE <order>      sent the moment a delineation finishes, so the
                        controller can time a cold run end to end
    @@RESULT <json>     one record: kind = setup | outlet | end

Usage (normally run by the controller):
    python worker.py --tool pysheds --config bench_config.json --megabasin 62 \
        --outlets job.json [--warmup]
"""

import time

T_FIRST_LINE = time.perf_counter()

import argparse  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import traceback  # noqa: E402
from pathlib import Path  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import winperf  # noqa: E402
from tools import TOOLS, geometry_stats  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool", required=True, choices=sorted(TOOLS))
    parser.add_argument("--config", required=True, help="bench_config.json")
    parser.add_argument("--megabasin", type=int, required=True)
    parser.add_argument("--outlets", required=True, help="JSON list of {id, lat, lng}")
    parser.add_argument("--warmup", action="store_true",
                        help="delineate the first outlet once more before the timed runs, and flag it")
    args = parser.parse_args()

    protocol = sys.stdout
    sys.stdout = sys.stderr  # keep tool output (e.g. delineator's prints) off the protocol channel

    def emit(tag: str, payload) -> None:
        protocol.write(f"{tag} {payload}\n")
        protocol.flush()

    acct = winperf.JobAccounting()
    cfg = json.loads(Path(args.config).read_text())
    outlets = json.loads(Path(args.outlets).read_text())
    tool = TOOLS[args.tool](cfg)

    before = acct.snapshot()
    t0 = time.perf_counter()
    try:
        tool.import_libraries()
        t1 = time.perf_counter()
        setup_phases = tool.setup(args.megabasin)
        t2 = time.perf_counter()
    except Exception as e:
        traceback.print_exc()
        emit("@@RESULT", json.dumps({"kind": "setup", "status": "error", "error": repr(e)[:500]}))
        return 1
    emit("@@RESULT", json.dumps({
        "kind": "setup", "status": "ok",
        "stdlib_import_s": t0 - T_FIRST_LINE, "import_s": t1 - t0, "setup_s": t2 - t1,
        "setup_phases": setup_phases, "covers_children": acct.covers_children,
        "available_memory_mb": winperf.memory_status_mb()["available_mb"],
        **acct.delta(before, acct.snapshot())}))

    plan = ([dict(outlets[0], warmup=True)] if args.warmup else []) + \
           [dict(o, warmup=False) for o in outlets]
    for order, outlet in enumerate(plan):
        record = {"kind": "outlet", "id": outlet["id"], "warmup": outlet["warmup"], "order": order}
        geom, phases = None, {}
        before = acct.snapshot()
        t = time.perf_counter()
        try:
            geom, phases = tool.delineate(float(outlet["lat"]), float(outlet["lng"]))
            record["status"] = "ok"
        except Exception as e:
            traceback.print_exc()
            record.update(status="error", error=repr(e)[:500])
        record["total_s"] = time.perf_counter() - t
        emit("@@DONE", order)
        record.update(acct.delta(before, acct.snapshot()))
        record["phases"] = phases
        if geom is not None:
            record.update(geometry_stats(geom))
        emit("@@RESULT", json.dumps(record))
        del geom

    emit("@@RESULT", json.dumps({"kind": "end", "peak_memory_mb": acct.peak_memory_mb(),
                                 "tool_version": getattr(tool, "version", None)}))
    if hasattr(tool, "teardown"):
        tool.teardown()
    return 0


if __name__ == "__main__":
    sys.exit(main())
