# Pilot: a fair speed comparison of delineator, pysheds and TauDEM

This pilot comes before re-running the speed experiments (Experiments 3 and 4).
It responds to a reviewer's comments on experimental hygiene: repetitions,
spread, cold versus warm file cache, stating the storage medium, and dropping
the polynomial fits. The pilot measures how noisy each tool's timings are on
10 outlets, so the main runs can use enough repetitions. Its summary also
previews the tables the paper will need.

The reusable parts live in [`../harness/`](../harness/), so that
Experiments 3 and 4 can be re-run the same way.

## What makes the comparison fair

1. **Same input, same watershed.** `prepare_outlets.py` snaps each outlet
   once and moves it to the exact centre of the snapped cell. All three
   tools read the same 3-arcsecond MERIT-Hydro grid (identical transform
   and size for megabasins 23, 27 and 62), so every tool gets the same cell
   with snapping off. `summarize.py` checks that the areas agree; in the
   smoke test they matched to 1 part in 10⁷. Snapping is left out of the
   timings for every tool, because each one snaps differently.
2. **Same end point.** Timing stops when a watershed polygon is in memory.
   No tool's timing includes writing output files.
3. **One-off costs reported separately.** These are Python imports and,
   for pysheds, reading the whole flow-direction raster. TauDEM starts
   `mpiexec` and reads the whole raster for every outlet, because that is
   how GageWatershed works, so that stays in its per-outlet time. TauDEM's
   own breakdown (read, compute, write) is recorded as well.
4. **Glue code is as fast as is reasonable.** The outlet shapefile is
   written with pyshp, which is lighter than geopandas. TauDEM's output
   raster is vectorized only within the watershed's bounding window, and
   GDAL decodes its LZW strips on all cores. In Experiment 3, polygonizing
   the full extent took 20–25 s per outlet in megabasins 23 and 62.
   pysheds reads its raster
   once: `Grid.from_raster(path)` followed by `read_raster(path)` would read
   it twice.
5. **Two scenarios, kept apart.**
   - *warm* (batch use): data files already in the OS file cache. Each
     block re-reads its files first, untimed, because a previous block's
     large allocations push them out of the cache. One worker process per
     tool and megabasin does setup, one discarded warm-up delineation, then
     every outlet in shuffled order.
   - *cold* (one-off use): before each run, the files that tool reads are
     purged from the Windows file cache, and the purge is checked. Then a
     fresh process delineates one outlet. This is timed from process
     launch to polygon, and also shown without Python start-up and
     imports.
6. **Repetition and order.** Every round (warm) or repetition (cold)
   visits each tool/outlet pair once, in an order shuffled with a fixed
   seed (`--seed`). Slow drift from heat or background tasks is then spread
   across all tools instead of hitting one.
7. **Measured the same way for every tool.** Timings use `time.perf_counter()`
   at full precision. Each worker puts itself in a Windows *job object*, so
   its CPU time, bytes read and peak memory include every process it
   starts, the TauDEM MPI ranks among them. delineator's internal phases
   are timed too (locate, upstream query, raster split, union). These show
   where its time goes, and whether it is waiting on disk or computing.
8. **Conditions recorded automatically.** `metadata.json` in each results
   folder records the CPU, RAM, the disk behind each data path (model, bus,
   SSD/HDD), the power plan and mode, and the antivirus state (real-time
   protection and, when run as administrator, which data folders were
   excluded). It also has every package version,
   the TauDEM and MS-MPI versions, the delineator commit, and the size and
   date of each data file. It leaves out the host name and user name.

### Purging the file cache without admin rights

Opening a file with `FILE_FLAG_NO_BUFFERING` makes Windows flush and drop
that file's pages from its cache (`harness/winperf.py`). Before each cold run,
the harness purges the tool's data files and samples 32 random 4 KiB reads.
It then purges again, and checks that those reads ran at disk speed rather
than RAM speed. A calibration at the start measures both speeds. On this
laptop's NVMe drive they are about 6 µs (cached) and 150 µs (not cached).
The purge fails if another program holds a file open (QGIS, a running
`delineator_serve`, a sync client), and the summary then reports
`purge_ok = False`.

Only the data files are purged. Python and its libraries stay cached, as
they would on any machine that has used them recently.

## Files

| File | Purpose |
|---|---|
| `prepare_outlets.py` | Picks the 10 pilot outlets from Experiment 3 and pre-snaps them; writes `pilot_outlets.csv` (already done) |
| `pilot_outlets.csv` | The outlets: requested and snapped coordinates, grid cell, upstream cell count, area |
| `bench_config.json` | Data paths, number of MPI processes for TauDEM, scratch folder |
| `run_pilot.py` | Runs the scenarios; writes `results/<timestamp>/` |
| `summarize.py` | Tables from a results folder: `summary.md` plus CSVs |
| `../harness/tools.py` | The three tools behind one interface |
| `../harness/worker.py` | The per-block worker process |
| `../harness/winperf.py` | Job-object accounting, cache purge, power and memory status |
| `../harness/machine_info.py` | Machine and software metadata |

Pilot outlets: two in Iceland (27), three in Western Europe (23) and five in
the Amazon basin (62), from about 100 km² to 5.9 million km² (the Amazon at
its mouth).

## Running it

```powershell
cd "benchmarking\Pilot - fair speed comparison"
..\..\.venv\Scripts\python.exe run_pilot.py --quick      # 2-3 min smoke test on two Iceland outlets
..\..\.venv\Scripts\python.exe run_pilot.py              # the pilot: 5 warm rounds, 3 cold repetitions
..\..\.venv\Scripts\python.exe summarize.py "results\<timestamp>"
```

Useful options: `--tools`, `--scenarios warm|cold`, `--ids`,
`--warm-rounds`, `--cold-reps`, `--seed`, and `--label`.

The full pilot should take about 3–3.5 hours: roughly 25 minutes per warm
round and 25 minutes per cold repetition. The Amazon at its mouth dominates:
on the smoke test it took about 17 s with delineator, 155 s with pysheds and
167 s with TauDEM.

### Before a run whose numbers you will use

- **Free up memory.** At the Amazon's mouth, pysheds peaks at about 24 GB of
  committed memory and TauDEM at about 16 GB. With other programs open, only
  18–20 GB of this laptop's 32 GB was free, so those runs page to disk.
  `run_pilot.py` warns when less than 26 GB is free. Each worker's setup
  record stores the free memory at its start.
- Plug in the laptop and set the power mode to *Best performance*. It was
  already set that way when this was written; `metadata.json` records it.
- Close QGIS, browsers, `delineator_serve`, and anything else that reads the
  data files or uses a lot of CPU. Pause OneDrive or Dropbox syncing.
- Don't use the laptop during the run, and stop it sleeping (for example,
  set *Sleep when plugged in* to *Never*).
- Microsoft Defender scans files as they are opened, which slows cold
  reads. Either add exclusions for the data folders and `_scratch`, or leave
  them out, but say which in the paper. `metadata.json` (under
  `antivirus`) always records whether real-time protection was on. Windows
  shows the exclusion list only to administrators, so run `run_pilot.py`
  from an administrator PowerShell if you want it to record which data
  folders were excluded; otherwise note it yourself, e.g. with
  `--label defender_excluded`.

## Output

`results/<timestamp>/`:

- `results.jsonl`: one record per worker setup (`kind = setup`), per
  delineation (`kind = outlet`), and per worker exit (`kind = end`, with
  peak memory). Outlet records hold `total_s`, `wall_since_launch_s` (cold
  runs), the phases, CPU time, bytes read and written, area and vertex
  count.
- `metadata.json`: machine, software, data files, settings, and the purge
  calibration.
- `logs/`: each worker's stderr, including TauDEM's output.
- After `summarize.py`: `summary.md`, `timings.csv`, `setup.csv`,
  `speedups.csv`, `areas.csv`, `phases.csv`, `resources.csv`,
  `peak_memory.csv`.

`summarize.py` ends with the number of repetitions each tool needs, in each
scenario, for the standard error of every outlet's median to be within 2%.
The noisiest outlet sets it.

## Limitations

- Windows only: the job-object accounting and the cache purge use the
  Windows API. On Linux, `posix_fadvise(POSIX_FADV_DONTNEED)` and
  `/proc/<pid>/io` would do the same jobs.
- Bytes read are *logical* reads, cache hits included. Reads through
  memory-mapped files are not counted; none of the tools use them for these
  files.
- TauDEM's CPU time is inflated, because idle MPI ranks spin while they
  wait.
