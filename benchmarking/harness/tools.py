"""
The three delineation tools behind one interface, so the harness times all of
them over the same span: from an outlet coordinate to a watershed polygon in
memory. Writing output files is not timed for any tool.

Each tool has:
    import_libraries()   timed separately as "import"
    setup(megabasin)     timed separately as "setup" (one-off cost per batch)
    delineate(lat, lng)  timed per outlet; returns (polygon, {phase: seconds})
    data_files(megabasin) the files it reads, so the harness can purge them

Outlets are pre-snapped to the centre of a stream cell (see
prepare_outlets.py), and every tool delineates from that cell without
snapping. All three read the same 3-arcsecond MERIT-Hydro grid, so they
delineate the same watershed. That is checked afterwards by comparing areas.

Heavy libraries are imported only in import_libraries(), so the controller
can import this module cheaply.
"""

import os
import re
import shutil
import subprocess
from collections import defaultdict
from importlib.util import find_spec
from pathlib import Path
from time import perf_counter

# ESRI/MERIT D8 codes in pysheds order (N, NE, E, SE, S, SW, W, NW); same as delineator.constants.DIRMAP
DIRMAP = (64, 128, 1, 2, 4, 8, 16, 32)

WGS84_PRJ = ('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],'
             'PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')


def geometry_stats(geom) -> dict:
    """Area on the WGS84 ellipsoid and vertex count, computed the same way for every tool."""
    import shapely
    from pyproj import Geod
    area_m2, _ = Geod(ellps="WGS84").geometry_area_perimeter(geom)
    return {"area_km2": abs(area_m2) / 1e6,
            "n_vertices": int(shapely.get_num_coordinates(geom)),
            "n_parts": len(getattr(geom, "geoms", [geom]))}


class DelineatorTool:
    """delineator's public delineate(), with snapping off and all post-processing off."""

    name = "delineator"

    # Functions in delineator.core whose run time is recorded as a phase.
    # Wrapping costs about a microsecond per call. Names missing from a
    # future version are skipped; their time then shows up under "other".
    PHASES = {
        "_point_in_polygon_analysis": "locate",
        "_get_upstream_catchments": "upstream_query",
        "get_upstream_area": "upstream_area",
        "get_home_unit_catchment_geom_and_area": "home_geometry",
        "_split_catchment": "raster_split",
        "_get_upstream_union": "upstream_union",
    }

    def __init__(self, cfg: dict):
        self.data_dir = Path(cfg["delineator_data_dir"])
        self._phase_times = defaultdict(float)

    def data_files(self, megabasin: int) -> list[Path]:
        package_dir = Path(find_spec("delineator").submodule_search_locations[0])
        return [package_dir / "data" / "megabasins.db"] + [
            self.data_dir / f"{stem}{megabasin}.{ext}"
            for stem, ext in (("basins", "db"), ("rivers", "db"), ("flowdir", "tif"), ("accum", "tif"))]

    def import_libraries(self):
        import delineator.core as core
        from delineator import DelineatorConfig, delineate
        self._delineate = delineate
        self.config = DelineatorConfig(
            rivers=False, outlets=False, simplify=False, fill=False, clean=False, smooth=False,
            calc_area=False, snapping=False, cache=False, auto_download=False,
            data_dir=str(self.data_dir))
        for func_name, phase in self.PHASES.items():
            func = getattr(core, func_name, None)
            if func is not None:
                setattr(core, func_name, self._timed(func, phase))

    def _timed(self, func, phase):
        def wrapper(*args, **kwargs):
            t = perf_counter()
            try:
                return func(*args, **kwargs)
            finally:
                self._phase_times[phase] += perf_counter() - t
        return wrapper

    def setup(self, megabasin: int) -> dict:
        return {}

    def delineate(self, lat: float, lng: float):
        self._phase_times.clear()
        t = perf_counter()
        watershed, _, _ = self._delineate(lat, lng, self.config)
        total = perf_counter() - t
        if watershed is None:
            raise RuntimeError("delineate() returned None; see the worker log for its warning")
        phases = dict(self._phase_times)
        phases["other"] = total - sum(phases.values())
        return watershed.geometry.iloc[0], phases


class PyshedsTool:
    """
    pysheds on the whole megabasin flow-direction raster, the way a pysheds
    user would do it: load the raster once, then delineate each outlet.
    The calls mirror delineator's own use of pysheds (merit_detailed.py).
    """

    name = "pysheds"

    def __init__(self, cfg: dict):
        self.fdir_pattern = cfg["pysheds_fdir"]

    def data_files(self, megabasin: int) -> list[Path]:
        return [Path(self.fdir_pattern.format(mb=megabasin))]

    def import_libraries(self):
        import numpy as np
        from pysheds.grid import Grid
        from shapely.geometry import MultiPolygon, shape
        self._np, self._Grid, self._MultiPolygon, self._shape = np, Grid, MultiPolygon, shape

    def setup(self, megabasin: int) -> dict:
        path = str(self.data_files(megabasin)[0])
        t = perf_counter()
        # Grid.from_raster(path) would read the whole raster just to define the
        # grid, and read_raster() would then read it again. Read it once.
        self.fdir = self._Grid().read_raster(path, nodata=0)
        self.grid = self._Grid.from_raster(self.fdir)
        return {"read_fdir": perf_counter() - t}

    def delineate(self, lat: float, lng: float):
        grid, fdir = self.grid, self.fdir
        t0 = perf_counter()
        try:
            catch = grid.catchment(x=lng, y=lat, fdir=fdir, dirmap=DIRMAP,
                                   xytype="coordinate", snap="center")
            t1 = perf_counter()
            grid.clip_to(catch)
            view = grid.view(catch, dtype=self._np.uint8)
            t2 = perf_counter()
            polygons = [self._shape(geom) for geom, _ in grid.polygonize(view)]
            geom = self._MultiPolygon(polygons)
        finally:
            grid.viewfinder = fdir.viewfinder  # undo clip_to for the next outlet
        t3 = perf_counter()
        return geom, {"catchment": t1 - t0, "clip_view": t2 - t1, "polygonize": t3 - t2}


class TaudemTool:
    """
    TauDEM GageWatershed run through mpiexec, then the output raster
    vectorized. Per outlet: write a one-point shapefile (TauDEM's input
    format), run GageWatershed, and vectorize only the watershed's bounding
    window. Polygonizing the whole megabasin extent instead took about
    20 s per outlet in Experiment 3; this keeps TauDEM's vectorizing as
    fast as is reasonable.
    """

    name = "taudem"
    OUTLET_ID = 1
    _TIME_LINE = re.compile(r"^\s*([A-Za-z][A-Za-z ]*?)\s+time\s*[:=]?\s*([-+0-9.eE]+)\s*$", re.I)

    def __init__(self, cfg: dict):
        self.fdir_pattern = cfg["taudem_fdir"]
        self.mpiexec = cfg.get("mpiexec", "mpiexec")
        self.gagewatershed = cfg.get("gagewatershed", "GageWatershed")
        self.processes = int(cfg.get("taudem_processes", 4))
        self.workdir = Path(cfg["scratch_dir"]) / f"taudem_{os.getpid()}"

    def data_files(self, megabasin: int) -> list[Path]:
        return [Path(self.fdir_pattern.format(mb=megabasin))]

    def import_libraries(self):
        import numpy as np
        import rasterio
        import shapefile  # pyshp: writes the one-point shapefile without importing geopandas
        from rasterio.features import shapes
        from rasterio.windows import Window
        from shapely.geometry import MultiPolygon, shape
        self._np, self._rasterio, self._shapefile = np, rasterio, shapefile
        self._shapes, self._Window, self._MultiPolygon, self._shape = shapes, Window, MultiPolygon, shape

    def setup(self, megabasin: int) -> dict:
        self.fdir = str(self.data_files(megabasin)[0])
        self.workdir.mkdir(parents=True, exist_ok=True)
        return {}

    def delineate(self, lat: float, lng: float):
        outlet_shp = self.workdir / "outlet.shp"
        watershed_tif = self.workdir / "watershed.tif"
        t0 = perf_counter()
        with self._shapefile.Writer(str(outlet_shp), shapeType=self._shapefile.POINT) as w:
            w.field("id", "N")
            w.point(lng, lat)
            w.record(self.OUTLET_ID)
        outlet_shp.with_suffix(".prj").write_text(WGS84_PRJ)
        t1 = perf_counter()
        proc = subprocess.run(
            [self.mpiexec, "-n", str(self.processes), self.gagewatershed,
             "-p", self.fdir, "-o", str(outlet_shp), "-gw", str(watershed_tif)],
            capture_output=True, text=True)
        t2 = perf_counter()
        print(proc.stdout.strip())  # the worker sends stdout to its log
        version = re.search(r"version\s+([\d.]+)", proc.stdout)
        if version:
            self.version = version.group(1)
        if proc.returncode != 0 or not watershed_tif.exists():
            raise RuntimeError(f"GageWatershed failed ({proc.returncode}): {proc.stderr.strip()[-400:]}")
        geom = self._vectorize(watershed_tif)
        t3 = perf_counter()
        phases = {"write_outlet": t1 - t0, "gagewatershed": t2 - t1, "vectorize": t3 - t2}
        # TauDEM's own breakdown of the GageWatershed run (read, compute, write, ...)
        for line in proc.stdout.splitlines():
            match = self._TIME_LINE.match(line)
            if match:
                phases["taudem_" + match.group(1).strip().lower().replace(" ", "_")] = float(match.group(2))
        self._cleanup()
        return geom, phases

    def _vectorize(self, path: Path):
        # GageWatershed writes a full-extent int32 LZW GeoTIFF in one-row
        # strips. Decompressing it dominates this step, so let GDAL decode
        # strips on all cores (about 3x faster than one thread) and read in
        # chunks of about 512 MB.
        np = self._np
        with self._rasterio.Env(GDAL_NUM_THREADS="ALL_CPUS"), self._rasterio.open(path) as src:
            rows_hit, cols_hit = [], np.zeros(src.width, dtype=bool)
            step = max(256, (512 << 20) // (src.width * 4))
            for r0 in range(0, src.height, step):
                h = min(step, src.height - r0)
                inside = src.read(1, window=self._Window(0, r0, src.width, h)) == self.OUTLET_ID
                hit = inside.any(axis=1)
                if hit.any():
                    rows_hit.append(r0 + np.flatnonzero(hit))
                    cols_hit |= inside.any(axis=0)
            if not rows_hit:
                raise RuntimeError("GageWatershed output contains no watershed cells")
            rows, cols = np.concatenate(rows_hit), np.flatnonzero(cols_hit)
            window = self._Window(cols[0], rows[0], cols[-1] - cols[0] + 1, rows[-1] - rows[0] + 1)
            data = src.read(1, window=window)
            transform = src.window_transform(window)
        inside = data == self.OUTLET_ID
        polygons = [self._shape(geom) for geom, _ in
                    self._shapes(inside.astype(np.uint8), mask=inside, transform=transform)]
        return self._MultiPolygon(polygons)

    def _cleanup(self):
        for p in self.workdir.glob("*"):
            try:
                p.unlink()
            except OSError:
                pass

    def teardown(self):
        shutil.rmtree(self.workdir, ignore_errors=True)


TOOLS = {cls.name: cls for cls in (DelineatorTool, PyshedsTool, TaudemTool)}
