"""Demo: hybrid GRIT watershed delineation with `pourpoint` + `delineator`.

`pourpoint` (pip install pourpoint) reads HFX datasets compiled from the GRIT
river network and handles the *vector* side of delineation: it resolves an
outlet to a GRIT drainage unit, traverses everything upstream, and returns the
whole ("pre-merge", unrefined) unit polygons. Its own docs warn that unioning
those whole units is not a refined watershed: the terminal unit is included in
full, even the part downstream of the outlet.

That is exactly the gap delineator's raster splitter fills. This demo:

1. asks pourpoint for the unrefined upstream unit polygons (streamed from the
   hosted public GRIT HFX dataset -- nothing to download in advance),
2. splits the terminal unit at the outlet with delineator.grit_detailed,
   using the local GRIT flow direction + drainage area VRT mosaics
   (see data_mgmt/build_grit_vrt.py),
3. splices the split into the union of the upstream units, and
4. writes a GeoPackage comparing the hybrid watershed against pourpoint's
   own merged result and the raw unrefined union.

The default outlet is the French Broad River at Asheville, NC (~2,450 km2).
The default raster mosaics cover North America; pass --fdir/--accum to use
others. Note the HFX unit polygons are simplified (~108 m median segment
length, smoothed boundaries) while the raster split follows exact 30 m pixel
edges, so the two never edge-match; the assembly step repairs the resulting
seam slivers.

Usage:
    python demo_grit_pourpoint.py [--point LAT LNG] [--output demo.gpkg]
"""
import argparse
import time

import geopandas as gpd
import pourpoint
from pyproj import Geod
from shapely import wkb
from shapely.geometry import Point
from shapely.ops import unary_union

from delineator.grit_detailed import _split_catchment
from delineator.settings import DelineatorConfig
from delineator.util import _close_holes

HFX_URL = "https://basin-delineations-public.upstream.tech/grit/hfx-v0.3.0/"
GEOD = Geod(ellps="WGS84")


def geod_area_km2(geom) -> float:
    return abs(GEOD.geometry_area_perimeter(geom)[0]) / 1e6


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--point", nargs=2, type=float, default=[35.61, -82.58],
                        metavar=("LAT", "LNG"),
                        help="outlet coordinates (default: French Broad R. at Asheville)")
    parser.add_argument("--dataset", default=HFX_URL, help="HFX dataset path or URL")
    parser.add_argument("--fdir", default=r"C:\Data\GIS\GRIT\fdir\grit_fdir_na.vrt",
                        help="GRIT drainage direction raster (VRT mosaic)")
    parser.add_argument("--accum", default=r"C:\Data\GIS\GRIT\accum_mainstem\grit_accum_mainstem_na.vrt",
                        help="GRIT drainage area raster (VRT mosaic)")
    parser.add_argument("--output", default="grit_pourpoint_demo.gpkg",
                        help="output GeoPackage path")
    args = parser.parse_args()
    lat, lng = args.point

    # ------------------------------------------------------------------
    # 1. pourpoint: unrefined vector units, streamed from the HFX dataset
    # ------------------------------------------------------------------
    print(f"Opening HFX dataset: {args.dataset}")
    t0 = time.time()
    engine = pourpoint.Engine(args.dataset)
    print(f"  engine ready in {time.time() - t0:.1f} s")

    t0 = time.time()
    level = engine.select_level(pourpoint.LevelSelection.FINEST)
    outlet = engine.resolve_outlet(level, lat=lat, lon=lng)
    upstream = engine.traverse(outlet)
    units = engine.pre_merge_units(upstream)
    print(f"  resolved + traversed {len(units.units)} upstream units "
          f"in {time.time() - t0:.1f} s (terminal unit {units.terminal_unit_id})")

    geoms = [wkb.loads(bytes(b)) for b in units.unit_geometry_wkb]
    recs = units.units

    terminal_idx = next(i for i, r in enumerate(recs) if r.id == units.terminal_unit_id)
    terminal_rec = recs[terminal_idx]
    terminal_poly = geoms[terminal_idx]
    print(f"  terminal unit area {terminal_rec.area_km2:.2f} km2, "
          f"upstream area {terminal_rec.up_area_km2:.0f} km2")

    # ------------------------------------------------------------------
    # 2. delineator: raster split of the terminal unit at the outlet
    # ------------------------------------------------------------------
    print("Splitting terminal unit with GRIT rasters (delineator.grit_detailed)")
    t0 = time.time()
    config = DelineatorConfig()
    split_poly, lat_snap, lng_snap = _split_catchment(
        lat, lng, terminal_poly, terminal_rec.area_km2, config,
        fdir_path=args.fdir, accum_path=args.accum,
    )
    if split_poly is None:
        raise SystemExit("Raster split failed; see log output above.")
    print(f"  split done in {time.time() - t0:.1f} s; outlet snapped to "
          f"({lat_snap:.5f}, {lng_snap:.5f})")

    # ------------------------------------------------------------------
    # 3. Assemble the hybrid watershed: upstream units + split terminal
    # ------------------------------------------------------------------
    upstream_geoms = [g for i, g in enumerate(geoms) if i != terminal_idx]
    hybrid = unary_union(upstream_geoms + [split_poly])
    if not hybrid.is_valid:
        hybrid = hybrid.buffer(0)

    # Seam hygiene. Tiny slivers open up along unit boundaries because the
    # split polygon follows exact 30 m pixel edges while the HFX unit
    # polygons are simplified (~108 m segments, no stair-step), so their
    # edges never match. Fill pinhole interior rings and drop crumb parts
    # well below one pixel.
    hybrid = _close_holes(hybrid, area_max_km2=0.01)
    if hybrid.geom_type == "MultiPolygon":
        parts = [p for p in hybrid.geoms if geod_area_km2(p) > 0.01]
        hybrid = unary_union(parts)

    # ------------------------------------------------------------------
    # 4. Compare against pourpoint's own merged result
    # ------------------------------------------------------------------
    result = engine.delineate(lat=lat, lon=lng)
    pourpoint_poly = wkb.loads(bytes(result.geometry_wkb))

    unrefined = unary_union(geoms)
    area_hybrid = geod_area_km2(hybrid)
    area_pourpoint = geod_area_km2(pourpoint_poly)
    area_unrefined = geod_area_km2(unrefined)
    print("\nWatershed areas:")
    print(f"  unrefined union of whole units:  {area_unrefined:10.2f} km2")
    print(f"  pourpoint merged result:         {area_pourpoint:10.2f} km2")
    print(f"  hybrid (delineator raster split):{area_hybrid:10.2f} km2")
    print(f"  terminal unit kept by the split: "
          f"{terminal_rec.area_km2 - (area_unrefined - area_hybrid):.2f} "
          f"of {terminal_rec.area_km2:.2f} km2")

    # ------------------------------------------------------------------
    # 5. Write everything to a GeoPackage for inspection
    # ------------------------------------------------------------------
    crs = "EPSG:4326"
    gpd.GeoDataFrame(
        {"id": [r.id for r in recs],
         "area_km2": [r.area_km2 for r in recs],
         "up_area_km2": [r.up_area_km2 for r in recs],
         "is_terminal": [i == terminal_idx for i in range(len(recs))]},
        geometry=geoms, crs=crs,
    ).to_file(args.output, layer="units_unrefined", driver="GPKG")
    gpd.GeoDataFrame({"source": ["hybrid"]}, geometry=[hybrid], crs=crs).to_file(
        args.output, layer="watershed_hybrid", driver="GPKG")
    gpd.GeoDataFrame({"source": ["pourpoint"]}, geometry=[pourpoint_poly], crs=crs).to_file(
        args.output, layer="watershed_pourpoint", driver="GPKG")
    gpd.GeoDataFrame(
        {"which": ["input", "pourpoint_resolved", "grit_snapped"]},
        geometry=[Point(lng, lat),
                  Point(*outlet.resolved_outlet),
                  Point(lng_snap, lat_snap)],
        crs=crs,
    ).to_file(args.output, layer="outlets", driver="GPKG")
    print(f"\nWrote {args.output} (layers: units_unrefined, watershed_hybrid, "
          f"watershed_pourpoint, outlets)")


if __name__ == "__main__":
    main()
