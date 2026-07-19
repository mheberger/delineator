"""
Benchmark for DelineatorConfig(cache=True).

Delineates a ring of nearby outlets around a point on the Amazon mainstem
near Manaus, Brazil, first without and then with caching. All the outlets
fall in the same unit catchment, so with cache=True the first call pays for
reading and dissolving the upstream geometries and the rest reuse it.

Rivers and outlets are turned off so the timings isolate the part of the
pipeline the cache affects. Note that the cached calls are not free: the
raster split of the home catchment and the upstream comid query still run
on every call.
"""
import time

from delineator import DelineatorConfig, delineate, clear_cache

# Centered on the COMID for the outlet of Rio Negro near Manaus
lat, lng = -3.157, -59.942

offsets = (
    (1, 0),
    (1, 1),
    (0, 1),
    (-1, 1),
    (-1, 0),
    (-1, -1),
    (0, -1),
    (1, -1),
)


def run(use_cache: bool) -> None:
    config = DelineatorConfig(cache=use_cache, rivers=False, outlets=False)
    print(f"\n--- cache={use_cache} ---")
    for offset in offsets:
        mylat = lat + offset[0] * 0.001
        mylng = lng + offset[1] * 0.001
        start = time.time()
        watershed, _, _ = delineate(mylat, mylng, config)
        status = "ok" if watershed is not None else "FAILED"
        print(f"({mylat:.3f}, {mylng:.3f}) {status}: {time.time() - start:.2f} s")


run(use_cache=False)
clear_cache()
run(use_cache=True)
