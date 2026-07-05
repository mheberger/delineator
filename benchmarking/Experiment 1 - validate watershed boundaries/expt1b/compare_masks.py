"""
Compare two raster masks using the intersection over union (IOU) metric.
Matthew Heberger, June 2026

This was part of `delineator` validation Experiment 0.
Since raster methods are deterministic when using the D-8 algorithm,
different methods should produce the same output.

So this was a basic check to see if my `delineator` package produces the same
output as `pysheds` and `TauDEM`.

"""

import numpy as np
import pandas as pd
import rasterio
from matplotlib import pyplot as plt


def load_bool(path):
    """
    Load a raster mask GeoTiff
    Returns a numpy array (TODO: check) of booleans where True is in the watershed
    and False is the background.
    """
    with rasterio.open(path) as src:
        data = src.read(1)
        # NaN == 1 is False, so this cleanly maps {1 -> True, NaN/other -> False}
        return (data == 1), src.transform, src.crs, src.res


def mask_iou(path_a: str, path_b: str) -> tuple[float, int, int]:
    """
    Calculate the Intersection over Union (IOU) for two raster masks
    
    This is the equivalent to the Coefficient of Areal Correspondance (CAC)
    but here we are using the number of pixels rather than areas.
    The masks used here are GeoTiff files containing values of 0 and 1, 
    where 1 represents a pixel in the watershed and 0 is the background.
    
    Args:
        path_a: path to a GeoTiff file 
        path_b: path to second GeoTiff
    
    Returns:
        iou: the intersection over union (floating point from 0 - 1.0)
        union: the number of pixels in the Union of masks A and B
        inter: the number of pixels in the Intersection of masks A and B
    """
    
    # Load the mask data
    a, ta, crs_a, res_a = load_bool(path_a)
    b, tb, crs_b, res_b = load_bool(path_b)

    # To do the comparison, the grids must have the same CRS and resolution
    # although they may have different extents
    assert crs_a == crs_b, "CRS mismatch"
    assert res_a == res_b, "resolution mismatch"
    
    # The following code maps the two grids onto the same shape array
    # so they can be directly compared.
    px, py = res_a  # pixel width, height (both positive)

    # Union footprint (assumes north-up, no rotation)
    left = min(ta.c, tb.c)
    top = max(ta.f, tb.f)
    right = max(ta.c + a.shape[1] * px, tb.c + b.shape[1] * px)
    bottom = min(ta.f - a.shape[0] * py, tb.f - b.shape[0] * py)

    # Round to kill floating point noise
    width = int(round((right - left) / px))
    height = int(round((top - bottom) / py))

    A = np.zeros((height, width), dtype=bool)
    B = np.zeros((height, width), dtype=bool)

    # Integer offsets into the union grid 
    for arr, t, dst in ((a, ta, A), (b, tb, B)):
        col = int(round((t.c - left) / px))
        row = int(round((top - t.f) / py))
        dst[row:row + arr.shape[0], col:col + arr.shape[1]] = arr

    inter = np.logical_and(A, B).sum()
    union = np.logical_or(A, B).sum()
    iou = inter / union if union else float("nan")
    return iou, int(inter), int(union)


# List of the watershed IDs (sample of 20 GRDC gages in Iceland)
ids = [6401090, 6401111, 6401120, 6401140, 6401150, 6401160, 6401200, 6401250, 6401300, 6401310, 6401400, 6401460,
       6401500, 6401601, 6401610, 6401702, 6401703, 6401800, 6401080, 6401070]

ids.sort()

results = {}

folder = r'C:\Users\mheberger\Dropbox\RESEARCH\Watershed Article\benchmarking\Experiment 1 - validate watershed boundaries\expt1b'

for id in ids:
    print(f"Comparing {id}...")
    raster_d = f"{folder}/delineator_out/{id}_del.tif"
    raster_p = f"{folder}/pysheds_out/{id}_py.tif"
    raster_t = f"{folder}/taudem_out/{id}_tau.tif"
    d_p, _, _ = mask_iou(raster_d, raster_p)
    d_t, _, _ = mask_iou(raster_d, raster_t)
    p_t, _, _ = mask_iou(raster_p, raster_t)

    results[id] = [d_p, d_t, p_t]
    print(results[id])

results_df = pd.DataFrame(results).T
results_df.columns = ["d_p", "d_t", "p_t"]
results_df.to_csv("mask_ious.csv")

# Make a boxplot of our results
results_df.boxplot()
plt.show()
