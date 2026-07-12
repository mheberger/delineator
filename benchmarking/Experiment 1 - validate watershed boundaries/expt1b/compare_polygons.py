""""
Compare watershed polygons created by `delineator` to those created by other software
Matthew Heberger, June 2026

Calculate the Coefficient of Areal Correspondence (CAC) or Intersection over Union (IOU)
or Jaccard index between the watershed polygons created by `delineator`, `pysheds` and `TauDEM`,
for the watershed polygons in GeoJSON format.
"""

import geopandas as gpd
import pandas as pd
from matplotlib import pyplot as plt

ICELAND_CRS = 3057  # ISN93 / Lambert 1993

def cac(file_a, file_b):
    a = gpd.read_file(file_a).to_crs(ICELAND_CRS)
    b = gpd.read_file(file_b).to_crs(ICELAND_CRS)

    # Repair any invalid rings/self-intersections before set operations
    a["geometry"] = a.geometry.make_valid()
    b["geometry"] = b.geometry.make_valid()

    # Dissolve each shapefile's features into a single (multi)polygon,
    # so the comparison is set-vs-set rather than feature-by-feature
    geom_a = a.union_all()
    geom_b = b.union_all()

    inter = geom_a.intersection(geom_b).area
    union = geom_a.union(geom_b).area
    cac = inter / union
    return cac, inter, union

ids = [6401070, 6401080, 6401090, 6401111, 6401120, 6401140, 6401150, 6401160, 6401200, 6401250, 6401300, 6401310,
       6401400, 6401460, 6401500, 6401601, 6401610, 6401702, 6401703, 6401800]

results = {}

folder = r'C:\Users\mheberger\Dropbox\RESEARCH\Watershed Article\benchmarking\Experiment 1 - validate watershed boundaries\expt1b'

for id in ids:
    print(f"Comparing {id}...")
    file_d = f"{folder}/delineator_out/{id}_del.geojson"
    file_p = f"{folder}/pysheds_out/{id}_py.geojson"
    file_t = f"{folder}/taudem_out/{id}_tau.geojson"

    d_p, _, _ = cac(file_d, file_p)
    d_t, _, _ = cac(file_d, file_t)
    p_t, _, _ = cac(file_p, file_t)

    results[id] = [d_p, d_t, p_t]
    print(results[id])


results_df = pd.DataFrame(results).T
results_df.columns = ["del_pys", "del_tau", "pys_tau"]
results_df.to_csv("polygons_cacs.csv")

# Make a boxplot of our results
results_df.boxplot()
plt.show()
