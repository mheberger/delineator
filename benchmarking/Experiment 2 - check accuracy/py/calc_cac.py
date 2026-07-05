"""
Calculate the Coefficient of Areal Correspondence or Jaccard Index
for the validation set of 120 watersheds.

"""

from shapely.geometry.base import BaseGeometry
from shapely.validation import make_valid
import geopandas as gpd



def areal_correspondence(a: BaseGeometry, b: BaseGeometry, fix_invalid: bool = True) -> float:
    """Coefficient of Areal Correspondence (Jaccard index) for two polygons.

    Returns the ratio of the intersection area to the union area, a value
    in [0, 1] where 1 means the polygons are coincident and 0 means they
    do not overlap at all.

    Parameters
    ----------
    a, b : shapely geometry
        Polygon or MultiPolygon geometries (in the same CRS / projected
        units if the result is to be meaningful).
    fix_invalid : bool
        If True, repair invalid geometries with make_valid() before
        computing areas.

    Returns
    -------
    float
        CAC / Jaccard index. Returns 0.0 if the union has zero area.
    """
    if fix_invalid:
        if not a.is_valid:
            a = make_valid(a)
        if not b.is_valid:
            b = make_valid(b)

    inter = a.intersection(b).area
    union = a.union(b).area

    return inter / union if union > 0 else 0.0


def compare_shapefiles(path_a, path_b, equal_area_crs="EPSG:5070"):
    """Read two shapefiles and compute the CAC between their dissolved extents.

    equal_area_crs defaults to EPSG:5070 (CONUS Albers). Pick the right
    equal-area projection for your region if you're not working in CONUS.
    """
    gdf_a = gpd.read_file(path_a)
    gdf_b = gpd.read_file(path_b)

    # Reproject both to a common equal-area CRS so .area is meaningful
    gdf_a = gdf_a.to_crs(equal_area_crs)
    gdf_b = gdf_b.to_crs(equal_area_crs)

    # Collapse each layer to a single (multi)polygon
    geom_a = gdf_a.union_all()   # use .unary_union on GeoPandas < 1.0
    geom_b = gdf_b.union_all()

    return areal_correspondence(geom_a, geom_b)


if __name__ == "__main__":
    gages = ["01015800","01092000","01094000","01170500","01209901","01318500","01403060","01431500","01464907","01480870","01488500","01520500","01555500","01567000","02070500","02110500","02131000","02225500","02226500","02299950","02365200","02374500","03050000","03072655","03086000","03102500","03155000","03171000","03180500","03208500","03252500","03335000","03410500","03424730","03479000","03588500","04099000","04100500","04172000","05053000","05060000","05122000","05127000","05275000","05433000","05463000","05479000","05481000","05495000","05515500","05534500","06050000","06120500","06178000","06278300","06336000","06386500","06478000","06481000","06600500","06620000","06715000","06719505","06856600","06870300","06876000","06904500","06906500","07026040","07063000","07086000","07156900","07195400","07238000","07242000","07274000","07282000","07288500","07312200","07326500","07328500","07331000","07331600","07375300","08026000","08044000","08111700","08176900","08279000","08287000","08390800","08407500","08449000","09050700","09115500","09217900","09306242","09379500","09512500","09519800"]

    usgs_folder = r"C:\Users\mheberger\Dropbox\RESEARCH\Watershed Article\benchmarking\Experiment 2 - check accuracy\watersheds\usgs"
    esri_folder = r"C:\Users\mheberger\Dropbox\RESEARCH\Watershed Article\benchmarking\Experiment 2 - check accuracy\watersheds\esri"
    mine_folder = r"C:\Users\mheberger\Dropbox\RESEARCH\Watershed Article\benchmarking\Experiment 2 - check accuracy\watersheds\mine"

    for gage in gages:
        usgs_shape = f"{usgs_folder}/gid_{gage}.shp"
        esri_shape = f"{esri_folder}/gid_{gage}.shp"
        mine_shape = f"{mine_folder}/gid_{gage}.shp"

        cac_esri = compare_shapefiles(usgs_shape, esri_shape)
        cac_mine = compare_shapefiles(usgs_shape, mine_shape)
        print(cac_esri, cac_mine)

