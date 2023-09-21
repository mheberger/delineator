"""
Faster Dissolve All with GeoPandas

For a layer with many polygons, it can be slow to dissolve to get the "outer boundary" or "outer perimeter"
using GeoPandas

I found a method that works a little bit more quickly.

(1) create a new rectangle out of the bounding box around all the features. 
(2) clip the rectangle using the input layer (containing polygons).


input: a geopandas dataframe with multiple polygons.
output: a geopandas dataseries with a single polygon
with no internal rings or "donut holes," which is what I was looking for
with my watershed boundaries. 

"""

import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon
gpd.options.use_pygeos = True


def buffer(poly: Polygon) -> Polygon:
    """
    Little trick that works wonders to remove slivers, dangles
    and other weird errors in a shapely polygon. We do a series of
    2 buffers, out and then in, and it magically fixes issues.

    """
    dist = 0.00001
    return poly.buffer(dist, join_style=2).buffer(-dist, join_style=2)


def close_holes(poly: Polygon or MultiPolygon, area_max: float) -> Polygon or MultiPolygon:
    """
    Close polygon holes by limitation to the exterior ring.
    Updated to accept a MultiPolygon as input
    Args:
        poly: Input shapely Polygon
        area_max: keep holes that are larger than this.
                  Fill any holes less than or equal to this.
                  We're working with unprojected lat, lng
                  so this needs to be in square decimal degrees...
    Example:
        df.geometry.apply(lambda p: close_holes(p))
    """

    if isinstance(poly, Polygon):
        # Handle Polygon case
        if area_max == 0:
            if poly.interiors:
                return Polygon(list(poly.exterior.coords))
            else:
                return poly

        else:
            list_interiors = []

            for interior in poly.interiors:
                p = Polygon(interior)
                if p.area > area_max:
                    list_interiors.append(interior)

            return Polygon(poly.exterior.coords, holes=list_interiors)

    elif isinstance(poly, MultiPolygon):
        # Handle MultiPolygon case
        result_polygons = []
        for sub_poly in poly:
            new_sub_poly = close_holes(sub_poly, area_max)
            result_polygons.append(new_sub_poly)
        return MultiPolygon(result_polygons)
    else:
        raise ValueError("Unsupported geometry type")



def dissolve_shp(shp: str) -> gpd.GeoDataFrame:
    """
    input is the path to a shapefile on disk. 
    
    Returns a GeoPandas dataframe containing the dissolved
    geometry
    """
    df = gpd.read_file(shp)
    return dissolve_geopandas(df)


def fill_geopandas(gdf: gpd.GeoDataFrame, area_max: float) -> gpd.GeoDataFrame:
    filled = gdf.geometry.apply(lambda p: close_holes(p, area_max))
    return filled


def dissolve_geopandas(df: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    input is a Geopandas dataframe with multiple polygons that you want 
      to merge and dissolve into a single polygon
      
    output is a Geopandas dataframe containing a single polygon

    This method is much faster than using GeoPandas dissolve()

    It creates a box around the polygons, then clips the box to
    that poly. The result is one feature instead of many.
    """
    
    [left, bottom, right, top] = df.total_bounds
    left -= 1
    right += 1
    top += 1
    bottom -= 1

    lat_point_list = [left, right, right, left, left]
    lon_point_list = [top, top, bottom, bottom, top]


    polygon_geom = Polygon(zip(lat_point_list, lon_point_list))
    rect = gpd.GeoDataFrame(index=[0], crs=df.crs, geometry=[polygon_geom])
    clipped = gpd.clip(rect, df)
    # This removes some weird artifacts that result from Merit-BASINS having lots
    # of little topology issues.

    clipped = clipped.geometry.apply(lambda p: buffer(p))

    return clipped
    