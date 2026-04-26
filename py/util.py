from shapely.geometry import MultiPolygon, Polygon


def get_largest(input_poly: MultiPolygon or Polygon) -> Polygon:
    """
    Converts a Shapely MultiPolygon to a Shapely Polygon
    For multipart polygons, will only keep the largest polygon
    in terms of area. In my testing, this was usually good enough

    Note: can also do this via PostGIS query... see myqueries_merit.py, query19a and query19b
          Not sure one approach is better than the other. They both seem to work well.
    Args:
        input_poly: A Shapely Polygon or MultiPolygon

    Returns:
        a shapely Polygon
    """
    if input_poly.geom_type == "MultiPolygon":
        areas = []
        polygons = list(input_poly.geoms)

        for poly in polygons:
            areas.append(poly.area)

        max_index = areas.index(max(areas))

        return polygons[max_index]
    else:
        return input_poly
