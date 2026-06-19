import math

from shapely.geometry import MultiPolygon, Polygon

from settings import DelineatorConfig


def get_overlap_lines(line, polygon):
    """Return the portion of a line that falls within a polygon, always as a list of LineStrings."""
    result = line.intersection(polygon)

    if result.is_empty:
        return []
    elif result.geom_type == "LineString":
        return [result]
    elif result.geom_type == "MultiLineString":
        return list(result.geoms)
    elif result.geom_type == "GeometryCollection":
        # Filter out any Point/MultiPoint components (tangent touches)
        return [g for g in result.geoms if g.geom_type in ("LineString", "MultiLineString")]
    else:
        return []  # Pure Point result — no real overlap


def validate_outlet_coordinates(lat: float, lng: float) -> bool:
    """
    Validates the geographical coordinates (latitude and longitude) of an outlet
     ensuring they fall within acceptable ranges.

    :param lat: Latitude of the outlet in decimal degrees. Must be a finite float value between -60 and 85 (exclusive).
    :type lat: float
    :param lng: Longitude of the outlet in decimal degrees. Must be a finite float value between -180 and 180 (exclusive).
    :type lng: float
    :return: True if the latitude and longitude are valid geographical coordinates.
    :rtype: bool
    :raises Exception: If latitude or longitude are not finite float values or fall outside the allowed range.
    """
    if not isinstance(lat, float):
        raise Warning(f"Latitude must be a float. Got {type(lat).__name__}: {lat}")

    if not isinstance(lng, float):
        raise Warning(f"Longitude must be a float. Got {type(lng).__name__}: {lng}")

    if not math.isfinite(lat):
        raise Warning(f"Latitude must be finite. Got {lat}")

    if not math.isfinite(lng):
        raise Warning(f"Longitude must be finite. Got {lng}")

    if lat <= -60:
        raise Warning("Latitude must be greater than -60°")

    if lat >= 85:
        raise Warning("Latitude must be less than 85°")

    if lng <= -180:
        raise Warning("Longitude must be greater than -180°")

    if lng >= 180:
        raise Warning("Longitude must be less than 180°")

    return True
