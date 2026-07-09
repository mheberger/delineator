import math
import logging

import numpy as np
from pyproj import Transformer, CRS
from shapely.ops import transform, unary_union
import shapely
from shapely.geometry import MultiPolygon, Polygon, LineString, MultiLineString
import geopandas as gpd

from delineator.settings import DelineatorConfig

logger = logging.getLogger(__name__)


def _get_overlap_lines(line: LineString | MultiLineString,
                       polygon: Polygon | MultiPolygon) \
        -> list[LineString | MultiLineString]:
    """
    Return the portion of a line that falls within a polygon,

    Always returns a list of LineStrings.
    Uses shapely library functions only, no spatialite queries.

    """
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


def _validate_outlet_coordinates(lat: float, lng: float) -> tuple[float, float]:
    """
    Validates the geographical coordinates (latitude and longitude) of an outlet
    ensuring they fall within acceptable ranges.

    Parameters
    ----------
    lat : float
        Latitude of the outlet in decimal degrees. Must be a finite float value
        between -60 and 85 (exclusive).
    lng : float
        Longitude of the outlet in decimal degrees. Must be a finite float
        value between -180 and 180 (exclusive).

    Returns
    -------
    tuple[float, float]
        Validated latitude and longitude.

    Raises
    ------
    TypeError
        If latitude or longitude are not float values.
    ValueError
        If latitude or longitude are not finite or fall outside the allowed
        range.
    """
    lat = float(lat)
    lng = float(lng)

    if not isinstance(lat, float):
        raise TypeError(f"Latitude must be a float. Got {type(lat).__name__}: {lat}")

    if not isinstance(lng, float):
        raise TypeError(f"Longitude must be a float. Got {type(lng).__name__}: {lng}")

    if not math.isfinite(lat):
        raise ValueError(f"Latitude must be finite. Got {lat}")

    if not math.isfinite(lng):
        raise ValueError(f"Longitude must be finite. Got {lng}")

    if lat <= -60:
        raise ValueError("Latitude must be greater than -60°")

    if lat >= 85:
        raise ValueError("Latitude must be less than 85°")

    if lng <= -180:
        raise ValueError("Longitude must be greater than -180°")

    if lng >= 180:
        raise ValueError("Longitude must be less than 180°")

    return lat, lng


# Kilometers per degree of great-circle arc on the authalic sphere
# (R = 6371.0088 km): pi * R / 180. Used for local degree <-> km conversions.
# TODO: move to delineator.constants alongside EARTH_RADIUS_KM (also defined
# in merit_detailed.py) so the package defines the Earth exactly once.
KM_PER_DEGREE = math.pi * 6371.0088 / 180  # ~111.195


def _ring_area_km2(ring) -> float:
    """
    Approximate area in km² of a shapely LinearRing (e.g. a polygon interior),
    using a local equirectangular conversion at the ring's own latitude:

        km² = deg² * KM_PER_DEGREE² * cos(latitude)

    Accurate to within about 1% (the residual is mostly the spherical-Earth
    approximation, not hole size); ample for threshold tests. For accurate
    areas of large geometries use _get_polygon_area instead.
    """
    p = Polygon(ring)
    lat = p.representative_point().y
    return p.area * KM_PER_DEGREE ** 2 * abs(math.cos(math.radians(lat)))


def _close_holes(poly: Polygon | MultiPolygon, area_max_km2: float) -> Polygon | MultiPolygon:
    """
    Close polygon holes by limitation to the exterior ring.
    Updated to accept a MultiPolygon as input.
    Example:
        df.geometry.apply(lambda p: close_holes(p))

    Parameters
    ----------
    poly: shapely Polygon or MultiPolygon
    area_max_km2: float
        Keep holes that are larger than this, in square kilometers.
        Fill any holes less than or equal to this.
        Enter 0 to fill all holes.
        The geometry is in unprojected lat/lng coordinates (EPSG:4326); each
        hole's area is converted to km² locally at the hole's own latitude.

    """

    if isinstance(poly, Polygon):
        # Handle Polygon case
        if area_max_km2 == 0:
            if poly.interiors:
                return Polygon(list(poly.exterior.coords))
            else:
                return poly

        else:
            list_interiors = []

            for interior in poly.interiors:
                if _ring_area_km2(interior) > area_max_km2:
                    list_interiors.append(interior)

            return Polygon(poly.exterior.coords, holes=list_interiors)

    elif isinstance(poly, MultiPolygon):
        # Handle MultiPolygon case
        result_polygons = []
        for sub_poly in poly.geoms:
            new_sub_poly = _close_holes(sub_poly, area_max_km2)
            result_polygons.append(new_sub_poly)
        return MultiPolygon(result_polygons)
    else:
        raise ValueError("Unsupported geometry type")


def _get_polygon_area(poly: Polygon | MultiPolygon) -> float:
    """
    Calculates the area of a Shapely Polygon (or MultiPolygon) in km²

    Assumes the input is in unprojected raw lat, lng coordinates (CRS 4326)
    and projects it first to  a locally-centered azimuthal equal-area projection.
    Works accurately anywhere on Earth.

    Parameters
    ----------
    poly: Shapely Polygon or MultiPolygon

    Returns
    -------
         area in km²
    """

    # Use the polygon centroid to define a custom LAEA projection
    lon, lat = poly.centroid.x, poly.centroid.y

    crs_laea = CRS.from_proj4(
        f"+proj=laea +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m"
    )
    crs_wgs84 = CRS.from_epsg(4326)

    transformer = Transformer.from_crs(crs_wgs84, crs_laea, always_xy=True)

    projected_poly = transform(transformer.transform, poly)

    # Get the area in m -> km²
    return projected_poly.area / 1e6


def safe_simplify(geom, tolerance):
    if geom is None or geom.is_empty:
        return geom
    # Flatten GeometryCollections to just polygonal parts
    if geom.geom_type == 'GeometryCollection':
        polys = [g for g in geom.geoms if g.geom_type in ('Polygon', 'MultiPolygon')]
        if not polys:
            return geom  # or return None
        geom = unary_union(polys)
    return geom.simplify(tolerance, preserve_topology=True)


def write_outputs(watershed_gdf: gpd.GeoDataFrame | None,
                  rivers_gdf: gpd.GeoDataFrame | None = None,
                  outlets_gdf: gpd.GeoDataFrame | None = None,
                  config: DelineatorConfig | None = None,
                  id: str | int | None = None) -> None:
    """
    Writes the watershed, rivers, and outlets to disk in the specified format.

    With geopackage, will write a single file with up to 3 layers.
    Other formats will result in up to 3 files.

    Supports common GeoPandas output formats, including:
    - shapefile: shp
    - geopackage: gpkg
    - geojson: geojson, json
    - parquet: parquet
    - kml: kml
    """
    if config is None:
        config = DelineatorConfig()

    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    output_format = config.output_format.lower().lstrip(".")

    drivers = {
        "shp": "ESRI Shapefile",
        "gpkg": "GPKG",
        "geojson": "GeoJSON",
        "json": "GeoJSON",
        "kml": "KML",
    }

    suffix = "" if id is None else f"_{id}"

    # Round coordinates to 6 decimal places
    grid_size = 10 ** -6
    watershed_gdf["geometry"] = shapely.set_precision(watershed_gdf.geometry, grid_size=grid_size)
    if rivers_gdf is not None:
        rivers_gdf["geometry"] = shapely.set_precision(rivers_gdf.geometry, grid_size=grid_size)

    layers = {
        "outlets": outlets_gdf,
        "watershed": watershed_gdf,
        "rivers": rivers_gdf,
    }

    if output_format in {"parquet", "geoparquet"}:
        for layer_name, gdf in layers.items():
            if gdf is not None:
                fname = output_dir / f"{layer_name}{suffix}.parquet"
                gdf.to_parquet(
                    fname,
                    index=False,
                )
                logger.info(f"Wrote {fname}")
        return

    if output_format == "gpkg":
        output_path = output_dir / f"watershed{suffix}.gpkg"

        for layer_name, gdf in layers.items():
            if gdf is not None:
                gdf.to_file(
                    output_path,
                    layer=layer_name,
                    driver="GPKG",
                    encoding="utf-8",
                    index=False,
                )
        logger.info(f"Wrote {output_path}")
        return

    driver = drivers.get(output_format)

    # All other formats (SHP, GeoJSON, KML)
    for layer_name, gdf in layers.items():
        if gdf is None:
            continue

        output_path = output_dir / f"{layer_name}{suffix}.{output_format}"

        write_kwargs = {
            "filename": output_path,
            "encoding": "utf-8",
            "index": False,
        }

        if driver is not None:
            write_kwargs["driver"] = driver

        if output_format == "shp":
            write_kwargs["layer"] = layer_name

        gdf.to_file(**write_kwargs)

        logger.info(f"Wrote {output_path}")


def earth_radius(lat=None, unit='m'):
    """
    Earth_radius gives the radius of the Earth.

    Python port from the Climate Data Toolbox for Matlab.

    Args:
    lat (float or array-like, optional): Latitude(s) at which to calculate the Earth's radius.
    unit (str, optional): Unit of the output, 'm' for meters or 'km' for kilometers. Default is 'm'.

    Returns:
    float or ndarray: The radius of the Earth. If lat is None, returns the nominal radius.
                      If lat is provided, returns the radius as a function of latitude.
    """
    # Nominal radius of the Earth in meters
    r = 6371000

    # If latitude is provided, calculate the radius as a function of latitude
    if lat is not None:
        a = 6378137  # equatorial radius in meters, WGS84
        b = 6356752  # polar radius in meters, WGS84
        lat = np.asarray(lat)
        lat_rad = np.deg2rad(lat)
        r = np.sqrt(((a ** 2 * np.cos(lat_rad)) ** 2 + (b ** 2 * np.sin(lat_rad)) ** 2) /
                    ((a * np.cos(lat_rad)) ** 2 + (b * np.sin(lat_rad)) ** 2))

    # Convert to kilometers if requested
    if unit.lower() == 'km':
        r = r / 1000

    return r


def pixel_area(latitude: float, lat_resolution: float, lon_resolution: float | None) -> float:
    """
    Calculate the area of a pixel on Earth's surface in CRS 4326
    for a given latitude and a given resolution

    Python port from the Climate Data Toolbox for Matlab.

    Args:
    latitude (float): Latitude of the pixel.
    lat_resolution (float): Latitude resolution in degrees.
    lon_resolution (float or None): Longitude resolution in degrees. If
      none, assume equal to lat_resolution (square pixels, most common)

    Returns:
    area: pixel area in square kilometers.
    """

    # If no longitude resolution given, assume square pixels
    if lon_resolution is None:
        lon_resolution = lat_resolution

    # Calculate the Earth's radius at each latitude
    R = earth_radius(lat=latitude)

    # Calculate dy and dx
    dy = lat_resolution * np.pi / 180 * R  # in meters
    dx = lon_resolution * np.pi / 180 * R * np.cos(np.deg2rad(latitude))  # in meters

    # Calculate the area of each pixel
    pixel_area = np.abs(dx * dy) / 1e6  # convert from square meters to square kilometers

    return pixel_area
