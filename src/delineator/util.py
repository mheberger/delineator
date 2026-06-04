import math
import logging
from pyproj import Transformer, CRS
from shapely.ops import transform
import shapely
from shapely.geometry import MultiPolygon, Polygon, LineString, MultiLineString
import geopandas as gpd

from delineator.settings import DelineatorConfig

logger = logging.getLogger(__name__)


def _get_overlap_lines(line: LineString | MultiLineString,
                       polygon: Polygon | MultiPolygon) \
        -> list[LineString | MultiLineString]:
    """Return the portion of a line that falls within a polygon, always as a list of LineStrings.
    Uses shapely library functions only, no spatialite queries.
    Params:

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

    :param lat: Latitude of the outlet in decimal degrees. Must be a finite
        float value between -60 and 85 (exclusive).
    :type lat: float
    :param lng: Longitude of the outlet in decimal degrees. Must be a finite
        float value between -180 and 180 (exclusive).
    :type lng: float
    :return: Validated latitude and longitude.
    :rtype: tuple[float, float]
    :raises TypeError: If latitude or longitude are not float values.
    :raises ValueError: If latitude or longitude are not finite or fall outside
        the allowed range.
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


def _close_holes(poly: Polygon | MultiPolygon, area_max: float) -> Polygon | MultiPolygon:
    """
    Close polygon holes by limitation to the exterior ring.
    Updated to accept a MultiPolygon as input
    Args:
        poly: Input shapely Polygon
        area_max: keep holes that are larger than this.
                  Fill any holes less than or equal to this.
                  Enter 0 to fill all holes.
                  We're working with unprojected lat, lng
                  so this needs to be in square decimal degrees.
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
        for sub_poly in poly.geoms:
            new_sub_poly = _close_holes(sub_poly, area_max)
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

    Args:
        poly: Shapely Polygon or MultiPolygon

    Returns:
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

    # Get the area in km²
    return projected_poly.area / 1e6


def write_outputs(watershed_gdf: gpd.GeoDataFrame | None,
                  rivers_gdf: gpd.GeoDataFrame | None,
                  outlets_gdf: gpd.GeoDataFrame | None,
                  config: DelineatorConfig,
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

    if config.simplify:
        watershed_gdf['geometry'] = watershed_gdf.geometry.simplify(
            tolerance=config.simplify_tolerance,
            preserve_topology=True
        )
        if rivers_gdf is not None:
            rivers_gdf['geometry'] = rivers_gdf.geometry.simplify(
                tolerance=config.simplify_tolerance,
                preserve_topology=True
            )

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
                gdf.to_parquet(
                    output_dir / f"{layer_name}{suffix}.parquet",
                    index=False,
                )
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
        return

    driver = drivers.get(output_format)

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

        logging.info(f"Wrote {output_path}")
