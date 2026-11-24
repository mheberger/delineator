import math
import logging

import numpy as np
from pyproj import Transformer, CRS
import shapely
from shapely.geometry import MultiPolygon, Polygon
import geopandas as gpd

from delineator.constants import KM_PER_DEGREE, ROUNDING_DECIMALS
from delineator.settings import DelineatorConfig

logger = logging.getLogger(__name__)


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
        If latitude or longitude cannot be converted to float
        (e.g. None or a non-numeric object).
    ValueError
        If latitude or longitude are unparseable strings, not finite,
        or fall outside the allowed range.
    """
    lat = float(lat)
    lng = float(lng)

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
    Fill the small interior holes ("donut holes") in a polygon.

    Parameters
    ----------
    poly : shapely Polygon or MultiPolygon
        The geometry, in unprojected lat/lng coordinates (EPSG:4326).
    area_max_km2 : float
        Keep holes that are larger than this, in square kilometers.
        Fill any holes less than or equal to this.
        Enter 0 to fill all holes.
        Each hole's area is converted to km² locally at the hole's own latitude.

    Returns
    -------
    shapely Polygon or MultiPolygon
        The geometry with small holes filled, of the same type as ``poly``.

    Raises
    ------
    ValueError
        If ``poly`` is not a Polygon or MultiPolygon.
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


def _reproject(geom, transformer: Transformer):
    """
    Reproject a shapely geometry with a pyproj Transformer (built with always_xy=True).

    Replaces the deprecated shapely.ops.transform(transformer.transform, geom).
    shapely.transform hands over all coordinates as one (N, 2) array, so pyproj
    transforms them in a single vectorized call.
    """
    return shapely.transform(
        geom, lambda xy: np.column_stack(transformer.transform(xy[:, 0], xy[:, 1]))
    )


def _get_polygon_area(poly: Polygon | MultiPolygon) -> float:
    """
    Calculates the area of a Shapely Polygon (or MultiPolygon) in km²

    Assumes the input is in unprojected raw lat, lng coordinates (CRS 4326)
    and projects it first to a locally-centered azimuthal equal-area projection.
    Works accurately anywhere on Earth.

    Parameters
    ----------
    poly : shapely Polygon or MultiPolygon
        The geometry, in EPSG:4326.

    Returns
    -------
    float
        The area, in km².
    """

    # Use the polygon centroid to define a custom LAEA projection
    lon, lat = poly.centroid.x, poly.centroid.y

    crs_laea = CRS.from_proj4(
        f"+proj=laea +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m"
    )
    crs_wgs84 = CRS.from_epsg(4326)

    transformer = Transformer.from_crs(crs_wgs84, crs_laea, always_xy=True)

    projected_poly = _reproject(poly, transformer)

    # Get the area in m -> km²
    return projected_poly.area / 1e6


def write_outputs(watershed_gdf: gpd.GeoDataFrame | None,
                  rivers_gdf: gpd.GeoDataFrame | None = None,
                  outlets_gdf: gpd.GeoDataFrame | None = None,
                  config: DelineatorConfig | None = None,
                  id: str | int | None = None) -> None:
    """
    Write the watershed, rivers, and outlets to disk.

    Files go in ``config.output_dir`` in the format ``config.output_format``.
    GeoPackage (``gpkg``) writes a single file, ``watershed.gpkg``, with up to
    three layers. Other formats write one file per layer: ``watershed``,
    ``rivers`` and ``outlets``. Supported formats include:

    - shapefile: ``shp``
    - GeoPackage: ``gpkg``
    - GeoJSON: ``geojson``, ``json``
    - GeoParquet: ``parquet``, ``geoparquet``
    - KML: ``kml``

    Parameters
    ----------
    watershed_gdf : GeoDataFrame or None
        The watershed, as returned by ``delineate``. Skipped if None.
    rivers_gdf : GeoDataFrame, optional
        The rivers. Skipped if None.
    outlets_gdf : GeoDataFrame, optional
        The outlet points. Skipped if None.
    config : DelineatorConfig, optional
        Supplies ``output_dir``, ``output_format`` and ``round_coordinates``.
        If not provided, one is created with default settings.
    id : str or int, optional
        Appended to each file name, as in ``watershed_{id}.gpkg``, so that
        outputs for several outlets do not overwrite each other.

    Notes
    -----
    When ``config.round_coordinates`` is True, watershed and river coordinates
    are rounded to 6 decimal places before writing. The rounding is applied to
    copies, so the GeoDataFrames passed in are not modified. If all three
    inputs are None, nothing is written and a warning is logged.
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

    layers = {
        "outlets": outlets_gdf,
        "watershed": watershed_gdf,
        "rivers": rivers_gdf,
    }

    if all(gdf is None for gdf in layers.values()):
        logger.warning("Nothing to write: the watershed, rivers and outlets are all None.")
        return

    # Round coordinates, if requested, on copies so the caller's data is untouched
    if config.round_coordinates:
        grid_size = 10 ** -ROUNDING_DECIMALS
        for layer_name in ("watershed", "rivers"):
            gdf = layers[layer_name]
            if gdf is not None:
                gdf = gdf.copy()
                gdf["geometry"] = shapely.set_precision(gdf.geometry, grid_size=grid_size)
                layers[layer_name] = gdf

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
    Return the radius of the Earth, optionally as a function of latitude.

    Python port from the Climate Data Toolbox for Matlab.

    Parameters
    ----------
    lat : float or array-like, optional
        Latitude(s) in decimal degrees at which to calculate the radius on the
        WGS84 ellipsoid. If None, returns the nominal radius of 6,371 km.
    unit : {'m', 'km'}, optional
        Unit of the output: 'm' for meters (default) or 'km' for kilometers.

    Returns
    -------
    float or ndarray
        The radius of the Earth, an array if ``lat`` is array-like.
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
    Calculate the area of a pixel on the Earth's surface in EPSG:4326.

    Python port from the Climate Data Toolbox for Matlab.

    Parameters
    ----------
    latitude : float
        Latitude of the pixel, in decimal degrees.
    lat_resolution : float
        Pixel height (latitude resolution), in degrees.
    lon_resolution : float or None
        Pixel width (longitude resolution), in degrees. If None, assumed equal
        to ``lat_resolution`` (square pixels, the most common case).

    Returns
    -------
    float
        The pixel area, in km².
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


