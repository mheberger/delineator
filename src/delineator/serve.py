"""
Spin up a simple interactive watershed
delineation web map with a local Flask server.

A lightweight web app inspired by Global Watersheds web app,
https://mghydro.com/watersheds, but able to run locally without
a server or even a web connection.

Lets the user open http://localhost:5000 in your browser, then
click anywhere on the map to delineate the upstream watershed.
Very fast, easy to use, and facilitates iteration and review.
"""

import copy
import json
import logging
import traceback
from importlib.resources import files

from flask import Flask, Response, jsonify, request

from delineator.core import delineate
from delineator.settings import DelineatorConfig

app = Flask(__name__)
logging.basicConfig(level=logging.WARNING)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# HTML / JS frontend (static single-page app; see static/index.html)
# ---------------------------------------------------------------------------

HTML = files("delineator").joinpath("static", "index.html").read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

# Flask config key under which serve() stashes the user-supplied DelineatorConfig.
_CONFIG_KEY = "DELINEATOR_CONFIG"

# The only DelineatorConfig fields the browser client is allowed to override per
# request. Deliberately excludes filesystem/download options (data_dir,
# output_dir, output_format, auto_download): those either don't apply here
# (results are streamed to the browser as GeoJSON, never written) or would let a
# page trigger arbitrary filesystem writes or large downloads.
_WEB_OPTIONS = frozenset({
    "calc_area",
    "clean",
    "fill",
    "fill_area_max",
    "high_res",
    "num_stream_orders",
    "outlets",
    "rivers",
    "round_coordinates",
    "search_dist",
    "simplify",
    "simplify_tolerance",
    "smooth",
    "snapping",
    "stream_threshold_km2",
})


def _apply_web_options(config: DelineatorConfig, options: dict) -> DelineatorConfig:
    """Apply whitelisted client options to a per-request config, in place.

    Only keys in _WEB_OPTIONS are considered; anything else is ignored. Values
    are validated (and coerced, e.g. int -> float) by constructing a throwaway
    DelineatorConfig, which raises ValueError/TypeError on bad input. The
    validated values are then copied onto the caller's config so its resolved
    data_dir/output_dir and custom_data_dir flag are preserved.
    """
    opts = {k: v for k, v in options.items() if k in _WEB_OPTIONS}
    if not opts:
        return config

    validated = DelineatorConfig(**opts)  # raises on invalid input
    for key in opts:
        setattr(config, key, getattr(validated, key))
    return config


def _request_config() -> DelineatorConfig:
    """Return the DelineatorConfig to use for one delineation request.

    Uses the config passed to serve() if there is one, otherwise a sensible
    default for the interactive app (rivers, outlets, and a cleaned boundary).

    The object is copied per request because _apply_web_options() writes the
    client's per-request option overrides onto it; without the copy those
    overrides would leak into every later request sharing the same config.
    """
    base = app.config.get(_CONFIG_KEY)
    if base is None:
        return DelineatorConfig(rivers=True, outlets=True, clean=True)
    return copy.deepcopy(base)


@app.route("/")
def index() -> Response:
    return Response(HTML, mimetype="text/html")


@app.route("/delineate", methods=["POST"])
def delineate_endpoint() -> Response:
    body = request.get_json(silent=True) or {}
    lat = body.get("lat")
    lng = body.get("lng")

    if lat is None or lng is None:
        return jsonify({"error": "Request body must include 'lat' and 'lng'."}), 400

    try:
        lat = float(lat)
        lng = float(lng)
    except (TypeError, ValueError):
        return jsonify({"error": "'lat' and 'lng' must be numeric."}), 400

    logger.info(f"Delineating watershed at lat={lat:.3f}, lng={lng:.3f}")

    config = _request_config()
    try:
        _apply_web_options(config, body.get("options") or {})
    except (ValueError, TypeError) as exc:
        return jsonify({"error": f"Invalid option: {exc}"}), 400

    try:
        watershed_gdf, rivers_gdf, outlets_gdf = delineate(lat, lng, config=config)
    except Exception as exc:
        logger.error(traceback.format_exc())
        return jsonify({"error": f"Delineation failed: {exc}"}), 500

    if watershed_gdf is None:
        return jsonify({
            "error": (
                "Could not delineate a watershed at that location. "
                "Make sure the point is over land (not ocean or a dry basin)."
            )
        }), 422

    def gdf_to_geojson(gdf):
        if gdf is None:
            return None
        return json.loads(gdf.to_json())

    return jsonify({
        "watershed": gdf_to_geojson(watershed_gdf),
        "rivers":    gdf_to_geojson(rivers_gdf),
        "outlets":   gdf_to_geojson(outlets_gdf),
    })


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def serve(
    config: DelineatorConfig | None = None,
    host: str = "127.0.0.1",
    port: int = 5000,
    debug: bool = False,
) -> None:
    """Start the local watershed delineation web app.

    Works both from the command line and when called from Python, e.g.::

    Usage:
    ------
        from delineator.serve import serve
        serve()

        # Customize how every click is delineated:
        from delineator import DelineatorConfig
        serve(DelineatorConfig(high_res=False, num_stream_orders=6, fill=False))

    Parameters:
    -----------
    config: DelineatorConfig applied to every delineation the app performs.
        If None, a default suited to interactive review is used
        (rivers, outlets, and a cleaned boundary). The output_format option
        is ignored here since results are streamed to the browser as GeoJSON;
        the map shows whichever layers the config produces.
    host: Interface to bind. Defaults to localhost only.
    port: Port to listen on.
    debug: Enable Flask's debug error pages. The auto-reloader is always
        left off (use_reloader=False) so that calling serve() from a
        Python session, notebook, or another script reliably starts the
        server. Leave debug off when distributing the package.

    """
    app.config[_CONFIG_KEY] = config

    print("\n  Global Watershed Delineator Local Web App")
    print("  -----------------------------------------")
    print(f"  Open http://{host}:{port} in your browser")
    print("  Click anywhere on the map to delineate\n")
    # use_reloader=False is essential: the reloader runs the real server in a
    # re-executed child process, which only works when launched as a command.
    app.run(host=host, port=port, debug=debug, use_reloader=False, threaded=True)


if __name__ == "__main__":
    serve()
