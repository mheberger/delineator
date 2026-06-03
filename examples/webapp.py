"""
webapp.py — Interactive watershed delineation web map (local Flask server)

Usage:
    pip install flask
    python app.py

Then open http://localhost:5000 in your browser.
Click anywhere on the map to delineate the upstream watershed.
"""

import json
import logging
import traceback

from flask import Flask, Response, jsonify, request

from delineator.core import delineate
from delineator.settings import DelineatorConfig

app = Flask(__name__)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# HTML / JS frontend (served inline so the app is a single file)
# ---------------------------------------------------------------------------

HTML = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>🌍 Global Watershed Delineator</title>

  <!-- Leaflet -->
  <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
  <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>

  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }

    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #1a1a2e;
      color: #e0e0e0;
      height: 100vh;
      display: flex;
      flex-direction: column;
    }

    /* ── Header ── */
    header {
      display: flex;
      align-items: center;
      gap: 12px;
      padding: 10px 16px;
      background: #16213e;
      border-bottom: 1px solid #0f3460;
      flex-shrink: 0;
      z-index: 1000;
    }
    header h1 {
      font-size: 1.1rem;
      font-weight: 600;
      color: #e0e0e0;
      letter-spacing: 0.02em;
    }
    header .subtitle {
      font-size: 0.78rem;
      color: #8899aa;
    }

    /* ── Status bar ── */
    #status-bar {
      padding: 6px 16px;
      background: #0f3460;
      font-size: 0.8rem;
      color: #a0c4ff;
      flex-shrink: 0;
      display: flex;
      align-items: center;
      gap: 8px;
      min-height: 30px;
    }
    #status-bar .spinner {
      display: none;
      width: 14px; height: 14px;
      border: 2px solid #a0c4ff44;
      border-top-color: #a0c4ff;
      border-radius: 50%;
      animation: spin 0.7s linear infinite;
      flex-shrink: 0;
    }
    #status-bar.loading .spinner { display: block; }
    @keyframes spin { to { transform: rotate(360deg); } }

    /* ── Map ── */
    #map {
      flex: 1;
      cursor: crosshair;
    }
    #map.loading { cursor: wait; }

    /* ── Info panel (bottom-right overlay) ── */
    #info-panel {
      position: absolute;
      bottom: 30px;
      right: 10px;
      z-index: 1000;
      background: rgba(22, 33, 62, 0.92);
      border: 1px solid #0f3460;
      border-radius: 8px;
      padding: 12px 16px;
      font-size: 0.78rem;
      max-width: 240px;
      display: none;
      backdrop-filter: blur(6px);
    }
    #info-panel h3 {
      font-size: 0.85rem;
      color: #a0c4ff;
      margin-bottom: 8px;
    }
    #info-panel .stat { margin-bottom: 4px; color: #c8d8e8; }
    #info-panel .stat span { color: #ffffff; font-weight: 600; }

    /* ── Legend ── */
    #legend {
      position: absolute;
      bottom: 30px;
      left: 10px;
      z-index: 1000;
      background: rgba(22, 33, 62, 0.88);
      border: 1px solid #0f3460;
      border-radius: 8px;
      padding: 10px 14px;
      font-size: 0.75rem;
      backdrop-filter: blur(4px);
    }
    #legend .legend-row {
      display: flex; align-items: center; gap: 8px; margin-bottom: 5px;
    }
    #legend .legend-row:last-child { margin-bottom: 0; }
    .swatch {
      width: 28px; height: 12px; border-radius: 2px; flex-shrink: 0;
    }
    .swatch.watershed {
      background: rgba(230, 92, 92, 0.3);
      border: 2px solid #f50a0a;
    }
    .swatch.river {
      background: #00bfff;
      height: 3px; border-radius: 2px;
    }
    .swatch.outlet-requested {
      width: 12px; height: 12px; border-radius: 50%;
      background: #fdce0b; border: 2px solid #fff;
      flex-shrink: 0;
    }
    .swatch.outlet-snapped {
      width: 12px; height: 12px; border-radius: 50%;
      background: #ff6b35; border: 2px solid #fff;
      flex-shrink: 0;
    }

    /* ── Error toast ── */
    #error-toast {
      position: absolute;
      top: 80px; left: 50%; transform: translateX(-50%);
      z-index: 2000;
      background: #c0392b;
      color: #fff;
      border-radius: 6px;
      padding: 8px 18px;
      font-size: 0.82rem;
      display: none;
      max-width: 90%;
      text-align: center;
      box-shadow: 0 4px 12px rgba(0,0,0,0.4);
    }
  </style>
</head>
<body>

<header>
  <div>
    <h1>🌍 Global Watershed Delineator</h1>
  </div>
  <div class="subtitle">Click anywhere on the map to delineate the upstream watershed</div>
</header>

<div id="status-bar">
  <div class="spinner"></div>
  <span id="status-text">Ready — click the map to delineate a watershed</span>
</div>

<div id="map"></div>

<!-- Legend -->
<div id="legend">
  <div class="legend-row"><div class="swatch watershed"></div> Watershed</div>
  <div class="legend-row"><div class="swatch river"></div> Rivers</div>
  <div class="legend-row"><div class="swatch outlet-requested"></div> Outlet, requested</div>
  <div class="legend-row"><div class="swatch outlet-snapped"></div> Outlet, snapped</div>
</div>

<!-- Info panel -->
<div id="info-panel">
  <h3>Watershed Info</h3>
  <div class="stat">Outlet Lat: <span id="info-lat">—</span></div>
  <div class="stat">Outlet Lng: <span id="info-lng">—</span></div>
  <div class="stat">Watershed area: <span id="info-area">—</span> km²</div>
  <div class="stat">River reaches returned: <span id="info-rivers">—</span></div>
</div>

<!-- Error toast -->
<div id="error-toast" id="error-toast"></div>

<script>
// ── Map init ──────────────────────────────────────────────────────────────
const map = L.map('map').setView([20, 0], 3);

L.tileLayer('https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png', {
  attribution: '© OpenStreetMap contributors © CARTO',
  subdomains: 'abcd',
  maxZoom: 20,
}).addTo(map);

// ── Layer groups ──────────────────────────────────────────────────────────
let watershedLayer = null;
let riversLayer    = null;
let outletLayer    = null;

function clearLayers() {
  if (watershedLayer) { map.removeLayer(watershedLayer); watershedLayer = null; }
  if (riversLayer)    { map.removeLayer(riversLayer);    riversLayer    = null; }
  if (outletLayer)    { map.removeLayer(outletLayer);    outletLayer    = null; }
}

// ── UI helpers ────────────────────────────────────────────────────────────
const statusBar  = document.getElementById('status-bar');
const statusText = document.getElementById('status-text');
const mapEl      = document.getElementById('map');
const infoPanel  = document.getElementById('info-panel');
const errorToast = document.getElementById('error-toast');
let errorTimer   = null;

function setLoading(msg) {
  statusBar.classList.add('loading');
  statusText.textContent = msg;
  mapEl.classList.add('loading');
}

function setReady(msg) {
  statusBar.classList.remove('loading');
  statusText.textContent = msg;
  mapEl.classList.remove('loading');
}

function showError(msg) {
  errorToast.textContent = msg;
  errorToast.style.display = 'block';
  clearTimeout(errorTimer);
  errorTimer = setTimeout(() => { errorToast.style.display = 'none'; }, 5000);
}

function updateInfoPanel(lat, lng, watershedArea, riverCount) {
  document.getElementById('info-lat').textContent    = lat.toFixed(5);
  document.getElementById('info-lng').textContent    = lng.toFixed(5);
  document.getElementById('info-area').textContent    = watershedArea;
  document.getElementById('info-rivers').textContent = riverCount !== null ? riverCount : 'none';
  infoPanel.style.display = 'block';
}

// ── Delineation ───────────────────────────────────────────────────────────
map.on('click', async (e) => {
  const { lat, lng } = e.latlng;

  clearLayers();
  setLoading(`Delineating watershed at (${lat.toFixed(4)}, ${lng.toFixed(4)})…`);
  infoPanel.style.display = 'none';

  try {
    const resp = await fetch('/delineate', {
      method:  'POST',
      headers: { 'Content-Type': 'application/json' },
      body:    JSON.stringify({ lat, lng }),
    });

    const data = await resp.json();

    if (!resp.ok || data.error) {
      showError(data.error || `Server error ${resp.status}`);
      setReady('Error — try a different location');
      return;
    }

    // ── Watershed polygon ──
    let watershedArea = null;
    if (data.watershed) {
      watershedArea = data.watershed.features[0].properties.area_km2;
      watershedLayer = L.geoJSON(data.watershed, {
        style: {
          color:       '#f50a0a',
          weight:      2,
          opacity:     0.9,
          fillColor:   '#e65c5c',
          fillOpacity: 0.22,
        },
      }).addTo(map);
      map.fitBounds(watershedLayer.getBounds(), { padding: [40, 40] });
    }

    // ── Rivers ──
    let riverCount = null;
    if (data.rivers) {
      const features = data.rivers.features || [];
      riverCount = features.length - 1;
      riversLayer = L.geoJSON(data.rivers, {
        style: {
          color:   '#00bfff',
          weight:  1.5,
          opacity: 0.85,
        },
      }).addTo(map);
    }

    // ── Outlets ──
    if (data.outlets) {
      outletLayer = L.geoJSON(data.outlets, {
        pointToLayer: (feature, latlng) => {
          const isSnapped = feature.properties?.type === 'snapped';
          return L.circleMarker(latlng, {
            radius:      isSnapped ? 6 : 5,
            fillColor:   isSnapped ? '#ff6b35' : '#ffcc00',
            color:       '#fff',
            weight:      2,
            opacity:     1,
            fillOpacity: 0.95,
          });
        },
        onEachFeature: (feature, layer) => {
          if (feature.properties?.type) {
            layer.bindTooltip(feature.properties.type, { permanent: false });
          }
        },
      }).addTo(map);
    }

    updateInfoPanel(lat, lng, watershedArea, riverCount);
    setReady(`Done — click elsewhere to delineate another watershed`);

  } catch (err) {
    console.error(err);
    showError('Network or server error. Is the Flask server running?');
    setReady('Error — see console for details');
  }
});
</script>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

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

    logger.info(f"Delineating watershed at lat={lat:.5f}, lng={lng:.5f}")

    config = DelineatorConfig(
        rivers=True,
        outlets=True,
        clean=True,
    )

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

if __name__ == "__main__":
    print("\n  Global Watershed Delineator")
    print("  ─────────────────────────────────────────")
    print("  Open http://localhost:5000 in your browser")
    print("  Click anywhere on the map to delineate\n")
    app.run(debug=True, port=5000)
