"""
batch_app.py — Batch watershed delineation from a CSV file, with a live results
table and Leaflet map viewer. Created with Claude.ai.

Usage:
    pip install flask
    python batch_app.py iceland_outlets.csv

    # Or to use a different CSV:
    python batch_app.py path/to/my_outlets.csv

Then open http://localhost:5001 in your browser.

Expected CSV columns: id, lat, lng, name, area
  - id   : unique identifier (will become a clickable link)
  - lat  : outlet latitude  (decimal degrees)
  - lng  : outlet longitude (decimal degrees)
  - name : human-readable label
  - area : a priori watershed area (km²), used to compute % difference

Results are saved to disk as each watershed completes:
  - output/<id>/watershed.geojson
  - output/<id>/rivers.geojson      (if rivers were returned)
  - output/<id>/outlets.geojson     (if outlets were returned)
  - output/results.csv              (summary table, appended row by row)

On restart, any outlet whose output directory already exists is skipped and its
previously computed results are loaded from results.csv instead of re-running
delineation. Delete the output directory (or a specific subdirectory) to force
a re-run.

GeoJSON is never held in memory or sent through the SSE stream. The browser
fetches geometry on demand via GET /geometry/<id> only when the user clicks a
row, so memory usage stays flat regardless of batch size.
"""

import csv
import json
import logging
import sys
import traceback
from pathlib import Path

from flask import Flask, Response, jsonify, stream_with_context

from delineator.core import delineate
from delineator.settings import DelineatorConfig

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)

# CSV is resolved once at startup (passed as CLI arg or defaults to iceland_outlets.csv)
_csv_path: Path = Path("iceland_outlets.csv")

# Output directory: one subdirectory per outlet, plus a summary CSV
_output_dir: Path = Path("output")

# Summary CSV written/appended as results arrive
_results_csv: Path = _output_dir / "results.csv"

# Ordered list of row dicts from the CSV (id, lat, lng, name, area)
_rows: list[dict] = []

# ---------------------------------------------------------------------------
# CSV loading
# ---------------------------------------------------------------------------

def _load_csv(path: Path) -> list[dict]:
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "id":   row["id"].strip(),
                "lat":  float(row["lat"]),
                "lng":  float(row["lng"]),
                "name": row.get("name", "").strip(),
                "area": float(row["area"]),
            })
    return rows


# ---------------------------------------------------------------------------
# Area computation (equal-area projection, accurate at high latitudes)
# ---------------------------------------------------------------------------

def _compute_area_km2(watershed_gdf) -> float:
    """Reproject to EPSG:6933 (WGS84 equal-area) and return area in km²."""
    projected = watershed_gdf.to_crs("EPSG:6933")
    return float(projected.geometry.area.sum()) / 1e6


def _pct_diff(apriori: float, computed: float) -> float:
    """Percent difference: (computed - apriori) / apriori * 100."""
    if apriori == 0:
        return float("nan")
    return (computed - apriori) / apriori * 100.0


# ---------------------------------------------------------------------------
# Disk I/O helpers
# ---------------------------------------------------------------------------

_RESULTS_CSV_FIELDS = ["id", "name", "lat", "lng", "apriori_area_km2", "computed_area_km2", "pct_diff"]


def _outlet_dir(row_id: str) -> Path:
    return _output_dir / row_id


def _is_cached(row_id: str) -> bool:
    """True if this outlet was already delineated (watershed file exists)."""
    return (_outlet_dir(row_id) / "watershed.geojson").exists()


def _save_geojson(gdf, path: Path) -> None:
    if gdf is not None:
        path.write_text(gdf.to_json(), encoding="utf-8")


def _save_result(row: dict, computed_area: float, pct_diff: float,
                 watershed_gdf, rivers_gdf, outlets_gdf) -> None:
    """Write GeoJSON files and append a row to results.csv."""
    out = _outlet_dir(row["id"])
    out.mkdir(parents=True, exist_ok=True)

    _save_geojson(watershed_gdf, out / "watershed.geojson")
    _save_geojson(rivers_gdf,    out / "rivers.geojson")
    _save_geojson(outlets_gdf,   out / "outlets.geojson")

    # Append to the summary CSV (write header only if file is new)
    write_header = not _results_csv.exists()
    with open(_results_csv, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=_RESULTS_CSV_FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow({
            "id":               row["id"],
            "name":             row.get("name", ""),
            "lat":              row["lat"],
            "lng":              row["lng"],
            "apriori_area_km2": row["area"],
            "computed_area_km2": round(computed_area, 2),
            "pct_diff":         round(pct_diff, 2),
        })


def _load_cached_result(row: dict) -> dict:
    """Read computed_area and pct_diff from results.csv for a cached outlet."""
    if not _results_csv.exists():
        return {}
    with open(_results_csv, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["id"] == row["id"]:
                return {
                    "computed_area_km2": float(r["computed_area_km2"]),
                    "apriori_area_km2":  float(r["apriori_area_km2"]),
                    "pct_diff":          float(r["pct_diff"]),
                }
    return {}


# ---------------------------------------------------------------------------
# HTML template (single-file, all CSS + JS inline)
# ---------------------------------------------------------------------------

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8"/>
  <meta name="viewport" content="width=device-width, initial-scale=1.0"/>
  <title>Batch Watershed Delineation</title>

  <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
  <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>

  <style>
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }

    :root {
      --bg:        #f5f6f8;
      --panel:     #ffffff;
      --border:    #dde1e7;
      --accent:    #2a7fd4;
      --accent2:   #0e6fad;
      --text:      #1a1f2e;
      --muted:     #6b7894;
      --good:      #2e8f5e;
      --warn:      #c47c00;
      --bad:       #c0392b;
      --row-hover: rgba(42,127,212,0.05);
      --row-active:rgba(42,127,212,0.11);
    }

    html, body { height: 100%; }

    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: var(--bg);
      color: var(--text);
      display: flex;
      flex-direction: column;
      height: 100vh;
      overflow: hidden;
    }

    /* ── Header ── */
    header {
      background: var(--panel);
      border-bottom: 1px solid var(--border);
      padding: 10px 20px;
      flex-shrink: 0;
      display: flex;
      align-items: baseline;
      gap: 16px;
    }
    header h1 { font-size: 1.05rem; font-weight: 600; }
    header .meta { font-size: 0.78rem; color: var(--muted); }

    /* ── Progress bar ── */
    #progress-wrap {
      flex-shrink: 0;
      background: var(--border);
      height: 3px;
    }
    #progress-bar {
      height: 3px;
      background: var(--accent);
      width: 0%;
      transition: width 0.3s ease;
    }

    /* ── Status strip ── */
    #status-strip {
      flex-shrink: 0;
      background: #e8edf4;
      padding: 5px 20px;
      font-size: 0.75rem;
      color: var(--accent);
      display: flex;
      align-items: center;
      gap: 8px;
      min-height: 26px;
    }
    .spinner {
      display: none;
      width: 12px; height: 12px;
      border: 2px solid #2a7fd433;
      border-top-color: var(--accent);
      border-radius: 50%;
      animation: spin 0.7s linear infinite;
      flex-shrink: 0;
    }
    .spinning { display: inline-block !important; }
    @keyframes spin { to { transform: rotate(360deg); } }

    /* ── Main split layout ── */
    #main {
      flex: 1;
      display: flex;
      overflow: hidden;
    }

    /* ── Left panel: table ── */
    #table-panel {
      width: 52%;
      min-width: 440px;
      display: flex;
      flex-direction: column;
      border-right: 1px solid var(--border);
      overflow: hidden;
    }

    #table-scroll {
      flex: 1;
      overflow-y: auto;
      overflow-x: auto;
    }

    table {
      width: 100%;
      border-collapse: collapse;
      font-size: 0.8rem;
    }
    thead {
      position: sticky;
      top: 0;
      z-index: 2;
      background: var(--panel);
    }
    thead th {
      padding: 9px 10px;
      text-align: left;
      font-weight: 600;
      color: var(--muted);
      border-bottom: 1px solid var(--border);
      white-space: nowrap;
      font-size: 0.72rem;
      text-transform: uppercase;
      letter-spacing: 0.05em;
    }
    thead th.num { text-align: right; }

    tbody tr {
      border-bottom: 1px solid rgba(15,52,96,0.5);
      cursor: pointer;
      transition: background 0.12s;
    }
    tbody tr:hover   { background: var(--row-hover); }
    tbody tr.active  { background: var(--row-active); }
    tbody tr.pending { opacity: 0.55; cursor: default; }
    tbody tr.error-row td { color: var(--bad); }

    tbody td {
      padding: 8px 10px;
      vertical-align: middle;
      white-space: nowrap;
    }
    tbody td.num { text-align: right; font-variant-numeric: tabular-nums; }

    /* ID link */
    .id-link {
      color: var(--accent);
      text-decoration: none;
      font-weight: 600;
      font-family: "SF Mono", "Fira Code", monospace;
      font-size: 0.78rem;
    }
    .id-link:hover { text-decoration: underline; }
    tbody tr.pending .id-link { pointer-events: none; color: var(--muted); }

    /* Status badge */
    .badge {
      display: inline-flex;
      align-items: center;
      gap: 5px;
      font-size: 0.7rem;
      padding: 2px 7px;
      border-radius: 10px;
      white-space: nowrap;
    }
    .badge.pending  { background: rgba(136,153,170,0.15); color: var(--muted); }
    .badge.running  { background: rgba(74,168,255,0.15);  color: var(--accent); }
    .badge.done     { background: rgba(76,175,130,0.15);  color: var(--good); }
    .badge.failed   { background: rgba(224,82,82,0.15);   color: var(--bad); }
    .badge .dot { width:6px; height:6px; border-radius:50%; background:currentColor; flex-shrink:0; }
    .badge.running .dot {
      animation: pulse 1s ease-in-out infinite;
    }
    @keyframes pulse {
      0%,100% { opacity:1; } 50% { opacity:0.3; }
    }

    /* Perc. Diff. colouring */
    .pbias-good { color: var(--good); }
    .pbias-warn { color: var(--warn); }
    .pbias-bad  { color: var(--bad);  }

    /* Summary row */
    #summary-row td {
      padding: 9px 10px;
      font-size: 0.75rem;
      color: var(--muted);
      border-top: 2px solid var(--border);
      background: var(--panel);
      position: sticky;
      bottom: 0;
    }

    /* ── Right panel: map ── */
    #map-panel {
      flex: 1;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }

    #map-header {
      flex-shrink: 0;
      padding: 8px 14px;
      background: var(--panel);
      border-bottom: 1px solid var(--border);
      font-size: 0.8rem;
      color: var(--muted);
      min-height: 36px;
      display: flex;
      align-items: center;
      gap: 8px;
    }
    #map-header strong { color: var(--text); }

    #map { flex: 1; }

    /* ── Legend overlay ── */
    #legend {
      position: absolute;
      bottom: 28px; 
      right: 10px;
      z-index: 1000;
      background: rgba(255,255,255,0.92);
      border: 1px solid var(--border);
      border-radius: 7px;
      padding: 9px 13px;
      font-size: 0.72rem;
      backdrop-filter: blur(4px);
    }
    .legend-row { display:flex; align-items:center; gap:7px; margin-bottom:4px; }
    .legend-row:last-child { margin-bottom:0; }
    .sw { width:26px; height:11px; border-radius:2px; flex-shrink:0; }
    .sw.watershed { background:rgba(245, 10, 10,0.18); border:2px solid #f50a0a; }
    .sw.river     { background:var(--accent2); height:3px; }
    .sw.outlet-requested    { width:10px;height:10px;border-radius:50%;background:#c7830c;border:2px solid #fff;flex-shrink:0; }
    .sw.outlet-snapped    { width:11px;height:11px;border-radius:50%;background:#e05c1a;border:2px solid #fff;flex-shrink:0; }

    /* ── Error toast ── */
    #toast {
      position: fixed;
      top: 70px; left: 50%; transform: translateX(-50%);
      z-index: 9999;
      background: #c0392b;
      color: #fff;
      border-radius: 6px;
      padding: 7px 16px;
      font-size: 0.8rem;
      display: none;
      max-width: 90%;
      text-align: center;
      box-shadow: 0 4px 12px rgba(0,0,0,0.5);
    }
  </style>
</head>
<body>

<header>
  <h1>📜 Batch Watershed Delineation</h1>
  <span class="meta" id="header-meta">Loading…</span>
</header>

<div id="progress-wrap"><div id="progress-bar"></div></div>

<div id="status-strip">
  <div class="spinner" id="spinner"></div>
  <span id="status-text">Starting…</span>
</div>

<div id="main">

  <!-- Left: table -->
  <div id="table-panel">
    <div id="table-scroll">
      <table>
        <thead>
          <tr>
            <th>ID</th>
            <th>Name</th>
            <th class="num">Lat</th>
            <th class="num">Lng</th>
            <th class="num">A priori (km²)</th>
            <th class="num">Computed (km²)</th>
            <th class="num">Perc. Diff. (%)</th>
            <th>Status</th>
          </tr>
        </thead>
        <tbody id="tbody"></tbody>
        <tfoot>
          <tr id="summary-row" style="display:none">
            <td colspan="8" id="summary-text"></td>
          </tr>
        </tfoot>
      </table>
    </div>
  </div>

  <!-- Right: map -->
  <div id="map-panel">
    <div id="map-header">
      <span id="map-label">← Click an ID to view its watershed</span>
    </div>
    <div id="map"></div>
    <div id="legend">
      <div class="legend-row"><div class="sw watershed"></div> Watershed</div>
      <div class="legend-row"><div class="sw river"></div> Rivers</div>
      <div class="legend-row"><div class="sw outlet-requested"></div> Outlet - requested</div>
      <div class="legend-row"><div class="sw outlet-snapped"></div> Outlet - snapped</div>
    </div>
  </div>

</div>

<div id="toast"></div>

<script>
// ── Map init ──────────────────────────────────────────────────────────────
const map = L.map('map').setView([65, -19], 6);   // Iceland default view

L.tileLayer('https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png', {
  attribution: '© OpenStreetMap contributors © CARTO',
  subdomains: 'abcd', maxZoom: 19,
}).addTo(map);

let watershedLayer = null, riversLayer = null, outletLayer = null;

function clearMapLayers() {
  if (watershedLayer) { map.removeLayer(watershedLayer); watershedLayer = null; }
  if (riversLayer)    { map.removeLayer(riversLayer);    riversLayer    = null; }
  if (outletLayer)    { map.removeLayer(outletLayer);    outletLayer    = null; }
}

// ── In-memory store: id -> lightweight metadata (no GeoJSON) ─────────────
const store = {};

// ── UI helpers ────────────────────────────────────────────────────────────
const spinner    = document.getElementById('spinner');
const statusText = document.getElementById('status-text');
const progressEl = document.getElementById('progress-bar');
const toastEl    = document.getElementById('toast');
let toastTimer   = null;

function setStatus(msg, loading=false) {
  statusText.textContent = msg;
  spinner.classList.toggle('spinning', loading);
}

function setProgress(done, total) {
  progressEl.style.width = total ? `${(done/total)*100}%` : '0%';
}

function showToast(msg) {
  toastEl.textContent = msg;
  toastEl.style.display = 'block';
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => { toastEl.style.display = 'none'; }, 5000);
}

function pbiasCls(v) {
  const a = Math.abs(v);
  if (a <= 10) return 'pbias-good';
  if (a <= 25) return 'pbias-warn';
  return 'pbias-bad';
}

function fmtNum(v, decimals=1) {
  if (v == null || isNaN(v)) return '—';
  return v.toLocaleString(undefined, {minimumFractionDigits: decimals, maximumFractionDigits: decimals});
}

function fmtPbias(v) {
  if (v == null || isNaN(v)) return '—';
  const sign = v > 0 ? '+' : '';
  return `${sign}${v.toFixed(1)}%`;
}

// ── Table row building ────────────────────────────────────────────────────
function makeBadge(status) {
  const labels = { pending:'Pending', running:'Running', done:'Done', failed:'Failed' };
  return `<span class="badge ${status}"><span class="dot"></span>${labels[status]||status}</span>`;
}

function makeRow(row) {
  const tr = document.createElement('tr');
  tr.id = `row-${row.id}`;
  tr.className = 'pending';
  tr.innerHTML = `
    <td><a href="#" class="id-link" data-id="${row.id}">${row.id}</a></td>
    <td>${row.name || '—'}</td>
    <td class="num">${fmtNum(row.lat, 4)}</td>
    <td class="num">${fmtNum(row.lng, 4)}</td>
    <td class="num">${fmtNum(row.area, 0)}</td>
    <td class="num" id="computed-${row.id}">—</td>
    <td class="num" id="pbias-${row.id}">—</td>
    <td id="status-${row.id}">${makeBadge('pending')}</td>
  `;
  tr.querySelector('.id-link').addEventListener('click', (e) => {
    e.preventDefault();
    showOnMap(row.id, row.name);
  });
  return tr;
}

function updateRow(id, payload) {
  const tr = document.getElementById(`row-${id}`);
  if (!tr) return;

  if (payload.error) {
    tr.classList.replace('pending', 'error-row');
    document.getElementById(`status-${id}`).innerHTML  = makeBadge('failed');
    document.getElementById(`computed-${id}`).textContent = '—';
    document.getElementById(`pbias-${id}`).textContent = '—';
    tr.title = payload.error;
    return;
  }

  tr.classList.remove('pending');
  document.getElementById(`status-${id}`).innerHTML = makeBadge('done');

  const computed = payload.computed_area_km2;
  const apriori  = payload.apriori_area_km2;
  const pbias    = payload.pct_diff;

  document.getElementById(`computed-${id}`).textContent = fmtNum(computed, 0);

  const pbCell = document.getElementById(`pbias-${id}`);
  pbCell.textContent = fmtPbias(pbias);
  pbCell.className   = `num ${pbiasCls(pbias)}`;
}

function markRunning(id) {
  const el = document.getElementById(`status-${id}`);
  if (el) el.innerHTML = makeBadge('running');
}

// ── Map display ───────────────────────────────────────────────────────────
async function showOnMap(id, name) {
  const meta = store[id];
  if (!meta || meta.error) {
    showToast(meta?.error || 'No data available for this watershed.');
    return;
  }

  // Highlight active row
  document.querySelectorAll('tbody tr').forEach(r => r.classList.remove('active'));
  const tr = document.getElementById(`row-${id}`);
  if (tr) tr.classList.add('active');

  clearMapLayers();

  const label = document.getElementById('map-label');
  label.innerHTML = `<strong>${id}</strong>${name ? ' — ' + name : ''} <span style="color:var(--muted);font-size:0.72rem">loading…</span>`;

  let data;
  try {
    const resp = await fetch(`/geometry/${encodeURIComponent(id)}`);
    data = await resp.json();
    if (!resp.ok) { showToast(data.error || `Error loading ${id}`); return; }
  } catch (err) {
    showToast(`Failed to load geometry for ${id}`);
    console.error(err);
    return;
  }

  label.innerHTML = `<strong>${id}</strong>${name ? ' — ' + name : ''}`;

  if (data.watershed) {
    watershedLayer = L.geoJSON(data.watershed, {
      style: {
        color: '#f50a0a', weight: 2, opacity: 0.9,
        fillColor: '#f50a0a', fillOpacity: 0.18,
      }
    }).addTo(map);
    map.fitBounds(watershedLayer.getBounds(), { padding: [40, 40] });
  }

  if (data.rivers) {
    riversLayer = L.geoJSON(data.rivers, {
      style: { color: '#0e6fad', weight: 1.5, opacity: 0.85 }
    }).addTo(map);
  }

  if (data.outlets) {
    outletLayer = L.geoJSON(data.outlets, {
      pointToLayer: (feature, latlng) => {
        const snapped = feature.properties?.type === 'snapped';
        return L.circleMarker(latlng, {
          radius: snapped ? 6 : 5,
          fillColor: snapped ? '#e05c1a' : '#c47c00',
          color: '#fff', weight: 2, opacity: 1, fillOpacity: 0.95,
        });
      },
      onEachFeature: (feature, layer) => {
        if (feature.properties?.type) layer.bindTooltip(feature.properties.type);
      }
    }).addTo(map);
  }
}

// ── Streaming SSE consumer ────────────────────────────────────────────────
let totalRows  = 0;
let doneCount  = 0;

async function startStream() {
  const resp = await fetch('/stream');
  const reader = resp.body.getReader();
  const decoder = new TextDecoder();
  let buffer = '';

  while (true) {
    const { value, done } = await reader.read();
    if (done) break;
    buffer += decoder.decode(value, { stream: true });

    // SSE: split on double-newline
    const events = buffer.split('\n\n');
    buffer = events.pop();           // last chunk may be incomplete

    for (const block of events) {
      if (!block.trim()) continue;
      // Each block: "data: {json}"
      const line = block.replace(/^data:\s*/, '');
      try {
        const msg = JSON.parse(line);
        handleMessage(msg);
      } catch (e) {
        console.warn('SSE parse error', e, line);
      }
    }
  }
}

function handleMessage(msg) {
  switch (msg.type) {

    case 'init': {
      totalRows = msg.total;
      document.getElementById('header-meta').textContent =
        `${msg.csv} · ${totalRows} outlets`;
      setProgress(0, totalRows);
      setStatus(`Processing 0 / ${totalRows}…`, true);

      const tbody = document.getElementById('tbody');
      tbody.innerHTML = '';
      for (const row of msg.rows) {
        tbody.appendChild(makeRow(row));
      }
      break;
    }

    case 'start': {
      markRunning(msg.id);
      setStatus(`Delineating ${msg.id} (${doneCount+1} / ${totalRows})…`, true);
      break;
    }

    case 'cached':
    case 'result': {
      doneCount++;
      // Normalize to pct_diff key regardless of source
      const meta = { ...msg, pct_diff: msg.pct_diff ?? msg.pbias };
      store[msg.id] = meta;
      updateRow(msg.id, meta);
      setProgress(doneCount, totalRows);
      const verb = msg.type === 'cached' ? 'loaded' : 'delineated';
      setStatus(`${verb} ${msg.id} (${doneCount} / ${totalRows})…`, doneCount < totalRows);

      if (doneCount === totalRows) {
        setStatus(`Done — ${totalRows} watersheds processed.`, false);
        showSummary();
      }
      break;
    }

    case 'error': {
      doneCount++;
      store[msg.id] = msg;
      updateRow(msg.id, msg);
      setProgress(doneCount, totalRows);
      break;
    }
  }
}

function showSummary() {
  const successes = Object.values(store).filter(r => !r.error);
  if (!successes.length) return;

  const pbVals = successes.map(r => r.pct_diff).filter(v => v != null && !isNaN(v));
  const mean   = pbVals.reduce((a,b) => a+b, 0) / pbVals.length;
  const errors = Object.values(store).filter(r => r.error).length;

  const row = document.getElementById('summary-row');
  row.style.display = '';
  document.getElementById('summary-text').textContent =
    `${successes.length} succeeded · ${errors} failed · ` +
    `Mean Perc. Diff. ${mean >= 0 ? '+' : ''}${mean.toFixed(1)}% `;
}

// ── Kick off ──────────────────────────────────────────────────────────────
startStream().catch(err => {
  setStatus('Stream error — see console', false);
  console.error(err);
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


@app.route("/geometry/<row_id>")
def geometry(row_id: str) -> Response:
    """
    Return the saved GeoJSON for a single outlet on demand.
    Called by the browser only when the user clicks a row.
    """
    out = _outlet_dir(row_id)
    if not out.exists():
        return jsonify({"error": f"No data on disk for '{row_id}'."}), 404

    def read(name):
        p = out / name
        return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None

    return jsonify({
        "watershed": read("watershed.geojson"),
        "rivers":    read("rivers.geojson"),
        "outlets":   read("outlets.geojson"),
    })


@app.route("/stream")
def stream() -> Response:
    """
    Server-Sent Events endpoint. Streams lightweight JSON messages as
    delineation progresses — GeoJSON is never included in the stream.

    Message types:
        init   — sent once at the start, includes the full row list for table init
        start  — sent just before delineating each watershed
        result — sent after a successful delineation (metadata only, no geometry)
        cached — sent for outlets already on disk from a previous run
        error  — sent if a watershed fails
    """

    def generate():
        rows = _rows

        # ── init message ────────────────────────────────────────────────
        init_payload = {
            "type":  "init",
            "csv":   _csv_path.name,
            "total": len(rows),
            "rows":  rows,
        }
        yield f"data: {json.dumps(init_payload)}\n\n"

        config = DelineatorConfig(
            rivers=True,
            outlets=True,
        )

        for row in rows:
            row_id  = row["id"]
            lat     = row["lat"]
            lng     = row["lng"]
            apriori = row["area"]

            # ── skip if already on disk ──────────────────────────────────
            if _is_cached(row_id):
                cached = _load_cached_result(row)
                payload = {
                    "type":             "cached",
                    "id":               row_id,
                    "computed_area_km2": cached.get("computed_area_km2"),
                    "apriori_area_km2":  apriori,
                    "pct_diff":         cached.get("pct_diff"),
                }
                logger.info(f"Skipping {row_id} (cached)")
                yield f"data: {json.dumps(payload)}\n\n"
                continue

            # ── start message ────────────────────────────────────────────
            yield f"data: {json.dumps({'type': 'start', 'id': row_id})}\n\n"

            try:
                watershed_gdf, rivers_gdf, outlets_gdf = delineate(lat, lng, config=config)

                if watershed_gdf is None:
                    raise ValueError("Delineation returned None — point may be over ocean or outside data coverage.")

                computed_area = _compute_area_km2(watershed_gdf)
                pct_diff      = _pct_diff(apriori, computed_area)

                _save_result(row, computed_area, pct_diff, watershed_gdf, rivers_gdf, outlets_gdf)

                payload = {
                    "type":             "result",
                    "id":               row_id,
                    "computed_area_km2": computed_area,
                    "apriori_area_km2":  apriori,
                    "pct_diff":         pct_diff,
                }

            except Exception as exc:
                logger.error(f"Failed to delineate {row_id}: {exc}")
                logger.debug(traceback.format_exc())
                payload = {
                    "type":  "error",
                    "id":    row_id,
                    "error": str(exc),
                }

            yield f"data: {json.dumps(payload)}\n\n"

    return Response(
        stream_with_context(generate()),
        mimetype="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "X-Accel-Buffering": "no",
        },
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    global _csv_path, _rows, _output_dir, _results_csv

    if len(sys.argv) > 1:
        _csv_path = Path(sys.argv[1])

    if not _csv_path.exists():
        print(f"Error: CSV file not found: {_csv_path}", file=sys.stderr)
        sys.exit(1)

    try:
        _rows = _load_csv(_csv_path)
    except Exception as e:
        print(f"Error reading CSV: {e}", file=sys.stderr)
        sys.exit(1)

    # Output directory sits next to the CSV
    _output_dir = _csv_path.parent / "output"
    _output_dir.mkdir(parents=True, exist_ok=True)
    _results_csv = _output_dir / "results.csv"

    # Report how many outlets are already cached
    cached = sum(1 for r in _rows if _is_cached(r["id"]))

    print("\n  Batch Watershed Delineator")
    print("  ─────────────────────────────────────────")
    print(f"  CSV:    {_csv_path} ({len(_rows)} outlets)")
    print(f"  Output: {_output_dir.resolve()}")
    if cached:
        print(f"  Cached: {cached} outlet(s) will be skipped")
    print("  Open http://localhost:5001 in your browser")
    print("  Delineation starts automatically.\n")

    app.run(debug=False, port=5001, threaded=True)


if __name__ == "__main__":
    main()