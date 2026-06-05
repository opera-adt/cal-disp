#!/usr/bin/env python
"""Build the OPERA DISP-S1 frame -> UNR grid-points parquet database.

For every North-America OPERA frame it finds the UNR gridded-timeseries points
that fall within the frame (plus a margin) and records the per-plate download URLs.

The output filename embeds the production date so each build is traceable, e.g.
    opera_disp_s1_frame_unr_points_20260604.parquet

Run inside the `my-cal-env` conda environment:
    conda activate my-cal-env

Examples
--------
# Build the database
    python build_unr_frame_db.py --outdir ../configs/data

# Build + write an interactive HTML map of all frames coloured by point count
    python build_unr_frame_db.py --html-map

# Inspect one frame (grid ids + URLs) and map just that frame + its UNR points,
# reusing an existing parquet (no rebuild / network)
    python build_unr_frame_db.py --from-parquet ../configs/data/opera_disp_s1_frame_unr_points_20260604.parquet \
        --frame-id 831 --plate IGS20 --html-map
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import pandas as pd
from opera_utils import get_frame_geojson
from shapely.geometry import mapping

from cal_disp.download._stage_unr import (
    GRID_BASE_URL,
    download_lookup_table,
    get_frame_grid_points,
    load_lookup_table,
)

# UNR plates to build download-URL columns for.
PLATES = ("IGS20", "NA", "PA")


def get_unr_url(grid_id: int, plate: str, version: str) -> str:
    """Build the UNR .tenv8 download URL for one grid point + plate."""
    filename = f"{plate}/{grid_id:06d}_{plate}.tenv8"
    return f"{GRID_BASE_URL.format(version=version)}/{filename}"


def load_grid_gdf(work_dir: Path, version: str) -> gpd.GeoDataFrame:
    """Download (if needed) and load the UNR grid lookup table as points.

    The GeoDataFrame is indexed by grid_id (so ``.loc[ids]`` selects points).
    """
    lookup_path = download_lookup_table(work_dir, version=version)
    lookup = load_lookup_table(lookup_path)
    return gpd.GeoDataFrame(
        lookup,
        geometry=gpd.points_from_xy(x=lookup.lon, y=lookup.lat),
        crs="EPSG:4326",
    )


def build_database(
    work_dir: Path,
    version: str,
    margin_deg: float,
    grid_gdf: gpd.GeoDataFrame | None = None,
) -> gpd.GeoDataFrame:
    """Compute the frame -> UNR grid-points GeoDataFrame (indexed by frame_id)."""
    if grid_gdf is None:
        grid_gdf = load_grid_gdf(work_dir, version)

    # OPERA frames restricted to North America
    frame_gdf = get_frame_geojson(as_geodataframe=True)
    opera_gdf = frame_gdf[frame_gdf.is_north_america].copy()

    # UNR grid points per frame
    def grid_points_for_row(row):
        temp_gdf = gpd.GeoDataFrame([row], geometry="geometry", crs=opera_gdf.crs)
        # get_frame_grid_points returns (grid_ids, grid_gdf_filtered); keep the ids.
        grid_ids, _ = get_frame_grid_points(temp_gdf, grid_gdf, margin_deg=margin_deg)
        return grid_ids

    opera_gdf["unr_grid_points"] = opera_gdf.apply(grid_points_for_row, axis=1)
    opera_gdf["unr_grid_count"] = opera_gdf["unr_grid_points"].apply(len)

    # Per-plate download URLs
    for plate in PLATES:
        opera_gdf[f"unr_urls_{plate}"] = opera_gdf["unr_grid_points"].apply(
            lambda ids, p=plate: [get_unr_url(g, p, version) for g in ids]
        )

    return opera_gdf


def report_frame(opera_gdf: gpd.GeoDataFrame, frame_id: int, plates: list[str]) -> None:
    """Print the UNR grid ids and download URLs for a single frame."""
    if frame_id not in opera_gdf.index:
        raise SystemExit(f"frame_id {frame_id} not found among North-America frames")
    row = opera_gdf.loc[frame_id]
    ids = list(row["unr_grid_points"])
    print(f"\nFrame {frame_id}: {len(ids)} UNR grid point(s)")
    print(f"  grid ids: {ids}")
    for plate in plates:
        urls = list(row[f"unr_urls_{plate}"])
        print(f"  [{plate}] {len(urls)} url(s):")
        for u in urls:
            print(f"    {u}")


def make_html_map(
    opera_gdf: gpd.GeoDataFrame,
    out_path: Path,
    frame_id: int | None = None,
    grid_gdf: gpd.GeoDataFrame | None = None,
) -> Path:
    """Write an interactive folium HTML map.

    With ``frame_id`` -> that single frame plus its UNR points.
    Otherwise -> all frames coloured by ``unr_grid_count``.
    """
    import folium  # required by GeoDataFrame.explore

    if frame_id is not None:
        if frame_id not in opera_gdf.index:
            raise SystemExit(f"frame_id {frame_id} not found among North-America frames")
        frame = opera_gdf.loc[[frame_id], ["geometry", "unr_grid_count"]]
        m = frame.explore(
            color="red",
            name=f"frame {frame_id}",
            tooltip=["unr_grid_count"],
            style_kwds={"fillOpacity": 0.1, "weight": 2},
        )
        ids = list(opera_gdf.loc[frame_id, "unr_grid_points"])
        if grid_gdf is not None and ids:
            pts = grid_gdf.loc[grid_gdf.index.isin(ids)]
            if len(pts):
                pts[["geometry"]].explore(
                    m=m, color="blue", name="UNR points", marker_kwds={"radius": 4}
                )
    else:
        # Drop heavy list columns (ids/urls) before serialising to keep the file small.
        slim = opera_gdf[["geometry", "unr_grid_count"]]
        m = slim.explore(
            column="unr_grid_count",
            cmap="viridis",
            legend=True,
            name="OPERA frames",
            tooltip=["unr_grid_count"],
            style_kwds={"fillOpacity": 0.5},
        )

    folium.LayerControl().add_to(m)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    m.save(str(out_path))
    return out_path


STANDALONE_TEMPLATE = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<title>OPERA DISP-S1 frame &rarr; UNR grid viewer</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { font-family: 'Segoe UI', Arial, sans-serif; background: #1a1a2e; color: #e0e0e0;
         height: 100vh; display: flex; flex-direction: column; overflow: hidden; }
  #top-bar { background: #16213e; border-bottom: 2px solid #0f3460; padding: 8px 16px;
             display: flex; align-items: center; gap: 16px; flex-shrink: 0; z-index: 1000; flex-wrap: wrap; }
  #top-bar h1 { font-size: 15px; color: #e94560; letter-spacing: 1px; white-space: nowrap; }
  #top-bar .stats { font-size: 12px; color: #888; }
  #search-box { padding: 5px 10px; border-radius: 4px; border: 1px solid #0f3460; background: #0d0d1a;
                color: #e0e0e0; font-size: 13px; width: 220px; }
  #search-box::placeholder { color: #555; }
  .nav-btn { background: #0f3460; color: #e0e0e0; border: 1px solid #1a4a80; padding: 4px 12px;
             border-radius: 4px; cursor: pointer; font-size: 13px; }
  .nav-btn:hover { background: #1a4a80; }
  #main { display: flex; flex: 1; overflow: hidden; }
  #map-panel { width: 50%; display: flex; flex-direction: column; border-right: 2px solid #0f3460; flex-shrink: 0; }
  #map { flex: 1; }
  #map-info { background: #16213e; padding: 6px 12px; font-size: 12px; color: #aaa;
              border-top: 1px solid #0f3460; height: 34px; display: flex; align-items: center; gap: 8px; }
  #map-info .hovered { color: #e94560; font-weight: bold; }
  #frame-list { background: #16213e; border-top: 1px solid #0f3460; overflow-y: auto; flex-shrink: 0; height: 150px; }
  #frame-list-inner { display: flex; flex-wrap: wrap; gap: 3px; padding: 6px; }
  .frame-pill { background: #0f3460; color: #aaa; padding: 2px 8px; border-radius: 10px; font-size: 11px;
                cursor: pointer; white-space: nowrap; border: 1px solid transparent; }
  .frame-pill:hover { background: #1a4a80; color: #fff; }
  .frame-pill.active { background: #e94560; color: #fff; border-color: #ff6b6b; }
  #detail-panel { flex: 1; display: flex; flex-direction: column; overflow: hidden; }
  #detail-header { background: #16213e; padding: 8px 14px; border-bottom: 1px solid #0f3460; flex-shrink: 0;
                   display: flex; align-items: center; gap: 10px; flex-wrap: wrap; }
  #frame-title { font-size: 14px; font-weight: bold; color: #e94560; min-width: 120px; }
  .tag { background: #0f3460; color: #7ec8e3; padding: 2px 8px; border-radius: 10px; font-size: 11px; }
  #pair-counter { font-size: 12px; color: #888; margin-left: auto; }
  #detail-content { flex: 1; overflow-y: auto; padding: 12px 14px; font-size: 13px; }
  .block { margin-bottom: 12px; }
  .block-h { font-weight: 600; margin-bottom: 4px; }
  .ids { font-family: ui-monospace, monospace; word-break: break-word; background: #0d0d1a;
         border: 1px solid #0f3460; padding: 6px; border-radius: 4px; }
  details { margin: 8px 0; border: 1px solid #0f3460; border-radius: 4px; padding: 4px 8px; }
  details summary { cursor: pointer; font-weight: 600; }
  .urls { list-style: none; margin: 6px 0 0; padding: 0; }
  .urls li { margin: 2px 0; }
  .urls a { color: #7ec8e3; font-family: ui-monospace, monospace; font-size: 12px; word-break: break-all; }
  .copybtn { background: #0f3460; color: #7ec8e3; border: 1px solid #1a4a80; margin-left: 8px;
             font-size: 11px; padding: 1px 6px; border-radius: 3px; cursor: pointer; }
  .legend { background:#16213e; padding:6px 8px; line-height:1.4; border-radius:4px; color:#ddd;
            box-shadow:0 1px 4px rgba(0,0,0,.4); font-size:12px; }
  .legend i { display:inline-block; width:14px; height:14px; margin-right:6px; vertical-align:-2px; }
</style>
</head>
<body>
<div id="top-bar">
  <h1>OPERA DISP-S1 &nbsp;&#8212;&nbsp; frame &rarr; UNR grid</h1>
  <input id="search-box" placeholder="Search frame_id (e.g. 831)" oninput="filterFrames(this.value)" />
  <span class="stats" id="global-stats"></span>
  <div style="margin-left:auto; display:flex; gap:6px; align-items:center;">
    <span style="font-size:11px;color:#666;">Frame:</span>
    <button class="nav-btn" onclick="goFramePrev()" title="Previous (,)">&#8249; Prev</button>
    <button class="nav-btn" onclick="goFrameNext()" title="Next (.)">Next &#8250;</button>
  </div>
</div>
<div id="main">
  <div id="map-panel">
    <div id="map"></div>
    <div id="map-info"><span>Click a frame to view its UNR grid ids &amp; URLs</span><span class="hovered" id="hover-label"></span></div>
    <div id="frame-list"><div id="frame-list-inner"></div></div>
  </div>
  <div id="detail-panel">
    <div id="detail-header">
      <span id="frame-title">No frame selected</span>
      <span class="tag" id="tag-pass"></span>
      <span class="tag" id="tag-count"></span>
      <span id="pair-counter"></span>
    </div>
    <div id="detail-content"><p style="color:#667;">Select a frame from the map or list.</p></div>
  </div>
</div>
<script>
const DATA = __DATA__;
const POINTS = __POINTS__;
const GRID_BASE = "__GRID_BASE__";
const PLATES = __PLATES__;
const VERSION = "__VERSION__";

const feats = DATA.features.slice().sort((a,b)=>a.properties.frame_id-b.properties.frame_id);
const byId = {}; feats.forEach(f=>byId[f.properties.frame_id]=f);
const orderedIds = feats.map(f=>f.properties.frame_id);
let maxCount = 0; feats.forEach(f=>maxCount=Math.max(maxCount,f.properties.count));
let selected = null;

function colorFor(c){ const s=['#440154','#3b528b','#21918c','#5ec962','#fde725'];
  if(!maxCount) return s[0]; return s[Math.min(s.length-1,Math.round(Math.min(1,c/maxCount)*(s.length-1)))]; }
function passColor(p){ return p==='ASCENDING' ? '#00e5ff' : p==='DESCENDING' ? '#ff9800' : '#aaaaaa'; }
function styleFor(p){ return {color:passColor(p.orbit_pass), weight:1.2, fillColor:colorFor(p.count), fillOpacity:0.45}; }
function urlFor(id,plate){ return GRID_BASE+'/'+plate+'/'+String(id).padStart(6,'0')+'_'+plate+'.tenv8'; }

// base layers (satellite default)
const sat = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
  {maxZoom:19, attribution:'Tiles &copy; Esri'});
const dark = L.tileLayer('https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}.png',
  {maxZoom:19, attribution:'&copy; OpenStreetMap &copy; CARTO'});
const streets = L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',
  {maxZoom:19, attribution:'&copy; OpenStreetMap'});

const map = L.map('map',{worldCopyJump:true, layers:[sat]}).setView([45,-110],3);

// frame layers split by orbit pass, toggleable + colour-coded
const layerById = {};
function mkPassLayer(passVal){
  return L.geoJSON({type:'FeatureCollection', features: feats.filter(f=>f.properties.orbit_pass===passVal)}, {
    style: f=>styleFor(f.properties),
    onEachFeature:(f,l)=>{ const id=f.properties.frame_id; layerById[id]=l;
      l.on('click', ()=>selectFrame(id));
      l.on('mouseover', ()=>{ document.getElementById('hover-label').textContent =
        'frame '+id+' · '+(f.properties.orbit_pass||'—')+' · '+f.properties.count+' pts'; }); }
  }); }
const ascLayer = mkPassLayer('ASCENDING').addTo(map);
const descLayer = mkPassLayer('DESCENDING').addTo(map);
const pointLayer = L.layerGroup().addTo(map);

L.control.layers(
  {Satellite:sat, Dark:dark, Streets:streets},
  {'Ascending frames':ascLayer, 'Descending frames':descLayer, 'UNR grid points':pointLayer},
  {collapsed:false}
).addTo(map);

try { map.fitBounds(L.featureGroup([ascLayer,descLayer]).getBounds(), {padding:[20,20]}); } catch(e) {}

const legend = L.control({position:'bottomright'});
legend.onAdd = () => { const d=L.DomUtil.create('div','legend');
  d.innerHTML = '<b>UNR pts / frame</b><br>' + [0,0.25,0.5,0.75,1].map(t=>{
      const c=Math.round(t*maxCount); return '<i style="background:'+colorFor(c)+'"></i>'+c; }).join('<br>')
    + '<br><b>Orbit pass</b><br><i style="background:#00e5ff"></i>Ascending'
    + '<br><i style="background:#ff9800"></i>Descending';
  return d; };
legend.addTo(map);

// grid-point markers for the selected frame
function showPoints(ids){ pointLayer.clearLayers();
  ids.forEach(i=>{ const c=POINTS[i]; if(!c) return;
    L.circleMarker([c[1],c[0]], {radius:3, color:'#ffd54f', weight:1, fillColor:'#ffd54f', fillOpacity:0.9})
      .bindTooltip('grid '+i).addTo(pointLayer); }); }

const listInner = document.getElementById('frame-list-inner');
function buildList(filter){ listInner.innerHTML=''; const q=(filter||'').toLowerCase();
  feats.forEach(f=>{ const id=f.properties.frame_id; if(q && !(''+id).includes(q)) return;
    const pill=document.createElement('div'); pill.className='frame-pill'+(id===selected?' active':'');
    pill.textContent=id; pill.onclick=()=>selectFrame(id); listInner.appendChild(pill); }); }
buildList();
function filterFrames(v){ buildList(v); }

document.getElementById('global-stats').textContent =
  feats.length+' frames · max '+maxCount+' pts/frame · UNR v'+VERSION;

function selectFrame(id){ if(byId[id]===undefined) return;
  if(selected!==null && layerById[selected]) layerById[selected].setStyle(styleFor(byId[selected].properties));
  selected=id; const l=layerById[id];
  if(l){ l.setStyle({color:'#ffffff', weight:3, fillOpacity:0.15}); l.bringToFront();
    try{ map.fitBounds(l.getBounds(),{maxZoom:8,padding:[40,40]}); }catch(e){} }
  showPoints(byId[id].properties.ids);
  renderDetails(id); buildList(document.getElementById('search-box').value); }

function renderDetails(id){ const p=byId[id].properties;
  document.getElementById('frame-title').textContent='Frame '+id;
  document.getElementById('tag-pass').textContent=p.orbit_pass||'—';
  document.getElementById('tag-count').textContent=p.count+' UNR pts';
  document.getElementById('pair-counter').textContent=(orderedIds.indexOf(id)+1)+' / '+orderedIds.length;
  let html='<div class="block"><div class="block-h">Grid ids '
    +'<button class="copybtn" data-copy="'+p.ids.join(',')+'">copy</button></div>'
    +'<div class="ids">'+(p.ids.length?p.ids.join(', '):'<i>none</i>')+'</div></div>';
  PLATES.forEach(plate=>{ const urls=p.ids.map(i=>urlFor(i,plate));
    html+='<details><summary>'+plate+' URLs ('+urls.length+')'
      +(urls.length?' <button class="copybtn" data-copy="'+urls.join('&#10;')+'">copy</button>':'')
      +'</summary><ul class="urls">'
      +urls.map(u=>'<li><a href="'+u+'" target="_blank" rel="noopener">'+u+'</a></li>').join('')
      +'</ul></details>'; });
  const c=document.getElementById('detail-content'); c.innerHTML=html;
  c.querySelectorAll('.copybtn').forEach(b=>b.onclick=()=>{
    navigator.clipboard.writeText(b.dataset.copy.replace(/&#10;/g,'\n'));
    b.textContent='copied'; setTimeout(()=>b.textContent='copy',1200); }); }

function goFrameNext(){ const i=orderedIds.indexOf(selected);
  selectFrame(orderedIds[i<0?0:Math.min(orderedIds.length-1,i+1)]); }
function goFramePrev(){ const i=orderedIds.indexOf(selected);
  selectFrame(orderedIds[i<0?0:Math.max(0,i-1)]); }
document.addEventListener('keydown',e=>{ if(e.target.id==='search-box') return;
  if(e.key===','||e.key==='ArrowLeft') goFramePrev();
  if(e.key==='.'||e.key==='ArrowRight') goFrameNext(); });
</script>
<!-- generated __GENERATED__ -->
</body>
</html>
'''


def _round_geom(geom, ndigits: int = 5) -> dict:
    """GeoJSON mapping of a geometry with coordinates rounded to shrink the file."""
    def rnd(c):
        if isinstance(c, (list, tuple)):
            if c and isinstance(c[0], (int, float)):
                return [round(float(c[0]), ndigits), round(float(c[1]), ndigits)]
            return [rnd(x) for x in c]
        return c

    gj = dict(mapping(geom))
    gj["coordinates"] = rnd(gj["coordinates"])
    return gj


def write_standalone_html(
    opera_gdf: gpd.GeoDataFrame,
    out_path: Path,
    grid_base: str,
    version: str,
    grid_gdf: gpd.GeoDataFrame | None = None,
) -> Path:
    """Write a self-contained interactive HTML viewer with the data embedded.

    URLs are not embedded; they are rebuilt in-browser from grid id + plate +
    ``grid_base`` to keep the file small. When ``grid_gdf`` is given, the lon/lat
    of every referenced grid point is embedded so the selected frame's points can
    be drawn on the map.
    """
    feats = []
    used: set[int] = set()
    for fid, row in opera_gdf.iterrows():
        pts = row["unr_grid_points"] if row["unr_grid_points"] is not None else []
        ids = [int(g) for g in pts]
        used.update(ids)
        op = row.get("orbit_pass")
        feats.append(
            {
                "type": "Feature",
                "properties": {
                    "frame_id": int(fid),
                    "orbit_pass": None if op is None or pd.isna(op) else str(op),
                    "count": len(ids),
                    "ids": ids,
                },
                "geometry": _round_geom(row["geometry"]),
            }
        )

    # id -> [lon, lat] for every referenced point (shared across frames -> compact)
    points: dict[int, list[float]] = {}
    if grid_gdf is not None and used:
        sub = grid_gdf.loc[grid_gdf.index.isin(used), ["lon", "lat"]]
        for i, lon, lat in sub.itertuples():
            points[int(i)] = [round(float(lon), 5), round(float(lat), 5)]

    data = json.dumps({"type": "FeatureCollection", "features": feats}, separators=(",", ":"))
    data = data.replace("</", "<\\/")  # keep the </script> tag safe

    html = (
        STANDALONE_TEMPLATE.replace("__DATA__", data)
        .replace("__POINTS__", json.dumps(points, separators=(",", ":")))
        .replace("__GRID_BASE__", grid_base)
        .replace("__PLATES__", json.dumps(list(PLATES)))
        .replace("__VERSION__", version)
        .replace("__GENERATED__", datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"))
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html)
    return out_path


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=Path("../configs/data"),
        help="Directory to write the parquet (and default HTML map) into.",
    )
    p.add_argument(
        "--work-dir",
        type=Path,
        default=Path.cwd(),
        help="Scratch dir for the downloaded UNR lookup table.",
    )
    p.add_argument("--version", default="0.3", help="UNR gridded data version (default: 0.3).")
    p.add_argument("--margin-deg", type=float, default=0.5, help="Frame buffer in degrees (default: 0.5).")
    p.add_argument(
        "--prod-date",
        default=datetime.now(timezone.utc).strftime("%Y%m%d"),
        help="Production date stamped into the filename (YYYYMMDD; default: today UTC).",
    )

    # Inspection / visualization (all optional)
    p.add_argument(
        "--from-parquet",
        type=Path,
        default=None,
        help="Load an existing parquet instead of rebuilding (fast for inspection/maps).",
    )
    p.add_argument(
        "--frame-id",
        type=int,
        default=None,
        help="Report UNR grid ids + URLs for this frame_id (and focus the HTML map on it).",
    )
    p.add_argument(
        "--plate",
        choices=[*PLATES, "all"],
        default="all",
        help="Which plate's URLs to report for --frame-id (default: all).",
    )
    p.add_argument(
        "--html-map",
        action="store_true",
        help="Also write an interactive HTML map (folium).",
    )
    p.add_argument(
        "--map-out",
        type=Path,
        default=None,
        help="HTML map output path (default: alongside the parquet).",
    )
    p.add_argument(
        "--standalone-html",
        action="store_true",
        help="Write a self-contained interactive viewer (data embedded; no server needed).",
    )
    p.add_argument(
        "--standalone-out",
        type=Path,
        default=None,
        help="Standalone viewer output path (default: alongside the parquet).",
    )
    args = p.parse_args()

    plates = list(PLATES) if args.plate == "all" else [args.plate]

    # 1. Get the database: load existing, or build (and persist) a new one.
    if args.from_parquet is not None:
        print(f"Loading existing DB: {args.from_parquet}")
        opera_gdf = gpd.read_parquet(args.from_parquet)
    else:
        print(f"Building UNR frame DB (version={args.version}, margin={args.margin_deg} deg) ...")
        opera_gdf = build_database(args.work_dir, args.version, args.margin_deg)
        args.outdir.mkdir(parents=True, exist_ok=True)
        out_path = args.outdir / f"opera_disp_s1_frame_unr_points_{args.prod_date}.parquet"
        opera_gdf.to_parquet(out_path)
        print(f"Wrote {len(opera_gdf)} frames -> {out_path}")
        print(opera_gdf["unr_grid_count"].describe().to_string())

    # 2. Optional: report a single frame's grid ids + URLs.
    if args.frame_id is not None:
        report_frame(opera_gdf, args.frame_id, plates)

    # 3. Optional: write an HTML map.
    if args.html_map:
        if args.map_out is not None:
            map_path = args.map_out
        else:
            base = (
                args.from_parquet.stem
                if args.from_parquet is not None
                else f"opera_disp_s1_frame_unr_points_{args.prod_date}"
            )
            suffix = f"_frame{args.frame_id}" if args.frame_id is not None else "_map"
            map_dir = args.from_parquet.parent if args.from_parquet is not None else args.outdir
            map_path = map_dir / f"{base}{suffix}.html"

        # Grid point geometries are only needed to plot a selected frame's points.
        grid_gdf = load_grid_gdf(args.work_dir, args.version) if args.frame_id is not None else None
        make_html_map(opera_gdf, map_path, frame_id=args.frame_id, grid_gdf=grid_gdf)
        print(f"Wrote HTML map -> {map_path}")

    # 4. Optional: write a self-contained interactive viewer (data embedded).
    if args.standalone_html:
        if args.standalone_out is not None:
            sa_path = args.standalone_out
        else:
            base = (
                args.from_parquet.stem
                if args.from_parquet is not None
                else f"opera_disp_s1_frame_unr_points_{args.prod_date}"
            )
            sa_dir = args.from_parquet.parent if args.from_parquet is not None else args.outdir
            sa_path = sa_dir / f"{base}_viewer.html"
        grid_base = GRID_BASE_URL.format(version=args.version)
        # Grid point lon/lat are needed to draw the selected frame's points.
        sa_grid_gdf = load_grid_gdf(args.work_dir, args.version)
        write_standalone_html(opera_gdf, sa_path, grid_base, args.version, sa_grid_gdf)
        size_mb = sa_path.stat().st_size / 1e6
        print(f"Wrote standalone viewer -> {sa_path}  ({size_mb:.1f} MB)")


if __name__ == "__main__":
    main()
