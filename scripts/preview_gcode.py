"""Previsualisation 3D du G-code genere par gen_l_overhang_gcode.py.

Parse le fichier G-code (format machine Apex : G00 Y<C> X<R> Z<Z> B<B> F
[E<e>]), reconvertit chaque point (R,C,Z) -> (X,Y,Z) cartesien, distingue
les segments d'EXTRUSION (E present et croissant) des trajets A VIDE, et
colore par PHASE (couches pleines / pilier / coude / surplomb, detectees
via les commentaires "; --- ... ---" ecrits par gen_l_overhang_gcode.py).

Sortie : un fichier HTML autonome avec un viewer 3D (Canvas, rotation a la
souris, molette = zoom) -- aucune dependance externe (pas de matplotlib,
pas de CDN), ouvrable direct dans un navigateur ou publiable en artefact.

Usage : python scripts/preview_gcode.py [chemin_gcode] [chemin_html_sortie]
Defauts : out/l_overhang.gcode -> out/l_overhang_preview.html
"""

from __future__ import annotations

import json
import math
import os
import re
import sys

LINE_RE = re.compile(
    r"^G00"
    r"(?: Y(?P<c>-?\d+\.\d+))?"
    r"(?: X(?P<r>-?\d+\.\d+))?"
    r"(?: Z(?P<z>-?\d+\.\d+))?"
    r"(?: B(?P<b>-?\d+\.\d+))?"
    r"(?: E(?P<e>-?\d+\.\d+))?"
)

PHASE_MARKERS = [
    ("couches de base", "pleines"),
    ("reste du pilier", "pilier"),
    ("coude", "coude"),
    ("surplomb", "surplomb"),
]


def parse_gcode(path: str) -> list[dict]:
    """Retourne une liste de points {x,y,z,extruding,phase}."""
    points: list[dict] = []
    r = c = z = b = 0.0
    e_prev = 0.0
    phase = "pleines"
    have_pos = False

    with open(path, "r", encoding="ascii") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(";"):
                for needle, label in PHASE_MARKERS:
                    if needle in line:
                        phase = label
                continue
            if not line.startswith("G00"):
                continue
            m = LINE_RE.match(line)
            if not m:
                continue
            if m.group("c") is not None:
                c = float(m.group("c"))
            if m.group("r") is not None:
                r = float(m.group("r"))
            if m.group("z") is not None:
                z = float(m.group("z"))
            if m.group("b") is not None:
                b = float(m.group("b"))
            e_str = m.group("e")

            c_rad = math.radians(c)
            x = r * math.cos(c_rad)
            y = r * math.sin(c_rad)

            extruding = False
            if e_str is not None:
                e_val = float(e_str)
                if e_val > e_prev + 1e-9:
                    extruding = True
                e_prev = e_val

            points.append({"x": x, "y": y, "z": z, "e": extruding, "p": phase})
            have_pos = True

    if not have_pos:
        raise ValueError(f"aucun point G00 trouve dans {path}")
    return points


def decimate(points: list[dict], target: int) -> list[dict]:
    """Sous-echantillonne uniformement a ~target points, sans casser les
    transitions de phase (chaque point de changement de phase est garde)."""
    n = len(points)
    if n <= target:
        return points
    stride = max(1, n // target)
    out = []
    prev_phase = None
    for i in range(n):
        if points[i]["p"] != prev_phase or i % stride == 0 or i == n - 1:
            out.append(points[i])
            prev_phase = points[i]["p"]
    return out


PHASE_LABELS = {
    "pleines": "Couches pleines (adherence)",
    "pilier": "Pilier (spirale continue)",
    "coude": "Coude (reorientation)",
    "surplomb": "Surplomb (B=45 constant)",
}
PHASE_COLORS = {
    "pleines": "#3ddc97",
    "pilier": "#8fa3bf",
    "coude": "#ffce54",
    "surplomb": "#ff8a3d",
}


def build_html(points: list[dict], src_name: str) -> str:
    n_total = len(points)
    n_extrude = sum(1 for p in points if p["e"])
    bounds = {
        "x": (min(p["x"] for p in points), max(p["x"] for p in points)),
        "y": (min(p["y"] for p in points), max(p["y"] for p in points)),
        "z": (min(p["z"] for p in points), max(p["z"] for p in points)),
    }
    height = bounds["z"][1] - bounds["z"][0]

    data_json = json.dumps(points, separators=(",", ":"))
    phase_order = ["pleines", "pilier", "coude", "surplomb"]
    legend_rows = "\n".join(
        f'<div class="legend-row"><span class="swatch" style="background:{PHASE_COLORS[p]}"></span>'
        f'<span class="legend-label">{PHASE_LABELS[p]}</span></div>'
        for p in phase_order
    )

    return f"""<title>Trajectoire L-Surplomb</title>
<meta name="viewport" content="width=device-width, initial-scale=1" />
<style>
  :root {{
    --bg: #eef1f6;
    --panel-bg: #ffffff;
    --panel-border: #d5dbe6;
    --text: #1b2230;
    --text-dim: #5c6779;
    --accent: #d9622b;
    --accent-soft: #f4dcc9;
    --mono-line: #c7cedb;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      --bg: #0d0f14;
      --panel-bg: #161a22;
      --panel-border: #262d3a;
      --text: #eef2f7;
      --text-dim: #8fa3bf;
      --accent: #ff8a3d;
      --accent-soft: #3a2a1a;
      --mono-line: #232936;
    }}
  }}
  :root[data-theme="dark"] {{
    --bg: #0d0f14;
    --panel-bg: #161a22;
    --panel-border: #262d3a;
    --text: #eef2f7;
    --text-dim: #8fa3bf;
    --accent: #ff8a3d;
    --accent-soft: #3a2a1a;
    --mono-line: #232936;
  }}

  * {{ box-sizing: border-box; }}
  html, body {{ margin: 0; padding: 0; height: 100%; overflow: hidden; }}
  body {{
    background: var(--bg);
    color: var(--text);
    font-family: 'Archivo', 'Inter', system-ui, sans-serif;
  }}

  #canvas-wrap {{ position: fixed; inset: 0; }}
  canvas {{ display: block; width: 100%; height: 100%; cursor: grab; touch-action: none; }}
  canvas.dragging {{ cursor: grabbing; }}

  header {{
    position: fixed; top: 0; left: 0; right: 0;
    display: flex; align-items: baseline; gap: 1.25rem;
    padding: 0.85rem 1.25rem;
    background: linear-gradient(to bottom, color-mix(in srgb, var(--bg) 88%, transparent), transparent);
    pointer-events: none;
    z-index: 10;
  }}
  h1 {{
    margin: 0;
    font-family: 'Archivo', sans-serif;
    font-weight: 700;
    font-size: 1.05rem;
    letter-spacing: 0.01em;
  }}
  h1 .unit {{ color: var(--accent); }}
  .stats {{
    display: flex; gap: 1rem; flex-wrap: wrap;
    font-family: 'JetBrains Mono', ui-monospace, monospace;
    font-variant-numeric: tabular-nums;
    font-size: 0.72rem;
    color: var(--text-dim);
  }}
  .stats b {{ color: var(--text); font-weight: 600; }}

  #legend {{
    position: fixed; left: 1.25rem; bottom: 1.25rem;
    background: var(--panel-bg);
    border: 1px solid var(--panel-border);
    border-radius: 10px;
    padding: 0.7rem 0.9rem;
    font-size: 0.74rem;
    box-shadow: 0 4px 18px rgba(0,0,0,0.18);
    z-index: 10;
  }}
  .legend-row {{ display: flex; align-items: center; gap: 0.5rem; padding: 0.18rem 0; }}
  .swatch {{ width: 0.7rem; height: 0.7rem; border-radius: 2px; flex: none; }}
  .legend-label {{ color: var(--text-dim); font-family: 'JetBrains Mono', monospace; }}
  #legend .hint {{
    margin-top: 0.5rem; padding-top: 0.45rem;
    border-top: 1px solid var(--panel-border);
    color: var(--text-dim); font-size: 0.68rem; line-height: 1.4;
  }}

  #travel-toggle {{
    position: fixed; right: 1.25rem; bottom: 1.25rem;
    background: var(--panel-bg);
    border: 1px solid var(--panel-border);
    color: var(--text-dim);
    border-radius: 8px;
    padding: 0.5rem 0.8rem;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.72rem;
    cursor: pointer;
    z-index: 10;
  }}
  #travel-toggle[aria-pressed="true"] {{ color: var(--accent); border-color: var(--accent); }}
  #travel-toggle:focus-visible {{ outline: 2px solid var(--accent); outline-offset: 2px; }}
</style>

<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Archivo:wght@500;700&family=JetBrains+Mono:wght@400;600&display=swap">

<div id="canvas-wrap"><canvas id="c"></canvas></div>

<header>
  <h1>Trajectoire <span class="unit">L-Surplomb</span></h1>
  <div class="stats">
    <span><b>{n_total:,}</b> points</span>
    <span><b>{n_extrude:,}</b> extrudes</span>
    <span>hauteur <b>{height:.1f}</b> mm</span>
    <span>source <b>{src_name}</b></span>
  </div>
</header>

<div id="legend">
{legend_rows}
  <div class="hint">glisser: tourner &middot; molette: zoom &middot; maj+glisser: deplacer</div>
</div>

<button id="travel-toggle" aria-pressed="false">trajets a vide : masques</button>

<script>
const POINTS = {data_json};
const PHASE_COLORS = {json.dumps(PHASE_COLORS)};
const TRAVEL_COLOR = "#e0455a";

const canvas = document.getElementById('c');
const ctx = canvas.getContext('2d');
let W, H, DPR;
function resize() {{
  DPR = window.devicePixelRatio || 1;
  W = canvas.clientWidth; H = canvas.clientHeight;
  canvas.width = W * DPR; canvas.height = H * DPR;
}}
window.addEventListener('resize', resize);
resize();

// centre du modele (pour tourner autour de son propre milieu)
let cx = 0, cy = 0, cz = 0;
for (const p of POINTS) {{ cx += p.x; cy += p.y; cz += p.z; }}
cx /= POINTS.length; cy /= POINTS.length; cz /= POINTS.length;

let maxExtent = 1;
for (const p of POINTS) {{
  maxExtent = Math.max(maxExtent, Math.abs(p.x - cx), Math.abs(p.y - cy), Math.abs(p.z - cz));
}}

let rotX = -0.55, rotY = 0.65;
let zoom = 1.0;
let panX = 0, panY = 0;
let showTravel = false;

function project(p) {{
  let x = p.x - cx, y = p.y - cy, z = p.z - cz;
  // rotation Y (autour de l'axe vertical ecran)
  let cosY = Math.cos(rotY), sinY = Math.sin(rotY);
  let x1 = x * cosY - y * sinY;
  let y1 = x * sinY + y * cosY;
  // rotation X (bascule)
  let cosX = Math.cos(rotX), sinX = Math.sin(rotX);
  let y2 = y1 * cosX - z * sinX;
  let z2 = y1 * sinX + z * cosX;
  const scale = (Math.min(W, H) * 0.38 / maxExtent) * zoom;
  return {{
    sx: W / 2 + x1 * scale + panX,
    sy: H / 2 - z2 * scale + panY,
    depth: y2,
  }};
}}

function draw() {{
  ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
  const bg = getComputedStyle(document.body).getPropertyValue('background-color') || '#0d0f14';
  ctx.clearRect(0, 0, W, H);

  // segments d'extrusion, groupes par run contigu de la meme phase pour
  // eviter de relier deux extrusions separees par un trajet a vide.
  ctx.lineWidth = 1.15;
  ctx.lineJoin = 'round';
  let prevProj = null, prevExtruding = null, prevPhase = null;

  for (let i = 0; i < POINTS.length; i++) {{
    const p = POINTS[i];
    const proj = project(p);
    if (prevProj) {{
      const isTravel = !(p.e && prevExtruding);
      if (isTravel) {{
        if (showTravel) {{
          ctx.strokeStyle = TRAVEL_COLOR;
          ctx.globalAlpha = 0.55;
          ctx.setLineDash([2, 3]);
          ctx.beginPath();
          ctx.moveTo(prevProj.sx, prevProj.sy);
          ctx.lineTo(proj.sx, proj.sy);
          ctx.stroke();
          ctx.setLineDash([]);
          ctx.globalAlpha = 1;
        }}
      }} else {{
        ctx.strokeStyle = PHASE_COLORS[p.p] || '#8fa3bf';
        ctx.globalAlpha = 0.92;
        ctx.beginPath();
        ctx.moveTo(prevProj.sx, prevProj.sy);
        ctx.lineTo(proj.sx, proj.sy);
        ctx.stroke();
        ctx.globalAlpha = 1;
      }}
    }}
    prevProj = proj;
    prevExtruding = p.e;
    prevPhase = p.p;
  }}
}}

// --- interactions : rotation (glisser), pan (maj+glisser), zoom (molette) ---
let dragging = false, shiftDrag = false, lastX = 0, lastY = 0;
canvas.addEventListener('pointerdown', (ev) => {{
  dragging = true;
  shiftDrag = ev.shiftKey;
  lastX = ev.clientX; lastY = ev.clientY;
  canvas.classList.add('dragging');
  canvas.setPointerCapture(ev.pointerId);
}});
canvas.addEventListener('pointermove', (ev) => {{
  if (!dragging) return;
  const dx = ev.clientX - lastX, dy = ev.clientY - lastY;
  lastX = ev.clientX; lastY = ev.clientY;
  if (shiftDrag) {{
    panX += dx; panY += dy;
  }} else {{
    rotY += dx * 0.006;
    rotX += dy * 0.006;
    rotX = Math.max(-1.5, Math.min(1.5, rotX));
  }}
  draw();
}});
window.addEventListener('pointerup', () => {{ dragging = false; canvas.classList.remove('dragging'); }});
canvas.addEventListener('wheel', (ev) => {{
  ev.preventDefault();
  zoom *= Math.pow(1.0015, -ev.deltaY);
  zoom = Math.max(0.15, Math.min(12, zoom));
  draw();
}}, {{ passive: false }});

document.getElementById('travel-toggle').addEventListener('click', (ev) => {{
  showTravel = !showTravel;
  ev.target.setAttribute('aria-pressed', String(showTravel));
  ev.target.textContent = 'trajets a vide : ' + (showTravel ? 'affiches' : 'masques');
  draw();
}});

draw();
</script>
"""


def main() -> None:
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join("out", "l_overhang.gcode")
    dst = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        "out", "l_overhang_preview.html")

    points = parse_gcode(src)
    points = decimate(points, target=20000)

    html = build_html(points, os.path.basename(src))
    out_dir = os.path.dirname(dst)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(dst, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"lu     : {src}")
    print(f"ecrit  : {dst}")
    print(f"points : {len(points)} (apres decimation)")


if __name__ == "__main__":
    main()
