"""Exporte le G-code Apex en page HTML ANIMEE (plateau + tete + depot).

Contrairement a scripts/preview_gcode.py (trace statique, repere PIECE), cette
vue anime le repere MACHINE : le plateau TOURNE (C), la tete se deplace dans
son demi-plan (R, Z) et s'incline (B). C'est le seul cadrage qui rend visible
l'oscillation de C et de R -- dans le repere piece, la rotation du plateau est
invisible par construction, puisqu'on la soustrait.

Le fichier de sortie est autonome (aucun CDN, aucune dependance).

Usage : python scripts/export_machine_animation.py [gcode] [sortie.html]
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

PHASES = [("couche 1", "base"), ("pilier", "pilier"),
          ("coude", "coude"), ("surplomb", "surplomb")]


def parse(path: str) -> list[dict]:
    pts: list[dict] = []
    r = c = z = b = 0.0
    e_prev = 0.0
    phase = "base"
    for line in open(path, "r", encoding="ascii"):
        line = line.strip()
        if line.startswith(";"):
            # Seules les lignes de SECTION ("; --- ... ---" et "; couche 1")
            # marquent une phase. L'en-tete du fichier contient les mots
            # "coude" et "surplomb" dans sa description : les lire comme des
            # marqueurs classait le premier deplacement en phase "coude".
            if line.startswith("; ---") or line.startswith("; couche"):
                for needle, label in PHASES:
                    if needle in line:
                        phase = label
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
        ext = False
        if m.group("e") is not None:
            ev = float(m.group("e"))
            ext = ev > e_prev + 1e-9
            e_prev = ev
        pts.append({"r": r, "c": c, "z": z, "b": b, "e": ext, "p": phase})
    if not pts:
        raise SystemExit(f"aucun point G00 dans {path}")
    return pts


def decimate(pts: list[dict], target: int) -> list[dict]:
    """Garde ~target points, SANS jamais fusionner deux points qui changent
    de phase ou d'etat d'extrusion : sinon un trajet a vide se retrouverait
    dessine comme un cordon depose (faux positif visuel)."""
    if len(pts) <= target:
        return pts
    stride = max(1, len(pts) // target)
    out = []
    prev = None
    for i, p in enumerate(pts):
        keep = (i % stride == 0 or i == len(pts) - 1 or prev is None
                or p["p"] != prev["p"] or p["e"] != prev["e"])
        if keep:
            out.append(p)
            prev = p
    return out


def main() -> None:
    src = sys.argv[1] if len(sys.argv) > 1 else "out/elbow_tube.gcode"
    dst = sys.argv[2] if len(sys.argv) > 2 else "out/elbow_tube_machine.html"

    # Decimation volontairement TRES haute : mesure sur cette piece, passer
    # de 20077 a 6730 points triple le pas angulaire median (0.41 -> 1.23
    # deg) et fait passer le p95 de 2.4 a 7.1 deg. Un saut de plateau de 7
    # deg entre deux images se voit comme un a-coup. 20077 points quantifies
    # a 2 decimales pesent ~500 Ko : largement sous la limite de 16 Mo, donc
    # rien ne justifie de les jeter.
    pts = decimate(parse(src), 40000)

    # quantification : 2 decimales suffisent a 0.01 mm / 0.01 deg, et divise
    # le poids du JSON par ~2 par rapport aux flottants bruts.
    data = [[round(p["r"], 2), round(p["c"], 2), round(p["z"], 2),
             round(p["b"], 2), 1 if p["e"] else 0,
             ["base", "pilier", "coude", "surplomb"].index(p["p"])]
            for p in pts]

    stats = {
        "n": len(pts),
        "cmin": min(p["c"] for p in pts), "cmax": max(p["c"] for p in pts),
        "rmin": min(p["r"] for p in pts), "rmax": max(p["r"] for p in pts),
        "zmax": max(p["z"] for p in pts),
        "bmax": max(p["b"] for p in pts),
        "src": os.path.basename(src),
    }

    html = TEMPLATE.replace("__DATA__", json.dumps(data, separators=(",", ":")))
    html = html.replace("__STATS__", json.dumps(stats, separators=(",", ":")))

    if os.path.dirname(dst):
        os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"lu    : {src}")
    print(f"ecrit : {dst}  ({os.path.getsize(dst)//1024} Ko)")
    print(f"points: {len(pts)}")
    print(f"C {stats['cmin']:.1f}..{stats['cmax']:.1f}  "
          f"R {stats['rmin']:.1f}..{stats['rmax']:.1f}  B max {stats['bmax']:.1f}")


TEMPLATE = r"""<title>Apex 4 Axes</title>
<style>
  :root {
    color-scheme: light;
    --bg: #ecebe6;
    --panel: #ffffff;
    --line: #cfccc2;
    --ink: #1a1c1a;
    --ink-dim: #6b6f68;
    --steel: #8d9491;
    --bed: #d8d5cc;
    --hot: #c8501e;
    --base: #4a7c6f;
    --pilier: #6f8ca8;
    --coude: #c99a2e;
    --surplomb: #c8501e;
    --travel: #b8332f;
  }
  @media (prefers-color-scheme: dark) {
    :root:not([data-theme="light"]) {
      color-scheme: dark;
      --bg: #14161a;
      --panel: #1d2026;
      --line: #2e333c;
      --ink: #eceef2;
      --ink-dim: #8b929c;
      --steel: #69717c;
      --bed: #272b32;
      --hot: #ff7a45;
      --base: #5fae9a;
      --pilier: #88a8c8;
      --coude: #e8b84b;
      --surplomb: #ff7a45;
      --travel: #e8555a;
    }
  }
  :root[data-theme="dark"] {
    color-scheme: dark;
    --bg: #14161a; --panel: #1d2026; --line: #2e333c;
    --ink: #eceef2; --ink-dim: #8b929c; --steel: #69717c; --bed: #272b32;
    --hot: #ff7a45; --base: #5fae9a; --pilier: #88a8c8;
    --coude: #e8b84b; --surplomb: #ff7a45; --travel: #e8555a;
  }

  * { box-sizing: border-box; }
  html, body { height: 100%; margin: 0; }
  body {
    background: var(--bg); color: var(--ink);
    font-family: ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif;
    display: flex; flex-direction: column;
  }
  canvas { display: block; width: 100%; height: 100%; touch-action: none; cursor: grab; }
  canvas.drag { cursor: grabbing; }
  #stage { flex: 1; min-height: 0; position: relative; }

  header {
    display: flex; align-items: baseline; gap: 1rem; flex-wrap: wrap;
    padding: 0.7rem 1rem 0.55rem;
    border-bottom: 1px solid var(--line);
  }
  h1 { margin: 0; font-size: 0.98rem; font-weight: 650; letter-spacing: -0.01em; }
  h1 em { font-style: normal; color: var(--hot); }
  .sub { font-size: 0.72rem; color: var(--ink-dim); }

  #hud {
    position: absolute; top: 0.75rem; left: 0.75rem;
    background: color-mix(in srgb, var(--panel) 92%, transparent);
    border: 1px solid var(--line); border-radius: 9px;
    padding: 0.6rem 0.75rem;
    font-family: ui-monospace, "SF Mono", Menlo, monospace;
    font-size: 0.72rem; font-variant-numeric: tabular-nums;
    line-height: 1.55; min-width: 12.5rem;
  }
  #hud .row { display: flex; justify-content: space-between; gap: 1.1rem; }
  #hud .k { color: var(--ink-dim); }
  #hud .phase { margin-top: 0.4rem; padding-top: 0.4rem; border-top: 1px solid var(--line); }
  #hud .dot { display: inline-block; width: 0.55rem; height: 0.55rem; border-radius: 2px; margin-right: 0.4rem; }

  #controls {
    display: flex; align-items: center; gap: 0.7rem; flex-wrap: wrap;
    padding: 0.6rem 1rem;
    border-top: 1px solid var(--line);
    background: var(--panel);
  }
  button {
    font: inherit; font-size: 0.78rem;
    background: var(--panel); color: var(--ink);
    border: 1px solid var(--line); border-radius: 7px;
    padding: 0.35rem 0.8rem; cursor: pointer;
  }
  button:hover { border-color: var(--steel); }
  button.primary { background: var(--hot); border-color: var(--hot); color: #fff; font-weight: 600; }
  button:focus-visible { outline: 2px solid var(--hot); outline-offset: 2px; }
  input[type=range] { accent-color: var(--hot); flex: 1; min-width: 8rem; }
  input[type=range].speedrange { flex: 0 1 7rem; min-width: 4.5rem; }
  label { font-size: 0.74rem; color: var(--ink-dim); display: flex; align-items: center; gap: 0.4rem; }
  input[type=number] {
    font: inherit; font-size: 0.78rem;
    font-family: ui-monospace, "SF Mono", Menlo, monospace;
    font-variant-numeric: tabular-nums;
    width: 5rem; padding: 0.28rem 0.4rem;
    background: var(--panel); color: var(--ink);
    border: 1px solid var(--line); border-radius: 6px;
  }
  input[type=number]:focus-visible { outline: 2px solid var(--hot); outline-offset: 1px; }
  .unit { font-size: 0.68rem; color: var(--ink-dim); }
  .legend { display: flex; gap: 0.8rem; flex-wrap: wrap; font-size: 0.7rem; color: var(--ink-dim); }
  .legend span.dot { width: 0.6rem; height: 0.6rem; border-radius: 2px; display: inline-block; margin-right: 0.3rem; }
  @media (max-width: 560px) { #hud { font-size: 0.66rem; min-width: 10rem; } }
</style>

<header>
  <h1>Apex <em>4 axes</em> &mdash; tube coud&eacute; en surplomb</h1>
  <div class="sub" id="sub"></div>
</header>

<div id="stage">
  <canvas id="c"></canvas>
  <div id="hud">
    <div class="row"><span class="k">R</span><span id="vr">&mdash;</span></div>
    <div class="row"><span class="k">C</span><span id="vc">&mdash;</span></div>
    <div class="row"><span class="k">Z</span><span id="vz">&mdash;</span></div>
    <div class="row"><span class="k">B</span><span id="vb">&mdash;</span></div>
    <div class="phase"><span class="dot" id="pdot"></span><span id="pname">&mdash;</span></div>
  </div>
</div>

<div id="controls">
  <button class="primary" id="play">Pause</button>
  <input type="range" id="scrub" min="0" max="100" value="0" step="0.01" aria-label="Avancement" />
  <label for="speed">vitesse
    <input type="number" id="speed" value="12" min="0.25" max="2000" step="0.25" />
    <span class="unit">pts/frame</span>
  </label>
  <input type="range" id="speedrange" min="0" max="100" value="38" step="1"
         aria-label="Vitesse" class="speedrange" />
  <label><input type="checkbox" id="showtravel" /> trajets &agrave; vide</label>
  <div class="legend">
    <span><span class="dot" style="background:var(--base)"></span>base</span>
    <span><span class="dot" style="background:var(--pilier)"></span>pilier</span>
    <span><span class="dot" style="background:var(--coude)"></span>coude</span>
    <span><span class="dot" style="background:var(--surplomb)"></span>surplomb</span>
  </div>
</div>

<script>
const PTS = __DATA__;
const S = __STATS__;
const PHASE_NAMES = ["couche 1 (base)", "pilier", "coude", "surplomb"];
const PHASE_VARS  = ["--base", "--pilier", "--coude", "--surplomb"];

document.getElementById('sub').textContent =
  `${S.n.toLocaleString('fr')} points · C ${S.cmin.toFixed(0)}° → ${S.cmax.toFixed(0)}° ` +
  `· R ${S.rmin.toFixed(0)}–${S.rmax.toFixed(0)} mm · B max ${S.bmax.toFixed(0)}°`;

const cv = document.getElementById('c'), ctx = cv.getContext('2d');
let W = 0, H = 0, DPR = 1;
function resize() {
  DPR = window.devicePixelRatio || 1;
  W = cv.clientWidth; H = cv.clientHeight;
  cv.width = Math.max(1, W * DPR); cv.height = Math.max(1, H * DPR);
  draw();
}
addEventListener('resize', resize);

function css(v) { return getComputedStyle(document.documentElement).getPropertyValue(v).trim(); }

// --- camera -------------------------------------------------------------
let rotZ = -0.62, tilt = 1.02, zoom = 1, panX = 0, panY = 0;
const BED_R = Math.max(60, S.rmax + 12);

// Projection : rotation autour de l'axe Z machine, puis bascule. Z va vers
// le HAUT de l'ecran (sy est negatif quand z croit).
function project(x, y, z) {
  const cr = Math.cos(rotZ), sr = Math.sin(rotZ);
  const x1 = x * cr - y * sr, y1 = x * sr + y * cr;
  const ct = Math.cos(tilt), st = Math.sin(tilt);
  const depth = y1 * ct + z * st;        // vers l'observateur
  const sy = -(z * ct - y1 * st);        // ecran : Z vers le haut
  const s = (Math.min(W, H) * 0.33 / BED_R) * zoom;
  return { x: W / 2 + x1 * s + panX, y: H / 2 + sy * s + panY + H * 0.12, d: depth };
}

// --- etat d'animation ---------------------------------------------------
let head = 0;          // index du point courant (partie entiere de headF)
let headF = 0;         // position FRACTIONNAIRE : interpolation entre points
let playing = true;
let speed = 12;        // points par frame
let showTravel = false;

// UNE SEULE rotation plateau, partagee par TOUT ce qui est solidaire du
// plateau (rayons repere ET matiere deposee). Deux formules separees ont
// deja diverge une fois : la matiere tournait de +C, les rayons de -C --
// plateau et piece tournaient en sens opposes, et la piece etait affichee
// en miroir (invisible sur un tube symetrique).
//
// Convention du generateur : C = atan2(y, x) du point PIECE a amener sous
// la tete (demi-plan machine y=0). Un point piece d'angle theta arrive a
// l'angle machine 0 quand le plateau a tourne de -theta : la rotation
// plateau vaut donc -C.
function plateToMachine(thetaPartDeg, r, z, cNow) {
  const a = (thetaPartDeg - cNow) * Math.PI / 180;
  return [r * Math.cos(a), r * Math.sin(a), z];
}

function pointXYZ(p, cNow) {
  // p = [r, c, z, b, e, phase]. Depose a (r, 0, z) machine quand le plateau
  // etait a C = p[1], donc a l'angle PIECE p[1] : c'est cet angle piece,
  // fixe, que le plateau emporte ensuite avec lui.
  return plateToMachine(p[1], p[0], p[2], cNow);
}

function drawBed(cNow) {
  const steps = 72;
  ctx.strokeStyle = css('--line'); ctx.lineWidth = 1;
  ctx.beginPath();
  for (let i = 0; i <= steps; i++) {
    const t = i / steps * Math.PI * 2;
    const q = project(BED_R * Math.cos(t), BED_R * Math.sin(t), 0);
    i ? ctx.lineTo(q.x, q.y) : ctx.moveTo(q.x, q.y);
  }
  ctx.closePath();
  ctx.fillStyle = css('--bed'); ctx.globalAlpha = 0.55; ctx.fill(); ctx.globalAlpha = 1;
  ctx.stroke();

  // rayons solidaires du plateau : c'est EUX qui rendent la rotation visible
  ctx.strokeStyle = css('--steel'); ctx.globalAlpha = 0.5;
  for (let k = 0; k < 12; k++) {
    const e = plateToMachine(k * 30, BED_R, 0, cNow);
    const a = project(0, 0, 0), b = project(e[0], e[1], e[2]);
    ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke();
  }
  ctx.globalAlpha = 1;
}

function drawHead(p) {
  // La tete vit dans le demi-plan machine y=0 : elle ne tourne JAMAIS avec
  // le plateau. C'est la contrainte cinematique que cette vue doit montrer.
  const r = p[0], z = p[2], b = p[3] * Math.PI / 180;
  const tip = project(r, 0, z);
  // buse : B>0 penche la pointe vers -R (vers le centre)
  const L = BED_R * 0.42;
  const bx = r + L * Math.sin(b), bz = z + L * Math.cos(b);
  const top = project(bx, 0, bz);

  ctx.strokeStyle = css('--steel'); ctx.lineWidth = 5; ctx.lineCap = 'round';
  ctx.beginPath(); ctx.moveTo(tip.x, tip.y); ctx.lineTo(top.x, top.y); ctx.stroke();

  // corps
  ctx.lineWidth = 11; ctx.globalAlpha = 0.34;
  const mid = project(r + L * 0.55 * Math.sin(b), 0, z + L * 0.55 * Math.cos(b));
  ctx.beginPath(); ctx.moveTo(mid.x, mid.y); ctx.lineTo(top.x, top.y); ctx.stroke();
  ctx.globalAlpha = 1;

  // pointe chaude
  ctx.fillStyle = css('--hot');
  ctx.beginPath(); ctx.arc(tip.x, tip.y, 4.2, 0, Math.PI * 2); ctx.fill();
  ctx.lineWidth = 1;
}

// Etat INTERPOLE entre deux points commandes. Sans ca, la tete saute de
// point en point : a vitesse lente c'est un defilement saccade, a vitesse
// rapide des bonds. Le firmware, lui, interpole ses axes lineairement entre
// deux consignes -- on affiche donc ce qu'il fera reellement.
function stateAt(f) {
  const i = Math.min(Math.floor(f), PTS.length - 1);
  const j = Math.min(i + 1, PTS.length - 1);
  const t = j > i ? f - i : 0;
  const a = PTS[i], b = PTS[j];
  return [a[0] + (b[0] - a[0]) * t,
          a[1] + (b[1] - a[1]) * t,
          a[2] + (b[2] - a[2]) * t,
          a[3] + (b[3] - a[3]) * t,
          b[4], b[5]];
}

function draw() {
  ctx.setTransform(DPR, 0, 0, DPR, 0, 0);
  ctx.clearRect(0, 0, W, H);
  if (!PTS.length) return;

  const cur = stateAt(headF);
  const cNow = cur[1];

  drawBed(cNow);

  // matiere deposee, exprimee a l'angle plateau COURANT
  ctx.lineWidth = 1.7; ctx.lineJoin = 'round'; ctx.lineCap = 'round';
  let prev = null, prevP = null;
  const last = Math.min(Math.floor(headF), PTS.length - 1);
  for (let i = 0; i <= last + 1; i++) {
    // le dernier element est la position INTERPOLEE : le cordon en cours de
    // depot s'arrete sous la buse, il n'apparait pas d'un bloc.
    const p = (i > last) ? cur : PTS[i];
    const v = pointXYZ(p, cNow);
    const q = project(v[0], v[1], v[2]);
    if (prev) {
      if (p[4]) {
        ctx.strokeStyle = css(PHASE_VARS[p[5]]);
        ctx.globalAlpha = 0.95;
        ctx.beginPath(); ctx.moveTo(prev.x, prev.y); ctx.lineTo(q.x, q.y); ctx.stroke();
      } else if (showTravel) {
        ctx.strokeStyle = css('--travel'); ctx.globalAlpha = 0.5;
        ctx.setLineDash([2, 3]);
        ctx.beginPath(); ctx.moveTo(prev.x, prev.y); ctx.lineTo(q.x, q.y); ctx.stroke();
        ctx.setLineDash([]);
      }
      ctx.globalAlpha = 1;
    }
    prev = q; prevP = p;
  }

  drawHead(cur);

  document.getElementById('vr').textContent = cur[0].toFixed(2) + ' mm';
  document.getElementById('vc').textContent = cur[1].toFixed(2) + '°';
  document.getElementById('vz').textContent = cur[2].toFixed(2) + ' mm';
  document.getElementById('vb').textContent = cur[3].toFixed(2) + '°';
  document.getElementById('pname').textContent = PHASE_NAMES[cur[5]];
  document.getElementById('pdot').style.background = css(PHASE_VARS[cur[5]]);
  const sc = document.getElementById('scrub');
  if (document.activeElement !== sc) sc.value = String(headF / (PTS.length - 1) * 100);
}

let raf = null;
// Avancement en temps REEL plutot qu'en points par image : a 0.25 pt/frame
// l'ancienne version ne bougeait qu'une image sur quatre (saccade), et le
// rythme dependait du taux de rafraichissement de l'ecran.
let lastT = 0;
function tick(now) {
  const dt = lastT ? Math.min(0.05, (now - lastT) / 1000) : 0;
  lastT = now;
  if (playing) {
    headF += speed * 60 * dt;          // `speed` reste exprime en pts/frame a 60 Hz
    if (headF >= PTS.length - 1) { headF = PTS.length - 1; playing = false; syncPlay(); }
    head = Math.floor(headF);
    draw();
  }
  raf = requestAnimationFrame(tick);
}

function syncPlay() {
  document.getElementById('play').textContent =
    playing ? 'Pause' : (headF >= PTS.length - 1 ? 'Rejouer' : 'Lecture');
}

document.getElementById('play').addEventListener('click', () => {
  if (!playing && headF >= PTS.length - 1) { headF = 0; head = 0; }
  playing = !playing; lastT = 0; syncPlay();
});
document.getElementById('scrub').addEventListener('input', (ev) => {
  headF = ev.target.value / 100 * (PTS.length - 1);
  head = Math.floor(headF);
  playing = false; syncPlay(); draw();
});
// Vitesse : champ libre + curseur, tenus synchronises. Le curseur est
// LOGARITHMIQUE -- de 0.25 a 2000 pts/frame en lineaire, les neuf premiers
// dixiemes de la course seraient au-dessus de 200, donc inutilisables pour
// les reglages lents (ceux qui servent a examiner le coude image par image).
const SPEED_MIN = 0.25, SPEED_MAX = 2000;
const speedBox = document.getElementById('speed');
const speedRange = document.getElementById('speedrange');

function speedToSlider(v) {
  const t = (Math.log(v) - Math.log(SPEED_MIN)) / (Math.log(SPEED_MAX) - Math.log(SPEED_MIN));
  return Math.round(Math.min(1, Math.max(0, t)) * 100);
}
function sliderToSpeed(s) {
  const v = Math.exp(Math.log(SPEED_MIN) + (s / 100) * (Math.log(SPEED_MAX) - Math.log(SPEED_MIN)));
  return v >= 10 ? Math.round(v) : Math.round(v * 4) / 4;
}
function setSpeed(v, from) {
  if (!isFinite(v)) return;
  speed = Math.min(SPEED_MAX, Math.max(SPEED_MIN, v));
  if (from !== 'box') speedBox.value = String(speed);
  if (from !== 'range') speedRange.value = String(speedToSlider(speed));
}

speedBox.addEventListener('input', () => setSpeed(parseFloat(speedBox.value), 'box'));
// Au blur, on renormalise le champ : une saisie hors bornes (ou vide) doit
// afficher la valeur reellement appliquee, jamais rester en desaccord.
speedBox.addEventListener('blur', () => setSpeed(parseFloat(speedBox.value) || speed, null));
speedRange.addEventListener('input', () => setSpeed(sliderToSpeed(+speedRange.value), 'range'));
setSpeed(12, null);
document.getElementById('showtravel').addEventListener('change', (ev) => {
  showTravel = ev.target.checked; draw();
});

// --- souris -------------------------------------------------------------
let drag = false, shift = false, lx = 0, ly = 0;
cv.addEventListener('pointerdown', (e) => {
  drag = true; shift = e.shiftKey; lx = e.clientX; ly = e.clientY;
  cv.classList.add('drag'); cv.setPointerCapture(e.pointerId);
});
cv.addEventListener('pointermove', (e) => {
  if (!drag) return;
  const dx = e.clientX - lx, dy = e.clientY - ly; lx = e.clientX; ly = e.clientY;
  if (shift) { panX += dx; panY += dy; }
  else { rotZ += dx * 0.007; tilt = Math.max(0.05, Math.min(1.5, tilt + dy * 0.006)); }
  draw();
});
addEventListener('pointerup', () => { drag = false; cv.classList.remove('drag'); });
cv.addEventListener('wheel', (e) => {
  e.preventDefault();
  zoom = Math.max(0.3, Math.min(6, zoom * Math.pow(1.0015, -e.deltaY)));
  draw();
}, { passive: false });

matchMedia('(prefers-color-scheme: dark)').addEventListener('change', draw);
resize();
syncPlay();
tick();
</script>
"""


if __name__ == "__main__":
    main()
