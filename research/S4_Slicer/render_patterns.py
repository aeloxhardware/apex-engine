"""
Planche des motifs de P4 : une vignette par (forme, motif) du banc
test_patterns.py (output_models/patterns_bench.pkl). Surface grise, passes en
degrade d'ordre de depot (violet = premieres, jaune = dernieres), trajets a
vide en fins traits rouges. Vue le long de la normale moyenne de la zone.

Usage (venv de la recherche) :
  .venv/Scripts/python.exe render_patterns.py [--patterns=a,b] [--shapes=mot1,mot2] [--out=nom.png] [--data=x.pkl]
"""
import os
import pickle
import sys

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
TILE = 420


def polylines(paths, scalars):
    pts, cells, vals = [], [], []
    start = 0
    for p, s in zip(paths, scalars):
        pts.append(p)
        vals.append(s)
        cells.append(np.r_[len(p), np.arange(start, start + len(p))])
        start += len(p)
    mesh = pv.PolyData(np.vstack(pts), lines=np.concatenate(cells))
    mesh.point_data["order"] = np.concatenate(vals)
    return mesh


def tile(V, F, paths, title):
    surf = pv.PolyData(V, np.column_stack([np.full(len(F), 3), F]).ravel())
    tri = V[F]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]).sum(0)
    n = n / np.linalg.norm(n) if np.linalg.norm(n) > 1e-9 else np.array([0, 0, 1.0])
    plotter = pv.Plotter(off_screen=True, window_size=(TILE, TILE))
    plotter.set_background("white")
    plotter.add_mesh(surf, color="#d8d8d8", lighting=True)
    paths = [p for p in paths if len(p) >= 2]
    if paths:
        total = sum(len(p) for p in paths)
        acc, scal = 0, []
        for p in paths:
            scal.append((acc + np.arange(len(p))) / max(total - 1, 1))
            acc += len(p)
        # leger decalage vers l'observateur : lignes jamais cachees par la surface
        lift = 0.02 * np.linalg.norm(np.ptp(V, axis=0)) / 50
        plotter.add_mesh(polylines([p + n * lift for p in paths], scal), scalars="order", cmap="plasma",
                         line_width=2, show_scalar_bar=False)
        if len(paths) > 1:
            a = np.array([p[-1] for p in paths[:-1]]) + n * lift
            b = np.array([p[0] for p in paths[1:]]) + n * lift
            seg = np.empty((2 * len(a), 3))
            seg[0::2], seg[1::2] = a, b
            cells = np.column_stack([np.full(len(a), 2), np.arange(0, 2 * len(a), 2),
                                     np.arange(1, 2 * len(a), 2)]).ravel()
            plotter.add_mesh(pv.PolyData(seg, lines=cells), color="red", line_width=1, opacity=0.5)
    center = V.mean(0)
    up = np.array([0, 0, 1.0]) if abs(n[2]) < 0.9 else np.array([0, 1.0, 0])
    plotter.camera_position = [tuple(center + n * 3 * np.linalg.norm(np.ptp(V, axis=0))), tuple(center), tuple(up)]
    plotter.reset_camera()
    plotter.camera.zoom(1.25)
    img = Image.fromarray(plotter.screenshot(return_img=True))
    plotter.close()
    ImageDraw.Draw(img).text((6, 6), title, fill=(0, 0, 0))
    return img


if __name__ == "__main__":
    flag = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if "=" in a}
    data = pickle.load(open(os.path.join(HERE, "output_models", flag.get("--data", "patterns_bench.pkl")), "rb"))
    shapes = list(dict.fromkeys(k[0] for k in data))
    patterns = list(dict.fromkeys(k[1] for k in data))
    if "--patterns" in flag:
        patterns = [p for p in flag["--patterns"].split(",") if p in patterns]
    if "--shapes" in flag:
        words = flag["--shapes"].split(",")
        shapes = [s for s in shapes if any(w in s for w in words)]
    sheet = Image.new("RGB", (TILE * len(patterns), TILE * len(shapes)), "white")
    for r, s in enumerate(shapes):
        for c, p in enumerate(patterns):
            if (s, p) in data:
                V, F, paths = data[(s, p)]
                sheet.paste(tile(V, F, paths, f"{s}\n{p}  ({len(paths)} passes)"), (c * TILE, r * TILE))
                print(s, p, flush=True)
    out = os.path.join(HERE, "output_models", flag.get("--out", "patterns_bench.png"))
    sheet.save(out)
    print(out)
