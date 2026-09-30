"""
Rendu des couches de lissage non planaire (skin_paths.py) : piece d'origine
transparente, trajectoires colorees par moment d'impression (violet = tot,
jaune = tard), trois vues.

Usage (venv de la recherche) :
  .venv/Scripts/python.exe render_skin_paths.py <modele> <passe> [--print-order]
"""
import os
import pickle
import sys

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
model, k = sys.argv[1], int(sys.argv[2])
prefix = os.path.join(HERE, "output_models", f"loop_{model}")
skin = np.load(f"{prefix}_skin_pass{k}.npz")
data = pickle.load(open(f"{prefix}_skin_paths_pass{k}.pkl", "rb"))
if "--print-order" in sys.argv:          # zones gardees par P5 seulement
    order = pickle.load(open(f"{prefix}_print_order_pass{k}.pkl", "rb"))
    kept = {ps["zone"] for ps in order["passes"] if ps["kind"] == "skin"}
    data["zones"] = {z: e for z, e in data["zones"].items() if z in kept}
V, F = skin["surface_vertices"], skin["surface_faces"]
part = pv.PolyData(V, np.column_stack([np.full(len(F), 3), F]).ravel())

pts, cells, vals, start = [], [], [], 0
moments = [e["moment"] for e in data["zones"].values()]
lo, hi = min(moments), max(moments)
for e in data["zones"].values():
    for p in e["paths"]:
        pts.append(p)
        vals.append(np.full(len(p), (e["moment"] - lo) / max(hi - lo, 1)))
        cells.append(np.r_[len(p), np.arange(start, start + len(p))])
        start += len(p)
lines = pv.PolyData(np.vstack(pts), lines=np.concatenate(cells))
lines.point_data["moment"] = np.concatenate(vals)

shots = []
for azimuth, elevation in ((-60, 25), (120, 25), (30, -30)):
    plotter = pv.Plotter(off_screen=True, window_size=(800, 800))
    plotter.set_background("white")
    plotter.add_mesh(part, color="#cccccc", opacity=0.25)
    plotter.add_mesh(lines, scalars="moment", cmap="plasma", line_width=1.5, show_scalar_bar=False)
    center = np.array(part.center)
    az, el = np.radians(azimuth), np.radians(elevation)
    eye = center + 3 * part.length * np.array([np.cos(az) * np.cos(el), np.sin(az) * np.cos(el), np.sin(el)])
    plotter.camera_position = [tuple(eye), tuple(center), (0, 0, 1)]
    plotter.reset_camera()
    plotter.camera.zoom(1.3)
    shots.append(plotter.screenshot(return_img=True))
    plotter.close()
image = Image.fromarray(np.hstack(shots))
n_zones = sum(1 for e in data["zones"].values() if e["paths"])
ImageDraw.Draw(image).text(
    (10, 10), f"{model} passe {k} : {n_zones} zones, motif {data['general_pattern']} ; "
              f"couleur = moment d'impression (violet tot, jaune tard)", fill=(0, 0, 0))
out = f"{prefix}_skin_paths{'_final' if '--print-order' in sys.argv else ''}_pass{k}.png"
image.save(out)
print(out)
