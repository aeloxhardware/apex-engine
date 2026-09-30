"""
Rendu des raccords de measure_seams.py : morceaux de lissage en couleurs
(gris clair = pas de lissage), raccords traces en noir (arete franche,
< 135 deg entre faces) ou en magenta (arete douce, 135 deg au seuil : ceux qui
peuvent venir des facettes d'un STL grossier).

Usage (venv de la recherche) :
  .venv/Scripts/python.exe render_seams.py <fichier _seams_pass<k>.npz>
"""
import sys

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw

path = sys.argv[1]
d = np.load(path)
faces = np.column_stack([np.full(len(d["faces"]), 3), d["faces"]]).ravel()
mesh = pv.PolyData(d["vertices"], faces)
mesh.cell_data["rgb"] = d["colors"]


def segments(mask):
    a, b = d["seam_a"][mask], d["seam_b"][mask]
    if not len(a):
        return None
    pts = np.empty((2 * len(a), 3))
    pts[0::2], pts[1::2] = a, b
    lines = np.column_stack([np.full(len(a), 2), np.arange(0, 2 * len(a), 2), np.arange(1, 2 * len(a), 2)]).ravel()
    return pv.PolyData(pts, lines=lines)


sharp, soft = segments(d["seam_angle"] < 135), segments(d["seam_angle"] >= 135)
shots = []
for azimuth, elevation in ((-60, 20), (120, 20), (-60, -35)):
    plotter = pv.Plotter(off_screen=True, window_size=(600, 600))
    plotter.set_background("white")
    plotter.add_mesh(mesh, scalars="rgb", rgb=True, lighting=True)
    for seg, color in ((sharp, "black"), (soft, "magenta")):
        if seg is not None:
            plotter.add_mesh(seg, color=color, line_width=4)
    center = np.array(mesh.center)
    az, el = np.radians(azimuth), np.radians(elevation)
    eye = center + 3 * mesh.length * np.array([np.cos(az) * np.cos(el), np.sin(az) * np.cos(el), np.sin(el)])
    plotter.camera_position = [tuple(eye), tuple(center), (0, 0, 1)]
    plotter.reset_camera()
    plotter.camera.zoom(1.3)
    shots.append(plotter.screenshot(return_img=True))
    plotter.close()
image = Image.fromarray(np.hstack(shots))
draw = ImageDraw.Draw(image)
for i, (label, color) in enumerate([("raccord sur arete franche (< 135 deg)", (0, 0, 0)),
                                    ("raccord sur arete douce (135 deg au seuil)", (255, 0, 255)),
                                    ("couleurs : un morceau de lissage chacun ; gris clair : pas de lissage", (200, 200, 200))]):
    draw.rectangle([10, 10 + 22 * i, 28, 28 + 22 * i], fill=color)
    draw.text((34, 12 + 22 * i), label, fill=(0, 0, 0))
out = path.replace(".npz", ".png")
image.save(out)
print(out)
