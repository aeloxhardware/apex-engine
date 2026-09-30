"""
Carte des zones de lissage non-planaire (surface_skin.py) : faces colorees
par statut, trois vues, legende.

  vert   retenue            violet  inatteignable (plage d'axes)
  orange collision plateau  rouge   collision objet
  gris   contact plateau    bleu    zone trop petite

Usage (venv de la recherche) :
  .venv/Scripts/python.exe render_skin_map.py <fichier _skin_map_pass<k>.ply>
"""
import sys

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw

path = sys.argv[1]
mesh = pv.read(path)
LEGEND = [("retenue", (40, 170, 70)), ("inatteignable", (150, 0, 150)), ("collision plateau", (230, 140, 0)),
          ("collision objet", (220, 30, 30)), ("contact plateau", (120, 120, 120)), ("zone trop petite", (90, 140, 230))]

shots = []
for azimuth, elevation in ((-60, 20), (120, 20), (-60, -35)):
    plotter = pv.Plotter(off_screen=True, window_size=(600, 600))
    plotter.set_background("white")
    plotter.add_mesh(mesh, scalars="RGBA" if "RGBA" in mesh.cell_data else None, rgb=True,
                     show_edges=False, lighting=True)
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
legend = LEGEND if "_zones" not in path else [("une couleur par zone de lissage (morceau imprime d'un seul tenant)", (200, 200, 200)),
                                              ("gris clair : pas de couche de lissage", (225, 225, 225))]
for i, (label, color) in enumerate(legend):
    draw.rectangle([10, 10 + 22 * i, 28, 28 + 22 * i], fill=color)
    draw.text((34, 12 + 22 * i), label, fill=(0, 0, 0))
out = path.replace(".ply", ".png")
image.save(out)
print(out)
