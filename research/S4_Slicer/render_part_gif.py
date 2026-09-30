"""
GIF d'ordre d'impression des trajectoires (parois + remplissage) DANS la
piece d'origine affichee en transparence : le rendu fil de fer seul
(render_layers_gif.py) ne permet pas de juger si un morceau de trajectoire
appartient a la piece.

Couleurs : couches deposees du bleu (premiere) au jaune (derniere), couche en
cours de depot en blanc. Pas IMPAIR entre images, pour que le remplissage
croise-alterne se voie changer de sens d'une image a l'autre.

Usage : .venv/Scripts/python.exe render_part_gif.py <couches.pkl> <piece.stl> [azimut] [elevation]
"""
import pickle
import sys

import numpy as np
import pyvista as pv
from PIL import Image

FIRST_COLOR = np.array([0x00, 0x00, 0xFF])
LAST_COLOR = np.array([0xFF, 0xF7, 0x00])

path, stl = sys.argv[1], sys.argv[2]
azimuth = float(sys.argv[3]) if len(sys.argv) > 3 else -60.0
elevation = float(sys.argv[4]) if len(sys.argv) > 4 else 15.0
layers = pickle.load(open(path, "rb"))
total = len(layers)
stride = max(1, total // 90)
if stride % 2 == 0:
    stride += 1


def layer_polydata(indices):
    """Lignes des couches `indices`, avec une couleur RGB par cellule."""
    points, cells, colors, offset = [], [], [], 0
    for i in indices:
        color = FIRST_COLOR + (LAST_COLOR - FIRST_COLOR) * (i / max(total - 1, 1))
        for t in layers[i]["toolpath_3d"]:
            t = np.asarray(t)[:, :3]
            if len(t) < 2:
                continue
            points.append(t)
            cells.append(np.concatenate([[len(t)], np.arange(offset, offset + len(t))]))
            colors.append(color)
            offset += len(t)
    if not points:
        return None
    mesh = pv.PolyData(np.vstack(points), lines=np.concatenate(cells))
    mesh.cell_data["rgb"] = np.array(colors, dtype=np.uint8)
    return mesh


part = pv.read(stl)
plotter = pv.Plotter(off_screen=True, window_size=(800, 800))
plotter.set_background("black")
plotter.add_mesh(part, color="lightgray", opacity=0.15)
# camera explicite, Z vers le haut : azimut autour de Z, elevation au-dessus
# de l'horizontale (memes conventions que render_layers_gif.py)
center = np.array(part.center)
az, el = np.radians(azimuth), np.radians(elevation)
eye = center + 3.0 * part.length * np.array([np.cos(az) * np.cos(el), np.sin(az) * np.cos(el), np.sin(el)])
plotter.camera_position = [tuple(eye), tuple(center), (0.0, 0.0, 1.0)]
plotter.reset_camera()
plotter.camera.zoom(1.4)

frames, placed, current_actor = [], 0, None
shown = list(range(0, total, stride))
if shown[-1] != total - 1:
    shown.append(total - 1)
for i in shown:
    batch = layer_polydata(range(placed, i + 1))
    if batch is not None:
        plotter.add_mesh(batch, scalars="rgb", rgb=True, line_width=1)
    placed = i + 1
    if current_actor is not None:
        plotter.remove_actor(current_actor)
        current_actor = None
    current = layer_polydata([i]) if i != total - 1 else None
    if current is not None:
        current_actor = plotter.add_mesh(current, color="white", line_width=3)
    plotter.add_text(f"Couche {i + 1}/{total}", color="white", font_size=12, name="title")
    frames.append(Image.fromarray(plotter.screenshot(return_img=True)))
plotter.close()

out = path.replace(".pkl", "_in_part.gif")
frames[0].save(out, save_all=True, append_images=frames[1:], duration=100, loop=0)
for frac in (0.25, 0.5, 0.75, 1.0):
    frames[int(frac * (len(frames) - 1))].save(path.replace(".pkl", f"_in_part_{int(frac * 100)}pct.png"))
print(f"GIF : {out}  ({len(frames)} images)")
