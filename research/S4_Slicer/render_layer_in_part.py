"""
Une couche de trajectoire dans la piece d'origine, affichee en transparence :
pour juger si un morceau de remplissage appartient bien a la piece (le rendu
fil de fer du GIF ne montre pas la surface).

Usage : .venv/Scripts/python.exe render_layer_in_part.py <couches.pkl> <piece.stl> <couche 1-based> [...]
"""
import pickle
import sys

import numpy as np
import pyvista as pv

path, stl = sys.argv[1], sys.argv[2]
indices = [int(a) for a in sys.argv[3:]] or [118]
layers = pickle.load(open(path, "rb"))
part = pv.read(stl)

for number in indices:
    layer = layers[number - 1]
    polylines = [np.asarray(t)[:, :3] for t in layer["toolpath_3d"] if len(t) >= 2]
    points = np.vstack(polylines)
    cells, offset = [], 0
    for t in polylines:
        cells.append(np.concatenate([[len(t)], np.arange(offset, offset + len(t))]))
        offset += len(t)
    lines = pv.PolyData(points, lines=np.concatenate(cells))

    shots = []
    for view in ("iso", "xy", "xz"):
        plotter = pv.Plotter(off_screen=True, window_size=(700, 700))
        plotter.set_background("black")
        plotter.add_mesh(part, color="lightgray", opacity=0.18)
        plotter.add_mesh(lines, color="orangered", line_width=2)
        if view == "iso":
            plotter.view_isometric()
        elif view == "xy":
            plotter.view_xy()
        else:
            plotter.view_xz()
        plotter.add_text(f"couche {number}/{len(layers)} ({view})", color="white", font_size=10)
        shots.append(plotter.screenshot(return_img=True))
        plotter.close()
    out = path.replace(".pkl", f"_layer{number}_in_part.png")
    from PIL import Image

    Image.fromarray(np.hstack(shots)).save(out)
    print(out)
