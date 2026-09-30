"""
Fenetre 3D des recouvrements entre tetraedres NON VOISINS, passe par passe.

Affiche la piece (surface grise translucide) et, en couleur, les tetraedres
qui recouvrent un tetraedre eloigne :
  - ROUGE  : recouvrement avec un tetraedre de rang >= 3 (vrai repli :
             aucun lien proche dans le maillage d'origine) ;
  - ORANGE : rang 2 seulement (voisin d'un voisin -- zone grise, souvent
             simple compression locale).

Curseur : numero de passe. Touche "o" : bascule entre l'espace DEFORME (la
ou le recouvrement a lieu) et la piece REELLE (pour voir OU il se trouve
sur l'objet imprime).

Relit les positions sauvees par detect_folds.py, ne relance aucune
deformation.

Usage: python view_overlaps.py <modele> [mult]
  ex : python view_overlaps.py "benchy upsidedown tilted" 1
"""
import glob
import re
import sys

import numpy as np
import pyvista as pv

import tet_overlap

model_name = sys.argv[1] if len(sys.argv) > 1 else "real_overhang"
rotation_multiplier = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0

prefix = f"output_models/folds_{model_name.replace(' ', '_')}_m{rotation_multiplier:g}"
cells = np.load(f"{prefix}_cells.npy")
passes = sorted(int(re.search(r"_pass(\d+)\.npy$", p).group(1)) for p in glob.glob(f"{prefix}_pass*.npy"))
points = {k: np.load(f"{prefix}_pass{k}.npy") for k in passes}

print(f"{model_name} : {len(cells)} tetraedres, passes disponibles {passes}")
results = {}
for k in passes:
    r = tet_overlap.non_neighbour_overlaps(points[k], cells)
    results[k] = r
    a, b = r["rang>=2"], r["rang>=3"]
    print(f"  passe {k} : rang>=3 {b['tets']:>6} tetraedres ({b['pairs']} paires, penetration max "
          f"{b['depth_max']:.3f} mm) | rang 2 inclus {a['tets']:>6} tetraedres")

cell_block = np.hstack([np.full((len(cells), 1), 4), cells]).ravel()
cell_types = np.full(len(cells), pv.CellType.TETRA)
view = {"pass": passes[-1], "real_space": False}

plotter = pv.Plotter(window_size=(1000, 750), title=f"Recouvrements non voisins - {model_name}")


def redraw():
    k = view["pass"]
    shown = points[0] if view["real_space"] else points[k]
    grid = pv.UnstructuredGrid(cell_block, cell_types, shown)
    plotter.add_mesh(grid.extract_surface(algorithm="dataset_surface"), color="lightgrey",
                     opacity=0.18, name="surface")

    r = results[k]
    far = r["rang>=3"]["tet_ids"]
    near_only = np.setdiff1d(r["rang>=2"]["tet_ids"], far)
    for name, ids, color in (("far", far, "red"), ("near", near_only, "orange")):
        if len(ids):
            plotter.add_mesh(grid.extract_cells(ids), color=color, show_edges=True, name=name)
        else:
            plotter.remove_actor(name)

    espace = "piece REELLE" if view["real_space"] else "espace DEFORME"
    plotter.add_text(
        f"Passe {k}  -  {espace}  (touche o)\n"
        f"rouge  rang>=3 : {r['rang>=3']['tets']} tetraedres, penetration max {r['rang>=3']['depth_max']:.3f} mm\n"
        f"orange rang 2  : {len(near_only)} tetraedres\n"
        f"total : {len(cells)} tetraedres",
        position="upper_left", font_size=10, name="info")


def on_slider(value):
    k = int(round(value))
    if k != view["pass"] and k in points:
        view["pass"] = k
        redraw()


def toggle_space():
    view["real_space"] = not view["real_space"]
    redraw()


redraw()
plotter.add_slider_widget(on_slider, [min(passes), max(passes)], value=view["pass"], title="passe",
                          fmt="%.0f", pointa=(0.35, 0.08), pointb=(0.9, 0.08))
plotter.add_key_event("o", toggle_space)
plotter.show()
