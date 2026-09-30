"""
Banc de mesure de mesh_flatten.py sur tous les modeles disponibles.

Avant / apres aplatissement : tetraedres (TetGen, reglages par defaut de
Bird), tetraedres microscopiques (< 1e-4 mm3), deplacement maximal des
sommets, etancheite et auto-intersections (step1_ingest d'Apex).

Usage (venv d'Apex) : APEX_USE_CPP=1 ../../.venv/Scripts/python.exe test_mesh_flatten.py
"""
import glob
import os
import sys

import numpy as np
import tetgen
import trimesh

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..")))
from pipeline.step1_ingest import _find_self_intersecting_faces  # noqa: E402

import mesh_flatten  # noqa: E402

LAYER = 0.2
MICRO = 1e-4


def tets(mesh):
    t = tetgen.TetGen(np.asarray(mesh.vertices), np.asarray(mesh.faces))
    t.tetrahedralize()
    p, c = np.asarray(t.node), np.asarray(t.elem)[:, :4]      # pas de pyvista dans ce venv
    a, b, cc, d = (p[c[:, i]] for i in range(4))
    vol = np.abs(np.einsum("ij,ij->i", np.cross(b - a, cc - a), d - a)) / 6
    return len(c), int((vol < MICRO).sum())


def main():
    paths = sorted(glob.glob(os.path.join(HERE, "input_models", "*.stl")))
    skip = {"death.stl", "tilted benchy undertriangulated.stl", "tilted benchy.stl"}   # versions reparees ou non tetraedrisables
    print(f"{'modele':<42} {'avant':>14} {'dessous aplati':>14}   {'depl. max':>8}"
          f"  (tetraedres / microscopiques ; '!' = auto-intersection ou non etanche)")
    for path in paths:
        name = os.path.basename(path)
        if name in skip:
            continue
        mesh = trimesh.load(path)
        try:
            before = tets(mesh)
        except RuntimeError:
            before = ("refus", "-")
        bottom, d_bottom = mesh_flatten.flatten_bottom(mesh, LAYER)
        cells = []
        for m in (bottom,):
            crossing = len(_find_self_intersecting_faces(m))
            try:
                n, micro = tets(m)
                cells.append(f"{n}/{micro}" + ("" if crossing == 0 and m.is_watertight else "!"))
            except RuntimeError:
                cells.append(f"refus({crossing}x)")
        print(f"{name[:-4]:<42} {before[0]:>8}/{before[1]:<5} {cells[0]:>14}   {d_bottom.max(initial=0):>7.3f}m")

if __name__ == "__main__":
    main()
