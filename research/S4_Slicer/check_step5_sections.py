"""
Defaut de l'etape 5 d'Apex : sur le tranchage de reference (corps non creuse),
compare, couche par couche, la region delimitee par les contours de step5
(parite) a l'UNION DES SECTIONS REELLES des tetraedres coupes par z = h, dans
l'espace deforme. Les deux doivent coincider.

Usage (venv d'Apex) :
  ../../.venv/Scripts/python.exe -u check_step5_sections.py <modele> <passe> <variante> [couche,...]
"""
import os
import sys

import numpy as np
import shapely
import trimesh
from shapely.geometry import Polygon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bird_slice_deformed as bsd  # noqa: E402
import skin_carve  # noqa: E402


def parity_region(loops):
    out = None
    for c in loops:
        c = np.asarray(c)
        if len(c) < 4:
            continue
        p = Polygon(c[:, :2]).buffer(0)
        out = p if out is None else out.symmetric_difference(p)
    return out


def main():
    model, k, variant = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    only = {int(x) for x in sys.argv[4].split(",")} if len(sys.argv) > 4 else None
    prefix = os.path.join(HERE, "output_models", "loop_" + model.replace(" ", "_") + variant)
    cells = np.load(f"{prefix}_cells.npy")
    p0 = np.load(f"{prefix}_pass0.npy")
    pk = np.load(f"{prefix}_pass{k}.npy")
    faces = bsd.boundary_faces(cells, p0)
    surface = trimesh.Trimesh(pk, faces, process=False)
    surface.remove_unreferenced_vertices()
    mesh_data = bsd._build_mesh_data(surface)
    tet_mesh = bsd.TetMesh(tet_vertices=pk, tetrahedra=cells)
    z = pk[:, 2].copy()
    field = bsd.ScalarField(values=z, gradient=np.tile([0.0, 0.0, 1.0], (len(cells), 1)))
    graph = bsd.RegionDependencyGraph(
        regions={0: bsd.Region(id=0, face_ids=np.arange(mesh_data.face_count),
                               target_direction=np.array([0.0, 0.0, 1.0]))}, dependencies=[])
    layers = bsd.layers_run(tet_mesh, mesh_data, field, graph, bsd.PipelineParams().layers)
    layers.sort(key=lambda L: L.isovalue)
    psi = np.ones(len(pk))      # inf donnerait NaN (inf*0) : sections perdues
    rows = []
    for i, L in enumerate(layers):
        if only and (i + 1) not in only:
            continue
        h = L.isovalue
        truth = shapely.union_all([Polygon(xy) for xy, _pv, _t in skin_carve._sections(cells, pk, psi, h)])
        mine = parity_region(L.contour_3d)
        ta = truth.area
        ma = 0.0 if mine is None else mine.area
        lost = truth.difference(mine).area if mine is not None else ta
        added = mine.difference(truth).area if mine is not None else 0.0
        on_vertex = np.abs(z - h) < 1e-9
        rows.append((i + 1, h, ta, ma, lost, added, len(L.contour_3d), int(on_vertex.sum())))
    bad = [r for r in rows if (r[4] / r[2] > 0.01 or r[5] / max(r[2], 1e-9) > 0.01)]
    print(f"{len(layers)} couches ; {len(bad)} avec > 1 % d'ecart")
    for r in (bad if not only else rows):
        print(f"  couche {r[0]} z={r[1]:.4f} : sections {r[2]:.2f} mm2, step5 {r[3]:.2f} mm2, "
              f"perdu {r[4]:.2f} ({r[4] / max(r[2], 1e-9):.0%}), ajoute {r[5]:.2f}, "
              f"{r[6]} boucles, {r[7]} noeuds exactement sur le plan")


if __name__ == "__main__":
    main()
