"""
Voie directe Bird -> Apex : la hauteur DEFORMEE devient le champ scalaire phi
d'Apex, et l'etape 5 en extrait les couches courbes. Pas de Cura, pas de
retriangulation, pas de reprojection de G-code (decide avec l'utilisateur le
2026-09-26).

Pourquoi ca marche : Cura tranche la piece deformee par des plans z = cste.
Dans la piece d'origine, ces plans sont exactement les surfaces ou la hauteur
deformee est constante -- les isosurfaces de phi. La deformation de Bird
remplace donc les etapes 2 a 4 d'Apex (surplombs, segmentation, champ).

Deux maillages differents, relies par interpolation :
  - celui de Bird (bird_loop.py, ~33 000 tetraedres sur le Benchy), sur
    lequel la deformation a ete calculee ;
  - celui d'Apex (etape 1, plus fin), sur lequel l'etape 5 travaille.
Chaque noeud d'Apex est localise dans un tetraedre de Bird (coordonnees
barycentriques) et recoit la hauteur deformee interpolee. Faire tourner Bird
directement sur le maillage d'Apex couterait des dizaines de minutes par
passe.

L'etape 5 ordonne par rang de region PUIS par isovaleur. On lui donne une
seule region couvrant toute la piece : seul phi decide de l'ordre, puisque
c'est la deformation qui a organise l'arborescence.

A lancer avec le venv d'APEX :
  APEX_USE_CPP=1 ../../.venv/Scripts/python.exe bird_to_apex.py <modele> <passe>
  ex : bird_to_apex.py "benchy upsidedown tilted" 4
"""
import os
import sys
import time

import numpy as np
import trimesh
from scipy.spatial import cKDTree

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from core.mesh_types import Region, RegionDependencyGraph, ScalarField  # noqa: E402
from main import PipelineParams  # noqa: E402
from pipeline.step1_ingest import run as ingest_run  # noqa: E402
from pipeline.step2_overhang import run as overhang_run  # noqa: E402
from pipeline.step5_layers import run as layers_run  # noqa: E402
from pipeline.step6_toolpath2d import _lift_to_3d, run as toolpath2d_run  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def bird_offset(stl_path):
    """Decalage applique par bird_core.load_and_tetrahedralize : milieu en XY,
    pose a z=0 (bornes du STL, que TetGen conserve)."""
    v = trimesh.load(stl_path).vertices
    lo, hi = v.min(axis=0), v.max(axis=0)
    return np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, lo[2]])


def frame_offset(stl_path, prefix):
    """Repere du plateau (celui de pass0) -> repere du STL. Tient compte du
    decalage sur le plateau choisi a la boucle Bird (bird_loop.py --offset,
    sauvegarde dans <prefix>_shift.npy) : pass0 = recentre + decalage."""
    import os
    shift = np.load(f"{prefix}_shift.npy") if os.path.exists(f"{prefix}_shift.npy") else np.zeros(3)
    return bird_offset(stl_path) - shift


def barycentric(points, tets_xyz):
    """Coordonnees barycentriques de chaque point dans son tetraedre, (K,4)."""
    a = tets_xyz[:, 0]
    m = np.stack([tets_xyz[:, 1] - a, tets_xyz[:, 2] - a, tets_xyz[:, 3] - a], axis=2)  # (K,3,3)
    lam = np.linalg.solve(m, (points - a)[:, :, None])[:, :, 0]
    return np.column_stack([1.0 - lam.sum(axis=1), lam])


def transport_height(apex_nodes, bird_p0, bird_cells, bird_height, k=24):
    """Hauteur deformee interpolee en chaque noeud d'Apex. Parmi les k
    tetraedres de Bird les plus proches (par centre), on garde celui qui
    contient le point, ou a defaut celui dont il est le moins exterieur
    (noeuds de surface, a une fraction de tetraedre pres)."""
    tets = bird_p0[bird_cells]
    # k voisins : core.spatial.KDTree ne couvre que k=1 (le seul cas du pipeline).
    _, cand = cKDTree(tets.mean(axis=1)).query(apex_nodes, k=k)
    best_score = np.full(len(apex_nodes), -np.inf)
    best_bary = np.zeros((len(apex_nodes), 4))
    best_tet = np.zeros(len(apex_nodes), dtype=np.int64)
    for j in range(k):
        t = cand[:, j]
        bary = barycentric(apex_nodes, tets[t])
        score = bary.min(axis=1)             # >= 0 : le point est dedans
        better = score > best_score
        best_score[better], best_bary[better], best_tet[better] = score[better], bary[better], t[better]
    height = np.einsum("ij,ij->i", best_bary, bird_height[bird_cells[best_tet]])
    return height, best_score


def tet_gradients(nodes, tetrahedra, values):
    """Gradient exact d'un champ lineaire par tetraedre (P1), (T,3)."""
    p = nodes[tetrahedra]
    m = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=1)  # lignes = aretes
    dv = values[tetrahedra[:, 1:]] - values[tetrahedra[:, :1]]
    return np.linalg.solve(m, dv[:, :, None])[:, :, 0]


def main():
    model_name = sys.argv[1] if len(sys.argv) > 1 else "benchy upsidedown tilted"
    pass_index = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    stl = os.path.join(HERE, "input_models", f"{model_name}.stl")
    prefix = os.path.join(HERE, "output_models", "loop_" + model_name.replace(" ", "_"))

    bird_cells = np.load(f"{prefix}_cells.npy")
    bird_p0 = np.load(f"{prefix}_pass0.npy")
    bird_pk = np.load(f"{prefix}_pass{pass_index}.npy")

    params = PipelineParams()
    t0 = time.time()
    mesh_data, tet_mesh, report = ingest_run(
        stl, params.resolution,
        near_platform_resolution=params.near_platform_resolution,
        near_platform_margin=params.near_platform_margin,
    )
    print(f"etape 1 : {tet_mesh.tet_count} tetraedres Apex ({time.time() - t0:.0f} s), "
          f"maillage Bird : {len(bird_cells)} tetraedres")

    nodes_in_bird = tet_mesh.tet_vertices - bird_offset(stl)
    phi, inside = transport_height(nodes_in_bird, bird_p0, bird_cells, bird_pk[:, 2])
    outside = inside < -1e-6
    print(f"transport : {outside.sum()} noeuds hors du maillage de Bird sur {len(phi)} "
          f"(pire : {inside.min():.3f} en coordonnee barycentrique)")

    gradient = tet_gradients(tet_mesh.tet_vertices, tet_mesh.tetrahedra, phi)
    norms = np.linalg.norm(gradient, axis=1)
    print(f"|grad phi| : mediane {np.median(norms):.3f} (1 = couches d'epaisseur exacte), "
          f"5e-95e centiles {np.percentile(norms, 5):.3f}-{np.percentile(norms, 95):.3f}")
    scalar_field = ScalarField(values=phi, gradient=gradient)

    graph = RegionDependencyGraph(
        regions={0: Region(id=0, face_ids=np.arange(mesh_data.face_count),
                           target_direction=np.array([0.0, 0.0, 1.0]))},
        dependencies=[],
    )

    t0 = time.time()
    layers = layers_run(tet_mesh, mesh_data, scalar_field, graph, params.layers)
    print(f"etape 5 : {len(layers)} couches ({time.time() - t0:.0f} s), "
          f"hauteur visee {params.layers.layer_height_base} mm")

    # Etape 6 : parois + remplissage dans chaque couche. L'etape 2 ne sert ici
    # qu'aux ponts (reorientation du remplissage), pas au champ.
    t0 = time.time()
    overhang_map = overhang_run(mesh_data, params.overhang)
    layers = toolpath2d_run(layers, mesh_data, overhang_map, params.toolpath2d)
    n_paths = sum(len(layer.toolpath_2d or []) for layer in layers)
    print(f"etape 6 : {n_paths} polylignes (parois + remplissage) ({time.time() - t0:.0f} s), "
          f"{params.toolpath2d.wall_count} parois, remplissage {params.toolpath2d.infill_density:.0%}")

    np.save(f"{prefix}_apex_phi_pass{pass_index}.npy", phi)
    # Couches sauvees pour render_layers_gif.py (venv de la recherche : le venv
    # d'Apex n'a ni matplotlib ni Pillow). Ordre de la liste = ordre d'impression.
    import pickle
    out = f"{prefix}_apex_layers_pass{pass_index}.pkl"
    with open(out, "wb") as f:
        pickle.dump([{"contour_3d": [np.asarray(c) for c in layer.contour_3d],
                      "isovalue": getattr(layer, "isovalue", None),
                      # trajectoire de l'etape 6 relevee en 3D ; ordre de
                      # l'etape 6 : remplissage, parois interieures, contour.
                      "toolpath_2d": [np.asarray(t) for t in (layer.toolpath_2d or [])],
                      "toolpath_3d": [_lift_to_3d(np.asarray(t), layer.local_frame)
                                      for t in (layer.toolpath_2d or [])]} for layer in layers], f)
    print(f"couches sauvees : {out}")


if __name__ == "__main__":
    main()
