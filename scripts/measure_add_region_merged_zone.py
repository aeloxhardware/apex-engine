"""
scripts/measure_add_region_merged_zone.py

Suite de measure_add_region_n_zones.py : le bug limite #42 se generalise a
N regions simultanees (mesure confirmee : aucune des 3 zones fines
n'atteint le facteur de densite attendu face a "rest", une zone (junction_y)
est meme MOINS dense que "rest"). Ce script teste le fallback prevu par le
plan : fusionner les 3 zones de jonction en UNE SEULE region (2 add_region
au total, comme le cas plateau deja fonctionnel) plutot que N regions
separees.

Usage :
    .venv/Scripts/python.exe scripts/measure_add_region_merged_zone.py
"""

from __future__ import annotations

import numpy as np
import tetgen
import trimesh


def main() -> None:
    foot = trimesh.creation.box(extents=[4, 4, 20])
    foot.apply_translation([0, 0, 10])
    arm_x = trimesh.creation.box(extents=[20, 4, 4])
    arm_x.apply_translation([8, 0, 20])
    arm_y = trimesh.creation.box(extents=[4, 20, 4])
    arm_y.apply_translation([0, 8, 20])
    arm_diag = trimesh.creation.box(extents=[14, 14, 4])
    arm_diag.apply_translation([-7, -7, 20])

    mesh = foot.union(arm_x, engine="manifold")
    mesh = mesh.union(arm_y, engine="manifold")
    mesh = mesh.union(arm_diag, engine="manifold")

    coarse_volume = 1.5 ** 3
    fine_volume = 0.3 ** 3

    rest_seed = np.array([8.0, 0.0, 20.0])
    # Une seule region fine, centree sur le centroide des 3 zones -
    # couvre la jonction physique (sommet du pied, z~=20) d'un seul bloc.
    merged_seed = np.array([0.667, 0.667, 21.0])

    for name, seed in [("rest", rest_seed), ("merged_junction", merged_seed)]:
        print(f"seed {name} = {seed}  inside={mesh.contains([seed])[0]}")

    tet = tetgen.TetGen(mesh.vertices, mesh.faces)
    tet.add_region(1, rest_seed, coarse_volume)
    tet.add_region(2, merged_seed, fine_volume)
    tet.tetrahedralize(switches="pq1.2Aa")

    tet_vertices = np.asarray(tet.node)
    tetrahedra = np.asarray(tet.elem)
    tet_centers = tet_vertices[tetrahedra].mean(axis=1)

    print(f"\nTotal tetraedres : {len(tetrahedra)}")

    for name, seed in [("rest", rest_seed), ("merged_junction", merged_seed)]:
        dist = np.linalg.norm(tet_centers - seed, axis=1)
        local_count = int((dist < 3.0).sum())
        print(f"  {name:18s} seed={seed}  tetraedres a <3mm : {local_count}")

    # Verifie AUSSI la densite loin de toute region fine (dans le bras +Y,
    # jamais couvert par merged_seed ni par le rayon de rest) pour confirmer
    # que le "rest" garde bien SON volume grossier propre, pas ecrase par
    # le volume fin de la region fusionnee.
    far_point = np.array([0.0, 15.0, 20.0])  # milieu du bras +Y, loin de la jonction
    print(f"far_point (bras +Y, loin de tout) = {far_point}  "
          f"inside={mesh.contains([far_point])[0]}")
    dist_far = np.linalg.norm(tet_centers - far_point, axis=1)
    print(f"  far_point        tetraedres a <3mm : {int((dist_far < 3.0).sum())}")


if __name__ == "__main__":
    main()
