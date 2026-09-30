"""
scripts/measure_add_region_n_zones.py

Mesure ad hoc AVANT d'ecrire le code definitif du Livrable 2 (limite #48) :
le bug de la limite #42 (TetGen.add_region ignore silencieusement le
volume cible d'une region dont le seed tombe hors de sa coque, appliquant
alors le volume d'une AUTRE region a tout le maillage) n'a jamais ete
mesure au-dela de 2 regions simultanees. Ce script teste directement 4
regions actives en meme temps (1 "rest" + 3 zones de jonction, sur le
tripode) et verifie si CHACUNE obtient bien son propre volume cible, ou si
le bug se generalise.

Methode de mesure : construire un PLC simple (le mesh tripode complet, pas
de decoupe de surface - les zones de jonction sont des seeds internes
purs), ajouter 1 region "rest" (volume grossier) + 3 seeds internes a
volume FIN tres different (facteur ~30x), tetraedriser, puis mesurer la
densite locale de tetraedres pres de chaque seed vs loin de tous les
seeds - si le bug se generalise, au moins un seed n'aura aucun effet
detectable sur la densite locale.

Usage :
    .venv/Scripts/python.exe scripts/measure_add_region_n_zones.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import tetgen
import trimesh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def main() -> None:
    # Meme geometrie que la fixture tripode (tests/test_step5_layers.py) -
    # jonction physique au sommet du pied, z~=20.
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

    coarse_volume = 1.5 ** 3   # "rest"
    fine_volume = 0.3 ** 3     # zones de jonction (~30x plus fin)

    # 3 seeds internes distincts, un par bras pres du sommet du pied (zone
    # de jonction reelle), + 1 seed "rest" loin de la jonction (ex: bout
    # d'un bras).
    seeds = {
        "rest": np.array([8.0, 0.0, 20.0]),      # milieu du bras +X (extent Z 18..22), loin de la jonction
        "junction_x": np.array([4.0, 0.0, 21.0]),
        "junction_y": np.array([0.0, 4.0, 21.0]),
        "junction_diag": np.array([-2.0, -2.0, 21.0]),
    }

    # Verifie que chaque seed est bien a l'interieur du mesh (sinon la
    # mesure ne teste rien du bug vise).
    for name, seed in seeds.items():
        inside = mesh.contains([seed])[0]
        print(f"seed {name} = {seed}  inside={inside}")

    tet = tetgen.TetGen(mesh.vertices, mesh.faces)
    tet.add_region(1, seeds["rest"], coarse_volume)
    tet.add_region(2, seeds["junction_x"], fine_volume)
    tet.add_region(3, seeds["junction_y"], fine_volume)
    tet.add_region(4, seeds["junction_diag"], fine_volume)
    tet.tetrahedralize(switches="pq1.2Aa")

    tet_vertices = np.asarray(tet.node)
    tetrahedra = np.asarray(tet.elem)
    tet_centers = tet_vertices[tetrahedra].mean(axis=1)

    print(f"\nTotal tetraedres : {len(tetrahedra)}")

    # Densite locale = nombre de tetraedres dont le centre est a moins de
    # 3mm de chaque seed - une region fine doit avoir BEAUCOUP plus de
    # tetraedres locaux qu'une region grossiere, a rayon egal.
    for name, seed in seeds.items():
        dist = np.linalg.norm(tet_centers - seed, axis=1)
        local_count = int((dist < 3.0).sum())
        print(f"  {name:15s} seed={seed}  tetraedres a <3mm : {local_count}")

    print(
        "\nSi les 3 zones 'junction_*' ont une densite locale nettement "
        "superieure a 'rest' (facteur proche de fine_volume/coarse_volume "
        f"= {coarse_volume / fine_volume:.1f}x en nombre de tetraedres), "
        "le bug limite #42 NE SE GENERALISE PAS a N regions - chaque seed "
        "garde son propre volume cible. Sinon, le bug se generalise et le "
        "plan doit fusionner les zones en une seule region plutot que N."
    )


if __name__ == "__main__":
    main()
