"""
Detection de REPLI (overfold) passe par passe -- candidat "juge" pour
automatiser le bareme manuel de la methode Bird.

POURQUOI : l'utilisateur arretait le bareme manuel quand la piece se
repliait sur elle-meme. Un repli n'est PAS un surplomb : aucune mesure
d'angle de face ne le voit. On mesure donc l'integrite du volume :

  - tetraedres RETOURNES : volume signe de signe oppose a l'original
    (la matiere passe a travers elle-meme) -- le repli au sens strict ;
  - tetraedres ECRASES : |volume| < 5 % de l'original (repli imminent) ;
  - aretes de surface en LAME DE COUTEAU : deux faces adjacentes dont les
    normales font plus de 150 deg -- le pli mesure sur real_overhang
    (angle diedre 177,9 deg) etait de ce type, sans forcement retourner
    de tetraedre.

Usage: python detect_folds.py <modele> [mult] [bareme separe par virgules]
  ex : python detect_folds.py "benchy upsidedown tilted" 1 30,20,10,5,0
"""
import sys
import time

import numpy as np
import trimesh

import bird_core

model_name = sys.argv[1] if len(sys.argv) > 1 else "real_overhang"
rotation_multiplier = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
bareme = [float(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [30, 20, 10, 5, 0]
# Densite du maillage (optionnelle) : plafond de volume par tetraedre, switch
# `-a` de TetGen. Incluse dans le nom des fichiers pour ne jamais melanger
# deux series de densites differentes.
max_tet_volume = float(sys.argv[4]) if len(sys.argv) > 4 else None

CRUSH_RATIO = 0.05
KNIFE_EDGE_DEG = 150.0


def signed_volumes(points, cells):
    p = points[cells]
    a, b, c, d = p[:, 0], p[:, 1], p[:, 2], p[:, 3]
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def knife_edges(tet):
    surf = tet.extract_surface(algorithm="dataset_surface").triangulate()
    faces = surf.faces.reshape(-1, 4)[:, 1:]
    mesh = trimesh.Trimesh(surf.points, faces, process=False)
    angles = np.degrees(mesh.face_adjacency_angles)
    return int((angles > KNIFE_EDGE_DEG).sum()), float(angles.max()) if len(angles) else 0.0


def main():
    safe_name = model_name.replace(" ", "_") + (f"_a{max_tet_volume:g}" if max_tet_volume else "")
    print(f"Modele : {model_name}   ROTATION_MULTIPLIER={rotation_multiplier}   bareme MAX_OVERHANG={bareme}")

    tet = bird_core.load_and_tetrahedralize(f"input_models/{model_name}.stl", max_tet_volume=max_tet_volume)
    tet, state = bird_core.calculate_tet_attributes(tet)
    cells = np.asarray(tet.field_data["cells"])
    # Positions sauvegardees a chaque passe : view_overlaps.py les relit sans
    # relancer la deformation (plusieurs minutes sur le Benchy).
    prefix = f"output_models/folds_{safe_name}_m{rotation_multiplier:g}"
    np.save(f"{prefix}_cells.npy", cells)
    np.save(f"{prefix}_pass0.npy", np.asarray(tet.points))
    vol0 = signed_volumes(np.asarray(tet.points), cells)
    ke0, kemax0 = knife_edges(tet)
    print(f"tetraedres : {len(cells)}   aretes >{KNIFE_EDGE_DEG:.0f} deg dans l'ORIGINAL : {ke0} (max {kemax0:.1f} deg)")
    print()
    print(f"{'passe':>5} {'MAX_OVH':>7} {'retournes':>10} {'ecrases':>8} {'lames':>6} {'diedre max':>11} "
          f"{'surplomb int.':>13} {'duree':>7}")

    for i, max_overhang in enumerate(bareme):
        params = {
            "NEIGHBOUR_LOSS_WEIGHT": 20,
            "MAX_OVERHANG": max_overhang,
            "ROTATION_MULTIPLIER": rotation_multiplier,
            "SET_INITIAL_ROTATION_TO_ZERO": False,
            "INITIAL_ROTATION_FIELD_SMOOTHING": 30,
            "MAX_POS_ROTATION_DEG": 360,
            "MAX_NEG_ROTATION_DEG": -360,
            "ROTATE_ITERATIONS": 100,
            "DEFORM_ITERATIONS": 1000,
            "STEEP_OVERHANG_COMPENSATION": True,
        }
        t0 = time.time()
        tet, info = bird_core.run_one_iteration(tet, state, params)
        dt = time.time() - t0

        vol = signed_volumes(np.asarray(tet.points), cells)
        flipped = int((np.sign(vol) != np.sign(vol0)).sum())
        crushed = int((np.abs(vol) < CRUSH_RATIO * np.abs(vol0)).sum())
        ke, kemax = knife_edges(tet)
        print(f"{i + 1:>5} {max_overhang:>7.0f} {flipped:>10} {crushed:>8} {ke:>6} {kemax:>10.1f}  "
              f"{info['overhang_fraction'] * 100:>11.2f} % {dt:>6.0f}s", flush=True)

        out = f"output_models/folds_{safe_name}_m{rotation_multiplier:g}_pass{i + 1}_ovh{max_overhang:g}.stl"
        bird_core.export_surface_stl(tet, out)
        np.save(f"{prefix}_pass{i + 1}.npy", np.asarray(tet.points))


if __name__ == "__main__":
    main()
