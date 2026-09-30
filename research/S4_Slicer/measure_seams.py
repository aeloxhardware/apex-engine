"""
P4, mesure prealable : RACCORDS crees par la coupe des zones de lissage aux
aretes vives (angle entre faces < seuil). Sur un STL grossier, la coupe peut
tomber entre deux facettes d'une surface voulue lisse (coque du Benchy).

Sortie : longueur des raccords, repartition de l'angle entre faces des aretes
coupees, et output_models/loop_<...>_seams_pass<k>.npz (couleur par morceau +
segments de raccord) pour render_seams.py (venv de la recherche).

Usage (venv d'Apex) :
  ../../.venv/Scripts/python.exe measure_seams.py <modele> <passe> <variante> [angle_entre_faces_deg]
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from measure_unfold import split_at_sharp_edges  # noqa: E402


def main():
    model, k, variant = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    limit = float(sys.argv[4]) if len(sys.argv) > 4 else 150.0
    prefix = os.path.join(HERE, "output_models", "loop_" + model.replace(" ", "_") + variant)
    skin = np.load(f"{prefix}_skin_pass{k}.npz")
    v, faces, normals, zone = skin["vertices"], skin["faces"], skin["normals"], skin["zone"]
    piece = np.full(len(faces), -1)
    n_pieces = 0
    for z in range(zone.max() + 1):
        for ids in split_at_sharp_edges(faces, normals, np.flatnonzero(zone == z), limit):
            piece[ids] = n_pieces
            n_pieces += 1
    # aretes partagees par deux faces
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(faces)), 3)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    e, o = edges[order], owner[order]
    same = np.all(e[1:] == e[:-1], axis=1)
    fa, fb, ed = o[:-1][same], o[1:][same], e[:-1][same]
    seam = (zone[fa] >= 0) & (zone[fa] == zone[fb]) & (piece[fa] != piece[fb])
    length = np.linalg.norm(v[ed[seam, 0]] - v[ed[seam, 1]], axis=1)
    cosn = np.clip(np.einsum("ij,ij->i", normals[fa[seam]], normals[fb[seam]]), -1, 1)
    face_angle = 180.0 - np.degrees(np.arccos(cosn))          # angle ENTRE LES FACES
    zone_border = (zone[fa] != zone[fb]) & ((zone[fa] >= 0) | (zone[fb] >= 0))
    border_len = np.linalg.norm(v[ed[zone_border, 0]] - v[ed[zone_border, 1]], axis=1).sum()
    print(f"{model}{variant} : {zone.max() + 1} zones -> {n_pieces} morceaux (seuil {limit:.0f} deg entre faces)")
    print(f"raccords crees par la coupe : {length.sum():.0f} mm (bords des zones d'origine : {border_len:.0f} mm)")
    for lo, hi in ((0, 120), (120, 135), (135, 145), (145, limit)):
        m = (face_angle >= lo) & (face_angle < hi)
        print(f"  aretes coupees a {lo}-{hi:.0f} deg entre faces : {length[m].sum():6.0f} mm ({length[m].sum() / max(length.sum(), 1e-9):.0%})")
    rng = np.random.default_rng(1)
    palette = rng.integers(60, 230, size=(n_pieces, 3))
    colors = np.where(piece[:, None] >= 0, palette[np.maximum(piece, 0)], 200)
    colors[zone < 0] = (235, 235, 235)
    np.savez(f"{prefix}_seams_pass{k}.npz", vertices=v, faces=faces, colors=colors.astype(np.uint8),
             seam_a=v[ed[seam, 0]], seam_b=v[ed[seam, 1]], seam_angle=face_angle)
    print(f"{prefix}_seams_pass{k}.npz")


if __name__ == "__main__":
    main()
