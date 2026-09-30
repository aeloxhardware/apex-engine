"""
scripts/diagnose_loop_extent.py

Diagnostic ad hoc : le diagnostic precedent (diagnose_contour_overflow.py)
mesure la distance de CHAQUE POINT a la surface la plus proche, et trouve
un ecart maximal infime (0.008mm) - mais l'utilisateur conteste ce resultat
en montrant des captures ou l'ecart visuel est clairement bien plus grand.

Hypothese non testee jusqu'ici : une boucle de contour peut avoir des points
INDIVIDUELLEMENT proches de LEUR face la plus proche (le point precedent
mesure une distance ponctuelle), tout en dessinant un trace GLOBAL qui
s'eloigne largement de la zone locale attendue - typiquement un
"aller-retour" en pointe qui longe presque la surface sur tout son
parcours (donc distance ponctuelle faible partout) mais dont l'extremite
de la pointe est en realite loin de la position 3D correcte pour cette
zone du modele. Un point peut etre "proche d'UNE face" sans etre proche de
LA bonne face pour cette partie du contour, si la face la plus proche
trouvee est sur une autre portion de la surface (repli de geometrie
complexe/festonnee).

Ce script mesure directement, pour chaque boucle de contour, son ETENDUE
GEOMETRIQUE PROPRE (bounding box, longueur totale du perimetre) comparee a
celle du "patch" local du mesh source le plus proche - pour detecter un
contour qui zigzague sur une distance totale bien plus grande que la
portion de surface locale qu'il est cense suivre.

Usage :
    .venv/Scripts/python.exe scripts/diagnose_loop_extent.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import trimesh

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

MODEL_PATH = "res/test_models/death.stl"


def loop_perimeter(loop: np.ndarray) -> float:
    diffs = np.diff(loop, axis=0)
    return float(np.linalg.norm(diffs, axis=1).sum())


def loop_bbox_diag(loop: np.ndarray) -> float:
    bmin = loop.min(axis=0)
    bmax = loop.max(axis=0)
    return float(np.linalg.norm(bmax - bmin))


def main() -> None:
    params = PipelineParams()
    result = run_slicing_only(MODEL_PATH, params)
    layers = result.colored_layers

    source_mesh = trimesh.load(MODEL_PATH, force="mesh")
    prox = trimesh.proximity.ProximityQuery(source_mesh)

    # Pour chaque boucle, mesure : distance MAX au maillage source (deja
    # fait avant, ici en complement), le perimetre de la boucle, sa
    # bbox-diagonale, et le ratio perimetre/bbox-diagonale (un contour
    # "normal" suit une courbe simple ; un contour qui zigzague en
    # aller-retours aura un perimetre bien plus grand que sa bbox ne le
    # laisserait supposer).
    results = []
    for layer in layers:
        for loop_idx, loop in enumerate(layer.contour_3d):
            pts = np.asarray(loop)
            if len(pts) >= 2 and np.allclose(pts[0], pts[-1]):
                pts_open = pts[:-1]
            else:
                pts_open = pts
            if len(pts_open) < 3:
                continue

            _, dist, _ = prox.on_surface(pts_open)
            perim = loop_perimeter(pts)
            bbox_diag = loop_bbox_diag(pts)
            ratio = perim / bbox_diag if bbox_diag > 1e-9 else 0.0

            results.append({
                "layer_id": layer.id,
                "loop_idx": loop_idx,
                "n_points": len(pts_open),
                "max_dist_to_surface": float(dist.max()),
                "perimeter": perim,
                "bbox_diag": bbox_diag,
                "ratio": ratio,
            })

    # Top boucles par ratio perimetre/bbox (signature d'un zigzag
    # aller-retour, independant de la distance ponctuelle a la surface)
    results_sorted = sorted(results, key=lambda r: -r["ratio"])
    print("Top 20 boucles par ratio perimetre/bbox_diagonale (zigzag aller-retour) :")
    for r in results_sorted[:20]:
        print(f"  layer={r['layer_id']:4d} loop={r['loop_idx']:2d} n_pts={r['n_points']:4d} "
              f"perim={r['perimeter']:8.2f}mm bbox_diag={r['bbox_diag']:8.3f}mm "
              f"ratio={r['ratio']:8.1f} max_dist_surface={r['max_dist_to_surface']:.4f}mm")

    # Top boucles par distance max a la surface, pour recouper avec le
    # diagnostic precedent.
    results_by_dist = sorted(results, key=lambda r: -r["max_dist_to_surface"])
    print("\nTop 20 boucles par distance max a la surface :")
    for r in results_by_dist[:20]:
        print(f"  layer={r['layer_id']:4d} loop={r['loop_idx']:2d} n_pts={r['n_points']:4d} "
              f"max_dist_surface={r['max_dist_to_surface']:.4f}mm perim={r['perimeter']:8.2f}mm "
              f"bbox_diag={r['bbox_diag']:8.3f}mm ratio={r['ratio']:8.1f}")

    all_ratios = np.array([r["ratio"] for r in results])
    print(f"\nDistribution des ratios perimetre/bbox : p50={np.median(all_ratios):.2f} "
          f"p90={np.percentile(all_ratios, 90):.2f} p99={np.percentile(all_ratios, 99):.2f} "
          f"max={all_ratios.max():.2f}")


if __name__ == "__main__":
    main()
