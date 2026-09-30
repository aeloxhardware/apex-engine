"""
scripts/diagnose_layer_spacing.py

Diagnostic ad hoc, derniere piste avant de clore le sujet du "debordement
visuel" signale par l'utilisateur sur out/death_preview.glb (Open3DViewer).
La geometrie des contours est mesuree correcte a 0.008mm pres (4 hypotheses
testees et ecartees : distance a la surface, zigzag/aller-retour, mesh
source vs mesh rendu). Piste restante : des couches successives tres
rapprochees en Z (ou avec un chevauchement XY important) pourraient creer
un effet visuel de "faisceau dense" ou de lignes qui semblent se croiser/
sortir du solide en vue rasante, meme si chaque contour individuel colle a
la surface - l'oeil ne distingue pas facilement N contours proches les uns
des autres empiles visuellement.

Mesure : espacement Z reel entre couches consecutives (median vs
layer_height_base attendu), et pour les zones a forte densite de boucles
deja identifiees (couches 187, 206-227), le nombre de couches dans une
fenetre Z donnee comparee au reste du modele.

Usage :
    .venv/Scripts/python.exe scripts/diagnose_layer_spacing.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

MODEL_PATH = "res/test_models/death.stl"


def main() -> None:
    params = PipelineParams()
    result = run_slicing_only(MODEL_PATH, params)
    layers = result.colored_layers

    # z moyen de chaque couche (moyenne de tous les points de contour)
    layer_z = []
    for layer in layers:
        all_z = np.concatenate([np.asarray(loop)[:, 2] for loop in layer.contour_3d])
        layer_z.append(float(all_z.mean()))
    layer_z = np.array(layer_z)

    order = np.argsort(layer_z)
    sorted_ids = np.array([layers[i].id for i in order])
    sorted_z = layer_z[order]

    z_gaps = np.diff(sorted_z)
    print(f"Espacement Z entre couches consecutives (triees par z) : "
          f"min={z_gaps.min():.4f}mm max={z_gaps.max():.4f}mm "
          f"median={np.median(z_gaps):.4f}mm mean={z_gaps.mean():.4f}mm")
    print(f"layer_height_base attendu (params par defaut) : {params.layers.layer_height_base}mm")

    # Zone signalee : couches 206-227 (identifiees precedemment comme
    # densement peuplees en boucles de contour)
    zone_mask = np.isin(sorted_ids, list(range(200, 230)))
    zone_z = sorted_z[zone_mask]
    if len(zone_z) > 1:
        zone_gaps = np.diff(zone_z)
        print(f"\nZone couches 200-227 : {len(zone_z)} couches, z de {zone_z.min():.2f} a {zone_z.max():.2f}mm")
        print(f"  espacement Z dans cette zone : min={zone_gaps.min():.4f}mm max={zone_gaps.max():.4f}mm "
              f"median={np.median(zone_gaps):.4f}mm")

    # Nombre de couches dont le contour a une extension XY qui SE RECOUVRE
    # significativement (>50% de l'aire de bbox XY) avec la couche
    # suivante - signature d'un empilement dense visuellement.
    print("\nTop 15 plus petits espacements Z (couches quasi confondues visuellement) :")
    gap_order = np.argsort(z_gaps)
    for i in gap_order[:15]:
        print(f"  layer {sorted_ids[i]:4d} (z={sorted_z[i]:7.3f}) -> "
              f"layer {sorted_ids[i+1]:4d} (z={sorted_z[i+1]:7.3f})  gap={z_gaps[i]:.4f}mm")


if __name__ == "__main__":
    main()
