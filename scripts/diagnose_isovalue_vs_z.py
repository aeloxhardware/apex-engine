"""
scripts/diagnose_isovalue_vs_z.py

Diagnostic non-monotonie du champ scalaire (CLAUDE.md limite #47) : mesure
les paires de couches consecutives (triees par Z) dont l'ecart Z est
anormalement petit ET l'ecart d'isovaleur anormalement grand - signature
d'une jonction multi-region ou phi "ondule" en Z plutot que d'etre
monotone.

Utilise desormais le diagnostic structure permanent
(LayerParams.detect_junction_variance, Layer.junction_variance_score,
cf. limite #48) plutot qu'un recalcul manuel - un seul endroit de verite
pour la metrique, partage entre ce script et pipeline/step5_layers.py.

Historique : ce script a d'abord servi a INFIRMER l'hypothese initiale (le
clustering spatial split_disjoint_clusters n'etait pas la cause - il est
desactive par defaut sur ce run) et a CONFIRMER la vraie cause (le champ
resolu par step4 n'est pas monotone en Z a la jonction). Voir CLAUDE.md
limite #47 pour le diagnostic complet.

Usage :
    .venv/Scripts/python.exe scripts/diagnose_isovalue_vs_z.py
    .venv/Scripts/python.exe scripts/diagnose_isovalue_vs_z.py --assert-max-anomalies 0
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import replace

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

MODEL_PATH = "res/test_models/death.stl"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--assert-max-anomalies", type=int, default=None,
        help="Echoue (exit 1) si le nombre de couches marquees anormales depasse N. "
             "Sert de non-regression exploitable en script/CI.",
    )
    args = parser.parse_args(argv)

    params = PipelineParams()
    params.layers = replace(params.layers, detect_junction_variance=True)

    result = run_slicing_only(MODEL_PATH, params)
    layers = result.colored_layers

    anomalous = [layer for layer in layers if layer.junction_variance_score is not None]
    print(f"\n=== RESULTAT : {len(anomalous)}/{len(layers)} couches marquees anormales "
          f"(detect_junction_variance=True, seuils par defaut) ===")

    if anomalous:
        by_id = {layer.id: layer for layer in layers}
        layer_z = {}
        for layer in layers:
            all_z = np.concatenate([np.asarray(loop)[:, 2] for loop in layer.contour_3d])
            layer_z[layer.id] = float(all_z.mean())

        anomalous_sorted = sorted(anomalous, key=lambda l: -l.junction_variance_score)
        print("\nTop 15 couches anormales (score decroissant) :")
        for layer in anomalous_sorted[:15]:
            print(f"  layer {layer.id:4d}  z={layer_z[layer.id]:8.4f}  "
                  f"isovalue={layer.isovalue:8.4f}  regions={layer.region_ids}  "
                  f"score={layer.junction_variance_score:10.2f}")

    if args.assert_max_anomalies is not None:
        if len(anomalous) > args.assert_max_anomalies:
            print(f"\nECHEC : {len(anomalous)} anomalies > seuil autorise "
                  f"{args.assert_max_anomalies}")
            return 1
        print(f"\nOK : {len(anomalous)} anomalies <= seuil autorise "
              f"{args.assert_max_anomalies}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
