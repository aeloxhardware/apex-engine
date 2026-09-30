"""
scripts/diagnose_gradient_variance.py

Diagnostic ad hoc : confirme la cause racine de l'espacement de couches
anormalement fin (jusqu'a 40x plus fin que layer_height_base) mesure dans
la zone festonnee de death.stl (couches 200-227, z~62-67mm).

Hypothese (lecture de code, pipeline/step5_layers.py::_representative_gradient_norm
et _generate_isovalues) : delta = layer_height_base * MEDIANE GLOBALE de
||grad(phi)|| (ponderee par volume, sur TOUT le TetMesh). L'espacement
PHYSIQUE reel entre deux isovaleurs consecutives est delta / ||grad(phi)||
LOCAL. Si le gradient local dans la zone festonnee est bien plus grand que
la mediane globale, l'espacement physique local devient bien plus fin que
layer_height_base - exactement le symptome mesure.

Ce script compare directement ||grad(phi)|| sur les tetraedres proches de
la zone festonnee (z entre 62 et 67mm, rayon XY > 150 - zone identifiee
precedemment) contre la mediane globale utilisee pour calibrer delta.

Usage :
    .venv/Scripts/python.exe scripts/diagnose_gradient_variance.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402

MODEL_PATH = "res/test_models/death.stl"


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    values_sorted = values[order]
    cum_weights = np.cumsum(weights[order])
    total = cum_weights[-1]
    if total <= 0.0:
        return float(np.median(values))
    idx = int(np.searchsorted(cum_weights, 0.5 * total))
    idx = min(idx, len(values_sorted) - 1)
    return float(values_sorted[idx])


def tet_volumes(tet_mesh) -> np.ndarray:
    verts = tet_mesh.tet_vertices[tet_mesh.tetrahedra]
    a, b, c, d = verts[:, 0], verts[:, 1], verts[:, 2], verts[:, 3]
    return np.abs(np.einsum("ij,ij->i", a - d, np.cross(b - d, c - d))) / 6.0


def main() -> None:
    params = PipelineParams()
    result = run_slicing_only(MODEL_PATH, params)

    tet_mesh = result.tet_mesh
    scalar_field = result.scalar_field

    norms = np.linalg.norm(scalar_field.gradient, axis=1)
    volumes = tet_volumes(tet_mesh)
    global_median = weighted_median(norms, volumes)
    print(f"Mediane globale ponderee de ||grad(phi)|| : {global_median:.4f}")
    print(f"Distribution globale : p10={np.percentile(norms, 10):.4f} "
          f"p50={np.percentile(norms, 50):.4f} p90={np.percentile(norms, 90):.4f} "
          f"p99={np.percentile(norms, 99):.4f} max={norms.max():.4f}")

    delta = params.layers.layer_height_base * global_median
    print(f"delta (espacement de champ) = layer_height_base({params.layers.layer_height_base}) "
          f"* mediane_globale({global_median:.4f}) = {delta:.4f}")

    tet_centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
    z = tet_centers[:, 2]
    r_xy = np.linalg.norm(tet_centers[:, :2], axis=1)

    # Zone festonnee identifiee precedemment : z entre 62 et 67mm
    zone_mask = (z >= 62.0) & (z <= 67.0)
    print(f"\nTetraedres dans la zone festonnee (z in [62,67]) : {zone_mask.sum()}/{len(z)}")
    if zone_mask.sum() > 0:
        zone_norms = norms[zone_mask]
        zone_volumes = volumes[zone_mask]
        zone_median = weighted_median(zone_norms, zone_volumes)
        print(f"  Mediane LOCALE ponderee de ||grad(phi)|| dans cette zone : {zone_median:.4f}")
        print(f"  Distribution locale : p50={np.percentile(zone_norms, 50):.4f} "
              f"p90={np.percentile(zone_norms, 90):.4f} p99={np.percentile(zone_norms, 99):.4f} "
              f"max={zone_norms.max():.4f}")
        print(f"  Ratio mediane_locale / mediane_globale = {zone_median / global_median:.2f}x")
        local_delta_z = delta / zone_median
        print(f"  Espacement Z physique attendu localement (delta/mediane_locale) = {local_delta_z:.4f}mm")
        print(f"  (compare a layer_height_base attendu = {params.layers.layer_height_base}mm, "
              f"et a l'espacement Z MINIMUM mesure precedemment = 0.0024mm)")

    # Top 1% des tetraedres a plus fort gradient, ou sont-ils ?
    top_mask = norms > np.percentile(norms, 99)
    print(f"\nTop 1% tetraedres a plus fort gradient ({top_mask.sum()} tets) : "
          f"z min={z[top_mask].min():.2f} max={z[top_mask].max():.2f} "
          f"r_xy min={r_xy[top_mask].min():.2f} max={r_xy[top_mask].max():.2f}")


if __name__ == "__main__":
    main()
