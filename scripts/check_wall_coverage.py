"""
Couverture de la BANDE DE PAROIS par les cordons reellement extrudes, couche
par couche, dans le repere local de step6. Bande = matiere - matiere retrecie
de wall_count * nozzle_width ; cordons = chaque polyligne de toolpath_2d
dilatee d'une demi-largeur (step8 extrude tous les segments d'une polyligne,
jamais le saut entre deux polylignes). Un trou dans une paroi (limite #65)
apparait comme une couverture < 100 %.
Usage : APEX_USE_CPP=1 .venv/Scripts/python.exe -u scripts/check_wall_coverage.py <stl> [--legacy]
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
from shapely.geometry import LineString, Polygon  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

from main import PipelineParams  # noqa: E402
from pipeline import step6_toolpath2d as s6  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402


def main():
    model = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    params = PipelineParams(resolution=0.8)
    if "--legacy" in sys.argv:
        params.toolpath2d.continuous_rings = False
    w, n_walls = params.toolpath2d.nozzle_width, params.toolpath2d.wall_count
    res = run_slicing_only(model, params)
    fractions, missing, worst, infill_cov = [], 0.0, [], []
    for layer in res.layers:
        c2d = [s6._project_to_local_2d(c, layer.local_frame) for c in layer.contour_3d]
        loops, groups = s6._classify_contour_indices(c2d)
        material = unary_union([Polygon(loops[o], [loops[h] for h in hs]).buffer(0) for o, hs in groups])
        # bande interieure de (wall_count - 0,5) largeurs : le dernier cordon de
        # paroi deborde d'une demi-largeur sur le remplissage
        band = material.difference(material.buffer(-(n_walls - 0.5) * w, join_style="mitre"))
        if band.area < 1e-6:
            continue
        beads = unary_union([LineString(p[:, :2]).buffer(w / 2, cap_style="flat")
                             for p in layer.toolpath_2d or [] if len(p) >= 2])
        covered = band.intersection(beads).area / band.area
        inner = material.buffer(-(n_walls + 0.5) * w, join_style="mitre")  # hors recouvrement des parois
        if inner.area > 4.0:
            infill_cov.append(inner.intersection(beads).area / inner.area)
        fractions.append(covered)
        worst.append((covered, layer.layer_id if hasattr(layer, 'layer_id') else len(worst), band.area))
        missing += band.area - band.intersection(beads).area
    f = np.array(fractions)
    print(f"{os.path.basename(model)}{' (historique)' if '--legacy' in sys.argv else ''} : bande de parois couverte "
          f"mediane {100 * np.median(f):.2f} %, 5e centile {100 * np.percentile(f, 5):.2f} %, "
          f"min {100 * f.min():.2f} % ; aire non couverte totale {missing:.1f} mm2 sur {len(f)} couches")
    if infill_cov:
        g = np.array(infill_cov)
        print(f"  remplissage : part de la zone interieure couverte mediane {100 * np.median(g):.1f} % "
              f"(visee {100 * params.toolpath2d.infill_density:.0f} %), 5e-95e {100 * np.percentile(g, 5):.1f}-"
              f"{100 * np.percentile(g, 95):.1f} %, couches a moins de la moitie de la visee : "
              f"{int((g < params.toolpath2d.infill_density / 2).sum())}/{len(g)}")
    print("  pires couches (couverture, couche, aire de bande) :",
          [(round(c, 3), i, round(a, 1)) for c, i, a in sorted(worst)[:6]])


if __name__ == "__main__":
    main()
