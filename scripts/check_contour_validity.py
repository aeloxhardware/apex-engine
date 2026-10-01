"""
Contours de couche AUTO-INTERSECTANTS dans le repere local 2D de step6
(entree de shapely.buffer) : combien de groupes (contour + trous) forment un
polygone invalide, et quelle part de l'aire. Pipeline FEM par defaut.
Usage : APEX_USE_CPP=1 .venv/Scripts/python.exe -u scripts/check_contour_validity.py <stl>
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from shapely.validation import explain_validity  # noqa: E402

from main import PipelineParams  # noqa: E402
from pipeline import step6_toolpath2d as s6  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402


def main():
    model = sys.argv[1]
    res = run_slicing_only(model, PipelineParams(resolution=0.8))
    layers_bad = 0
    groups = bad = 0
    area_all = area_bad = 0.0
    kinds = {}
    for layer in res.layers:
        c2d = [s6._project_to_local_2d(c, layer.local_frame) for c in layer.contour_3d]
        loops, idx = s6._classify_contour_indices(c2d)
        hit = False
        for outer, holes in idx:
            groups += 1
            poly = s6._to_shapely_polygon(loops[outer], [loops[h] for h in holes])
            a = abs(poly.area)
            area_all += a
            if not poly.is_valid:
                bad += 1
                area_bad += a
                hit = True
                k = explain_validity(poly).split("[")[0].strip()
                kinds[k] = kinds.get(k, 0) + 1
        layers_bad += hit
    print(f"{os.path.basename(model)} : couches touchees {layers_bad}/{len(res.layers)}, "
          f"groupes invalides {bad}/{groups}, part de l'aire {100 * area_bad / max(area_all, 1e-12):.2f} %")
    print("motifs :", kinds)


if __name__ == "__main__":
    main()
