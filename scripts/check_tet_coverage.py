"""
Part de la surface du modele hors du volume tetraedrise (core/tet_coverage.py).

Usage :
  .venv/Scripts/python.exe scripts/check_tet_coverage.py <stl> [--resolution=0.8]   # maillage de step1
Affiche la fraction, l'aire et les plus grands amas (aire mm2, centre x y z).
"""
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import trimesh  # noqa: E402

from core.tet_coverage import surface_outside_volume  # noqa: E402


def show(label, rep, seconds):
    print(f"{label} : {rep.far_fraction:.3%} de la surface a plus de {rep.tolerance:.2f} mm de la frontiere du maillage "
          f"(~{rep.far_area:.1f} mm2 sur {rep.surface_area:.0f}) ; ecart max {rep.max_distance:.3f} mm ; "
          f"aire frontiere {rep.boundary_area:.0f} ; rapport de volumes {rep.volume_ratio:.4f} ; "
          f"{len(rep.clusters)} amas ; {seconds:.1f} s")
    for c in rep.clusters[:6]:
        print(f"   {c.area:7.1f} mm2  centre {np.round(c.center, 1)}")


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    from pipeline import step1_ingest
    stl = args[0]
    t0 = time.time()
    _, _, report = step1_ingest.run(stl, float(flags.get("--resolution", 0.8)), check_surface_coverage=True)
    print(f"{os.path.basename(stl)} : {report.surface_outside_volume_fraction:.3%} de la peau loin de la frontiere "
          f"(~{report.surface_outside_volume_area:.1f} mm2) ; ecart max {report.surface_max_distance_to_tets:.3f} mm ; rapport de volumes {report.tet_volume_ratio:.4f} ; {report.tetrahedra_generated} tetraedres ; {time.time() - t0:.1f} s (step1 compris)")
    for area, center in report.surface_outside_volume_clusters[:6]:
        print(f"   {area:7.1f} mm2  centre {np.round(center, 1)}")


if __name__ == "__main__":
    main()
