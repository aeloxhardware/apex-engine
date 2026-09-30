"""
scripts/diagnose_junction_zones_on_death.py

Diagnostic : le fix du Livrable 2 (junction_refinement_resolution) n'a
mesure AUCUN effet sur death.stl (27/228 anomalies avant ET apres,
tetraedres identiques) - verifie si detect_junction_zones() detecte ne
serait-ce qu'une seule zone sur ce modele, avec les parametres par defaut
(angle_threshold_deg=30.0).

Usage :
    .venv/Scripts/python.exe scripts/diagnose_junction_zones_on_death.py
"""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dataclasses import replace  # noqa: E402

from pipeline.step1_ingest import run as ingest_run  # noqa: E402
from pipeline.step1b_junction_detection import detect_junction_zones  # noqa: E402
from pipeline.step2_overhang import OverhangParams, run as overhang_run  # noqa: E402
from pipeline.step3_segmentation import SegmentationParams, run as segmentation_run  # noqa: E402

MODEL_PATH = "res/test_models/death.stl"


def main() -> None:
    mesh_data, _tet_mesh, ingest_report = ingest_run(MODEL_PATH, resolution=1.5)
    overhang_map = overhang_run(mesh_data, OverhangParams())
    seg_params = replace(SegmentationParams(), platform_z=ingest_report.platform_z_used)
    region_graph = segmentation_run(mesh_data, overhang_map, seg_params)

    print(f"platform_z_used = {ingest_report.platform_z_used}")
    print(f"n_regions = {len(region_graph.regions)}")

    for angle in [10.0, 20.0, 30.0, 45.0, 60.0, 90.0]:
        zones = detect_junction_zones(mesh_data, region_graph, angle_threshold_deg=angle, margin=3.0)
        print(f"angle_threshold_deg={angle:5.1f}  n_zones={len(zones)}")
        for zone in zones[:5]:
            print(f"    centroid={zone.centroid}  radius={zone.radius:.2f}")


if __name__ == "__main__":
    main()
