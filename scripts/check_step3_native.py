"""Pas 21b : step3 natif (chemin par defaut) contre Python : regions,
directions, profondeurs et dependances (ordre compris) identiques."""
import os
import sys
import time
from dataclasses import replace

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from core import _cpp  # noqa: E402
from main import PipelineParams  # noqa: E402
from pipeline import step1_ingest, step2_overhang, step3_segmentation as s3  # noqa: E402


def same(a, b) -> bool:
    if a.regions.keys() != b.regions.keys() or a.dependencies != b.dependencies:
        return False
    if a.region_print_depth != b.region_print_depth:
        return False
    return all(np.array_equal(a.regions[k].face_ids, b.regions[k].face_ids)
               and np.array_equal(a.regions[k].target_direction, b.regions[k].target_direction)
               for k in a.regions)


def main():
    p = PipelineParams(resolution=0.8)
    md, tm, rep = step1_ingest.run(sys.argv[1], p.resolution)
    om = step2_overhang.run(md, p.overhang)
    sp = replace(p.segmentation, platform_z=rep.platform_z_used)
    for merge in (None, 30.0):
        spm = replace(sp, max_region_merge_angle_deg=merge)
        _cpp.set_enabled(False); t0 = time.time(); a = s3.run(md, om, spm, tm); ta = time.time() - t0
        _cpp.set_enabled(True); t0 = time.time(); b = s3.run(md, om, spm, tm); tb = time.time() - t0
        print(f"{os.path.basename(sys.argv[1])} fusion={merge} : {len(a.regions)} regions, "
              f"{len(a.dependencies)} dependances, {'IDENTIQUE' if same(a, b) else 'DIFFERENT'} "
              f"(python {ta:.2f} s, natif {tb:.2f} s)")


if __name__ == "__main__":
    main()
