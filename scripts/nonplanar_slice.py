"""Tranche un maillage « surface en hauteur » en couches non planaires : STL -> G-code.

    python scripts/nonplanar_slice.py modele.stl sortie.gcode [--layer-height 0.2] [--center 155,155]
"""
import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import trimesh  # noqa: E402

from nonplanar import SliceParams, slice_mesh  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stl")
    ap.add_argument("gcode")
    ap.add_argument("--layer-height", type=float, default=None)
    ap.add_argument("--min-layer", type=float, default=0.1)
    ap.add_argument("--max-layer", type=float, default=0.25)
    ap.add_argument("--no-intermediate", action="store_true")
    ap.add_argument("--grid-cell", type=float, default=0.25)
    ap.add_argument("--center", help="x,y : centre de l'empreinte sur le plateau")
    a = ap.parse_args()
    p = SliceParams(layer_height=a.layer_height, grid_cell=a.grid_cell, min_layer_height=a.min_layer, max_layer_height=a.max_layer, intermediate_layers=not a.no_intermediate)
    if a.center:
        p.center_xy = tuple(float(v) for v in a.center.split(","))
    t0 = time.time()
    gcode, rep = slice_mesh(trimesh.load(a.stl), p)
    Path(a.gcode).write_text(gcode)
    print(f"{rep.layer_count} couches principales + {rep.intermediate_layer_count} inter-couches, épaisseur {rep.thickness_min:.3f}-{rep.thickness_max:.3f} mm, "
          f"pente max du dessus {rep.slope_max_deg:.1f}°, filament {rep.filament_length:.0f} mm, "
          f"{time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
