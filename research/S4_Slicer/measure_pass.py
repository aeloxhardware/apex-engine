"""
Mesures d'une passe DEJA calculee par bird_loop.py (sans rien recalculer) :
  - base : ecart de hauteur DEFORMEE des noeuds poses sur le plateau dans la
    piece d'origine (0 = base restee plate), et nombre de couches qui coupent
    donc la face du dessous ;
  - replis : plus grande zone de croisement de la passe (rapport de la boucle) ;
  - surplomb reel : piece d'ORIGINE, tete inclinee (Kabsch par tetraedre),
    memes fonctions et meme seuil que measure_true_overhang.py.

Usage (venv de la recherche) :
  .venv/Scripts/python.exe measure_pass.py <prefixe> [passe]
  <prefixe> : output_models/loop_<modele>[_pin]   (passe : derniere acceptee)
"""
import json
import sys

import numpy as np
import pyvista as pv

import bird_core

prefix = sys.argv[1]
_argv, sys.argv = sys.argv, sys.argv[:1]   # measure_true_overhang lit sys.argv a l'import
import measure_true_overhang as mto        # noqa: E402
sys.argv = _argv

LAYER_HEIGHT = 0.2
OVERHANG_THRESHOLD_DEG = 45.0              # celui de measure_true_overhang (6,78 % sur le Benchy)

report = json.load(open(f"{prefix}_report.json", encoding="utf-8"))
k = int(sys.argv[2]) if len(sys.argv) > 2 else report["accepted_passes"]
cells = np.load(f"{prefix}_cells.npy")
p0 = np.load(f"{prefix}_pass0.npy")
pk = np.load(f"{prefix}_pass{k}.npy") if k > 0 else p0

# --- base -----------------------------------------------------------------
# noeuds de la FACE posee sur le plateau (meme tolerance que
# measure_true_overhang.platform_contact_mask et bird_core.planar_base_mask)
on_plate = p0[:, 2] < p0[:, 2].min() + 1e-3
spread = float(np.ptp(pk[on_plate, 2]))
layers_on_bottom = int(np.floor(spread / LAYER_HEIGHT)) + 1

# --- replis ---------------------------------------------------------------
zones = [a["largest_zone_mm2"] for a in report["attempts"] if a.get("pass") == k]
fold = zones[0] if zones else float("nan")

# --- surplomb reel ----------------------------------------------------------
grid = pv.UnstructuredGrid(np.hstack([np.full((len(cells), 1), 4), cells]).ravel(),
                           np.full(len(cells), pv.CellType.TETRA), p0.copy())
original, _state = bird_core.calculate_tet_attributes(grid)
deformed = pv.UnstructuredGrid(original.cells, np.full(len(cells), pv.CellType.TETRA), pk.copy())
deformed.field_data["cells"] = cells
deformed.field_data["cell_vertices"] = pk
deformed.cell_data["cell_center"] = pk[cells].mean(axis=1)
rotations = mto.measure_cell_rotations(original, deformed)
face_rot, count = mto.face_rotations_from_cells(original, rotations)
surface = original.extract_surface()
normals = surface.face_normals
areas = surface.compute_cell_sizes(length=False, area=True, volume=False)["Area"]
centers = surface.cell_centers().points
exclude = mto.platform_contact_mask(centers)
compensated = mto.rotate_about_tangential(normals, centers, np.where(count > 0, face_rot, 0.0))
compensated /= np.maximum(np.linalg.norm(compensated, axis=1), 1e-12)[:, None]
threshold = 90.0 + OVERHANG_THRESHOLD_DEG


def overhang_percent(n):
    angle = np.degrees(np.arccos(np.clip(n[:, 2], -1, 1)))
    return 100.0 * areas[(angle > threshold) & ~exclude].sum() / areas.sum()


result = {
    "prefix": prefix, "pass": k, "planar_base": report.get("planar_base", False),
    "base_spread_mm": round(spread, 3), "layers_on_bottom_face": layers_on_bottom,
    "largest_fold_zone_mm2": fold,
    "overhang_original_pct": round(overhang_percent(normals), 2),
    "overhang_real_tilted_pct": round(overhang_percent(compensated), 2),
    "stop_reason": report["stop_reason"],
}
print(json.dumps(result, ensure_ascii=False))
