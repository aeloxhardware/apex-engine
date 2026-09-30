"""
Export G-code d'une trajectoire machine corrigee (full_step89.py --apply) via
file_io/gcode_export.py d'Apex.

Usage (venv d'Apex) :
  APEX_USE_CPP=1 ../../.venv/Scripts/python.exe -u export_gcode.py <modele> <passe> [variante]
      [--print-speed=25] [--travel-speed=80] [--retract=0.8] [--retract-speed=35] [--clearance=1.0]
Sortie : output_models/loop_<modele><variante>_pass<k>.gcode + rapport.
Degagements controles : points de l'axe de la tete (pointe -> pivot) contre
la matiere deja deposee (arbre k-d, contact franc = rayon de cordon) ;
hauteur mini de la VRAIE tete (aetheris-4) pour refuser un degagement vers le
plateau.
"""
import os
import pickle
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

import main as apex_main  # noqa: E402
from file_io.gcode_export import GcodeParams, export_gcode  # noqa: E402
from pipeline import step9_collision_sim as s9  # noqa: E402


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    model, k, variant = args[0], int(args[1]), args[2] if len(args) > 2 else ""
    prefix = os.path.join(HERE, "output_models", "loop_" + model.replace(" ", "_") + variant)
    tp = pickle.load(open(f"{prefix}_machine_toolpath_pass{k}.pkl", "rb"))
    machine = apex_main._default_machine_profile()
    params = GcodeParams(print_speed=float(flags.get("--print-speed", 25)),
                         travel_speed=float(flags.get("--travel-speed", 80)),
                         retract_length=float(flags.get("--retract", 0.8)),
                         retract_speed=float(flags.get("--retract-speed", 35)),
                         travel_clearance=float(flags.get("--clearance", 1.0)))
    t0 = time.time()
    sim = s9.CollisionSimParams(path_sampling_step=2.0, safety_margin=params.platform_margin)
    head = s9._filter_nozzle_contact_zone(s9._surface_sample_points(machine.nozzle_model, sim.toolhead_sample_spacing),
                                          sim.nozzle_contact_clearance)
    grid_b = np.arange(0.0, 185.0 + 1e-9, 0.5)
    zlow = np.array([s9._positioned_toolhead_vertices(head, np.zeros(3), [0.0, float(b)], machine.rotation_axes)
                     [:, 2].min() for b in grid_b])

    def head_zlow(b):
        return float(zlow[np.abs(grid_b - b).argmin()])

    pos = np.asarray(tp["position"], dtype=float)
    ext = np.asarray(tp["extrusion"], dtype=float)
    # matiere deposee densifiee, avec l'indice du point qui la depose
    dense, order = [], []
    for i in range(1, len(pos)):
        if ext[i] <= 0:
            continue
        seg = pos[i] - pos[i - 1]
        m = max(1, int(np.ceil(np.linalg.norm(seg) / 0.25)))
        t = np.arange(1, m + 1)[:, None] / m
        dense.append(pos[i - 1] + t * seg)
        order.append(np.full(m, i))
    dense, order = np.vstack(dense), np.concatenate(order)
    tree = cKDTree(dense)
    radius = machine.nozzle_diameter / 2
    axis_t = np.linspace(0.0, params.pivot_to_tip, 16)
    hits = []

    def travel_check(i, point, b):
        c = float(np.degrees(np.arctan2(point[0], -point[1])))
        axis_pts = s9._positioned_toolhead_vertices(np.column_stack([0 * axis_t, 0 * axis_t, axis_t]), point,
                                                    [c, b], machine.rotation_axes)
        for lst in tree.query_ball_point(axis_pts, radius):
            if lst and (order[lst] < i).any():
                hits.append((i, point.tolist(), b))
                return True
        return False

    out = f"{prefix}_pass{k}.gcode"
    lines, rep = export_gcode(tp, machine, params, out, head_zlow=head_zlow, travel_check=travel_check)
    print(f"{out}\n  {rep.lines} lignes, {rep.moves} mouvements, {rep.travels} trajets a vide, "
          f"{rep.retractions} retractions, degagements {rep.clearance_axes}")
    print(f"  filament {rep.filament_mm / 1000:.2f} m ; temps nominal {rep.print_time_nominal_s / 3600:.2f} h, "
          f"avec limites de B {rep.print_time_s / 3600:.2f} h")
    print(f"  controle de repere |Y - (C - 90)| max {rep.c_mismatch_deg:.2e} deg ; "
          f"degagements en contact franc : {rep.travel_collisions}")
    print(f"  duree de l'export {time.time() - t0:.0f} s")
    if hits:
        pickle.dump(hits, open(f"{prefix}_gcode_travel_hits_pass{k}.pkl", "wb"))


if __name__ == "__main__":
    main()
