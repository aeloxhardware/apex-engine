"""
STL/3MF -> G-code pour la machine actuelle, par le pipeline d'Apex
(step1 -> step8, trajectoire machine, export). Voie d'essai rapide : un petit
modele.

Controles avant d'ecrire : refus d'un fichier dont la pointe passe sous le
plateau (`strict`) et collisions de la VRAIE tete (step9 rapide) ; les collisions
sont comptees et affichees, pas corrigees (c'est le role de step10).

Usage : APEX_USE_CPP=1 .venv/Scripts/python.exe -u scripts/export_model_gcode.py <modele.stl>
        [--out=out/<nom>.gcode] [--resolution=0.8] [--b-sign=1] [--pivot-to-tip=<mm>]
L (pivot -> pointe) : derniere valeur de config/pivot_to_tip.json (file_io/pivot_length.py),
sinon 46 mm. A regenerer apres tout changement de L.
"""
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from file_io.gcode_export import GcodeParams, export_gcode, head_zlow_from_table, toolpath_arrays  # noqa: E402
from file_io.pivot_length import current_pivot_to_tip  # noqa: E402
from main import PipelineParams  # noqa: E402
from pipeline import step9_collision_sim as s9  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    model = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    name = os.path.splitext(os.path.basename(model))[0]
    out = flags.get("--out", os.path.join(ROOT, "out", f"{name}.gcode"))
    params = PipelineParams(resolution=float(flags.get("--resolution", 0.8)))
    t0 = time.time()
    result = run_slicing_only(model, params)
    machine = params.machine_profile
    tp = toolpath_arrays(result.toolpath)
    # collisions de la vraie tete (step9 rapide), comptees sans etre corrigees
    sim = s9.CollisionSimParams(**{**params.collision.__dict__, "fast_kdtree": True})
    rep9 = s9.run(result.toolpath, machine, sim)
    kinds = {}
    for e in rep9.events:
        kinds[e.obstacle_type] = kinds.get(e.obstacle_type, 0) + 1
    # point le plus bas de la tete selon B (marge plateau de l'export)
    head = s9._filter_nozzle_contact_zone(s9._surface_sample_points(machine.nozzle_model, sim.toolhead_sample_spacing),
                                          sim.nozzle_contact_clearance)
    lo_b, hi_b = machine.rotation_axes[1].angle_range_deg
    grid_b = np.arange(lo_b, hi_b + 1e-9, 0.5)
    zlow = np.array([s9._positioned_toolhead_vertices(head, np.zeros(3), [0.0, float(b)], machine.rotation_axes)
                     [:, 2].min() for b in grid_b])
    if "--pivot-to-tip" in flags:
        pivot_to_tip, origin = float(flags["--pivot-to-tip"]), "option --pivot-to-tip"
    else:
        pivot_to_tip, origin = current_pivot_to_tip()
    gp = GcodeParams(b_sign=float(flags.get("--b-sign", 1.0)), pivot_to_tip=pivot_to_tip, strict=True)
    lines, rep = export_gcode(tp, machine, gp, out, head_zlow=head_zlow_from_table(grid_b, zlow))
    b_fw = np.array([float(l.split(" B")[1].split()[0]) for l in lines if l.startswith("G00") and " B" in l])
    print(f"{out}\n  pivot -> pointe {pivot_to_tip:.3f} mm ({origin})")
    print(f"  {rep.lines} lignes, {rep.moves} mouvements, {rep.travels} trajets a vide, {rep.retractions} retractions")
    print(f"  filament {rep.filament_mm / 1000:.2f} m ; temps nominal {rep.print_time_nominal_s / 3600:.2f} h, "
          f"avec limites de B {rep.print_time_s / 3600:.2f} h ; tous axes (config.g recommande) {rep.print_time_fw_s / 3600:.2f} h")
    print(f"  B firmware {b_fw.min():g} a {b_fw.max():g} ; pointe la plus basse {rep.tip_z_min:.3f} mm ; "
          f"tete sous la marge : {rep.head_below_margin} ; degagements {rep.clearance_axes} ; Y cumule {rep.y_turns:.1f} tours")
    print(f"  collisions de la vraie tete (step9 rapide) : {kinds or 'aucune'} ; duree {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
