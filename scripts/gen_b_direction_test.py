"""
G-code d'essai A SEC (sans extrusion, sans chauffe) du sens de rotation de B.

La POINTE reste immobile en (R, 0, Z) pendant que B balaie 0 -> B_max -> 0 ; le
pivot, lui, se deplace (le firmware commande le pivot, limite #61). Deux
verifications visuelles, a faire buse en l'air, a vitesse lente :
  1. la pointe de la buse ne bouge pas (sinon `tilt_sign` / `pivot_to_tip` sont faux) ;
  2. le SENS : B firmware POSITIF (`b_sign = +1`) doit pencher le CORPS de la tete vers +R
     (l'exterieur du plateau), donc la buse vers le centre ; B NEGATIF envoie la pointe vers
     l'exterieur (croquis de l'utilisateur, 2026-10-06). Si c'est l'inverse, passer `--b-sign=-1`.

Variante `--preview` (pour un previsualiseur qui ne dessine que la matiere EXTRUDEE) : la
pointe avance en LIGNE DROITE de R=50 a R=90 mm, a Z=40 mm, pendant que B balaie
0 -> B_max -> 0, avec de l'extrusion (meme format que les G-code d'impression :
a ne PAS lancer sur la machine, l'extrudeur est froid). Lecture correcte
(pivot compense) : une droite ; lecture « X, Z = pointe » : une courbe.

Usage : .venv/Scripts/python.exe scripts/gen_b_direction_test.py [--out=out/test_b_direction.gcode]
        [--radius=60] [--height=40] [--b-max=60] [--step=20] [--b-sign=1] [--feed=1500] [--preview]
        [--pivot-to-tip=<mm>]   (defaut : derniere valeur de config/pivot_to_tip.json, sinon 46)

Avec le bon L, la pointe reste fixe ; une erreur dL la fait deriver d'environ dL sin(B) en XY.
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from file_io.gcode_export import GcodeParams, machine_to_firmware  # noqa: E402
from file_io.pivot_length import current_pivot_to_tip  # noqa: E402


def preview(out, b_max, b_sign, height, pivot_to_tip):
    """G-code d'apercu : pointe en ligne droite, B balaye, extrusion non nulle."""
    import main as apex_main
    from file_io.gcode_export import export_gcode
    from pipeline.step8_machine_toolpath import coupled_axis_angle
    machine = apex_main._default_machine_profile()
    machine.nozzle_model = None
    n = 61
    radius = np.linspace(50.0, 90.0, n)
    b = np.concatenate([np.linspace(0.0, b_max, n // 2 + 1), np.linspace(b_max, 0.0, n - n // 2 - 1)[1:]])
    b = np.interp(np.arange(n), np.linspace(0, n - 1, len(b)), b)
    pos = np.column_stack([radius, np.zeros(n), np.full(n, height)])
    c = [coupled_axis_angle(p, machine.platform_origin, machine.rotation_axes[0].axis, machine.radial_direction, 1e-6)
         for p in pos]
    tp = {"position": pos, "axis_angles": np.column_stack([c, b]),
          "extrusion": np.r_[0.0, np.ones(n - 1)], "is_skin": np.zeros(n, dtype=bool)}
    lines, rep = export_gcode(tp, machine, GcodeParams(b_sign=b_sign, pivot_to_tip=pivot_to_tip, tip_at_bead_top=False,
                                                       max_tip_deviation=None), out)
    print(f"{out} : apercu, {rep.moves} mouvements, B firmware de {min(l for l in [b_sign * x for x in b]):g} a "
          f"{max(b_sign * x for x in b) + 0.0:g}")


def ring_preview(out, b_sign, pivot_to_tip):
    """Deux tours de plateau (Y deroule de 0 a 720 deg) : le tour 1 a B = 0 (R=60, Z=40),
    le tour 2 a B = 60 Apex (B firmware -60) a Z = 45. Lecture correcte : deux cercles de
    rayon 60 mm. Teste la convention de Y (sens, repere) et son deroulement au-dela de 360."""
    import main as apex_main
    from file_io.gcode_export import export_gcode
    from pipeline.step8_machine_toolpath import coupled_axis_angle
    machine = apex_main._default_machine_profile()
    machine.nozzle_model = None
    m = 180
    t = np.linspace(0.0, 4.0 * np.pi, 2 * m + 1)
    pos = np.column_stack([60.0 * np.cos(t), 60.0 * np.sin(t), np.where(np.arange(len(t)) <= m, 40.0, 45.0)])
    b = np.where(np.arange(len(t)) <= m, 0.0, 60.0)
    c = [coupled_axis_angle(p, machine.platform_origin, machine.rotation_axes[0].axis, machine.radial_direction, 1e-6)
         for p in pos]
    tp = {"position": pos, "axis_angles": np.column_stack([c, b]),
          "extrusion": np.r_[0.0, np.ones(len(t) - 1)], "is_skin": np.zeros(len(t), dtype=bool)}
    lines, rep = export_gcode(tp, machine, GcodeParams(b_sign=b_sign, pivot_to_tip=pivot_to_tip, tip_at_bead_top=False,
                                                       max_tip_deviation=None), out)
    ys = [float(w[1:]) for l in lines if l.startswith("G00") for w in l.split() if w.startswith("Y")]
    print(f"{out} : anneaux, {rep.moves} mouvements, Y de {min(ys):g} a {max(ys):g}")


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    if "--pivot-to-tip" in flags:
        pivot_to_tip, origin = float(flags["--pivot-to-tip"]), "option --pivot-to-tip"
    else:
        pivot_to_tip, origin = current_pivot_to_tip()
    print(f"pivot -> pointe : {pivot_to_tip:.3f} mm ({origin})")
    if "--ring" in sys.argv:
        ring_preview(flags.get("--out", os.path.join(ROOT, "out", "test_y_ring_preview.gcode")),
                     float(flags.get("--b-sign", 1.0)), pivot_to_tip)
        return
    if "--preview" in sys.argv:
        preview(flags.get("--out", os.path.join(ROOT, "out", "test_b_preview.gcode")),
                float(flags.get("--b-max", 60)), float(flags.get("--b-sign", 1.0)), float(flags.get("--height", 40)),
                pivot_to_tip)
        return
    out = flags.get("--out", os.path.join(ROOT, "out", "test_b_direction.gcode"))
    radius, height = float(flags.get("--radius", 60)), float(flags.get("--height", 40))
    b_max, step = float(flags.get("--b-max", 60)), float(flags.get("--step", 20))
    feed = float(flags.get("--feed", 1500))
    params = GcodeParams(b_sign=float(flags.get("--b-sign", 1.0)), pivot_to_tip=pivot_to_tip)
    tip = np.array([radius, 0.0, height])
    up = np.arange(0.0, b_max + 1e-9, step)
    sequence = list(up) + list(up[-2::-1])
    lines = [
        "; Apex -- essai A SEC du sens de B (sans extrusion, sans chauffe)",
        f"; pointe immobile en R={radius:g} mm, Z={height:g} mm ; pivot a {params.pivot_to_tip:g} mm de la pointe",
        f"; b_sign = {params.b_sign:g} : B firmware = {params.b_sign:g} x B d'Apex",
        "; A VERIFIER : 1) la pointe ne bouge pas ; 2) B positif penche le corps de la tete vers +R (buse vers le centre)",
        "G21", "G90",
    ]
    prev_y = None
    for b in sequence:
        x, y, z, bf = machine_to_firmware(tip, float(b), params.pivot_to_tip, prev_y, params.singularity_radius,
                                          params.tilt_sign, params.z_reference, params.b_sign)
        prev_y = y
        lines.append(f"G00 X{x:.4f} Y{y:.4f} Z{z:.4f} B{(bf if bf else 0.0):g} F{feed:g} ; B Apex {b:g}")
    lines.append("; fin")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print(f"{out} : {len(sequence)} mouvements, B firmware de {min(params.b_sign * b for b in sequence):g} a "
          f"{max(params.b_sign * b for b in sequence) + 0.0:g}")


if __name__ == "__main__":
    main()
