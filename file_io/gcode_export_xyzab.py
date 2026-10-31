"""Export G-code d'Aetheris-4 Gen 1 : `G00 X Y Z A B C E F`, pivot compense par le slicer.

Le firmware ne fait aucune cinematique : chaque mouvement porte les positions des axes physiques
(`file_io/xyzab_kinematics.py`). Les six lettres XYZABC sont toujours ecrites (`file_io/gcode_axes.py`),
C reste a 0.

Comme l'export polaire (`gcode_export.py`) :
  - la pointe est relevee de h/2 le long de l'axe de la tete (les points sont au milieu de la couche) ;
  - le firmware interpole X Y Z A B LINEAIREMENT, donc un changement d'angle entre deux points proches
    ferait sortir la pointe de sa trajectoire : chaque mouvement est decoupe par bisection jusqu'a ce que
    la pointe, a mi-mouvement, reste a `max_tip_deviation` de la droite visee ;
  - trajet a vide : retraction, degagement le long de l'axe du corps de la tete (la tete s'ecarte de la
    matiere qu'elle vient de deposer), trajet, descente, deretraction.

Ce qui N'EST PAS fait (honnetement) : ralentissement aux limites de vitesse / acceleration de A et B
(`feed_plan` est specifique a B), controle de collision des degagements (`travel_check`), contour de la
vraie tete sous le plateau (`head_zlow`), estimation de temps par les vitesses du firmware. Le plateau est
verifie sur la pointe seulement.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from core import progress as _progress
from file_io import kinematics_description as _kin
from file_io.apex_signature import SIGNATURE, common_meta
from file_io.gcode_axes import axes_words
from file_io.gcode_export import GcodeReport
from file_io.xyzab_kinematics import HeadChain, body_axis, firmware_to_tip, tip_to_firmware


@dataclass
class GcodeParamsXYZAB:
    nozzle_temp: float = 195.0
    filament_diameter: float = 1.75
    layer_height: float = 0.2
    skin_thickness: float = 0.2
    line_width: float = 0.4
    flow: float = 1.0
    print_speed: float = 25.0
    travel_speed: float = 80.0
    rotation_feed: float = 3000.0           # F d'un mouvement sans deplacement de pointe
    max_feed: float = 12000.0
    retract_length: float = 0.8
    retract_speed: float = 35.0
    retract_min_travel: float = 2.0
    travel_clearance: float = 1.0           # degagement le long de l'axe du corps (mm)
    pivot_to_tip: float = 46.0              # A CONFIRMER (valeur de la machine polaire)
    pivot_offsets: Optional[tuple] = None   # decalages entre centres (voir HeadChain.offsets)
    a_sign: float = 1.0                     # A firmware = a_sign * A d'Apex
    b_sign: float = 1.0
    tip_at_bead_top: bool = True
    strict: bool = False
    platform_margin: float = 1.0
    decimals: int = 4
    max_tip_deviation: Optional[float] = 0.05
    # Rotation sans limite (plage >= 360 deg) : angle deroule le long du fichier, jamais de saut de 360.
    unwrap_full_turn_axes: bool = True


def params_from_gcode_params(gp) -> GcodeParamsXYZAB:
    """Reprend les champs communs d'un `GcodeParams` (requete de tranchage, interface)."""
    names = GcodeParamsXYZAB.__dataclass_fields__
    return GcodeParamsXYZAB(**{k: getattr(gp, k) for k in names if hasattr(gp, k)})


def chain_for(machine, params: GcodeParamsXYZAB) -> HeadChain:
    axes = [tuple(float(v) for v in spec.axis) for spec in machine.rotation_axes]
    letters = ["A" if abs(a[0]) > 0.5 else "B" for a in axes]
    return HeadChain(axes=axes, letters=letters, pivot_to_tip=params.pivot_to_tip, offsets=params.pivot_offsets)


def _encrypted_header(params, machine):
    from file_io.apex_header import header_block_lines
    desc = _kin.apex_gantry_ab(params, machine)
    meta = {**common_meta("Apex engine"), "machine": desc["machine"], "kinematics_description": desc,
            "kinematics": {"pivot_to_tip": float(params.pivot_to_tip), "a_sign": float(params.a_sign),
                           "b_sign": float(params.b_sign), "z_reference": "tip"},
            "layer_height": float(params.layer_height), "line_width": float(params.line_width),
            "print_speed": float(params.print_speed), "travel_speed": float(params.travel_speed),
            "nozzle_temp": float(params.nozzle_temp)}
    return header_block_lines(meta)


def _signed(angles, params):
    out = np.array(angles, dtype=float)
    out[:, 0] *= params.a_sign
    out[:, 1] *= params.b_sign
    return out


def _unwrapped(ang, machine, params):
    ang = np.array(ang, dtype=float)
    if params.unwrap_full_turn_axes:
        for k, spec in enumerate(machine.rotation_axes):
            lo, hi = spec.angle_range_deg
            if hi - lo >= 360.0 and len(ang) > 1:
                ang[:, k] = np.degrees(np.unwrap(np.radians(ang[:, k])))
    return ang


def export_gcode_xyzab(toolpath: dict, machine, params: GcodeParamsXYZAB, path: Optional[str] = None):
    """Ecrit le G-code (si `path`) et rend (lignes, GcodeReport). `toolpath` : `toolpath_arrays`."""
    pos = np.asarray(toolpath["position"], dtype=float)
    ang = np.asarray(toolpath["axis_angles"], dtype=float)
    ext = np.asarray(toolpath["extrusion"], dtype=float)
    skin = np.asarray(toolpath.get("is_skin", np.zeros(len(pos), bool)), dtype=bool)
    n = len(pos)
    if len(machine.rotation_axes) != 2:
        raise ValueError("l'export XYZAB attend exactement deux axes de rotation (A et B)")
    chain = chain_for(machine, params)
    ang = _unwrapped(ang, machine, params)
    rep = GcodeReport()
    dec = params.decimals
    fil_area = np.pi * (params.filament_diameter / 2) ** 2

    if params.tip_at_bead_top and n:
        half = np.where(skin, params.skin_thickness, params.layer_height) * 0.5
        pos = pos + body_axis(chain, ang) * half[:, None]

    def fmt(v):
        s = f"{v:.{dec}f}".rstrip("0").rstrip(".")
        return "0" if s in ("-0", "") else s

    out = [SIGNATURE, "; Apex -- G-code machine X Y Z A B (file_io/gcode_export_xyzab.py)",
           f"; pivot->pointe {params.pivot_to_tip} mm compense ; retraction {params.retract_length} mm a "
           f"{params.retract_speed} mm/s",
           "G21", "G90", "M82", "G92 E0", f"M109 S{params.nozzle_temp:g}", "; apex_sliced = 1"]
    out += _kin.header_lines(_kin.apex_gantry_ab(params, machine))
    out.append(f"; apex_kinematics pivot_to_tip={params.pivot_to_tip:g} a_sign={params.a_sign:g} "
               f"b_sign={params.b_sign:g} z_reference=tip order={''.join(chain.letters)}")
    out += _encrypted_header(params, machine)

    state = {"e": 0.0, "fw": None, "ang": None, "tip": None}
    fw_all = {"min": np.full(3, np.inf), "max": np.full(3, -np.inf)}

    def fw_of(point, a):
        return tip_to_firmware(chain, np.asarray(point, float)[None], np.asarray(a, float)[None])[0]

    def deviation(p0, a0, p1, a1):
        """Ecart de la pointe a mi-mouvement quand le firmware interpole X Y Z A B lineairement."""
        f0, f1 = fw_of(p0, a0), fw_of(p1, a1)
        mid = firmware_to_tip(chain, ((f0 + f1) / 2)[None], ((np.asarray(a0) + np.asarray(a1)) / 2)[None])[0]
        return float(np.linalg.norm(mid - (np.asarray(p0) + np.asarray(p1)) / 2.0))

    def emit(point, a, tip_speed, e_add=0.0, comment="", fw=None):
        point = np.asarray(point, float)
        a = np.asarray(a, float)
        fw = fw_of(point, a) if fw is None else fw
        prev = state["fw"]
        rep.tip_z_min = min(rep.tip_z_min, float(point[2]))
        rep.tip_below_plate += int(point[2] < 0.0)
        fw_all["min"] = np.minimum(fw_all["min"], fw)
        fw_all["max"] = np.maximum(fw_all["max"], fw)
        tip_d = 0.0 if prev is None else float(np.linalg.norm(point - state["tip"]))
        fw_d = 0.0 if prev is None else float(np.linalg.norm(fw - prev))
        feed = tip_speed * 60.0 * fw_d / tip_d if tip_d > 1e-9 else params.rotation_feed
        feed = min(max(feed, 1.0), params.max_feed)
        state["e"] += e_add
        sa = a * np.array([params.a_sign, params.b_sign])
        words = axes_words({"X": fw[0], "Y": fw[1], "Z": fw[2], "A": sa[0], "B": sa[1]}, fmt)
        if e_add:
            words.append(f"E{fmt(state['e'])}")
        out.append("G00 " + " ".join(words) + f" F{fmt(feed)}" + (f" ; {comment}" if comment else ""))
        state["fw"], state["tip"], state["ang"] = fw, point, a
        rep.moves += 1

    def move(point, a, tip_speed, e_add=0.0, comment=""):
        tol = params.max_tip_deviation
        if tol is None or state["fw"] is None:
            emit(point, a, tip_speed, e_add, comment)
            return
        p_start, a_start = state["tip"], state["ang"]
        p1, a1 = np.asarray(point, float), np.asarray(a, float)
        pieces = []
        # Subdivision UNIFORME (pointe et angles interpoles lineairement) : le plus petit k = 1, 2, 4, ...
        # dont tous les morceaux gardent la pointe a `tol` de leur droite a mi-chemin. Un seul calcul
        # vectorise par k (la bisection point par point coutait ~1 ms par appel).
        for k in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024):
            t = np.linspace(0.0, 1.0, k + 1)[:, None]
            nodes_p = p_start + t * (p1 - p_start)
            nodes_a = a_start + t * (a1 - a_start)
            fw = tip_to_firmware(chain, nodes_p, nodes_a)
            mid = firmware_to_tip(chain, (fw[1:] + fw[:-1]) / 2.0, (nodes_a[1:] + nodes_a[:-1]) / 2.0)
            dev = np.linalg.norm(mid - (nodes_p[1:] + nodes_p[:-1]) / 2.0, axis=1)
            if k == 1:
                rep.max_tip_deviation_before = max(rep.max_tip_deviation_before, float(dev[0]))
            if dev.max() <= tol or k == 1024:
                pieces = [(nodes_p[m], nodes_a[m], fw[m]) for m in range(1, k + 1)]
                break
        if len(pieces) > 1:
            rep.subdivided_moves += 1
        total = float(np.linalg.norm(np.asarray(point, float) - p_start))
        start = p_start
        for k, (p, aa, fwk) in enumerate(pieces):
            share = float(np.linalg.norm(p - start)) / total if total > 1e-12 else 1.0 / len(pieces)
            emit(p, aa, tip_speed, e_add * share, comment if k == len(pieces) - 1 else "", fw=fwk)
            start = p

    def retract(sign):
        state["e"] += sign * params.retract_length
        out.append(f"G01 E{fmt(state['e'])} F{fmt(params.retract_speed * 60.0)}")

    def clearance(point, a):
        """Point degage le long de l'axe du corps ; remonte seulement (jamais sous le plateau)."""
        d = body_axis(chain, np.asarray(a, float)[None])[0]
        up = np.asarray(point, float) + d * params.travel_clearance
        if up[2] < params.platform_margin and up[2] < point[2]:
            up = np.asarray(point, float) + np.array([0.0, 0.0, params.travel_clearance])
        return up

    # Pre-calcul VECTORISE : point commande de chaque point et ecart de la pointe a mi-mouvement de
    # chaque segment i-1 -> i (le chemin scalaire, ~1 ms par appel, ne tenait pas sur 670 000 points).
    if n:
        fw_all_pts = tip_to_firmware(chain, pos, ang)
    if n > 1:
        mid_tip = firmware_to_tip(chain, (fw_all_pts[1:] + fw_all_pts[:-1]) / 2.0, (ang[1:] + ang[:-1]) / 2.0)
        dev_all = np.concatenate([[0.0], np.linalg.norm(mid_tip - (pos[1:] + pos[:-1]) / 2.0, axis=1)])
    else:
        dev_all = np.zeros(n)

    last_kind = None
    for i in range(n):
        _progress.report(i, n, span=(0.4, 1.0))
        kind = "lissage" if skin[i] else "corps"
        if kind != last_kind:
            out.append(f"; --- {kind} ---")
            last_kind = kind
        a = ang[i]
        if i == 0 or ext[i] <= 0:
            rep.travels += 1
            if i == 0:
                move(clearance(pos[0], a), a, params.travel_speed, comment="approche")
                move(pos[0], a, params.travel_speed)
                continue
            if float(np.linalg.norm(pos[i] - pos[i - 1])) >= params.retract_min_travel:
                retract(-1)
                rep.retractions += 1
                move(clearance(pos[i - 1], ang[i - 1]), ang[i - 1], params.travel_speed, comment="degagement")
                move(clearance(pos[i], a), a, params.travel_speed)
                move(pos[i], a, params.travel_speed)
                retract(+1)
            else:
                move(pos[i], a, params.travel_speed)
            continue
        seg = float(np.linalg.norm(pos[i] - pos[i - 1]))
        h = params.skin_thickness if skin[i] else params.layer_height
        e_add = seg * h * params.line_width * params.flow / fil_area
        rep.filament_mm += e_add
        if params.max_tip_deviation is None or dev_all[i] <= params.max_tip_deviation:
            rep.max_tip_deviation_before = max(rep.max_tip_deviation_before, float(dev_all[i]))
            emit(pos[i], a, params.print_speed, e_add, fw=fw_all_pts[i])
        else:
            move(pos[i], a, params.print_speed, e_add=e_add)
    if n:
        retract(-1)
        move(clearance(pos[-1], ang[-1]), ang[-1], params.travel_speed, comment="degagement final")
    out += ["M104 S0", "M84"]
    if n > 1:
        seg = np.linalg.norm(np.diff(pos, axis=0), axis=1)
        pr = (ext > 0)[1:]
        rep.print_time_nominal_s = float(seg[pr].sum() / params.print_speed + seg[~pr].sum() / params.travel_speed)
        rep.print_time_s = rep.print_time_nominal_s       # pas de ralentissement aux limites de A et B (non fait)
    rep.lines = len(out)
    if np.isfinite(fw_all["min"]).all():
        rep.fw_x_range = (float(fw_all["min"][0]), float(fw_all["max"][0]))
        rep.fw_y_range = (float(fw_all["min"][1]), float(fw_all["max"][1]))
        rep.fw_z_range = (float(fw_all["min"][2]), float(fw_all["max"][2]))
    if params.strict and rep.tip_below_plate:
        return out, rep
    if path:
        rep.written = True
        with open(path, "w", newline="\n") as f:
            f.write("\n".join(out) + "\n")
    return out, rep
