"""
G-code de la machine XYZAB (`MachineProfile.kinematics == "head-ab"`, docs/CINEMATIQUE_XYZAB.md).

Architecture (decision de l'utilisateur, 2026-11-27) : plateau carre qui se deplace en Y (bedslinger), tete en X et Z ;
A tourne autour de X et porte B, qui tourne autour de Y au repos ; buse vers -Z au repos. Les axes de rotation sont
ceux du profil (`rotation_axes`, du plus externe au plus interne), jamais des constantes du code.

Le firmware ne fait aucune cinematique : le G-code commande le PIVOT des axes de la tete, a `pivot_to_tip` de la
pointe le long de l'axe de la buse. Avec d(A, B) la direction pointe -> pivot (R_A R_B (0, 0, 1), sens d'impression
de l'etape 8), pivot = pointe + L d. Le firmware interpolant X, Y, Z, A, B lineairement, un changement d'angle entre
deux points fait sortir la pointe de sa droite : chaque mouvement est decoupe jusqu'a ce qu'elle reste a moins de
`GcodeParams.max_tip_deviation` (meme regle que la machine polaire).

Limites de cette premiere version : chemin Python seul (le noyau natif ne connait que la machine polaire) ; pas de
controle de collision de la tete (aucun modele 3D de la nouvelle tete) ; pas de bloc d'en-tete chiffre.
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np

from file_io import kinematics_description as _kin
from core import line_types
from file_io.gcode_axes import axes_words

REFERENCE = np.array([0.0, 0.0, 1.0])     # sens d'impression a angles nuls (etape 8)


def _rotate(v: np.ndarray, axis: np.ndarray, angle_rad: np.ndarray) -> np.ndarray:
    """Rodrigues, un angle par ligne de `v`."""
    k = np.asarray(axis, dtype=float)
    k = k / np.linalg.norm(k)
    c, s = np.cos(angle_rad)[:, None], np.sin(angle_rad)[:, None]
    return v * c + np.cross(k, v) * s + np.outer(v @ k, k) * (1.0 - c)


def head_directions(angles_deg: np.ndarray, axes) -> np.ndarray:
    """Direction pointe -> pivot pour chaque jeu d'angles (N, len(axes)) : R_1 ... R_n (0, 0, 1)."""
    ang = np.atleast_2d(np.asarray(angles_deg, dtype=float))
    d = np.tile(REFERENCE, (len(ang), 1))
    for j in range(len(axes) - 1, -1, -1):         # le plus interne d'abord
        d = _rotate(d, axes[j].axis, np.radians(ang[:, j]))
    return d


def pivot_positions(tips: np.ndarray, angles_deg: np.ndarray, axes, pivot_to_tip: float) -> np.ndarray:
    return np.asarray(tips, dtype=float) + pivot_to_tip * head_directions(angles_deg, axes)


def description(machine, pivot_to_tip: float) -> dict:
    """`apex-kinematics/1` de la machine XYZAB."""
    a_axis, b_axis = (list(map(float, ax.axis)) for ax in machine.rotation_axes[:2])
    a_rng, b_rng = (list(map(float, ax.angle_range_deg)) for ax in machine.rotation_axes[:2])
    size = getattr(machine, "platform_size", None)
    axes = [
        {"letter": "X", "type": "linear", "side": "tool", "axis": [1, 0, 0],
         "range": None if size is None else [-size[0] / 2, size[0] / 2], "positive": "la tête va vers +X"},
        {"letter": "Y", "type": "linear", "side": "bed", "axis": [0, 1, 0],
         "range": None if size is None else [-size[1] / 2, size[1] / 2],
         "positive": "le point commandé va vers +Y de la pièce (le plateau recule)"},
        {"letter": "Z", "type": "linear", "side": "tool", "axis": [0, 0, 1], "range": None, "positive": "la tête monte"},
        {"letter": "A", "type": "rotary", "side": "tool", "axis": a_axis, "origin": "pivot de la tête", "range": a_rng,
         "positive": "main droite autour de axis"},
        {"letter": "B", "type": "rotary", "side": "tool", "axis": b_axis, "origin": "pivot de la tête", "range": b_rng,
         "positive": "main droite autour de axis (repère déjà tourné par A)"},
        _kin._unused("C"),
    ]
    origin = [float(v) for v in (machine.platform_origin if machine.platform_origin is not None else (0, 0, 0))]
    return {"schema": _kin.SCHEMA, "machine": "apex-xyzab", "units": {"linear": "mm", "angular": "deg"},
            "word_order": "XYZABC", "frame": "repère de la pièce, Z vers le haut",
            "axes": axes, "chain": {"bed": ["Y"], "tool": ["X", "Z", "A", "B"]},
            "tool_offset": [0, 0, -float(pivot_to_tip)], "part_origin": origin, "coupling": None,
            "conventions": {"rotation": _kin.ROTATION_RULE, "tool_offset_from": "pivot de A et B",
                            "pivot_compensation": "le G-code commande le PIVOT : (X, Y, Z) = pointe + L·R_A·R_B·(0, 0, 1)",
                            "pivot_to_tip": float(pivot_to_tip)}}


class _Writer:
    def __init__(self, params, axes, L: float, rep, out: list, types):
        self.params, self.axes, self.L, self.rep, self.out, self.types = params, axes, L, rep, out, types
        self.tip, self.ang, self.fw, self.e = None, None, None, 0.0
        self.guard = {"x": [], "y": [], "z": []}

    def fmt(self, v: float) -> str:
        return f"{v:.{self.params.decimals}f}"

    def _fw(self, tip, ang) -> np.ndarray:
        return np.concatenate([pivot_positions(np.atleast_2d(tip), np.atleast_2d(ang), self.axes, self.L)[0],
                               np.asarray(ang, dtype=float)[:2]])

    def _deviation(self, t0, a0, t1, a1) -> float:
        """Ecart de la pointe a la droite t0 -> t1 au milieu du mouvement interpole par le firmware."""
        fw0, fw1 = self._fw(t0, a0), self._fw(t1, a1)
        mid = 0.5 * (fw0 + fw1)
        tip_mid = mid[:3] - self.L * head_directions(mid[3:][None, :], self.axes)[0]
        return float(np.linalg.norm(tip_mid - 0.5 * (t0 + t1)))

    def _split(self, t0, a0, t1, a1, depth, pieces):
        limit = self.params.max_tip_deviation
        if limit is None or depth >= 12 or self._deviation(t0, a0, t1, a1) <= limit:
            pieces.append((t1, a1))
            return
        tm, am = 0.5 * (t0 + t1), 0.5 * (a0 + a1)
        self._split(t0, a0, tm, am, depth + 1, pieces)
        self._split(tm, am, t1, a1, depth + 1, pieces)

    def move(self, tip, ang, tip_speed, e_add=0.0, comment=""):
        tip, ang = np.asarray(tip, dtype=float), np.asarray(ang, dtype=float)
        if self.tip is None:
            self.emit(tip, ang, tip_speed, e_add, comment)
            return
        pieces = []
        self._split(self.tip, self.ang, tip, ang, 0, pieces)
        self.rep.subdivided_moves += int(len(pieces) > 1)
        total = float(np.linalg.norm(tip - self.tip))
        start = self.tip
        for k, (p, a) in enumerate(pieces):
            share = float(np.linalg.norm(p - start)) / total if total > 1e-12 else 1.0 / len(pieces)
            self.emit(p, a, tip_speed, e_add * share, comment if k == len(pieces) - 1 else "")
            start = p

    def emit(self, tip, ang, tip_speed, e_add=0.0, comment=""):
        params, rep = self.params, self.rep
        fw = self._fw(tip, ang)
        rep.tip_z_min = min(rep.tip_z_min, float(tip[2]))
        rep.tip_below_plate += int(tip[2] < 0.0)
        for key, v in zip("xyz", fw[:3]):
            self.guard[key].append(float(v))
        tip_d = 0.0 if self.tip is None else float(np.linalg.norm(tip - self.tip))
        fw_d = 0.0 if self.fw is None else float(np.linalg.norm(fw[:3] - self.fw[:3]))
        feed = tip_speed * 60.0 * fw_d / tip_d if tip_d > 1e-9 else params.rotation_feed
        feed = min(max(feed, 1.0), params.max_feed)
        self.e += e_add
        words = axes_words({"X": fw[0], "Y": fw[1], "Z": fw[2], "A": fw[3], "B": fw[4]}, self.fmt)
        if e_add:
            words.append(f"E{self.fmt(self.e)}")
        self.out.append("G01 " + " ".join(words) + f" F{self.fmt(feed)}" + (f" ; {comment}" if comment else ""))
        self.tip, self.ang, self.fw = tip, ang, fw
        rep.moves += 1

    def retract(self, sign: float):
        self.e += sign * self.params.retract_length
        self.out.append(f"G01 E{self.fmt(self.e)} F{self.fmt(self.params.retract_speed * 60.0)}")

    def lifted(self, tip, ang) -> np.ndarray:
        """Degagement : la pointe recule de `travel_clearance` le long de l'axe de la buse (vers le pivot)."""
        return np.asarray(tip, dtype=float) + self.params.travel_clearance * head_directions(
            np.atleast_2d(ang), self.axes)[0]

    def travel(self, tip, ang):
        params = self.params
        if float(np.linalg.norm(tip - self.tip)) < params.retract_min_travel:
            self.move(tip, ang, params.travel_speed)
            return
        self.retract(-1)
        self.rep.retractions += 1
        self.move(self.lifted(self.tip, self.ang), self.ang, params.travel_speed, comment="degagement")
        self.move(self.lifted(tip, ang), ang, params.travel_speed)
        self.move(tip, ang, params.travel_speed)
        self.retract(+1)


def _header(params, machine, L: float) -> list:
    from file_io.apex_signature import SIGNATURE
    out = [SIGNATURE, "; Apex -- G-code machine XYZAB (file_io/gcode_xyzab.py)",
           f"; pivot->pointe {L:g} mm compense ; retraction {params.retract_length} mm a {params.retract_speed} mm/s",
           "G21", "G90", "M82", "G92 E0", f"M109 S{params.nozzle_temp:g}", "; apex_sliced = 1"]
    return out + _kin.header_lines(description(machine, L))


def export(toolpath: dict, machine, params, path: Optional[str] = None):
    """(lignes, GcodeReport) de la machine XYZAB ; ecrit `path` si donne."""
    from file_io.gcode_export import GcodeReport
    pos = np.asarray(toolpath["position"], dtype=float)
    ang = np.asarray(toolpath["axis_angles"], dtype=float)
    ext = np.asarray(toolpath["extrusion"], dtype=float)
    skin = np.asarray(toolpath.get("is_skin", np.zeros(len(pos), bool)), dtype=bool)
    types = toolpath.get("line_type")
    types = None if types is None else np.asarray(types, dtype=np.uint8)
    axes = machine.rotation_axes
    if len(axes) != 2 or (len(pos) and ang.shape[1] != 2):
        raise ValueError(f"machine XYZAB : deux axes de tete attendus, {len(axes)} dans le profil")
    L = float(machine.pivot_to_tip if getattr(machine, "pivot_to_tip", None) is not None else params.pivot_to_tip)
    if params.tip_at_bead_top and len(pos):
        half = np.where(skin, params.skin_thickness, params.layer_height) * 0.5
        pos = pos + head_directions(ang, axes) * half[:, None]
    rep = GcodeReport()
    out = _header(params, machine, L)
    w = _Writer(params, axes, L, rep, out, types)
    fil_area = math.pi * (params.filament_diameter / 2) ** 2
    last_type = None
    for i in range(len(pos)):
        if i == 0:
            rep.travels += 1
            w.move(w.lifted(pos[0], ang[0]), ang[0], params.travel_speed, comment="approche")
            w.move(pos[0], ang[0], params.travel_speed)
            continue
        if ext[i] <= 0:
            rep.travels += 1
            w.travel(pos[i], ang[i])
            continue
        if types is not None and int(types[i]) != last_type:
            last_type = int(types[i])
            out.append(f";TYPE:{line_types.LABELS.get(last_type, 'CUSTOM')}")
        seg = float(np.linalg.norm(pos[i] - pos[i - 1]))
        h = params.skin_thickness if skin[i] else params.layer_height
        e_add = seg * h * params.line_width * params.flow / fil_area
        rep.filament_mm += e_add
        w.move(pos[i], ang[i], params.print_speed, e_add=e_add)
    if len(pos):
        w.retract(-1)
        w.move(w.lifted(pos[-1], ang[-1]), ang[-1], params.travel_speed, comment="degagement final")
    out += ["M104 S0", "M84"]
    _report(rep, out, pos, ext, params, w.guard)
    if params.strict and rep.tip_below_plate:
        return out, rep
    if path:
        rep.written = True
        with open(path, "w", newline="\n") as f:
            f.write("\n".join(out) + "\n")
    return out, rep


def _report(rep, out: list, pos: np.ndarray, ext: np.ndarray, params, guard: dict) -> None:
    if len(pos) > 1:
        seg = np.linalg.norm(np.diff(pos, axis=0), axis=1)
        pr = ext[1:] > 0
        rep.print_time_nominal_s = float(seg[pr].sum() / params.print_speed + seg[~pr].sum() / params.travel_speed)
        rep.print_time_s = rep.print_time_nominal_s
    rep.lines = len(out)
    if guard["x"]:
        rep.fw_x_range = (min(guard["x"]), max(guard["x"]))
        rep.fw_y_range = (min(guard["y"]), max(guard["y"]))
        rep.fw_z_range = (min(guard["z"]), max(guard["z"]))
