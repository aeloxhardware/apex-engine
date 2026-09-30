"""
Export G-code pour la machine ACTUELLE (R C Z B, limite #59 de CLAUDE.md).

Ecrit de zero le 2026-09-30, sans reprendre les generateurs de formes fixes
`scripts/gen_*_gcode.py` (decision de l'utilisateur). Seuls les FAITS machine
releves dans ces scripts, verifies avec l'utilisateur, sont repris :

- commande `G00 Y<C> X<R> Z<Z> B<B> E<e> F<f>` ; en-tete `G21 G90 M82 G92 E0
  M109 S<t>` ;
- firmware Y = angle polaire du point, `atan2(y, x)` en degres. Avec la
  convention d'Apex (`coupled_axis_angle`, +R = -Y a angles nuls), c'est
  exactement `C - 90` : l'export le recalcule depuis la position et controle
  l'accord avec `axis_angles[0]` (`GcodeReport.c_mismatch_deg`) ;
- firmware B = B d'Apex (B > 0 penche le CORPS de la tete vers +R, donc la
  buse vers le centre, meme convention que step9._positioned_toolhead_vertices) ;
- le firmware commande le PIVOT de B, a `pivot_to_tip` (45 mm) de la pointe,
  sans compensation : l'export la fait. Reference : a B = 0 les coordonnees
  firmware sont celles de la pointe (calibration buse verticale), donc
  X = R + L sin B et Z = z + L (cos B - 1) ;
- le firmware traite Y (degres) et B comme des axes lineaires : F est
  recalcule par segment pour garder la vitesse de la POINTE constante.

Les limites de vitesse et d'acceleration de B (MachineProfile) ralentissent
tout le mouvement, le firmware coordonnant les axes (limite #59) :
`feed_plan` donne ce facteur par segment.

Trajets a vide (decision de l'utilisateur, 2026-09-30) : retraction, puis
degagement DYNAMIQUE sur un seul axe a la fois, choisi selon l'orientation de
la tete : Z si la tete est plutot verticale (+Z tete vers le bas, -Z tete
retournee), sinon R dans le sens du corps de la tete. Un degagement qui ferait
passer la tete sous la marge du plateau bascule sur l'autre axe.

Entree : le dictionnaire des trajectoires machine (`position` (N,3) repere
machine, `axis_angles` (N,2) degres, `extrusion` (N,) -- 0 au premier point
de chaque polyligne, donc le segment qui y mene est un trajet a vide --,
`is_skin` (N,) optionnel), tel qu'ecrit par la recherche
(`full_step89.py --apply`) ; la conversion depuis `ToolpathMachine` est
`toolpath_arrays`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np


@dataclass
class GcodeParams:
    nozzle_temp: float = 195.0
    filament_diameter: float = 1.75
    layer_height: float = 0.2            # epaisseur du cordon du corps (mm)
    skin_thickness: float = 0.2          # epaisseur des couches de lissage (mm)
    line_width: float = 0.4
    flow: float = 1.0
    print_speed: float = 25.0            # mm/s a la pointe
    travel_speed: float = 80.0           # mm/s a la pointe
    rotation_feed: float = 3000.0        # F (unites firmware / min) d'un mouvement sans deplacement de pointe
    max_feed: float = 12000.0            # plafond de F
    retract_length: float = 0.8          # mm de filament (utilisateur, 2026-09-30)
    retract_speed: float = 35.0          # mm/s (utilisateur, 2026-09-30)
    retract_min_travel: float = 2.0      # pas de retraction pour un trajet plus court (mm)
    travel_clearance: float = 1.0        # degagement le long d'un seul axe (mm)
    pivot_to_tip: float = 45.0           # pivot de B -> pointe (mm)
    platform_margin: float = 1.0         # marge tete / plateau (mm)
    singularity_radius: float = 1e-6     # R sous lequel Y est herite du point precedent
    decimals: int = 4


@dataclass
class GcodeReport:
    lines: int = 0
    moves: int = 0
    travels: int = 0
    retractions: int = 0
    clearance_axes: dict = field(default_factory=lambda: {"Z": 0, "R": 0, "none": 0})
    filament_mm: float = 0.0
    print_time_nominal_s: float = 0.0
    print_time_s: float = 0.0            # avec le ralentissement aux limites de B
    c_mismatch_deg: float = 0.0          # ecart max |Y - (C - 90)|, controle de repere
    travel_collisions: int = 0           # degagements signales par `travel_check`


def toolpath_arrays(toolpath) -> dict:
    """ToolpathMachine (Apex) -> dictionnaire d'entree de `export_gcode`."""
    pts = toolpath.points
    return {"position": np.array([p.position for p in pts], dtype=float),
            "axis_angles": np.array([p.axis_angles for p in pts], dtype=float),
            "extrusion": np.array([p.extrusion for p in pts], dtype=float),
            "is_skin": np.array([bool(p.color_id) for p in pts])}


def head_body_direction_rz(b_deg):
    """Direction du corps de la tete (pointe -> pivot) dans le plan (R, Z)."""
    b = np.radians(b_deg)
    return np.sin(b), np.cos(b)


def machine_to_firmware(position, b_deg, pivot_to_tip, previous_y=None, singularity_radius=1e-6):
    """Pointe (x, y, z) + B -> (X, Y, Z, B) firmware. Y deroule par rapport a
    `previous_y` (jamais de saut de 360 deg) ; herite sous `singularity_radius`."""
    x, y, z = (float(v) for v in position)
    r = float(np.hypot(x, y))
    if r < singularity_radius and previous_y is not None:
        fy = previous_y
    else:
        fy = float(np.degrees(np.arctan2(y, x)))
        if previous_y is not None:
            fy = previous_y + ((fy - previous_y + 180.0) % 360.0 - 180.0)
    sr, cz = head_body_direction_rz(b_deg)
    return (r + pivot_to_tip * sr, fy, z + pivot_to_tip * (cz - 1.0), float(b_deg))


def firmware_to_machine(fw_x, fw_y, fw_z, fw_b, pivot_to_tip):
    """Inverse de `machine_to_firmware` (aperçu et controle aller-retour)."""
    sr, cz = head_body_direction_rz(fw_b)
    r = fw_x - pivot_to_tip * sr
    z = fw_z - pivot_to_tip * (cz - 1.0)
    a = np.radians(fw_y)
    return np.array([r * np.cos(a), r * np.sin(a), z])


def feed_plan(points, angles, max_velocity_deg, max_accel_deg, coupled_axes, speed):
    """Temps d'une suite de points quand le firmware COORDONNE les axes : les
    angles restent exacts, le mouvement entier ralentit la ou un axe
    d'orientation (non couple a la position) atteint sa vitesse ou son
    acceleration maximale. Mouvement 1D de l'angle par segment, vitesse
    plafonnee par passes avant/arriere (v^2 <= v0^2 + 2 a |dtheta|, remise a
    zero a chaque changement de sens). Rend (temps nominal, temps reel,
    facteur de vitesse par segment, <= 1)."""
    seg = np.linalg.norm(np.diff(points, axis=0), axis=1)
    t_nom = seg / max(speed, 1e-9)
    t_real = t_nom.copy()
    if len(seg) and max_velocity_deg > 0 and max_accel_deg > 0:
        for ax in range(angles.shape[1]):
            if ax in coupled_axes:
                continue
            d = np.degrees(np.diff(np.unwrap(np.radians(angles[:, ax]))))
            mag = np.abs(d)
            cap = np.full(len(d), float(max_velocity_deg))
            prev = 0.0
            for i in range(len(d)):
                if i and np.sign(d[i]) != np.sign(d[i - 1]):
                    prev = 0.0
                cap[i] = min(cap[i], np.sqrt(prev ** 2 + 2 * max_accel_deg * mag[i]))
                prev = cap[i]
            nxt = 0.0
            for i in range(len(d) - 1, -1, -1):
                if i < len(d) - 1 and np.sign(d[i]) != np.sign(d[i + 1]):
                    nxt = 0.0
                cap[i] = min(cap[i], np.sqrt(nxt ** 2 + 2 * max_accel_deg * mag[i]))
                nxt = cap[i]
            t_real = np.maximum(t_real, mag / np.maximum(cap, 1e-9))
    factor = np.where(t_real > 0, t_nom / np.maximum(t_real, 1e-12), 1.0)
    return float(t_nom.sum()), float(t_real.sum()), factor


def clearance_move(position, b_deg, params, head_zlow: Optional[Callable[[float], float]] = None):
    """Degagement sur UN seul axe (Z ou R) selon l'orientation de la tete.
    Rend (axe, point degage) ; axe 'none' si aucun des deux n'est sur."""
    sr, cz = head_body_direction_rz(b_deg)
    x, y, z = (float(v) for v in position)
    r = float(np.hypot(x, y))
    radial = np.array([x, y, 0.0]) / r if r > 1e-9 else np.array([1.0, 0.0, 0.0])
    order = ["Z", "R"] if abs(cz) >= abs(sr) else ["R", "Z"]
    for axis in order:
        if axis == "Z":
            step = np.array([0.0, 0.0, params.travel_clearance * (1.0 if cz >= 0 else -1.0)])
        else:
            if r < 1e-9 and sr < 0:
                continue
            step = radial * params.travel_clearance * (1.0 if sr >= 0 else -1.0)
        cand = np.asarray(position, dtype=float) + step
        lowest = cand[2] + (head_zlow(b_deg) if head_zlow is not None else 0.0)
        if lowest >= params.platform_margin or cand[2] >= position[2]:
            return axis, cand
    return "none", np.asarray(position, dtype=float)


def export_gcode(toolpath: dict, machine, params: GcodeParams, path: Optional[str] = None,
                 head_zlow: Optional[Callable[[float], float]] = None,
                 travel_check: Optional[Callable[[int, np.ndarray, float], bool]] = None):
    """Ecrit le G-code (si `path`) et rend (lignes, GcodeReport).

    `head_zlow(B)` : hauteur du point le plus bas de la vraie tete sous la
    pointe pour cet angle (negative), pour refuser un degagement vers le
    plateau. `travel_check(i, point, B) -> bool` : controle de collision d'un
    point de degagement (vrai = contact), compte dans le rapport."""
    pos = np.asarray(toolpath["position"], dtype=float)
    ang = np.asarray(toolpath["axis_angles"], dtype=float)
    ext = np.asarray(toolpath["extrusion"], dtype=float)
    skin = np.asarray(toolpath.get("is_skin", np.zeros(len(pos), bool)), dtype=bool)
    n = len(pos)
    rep = GcodeReport()
    fil_area = np.pi * (params.filament_diameter / 2) ** 2
    L, dec = params.pivot_to_tip, params.decimals
    coupled = list(getattr(machine, "position_coupled_axes", None) or [])
    vmax = float(np.degrees(getattr(machine, "max_angular_velocity", 0.0) or 0.0))
    amax = float(np.degrees(getattr(machine, "max_angular_acceleration", 0.0) or 0.0))
    _, _, slow = feed_plan(pos, ang, vmax, amax, coupled, params.print_speed) if n > 1 else (0, 0, np.ones(0))

    out = ["; Apex -- G-code machine R C Z B (file_io/gcode_export.py)",
           f"; pivot->pointe {L} mm compense ; retraction {params.retract_length} mm a {params.retract_speed} mm/s",
           "G21", "G90", "M82", "G92 E0", f"M109 S{params.nozzle_temp:g}"]
    state = {"fw": None, "e": 0.0, "retracted": False}

    def fmt(v):
        s = f"{v:.{dec}f}".rstrip("0").rstrip(".")
        return "0" if s in ("-0", "") else s

    def move(point, b, tip_speed, e_add=0.0, comment=""):
        prev = state["fw"]
        fw = machine_to_firmware(point, b, L, None if prev is None else prev[1], params.singularity_radius)
        tip_d = 0.0 if prev is None else float(np.linalg.norm(point - state["tip"]))
        fw_d = 0.0 if prev is None else float(np.linalg.norm(np.subtract(fw, prev)))
        if tip_d > 1e-9:
            feed = tip_speed * 60.0 * fw_d / tip_d
        else:
            feed = params.rotation_feed
        feed = min(max(feed, 1.0), params.max_feed)
        state["e"] += e_add
        words = [f"X{fmt(fw[0])}", f"Y{fmt(fw[1])}", f"Z{fmt(fw[2])}", f"B{fmt(fw[3])}"]
        if e_add:
            words.append(f"E{fmt(state['e'])}")
        out.append("G00 " + " ".join(words) + f" F{fmt(feed)}" + (f" ; {comment}" if comment else ""))
        state["fw"], state["tip"] = fw, np.asarray(point, dtype=float)
        rep.moves += 1

    def retract(sign):
        state["e"] += sign * params.retract_length
        out.append(f"G01 E{fmt(state['e'])} F{fmt(params.retract_speed * 60.0)}")
        state["retracted"] = sign < 0

    last_kind = None
    for i in range(n):
        kind = "lissage" if skin[i] else "corps"
        b = float(ang[i, 1])
        if i > 0:
            y_fw = machine_to_firmware(pos[i], b, L, None, params.singularity_radius)[1]
            if np.hypot(*pos[i, :2]) > 1e-3:
                rep.c_mismatch_deg = max(rep.c_mismatch_deg,
                                         abs((y_fw - (ang[i, 0] - 90.0) + 180.0) % 360.0 - 180.0))
        if kind != last_kind:
            out.append(f"; --- {kind} ---")
            last_kind = kind
        if i == 0 or ext[i] <= 0:
            # trajet a vide vers pos[i]
            rep.travels += 1
            if i == 0:
                axis, up = clearance_move(pos[0], b, params, head_zlow)
                move(up, b, params.travel_speed, comment="approche")
                rep.clearance_axes[axis] += 1
                move(pos[0], b, params.travel_speed)
                continue
            length = float(np.linalg.norm(pos[i] - pos[i - 1]))
            if length >= params.retract_min_travel:
                retract(-1)
                rep.retractions += 1
                b0 = float(ang[i - 1, 1])
                axis0, up0 = clearance_move(pos[i - 1], b0, params, head_zlow)
                axis1, up1 = clearance_move(pos[i], b, params, head_zlow)
                rep.clearance_axes[axis0] += 1
                if travel_check is not None:
                    rep.travel_collisions += int(travel_check(i, up0, b0)) + int(travel_check(i, up1, b))
                move(up0, b0, params.travel_speed, comment=f"degagement {axis0}")
                move(up1, b, params.travel_speed)
                move(pos[i], b, params.travel_speed)
                retract(+1)
            else:
                move(pos[i], b, params.travel_speed)
            continue
        seg = float(np.linalg.norm(pos[i] - pos[i - 1]))
        h = params.skin_thickness if skin[i] else params.layer_height
        e_add = seg * h * params.line_width * params.flow / fil_area
        rep.filament_mm += e_add
        speed = params.print_speed * (float(slow[i - 1]) if len(slow) >= i else 1.0)
        move(pos[i], b, speed, e_add=e_add)
    if n:
        retract(-1)
        axis, up = clearance_move(pos[-1], float(ang[-1, 1]), params, head_zlow)
        move(up, float(ang[-1, 1]), params.travel_speed, comment=f"degagement final {axis}")
    out += ["M104 S0", "M84"]
    printing = ext > 0
    if n > 1:
        seg = np.linalg.norm(np.diff(pos, axis=0), axis=1)
        pr = printing[1:]
        rep.print_time_nominal_s = float(seg[pr].sum() / params.print_speed + seg[~pr].sum() / params.travel_speed)
        rep.print_time_s = float((seg[pr] / (params.print_speed * np.maximum(slow[pr], 1e-9))).sum()
                                 + seg[~pr].sum() / params.travel_speed)
    rep.lines = len(out)
    if path:
        with open(path, "w", newline="\n") as f:
            f.write("\n".join(out) + "\n")
    return out, rep


def parse_gcode(lines):
    """Relit les mouvements G00 d'un fichier produit par `export_gcode` :
    liste de (X, Y, Z, B, extrude?) firmware, dans l'ordre."""
    moves, e_prev = [], 0.0
    cur = {}
    for line in lines:
        code = line.split(";", 1)[0].strip()
        if not code:
            continue
        words = code.split()
        if words[0] not in ("G00", "G0", "G01", "G1"):
            continue
        vals = {w[0]: float(w[1:]) for w in words[1:]}
        e = vals.get("E")
        extrude = e is not None and e > e_prev + 1e-12 and any(k in vals for k in "XYZB")
        if e is not None:
            e_prev = e
        if any(k in vals for k in "XYZB"):
            cur.update({k: vals[k] for k in "XYZB" if k in vals})
            moves.append((cur["X"], cur["Y"], cur["Z"], cur["B"], extrude))
    return moves
