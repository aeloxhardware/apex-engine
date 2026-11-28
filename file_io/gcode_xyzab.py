"""
G-code de la machine XYZAB (`MachineProfile.kinematics == "head-ab"`, docs/CINEMATIQUE_XYZAB.md, limite #97).

Architecture (decision de l'utilisateur, 2026-11-27) : plateau carre qui se deplace en Y (bedslinger), tete en X et Z ;
A tourne autour de X et porte B, qui tourne autour de Y au repos ; buse vers -Z au repos. Les axes de rotation sont
ceux du profil (`rotation_axes`, du plus externe au plus interne), jamais des constantes du code.

Le firmware ne fait aucune cinematique : le G-code commande le PIVOT des axes de la tete, a `pivot_to_tip` de la
pointe le long de l'axe de la buse. Avec d(A, B) la direction pointe -> pivot (R_A R_B (0, 0, 1), sens d'impression
de l'etape 8), pivot = pointe + L d. Le firmware interpolant X, Y, Z, A, B lineairement, la pointe s'ecarte de sa
droite pendant un changement d'angle, de L |(d0 + d1) / 2 - d(milieu)| au milieu de chaque morceau : chaque mouvement
est coupe en 2^k morceaux egaux, k le plus petit qui tient cet ecart sous `GcodeParams.max_tip_deviation`.

Ecriture en deux temps (2026-11-28, demande de l'utilisateur : « l'etape d'ecriture du G-code doit absolument etre
parallelisee ») :
  1. TOUT le calcul est vectoriel, sur tous les points a la fois : plan des mouvements (approche, impression,
     degagements), decoupage, pivots, avances, extrusion cumulee (somme cumulee) ;
  2. la mise en texte, independante d'un morceau de trajectoire a l'autre, est repartie entre des processus
     (`RuntimeOptions.gcode_workers`). Le texte ne depend pas du nombre de processus.
Avant : une boucle Python point par point, environ 9 min pour un Benchy.
"""
from __future__ import annotations

import io
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from typing import Optional

import numpy as np

from core import line_types
from core.geometry_utils import head_directions
from file_io import kinematics_description as _kin

_MAX_SPLIT_LEVEL = 12                     # au plus 4 096 morceaux par mouvement
_PARALLEL_MIN_LINES = 200_000             # en dessous, un seul processus (lancer les autres coute plus cher)
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COLUMNS = ("x", "y", "z", "a", "b", "feed", "e", "extrude", "comment")


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


# --------------------------------------------------------------------------------------------------------------------
# 1. Plan des mouvements (vectoriel)
# --------------------------------------------------------------------------------------------------------------------

def _plan(pos, ang, ext, skin, axes, params) -> dict:
    """Mouvements cibles, dans l'ordre : approche (degagement puis premier point), impression, trajets a vide (longs :
    retraction, degagement au-dessus du depart, au-dessus de l'arrivee, arrivee, remise), degagement final."""
    n = len(pos)
    lift = params.travel_clearance * head_directions(ang, axes)
    seg = np.r_[0.0, np.linalg.norm(np.diff(pos, axis=0), axis=1)]
    printing = (ext > 0) & (np.arange(n) > 0)
    long_travel = ~printing & (seg >= params.retract_min_travel) & (np.arange(n) > 0)
    count = np.where(long_travel, 3, 1)
    count[0] = 2
    owner = np.repeat(np.arange(n), count)
    step = np.arange(count.sum()) - np.repeat(np.cumsum(count) - count, count)
    last = step == count[owner] - 1
    src = np.where((step == 0) & long_travel[owner], owner - 1, owner)          # degagement au-dessus du depart
    tip = pos[src] + np.where(last[:, None], 0.0, 1.0) * lift[src]
    angles = ang[src]
    h = np.where(skin, params.skin_thickness, params.layer_height)
    fil_area = math.pi * (params.filament_diameter / 2) ** 2
    e = np.where(printing, seg * h * params.line_width * params.flow / fil_area, 0.0)
    move_e = np.where(last & printing[owner], e[owner], 0.0)
    speed = np.where(last & printing[owner], params.print_speed, params.travel_speed)
    pre = (step == 0) & long_travel[owner]
    post = last & long_travel[owner]
    comment = np.where(step == 0, np.where(owner == 0, 1, np.where(long_travel[owner], 2, 0)), 0)
    # degagement final, apres une retraction
    tip = np.vstack([tip, pos[-1:] + lift[-1:]])
    angles = np.vstack([angles, ang[-1:]])
    return {"tip": tip, "ang": angles, "e": np.r_[move_e, 0.0], "speed": np.r_[speed, params.travel_speed],
            "pre": np.r_[pre, True], "post": np.r_[post, False], "comment": np.r_[comment, 3],
            "owner": np.r_[owner, n - 1], "printing": np.r_[last & printing[owner], False]}


def _pieces(plan: dict, axes, L: float, limit: Optional[float]) -> np.ndarray:
    """Nombre de morceaux (puissance de 2) de chaque mouvement, le premier (depuis nulle part) entier."""
    m = len(plan["tip"])
    count = np.ones(m, dtype=np.int64)
    if limit is None or m < 2:
        return count
    a0, a1 = plan["ang"][:-1], plan["ang"][1:]
    todo = np.arange(1, m)
    for _level in range(_MAX_SPLIT_LEVEL):
        if not len(todo):
            break
        k = count[todo]
        row = np.repeat(np.arange(len(todo)), k)
        u = (np.arange(k.sum()) - np.repeat(np.cumsum(k) - k, k)) / k[row]
        du = 1.0 / k[row]
        lo, hi = a0[todo - 1][row], a1[todo - 1][row]
        da = hi - lo
        d0 = head_directions(lo + da * u[:, None], axes)
        d1 = head_directions(lo + da * (u + du)[:, None], axes)
        dm = head_directions(lo + da * (u + 0.5 * du)[:, None], axes)
        dev = L * np.linalg.norm(0.5 * (d0 + d1) - dm, axis=1)
        worst = np.zeros(len(todo))
        np.maximum.at(worst, row, dev)
        over = worst > limit
        count[todo[over]] *= 2
        todo = todo[over]
    return count


# --------------------------------------------------------------------------------------------------------------------
# 2. Lignes (vectoriel), puis texte (parallele)
# --------------------------------------------------------------------------------------------------------------------

def _lines(plan: dict, count: np.ndarray, axes, L: float, params) -> dict:
    """Colonnes de chaque ligne de mouvement (morceaux compris) et E cumule, retractions comprises."""
    m = len(plan["tip"])
    move = np.repeat(np.arange(m), count)
    j = np.arange(count.sum()) - np.repeat(np.cumsum(count) - count, count) + 1     # 1..count
    u = j / count[move]
    start = np.maximum(move - 1, 0)
    tip = plan["tip"][start] + (plan["tip"][move] - plan["tip"][start]) * u[:, None]
    angles = plan["ang"][start] + (plan["ang"][move] - plan["ang"][start]) * u[:, None]
    fw = pivot_positions(tip, angles, axes, L)
    e_piece = plan["e"][move] / count[move]
    r = params.retract_length
    # E cumule : chaque mouvement apporte -r (retraction avant), son extrusion, +r (remise apres)
    delta = np.where(plan["pre"], -r, 0.0) + plan["e"] + np.where(plan["post"], r, 0.0)
    before = np.r_[0.0, np.cumsum(delta)[:-1]] + np.where(plan["pre"], -r, 0.0)   # E apres la retraction du mouvement
    e_cum = before[move] + e_piece * j
    tip_prev = np.vstack([tip[:1], tip[:-1]])
    fw_prev = np.vstack([fw[:1], fw[:-1]])
    tip_d = np.linalg.norm(tip - tip_prev, axis=1)
    fw_d = np.linalg.norm(fw - fw_prev, axis=1)
    speed = plan["speed"][move]
    feed = np.where(tip_d > 1e-9, speed * 60.0 * fw_d / np.maximum(tip_d, 1e-12), params.rotation_feed)
    return {"fw": fw, "ang": angles, "feed": np.clip(feed, 1.0, params.max_feed), "e": e_cum,
            "extrude": plan["printing"][move], "move": move, "first": j == 1, "last": j == count[move],
            "tip_z": tip[:, 2], "before": before}


_COMMENTS = {1: " ; approche", 2: " ; degagement", 3: " ; degagement final"}


def _format_chunk(args) -> str:
    """Texte des lignes de mouvement d'un morceau (fonction de module : appelable dans un autre processus)."""
    cols, decimals = args
    x, y, z, a, b, feed, e, extrude, comment = cols
    f = f".{decimals}f"
    out = []
    for i in range(len(x)):
        line = (f"G01 X{x[i]:{f}} Y{y[i]:{f}} Z{z[i]:{f}} A{a[i]:{f}} B{b[i]:{f}} C{0.0:{f}}"
                + (f" E{e[i]:{f}}" if extrude[i] else "") + f" F{feed[i]:{f}}" + _COMMENTS.get(int(comment[i]), ""))
        out.append(line)
    return "\n".join(out)


def _workers() -> int:
    from core import runtime
    chosen = runtime.state().options.gcode_workers
    return max(1, int(chosen) if chosen is not None else min(8, os.cpu_count() or 1))


def _text(rows: dict, comment: np.ndarray, decimals: int) -> list[str]:
    """Lignes de mouvement, un morceau de trajectoire par processus ; meme texte quel que soit leur nombre."""
    cols = (rows["fw"][:, 0], rows["fw"][:, 1], rows["fw"][:, 2], rows["ang"][:, 0], rows["ang"][:, 1],
            rows["feed"], rows["e"], rows["extrude"], comment)
    n = len(cols[0])
    workers = _workers() if n >= _PARALLEL_MIN_LINES else 1
    bounds = np.linspace(0, n, workers + 1).astype(int)
    chunks = [tuple(c[lo:hi] for c in cols) for lo, hi in zip(bounds[:-1], bounds[1:])]
    if workers == 1:
        texts = [_format_chunk((c, decimals)) for c in chunks]
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:      # chaque fil attend SON processus
            texts = list(pool.map(lambda c: _format_in_process(c, decimals), chunks))
    return "\n".join(t for t in texts if t).split("\n") if n else []


def _format_in_process(chunk: tuple, decimals: int) -> str:
    """Texte d'un morceau par `file_io.gcode_text_worker`. Si le processus echoue, le meme texte est produit ici
    (fonction identique) : le fichier ne depend jamais de la reussite d'un processus."""
    buf = io.BytesIO()
    np.savez(buf, decimals=decimals, **dict(zip(COLUMNS, chunk)))
    try:
        done = subprocess.run([sys.executable, "-m", "file_io.gcode_text_worker"], input=buf.getvalue(),
                              capture_output=True, cwd=_ROOT, timeout=600)
        if done.returncode == 0:
            return done.stdout.decode("ascii")
        print(f"[gcode] processus d'ecriture en echec ({done.returncode}) : morceau ecrit ici", file=sys.stderr)
    except (OSError, subprocess.TimeoutExpired) as exc:
        print(f"[gcode] processus d'ecriture indisponible ({exc}) : morceau ecrit ici", file=sys.stderr)
    return _format_chunk((chunk, decimals))


def _assemble(moves: list[str], rows: dict, plan: dict, types, params) -> list[str]:
    """Lignes de mouvement + lignes intercalees (retractions, etiquettes `;TYPE:`), dans l'ordre."""
    fmt = lambda v: f"{v:.{params.decimals}f}"   # noqa: E731
    feed_r = fmt(params.retract_speed * 60.0)
    before, after = {}, {}
    first_rows = np.flatnonzero(rows["first"])
    last_rows = np.flatnonzero(rows["last"])
    r = params.retract_length
    for k in np.flatnonzero(plan["pre"]):
        before.setdefault(int(first_rows[k]), []).append(f"G01 E{fmt(rows['before'][k])} F{feed_r}")
    for k in np.flatnonzero(plan["post"]):
        row = int(last_rows[k])
        after.setdefault(row, []).append(f"G01 E{fmt(rows['e'][row] + r)} F{feed_r}")
    if types is not None:
        owner = plan["owner"]
        last_type = None
        for k, row in zip(np.flatnonzero(plan["printing"]), first_rows[np.flatnonzero(plan["printing"])]):
            t = int(types[owner[k]])
            if t != last_type:
                before.setdefault(int(row), []).insert(0, f";TYPE:{line_types.LABELS.get(t, 'CUSTOM')}")
                last_type = t
    out = []
    marks = sorted(set(before) | set(after))
    cursor = 0
    for row in marks:
        out += moves[cursor:row]
        out += before.get(row, [])
        if row in after:
            out.append(moves[row])
            out += after[row]
            cursor = row + 1
        else:
            cursor = row
    out += moves[cursor:]
    return out


def _header(params, machine, L: float) -> list:
    from file_io.apex_signature import SIGNATURE
    out = [SIGNATURE, "; Apex -- G-code machine XYZAB (file_io/gcode_xyzab.py)",
           f"; pivot->pointe {L:g} mm compense ; retraction {params.retract_length} mm a {params.retract_speed} mm/s",
           "G21", "G90", "M82", "G92 E0", f"M109 S{params.nozzle_temp:g}", "; apex_sliced = 1"]
    out += _kin.header_lines(description(machine, L))
    from file_io.gcode_export import _encrypted_header
    return out + _encrypted_header(params, machine, description(machine, L), {"pivot_to_tip": float(L)})


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
    rep = GcodeReport()
    out = _header(params, machine, L)
    if len(pos):
        if params.tip_at_bead_top:
            half = np.where(skin, params.skin_thickness, params.layer_height) * 0.5
            pos = pos + head_directions(ang, axes) * half[:, None]
        plan = _plan(pos, ang, ext, skin, axes, params)
        count = _pieces(plan, axes, L, params.max_tip_deviation)
        rows = _lines(plan, count, axes, L, params)
        comment = np.where(rows["last"], plan["comment"][rows["move"]], 0)
        moves = _text(rows, comment, params.decimals)
        out += _assemble(moves, rows, plan, types, params)
        _report(rep, rows, plan, count, pos, ext, params)
    out += ["M104 S0", "M84"]
    rep.lines = len(out)
    if params.strict and rep.tip_below_plate:
        return out, rep
    if path:
        rep.written = True
        with open(path, "w", newline="\n") as f:
            f.write("\n".join(out) + "\n")
    return out, rep


def _report(rep, rows: dict, plan: dict, count: np.ndarray, pos: np.ndarray, ext: np.ndarray, params) -> None:
    rep.moves = int(len(rows["move"]))
    rep.subdivided_moves = int((count > 1).sum())
    rep.retractions = int(plan["post"].sum())
    rep.travels = int((ext[1:] <= 0).sum() + 1)
    rep.filament_mm = float(plan["e"].sum())
    rep.tip_z_min = float(rows["tip_z"].min())
    rep.tip_below_plate = int((rows["tip_z"] < 0.0).sum())
    fw = rows["fw"]
    rep.fw_x_range = (float(fw[:, 0].min()), float(fw[:, 0].max()))
    rep.fw_y_range = (float(fw[:, 1].min()), float(fw[:, 1].max()))
    rep.fw_z_range = (float(fw[:, 2].min()), float(fw[:, 2].max()))
    if len(pos) > 1:
        seg = np.linalg.norm(np.diff(pos, axis=0), axis=1)
        pr = ext[1:] > 0
        rep.print_time_nominal_s = float(seg[pr].sum() / params.print_speed + seg[~pr].sum() / params.travel_speed)
        rep.print_time_s = rep.print_time_nominal_s
