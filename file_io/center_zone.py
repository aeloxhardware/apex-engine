"""
file_io/center_zone.py

Zone centrale SANS EXTRUSION autour de l'axe du plateau (2026-10-24), appliquee
A L'EXPORT G-CODE (demande de l'utilisateur : « quand on est immobile au centre,
alors retracter »).

Sur la machine actuelle l'angle du plateau est impose par la position du point
(limite #59) : une ligne qui passe pres de l'axe oblige le plateau a faire un
demi-tour presque instantane. Les axes ayant des vitesses et accelerations
limitees, le firmware ralentit tout le mouvement et la pointe quasi immobile
continue de recevoir son filament : suintement, amas (la « montagne » du
centre). Aucune retractation n'etait emise, la ligne etant une polyligne
continue.

Cette etape traite le centre comme une ZONE DE RETRACTATION : tout segment
extrude est decoupe au bord du disque de rayon `radius` centre sur l'axe, et la
partie interieure devient un trajet a vide. L'export G-code (`file_io.
gcode_export`) emet alors la retractation, le degagement et la
de-retractation tout seul (il retracte a chaque point sans extrusion), et la
simulation de collision n'ajoute aucune matiere dans la zone.

Appliquee JUSTE AVANT l'ecriture du fichier (`main.export_gcode_checked`) : la
trajectoire que voient step9 et step9b n'est pas modifiee (la simulation de
collision reste prudente : elle suppose la matiere presente).

Contrat avec step9 : une capsule n'existe qu'entre deux points consecutifs
extrudes. Apres la zone, on insere donc un point de trajet (extrusion 0) puis
un point extrude au meme endroit, comme step8 le fait au debut de chaque
polyligne.

Hypothese : l'axe du plateau est vertical (`rotation_axes[0].axis == +Z`).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np

from core.mesh_types import MachineProfile, ToolpathMachine, ToolpathPoint
from pipeline.step8_machine_toolpath import coupled_axis_angle

_T_EPS = 1e-9


@dataclass
class CenterZoneReport:
    """Ce que la zone a change - jamais un silence sur de la matiere retiree."""
    applied: bool = False
    skipped_reason: str = ""
    radius: float = 0.0
    crossings: int = 0            # segments extrudes qui traversaient le disque
    removed_length: float = 0.0   # mm de cordon NON depose (parties interieures)
    points_added: int = 0


def _segment_disc_params(a_xy: np.ndarray, b_xy: np.ndarray, centre: np.ndarray, radius: float) -> list[float]:
    """Parametres t dans (0, 1) ou le segment coupe le cercle de rayon `radius` (0, 1 ou 2)."""
    d = b_xy - a_xy
    f = a_xy - centre
    a = float(d @ d)
    if a < 1e-18:
        return []
    b = 2.0 * float(f @ d)
    c = float(f @ f) - radius * radius
    disc = b * b - 4.0 * a * c
    if disc <= 0.0:
        return []
    root = float(np.sqrt(disc))
    ts = [(-b - root) / (2.0 * a), (-b + root) / (2.0 * a)]
    return [t for t in ts if _T_EPS < t < 1.0 - _T_EPS]


def _interpolated_angles(
    a: ToolpathPoint, b: ToolpathPoint, t: float, machine_profile: MachineProfile, position: np.ndarray
) -> list[float]:
    """Angles au parametre t : interpolation (arc le plus court pour un axe a tour complet),
    puis l'axe COUPLE a la position est recalcule depuis la position elle-meme."""
    axes = machine_profile.rotation_axes
    angles: list[float] = []
    for k, (ang_a, ang_b) in enumerate(zip(a.axis_angles, b.axis_angles)):
        full_turn = False
        if k < len(axes):
            lo, hi = axes[k].angle_range_deg
            full_turn = (hi - lo) >= 360.0 - 1e-9
        delta = ang_b - ang_a
        if full_turn:
            delta = (delta + 180.0) % 360.0 - 180.0
        value = ang_a + t * delta
        angles.append(float(value % 360.0) if full_turn else float(value))
    coupled = list(getattr(machine_profile, "position_coupled_axes", None) or [])
    if coupled and machine_profile.platform_origin is not None and machine_profile.radial_direction is not None:
        for k in coupled:
            c = coupled_axis_angle(
                position, np.asarray(machine_profile.platform_origin, dtype=float),
                np.asarray(machine_profile.rotation_axes[k].axis, dtype=float),
                np.asarray(machine_profile.radial_direction, dtype=float), 1e-6,
            )
            if c is not None:
                angles[k] = c
    return angles


def auto_center_radius(print_speed_mm_s: float, max_plate_speed_deg_s: float) -> float:
    """Rayon sous lequel la pointe doit ralentir sous `print_speed` parce que le plateau
    atteint sa vitesse angulaire max : v = omega * r  =>  r = v / omega."""
    return float(print_speed_mm_s / np.radians(max_plate_speed_deg_s))


def apply_center_retract_zone(
    toolpath: ToolpathMachine, machine_profile: MachineProfile, radius: Optional[float]
) -> tuple[ToolpathMachine, CenterZoneReport]:
    """
    Decoupe tout segment extrude au bord du disque de rayon `radius` (mm) centre
    sur l'axe du plateau et rend la partie interieure sans extrusion.

    Fonction pure : `toolpath` n'est jamais mute. `radius` None ou <= 0, ou
    machine sans axe de plateau vertical connu : trajectoire rendue telle
    quelle, avec la raison dans le rapport.
    """
    report = CenterZoneReport(radius=float(radius or 0.0))
    if radius is None or radius <= 0.0:
        report.skipped_reason = "desactive"
        return toolpath, report
    origin = machine_profile.platform_origin
    axes = machine_profile.rotation_axes
    if origin is None or not axes:
        report.skipped_reason = "pas d'axe de plateau dans le profil machine"
        return toolpath, report
    axis = np.asarray(axes[0].axis, dtype=float)
    if np.linalg.norm(axis) < 1e-12 or abs(axis[2] / np.linalg.norm(axis)) < 1.0 - 1e-6:
        report.skipped_reason = "l'axe du plateau n'est pas vertical"
        return toolpath, report

    centre = np.asarray(origin, dtype=float)[:2]
    pts = toolpath.points
    out: list[ToolpathPoint] = []
    if pts:
        out.append(pts[0])
    for i in range(1, len(pts)):
        a, b = pts[i - 1], pts[i]
        pa = np.asarray(a.position, dtype=float)
        pb = np.asarray(b.position, dtype=float)
        ts = _segment_disc_params(pa[:2], pb[:2], centre, radius) if b.extrusion > 0.0 else []
        # un segment entierement a l'interieur n'a pas de coupe : on le teste aussi
        starts_inside = float(np.linalg.norm(pa[:2] - centre)) < radius
        ends_inside = float(np.linalg.norm(pb[:2] - centre)) < radius
        if b.extrusion <= 0.0 or not (ts or (starts_inside and ends_inside)):
            out.append(b)
            continue

        bounds = [0.0] + ts + [1.0]
        crossed = False
        for s, e in zip(bounds[:-1], bounds[1:]):
            mid = pa + 0.5 * (s + e) * (pb - pa)
            inside = float(np.linalg.norm(mid[:2] - centre)) < radius
            if e >= 1.0 - _T_EPS:
                end_point = b
            else:
                position = pa + e * (pb - pa)
                end_point = replace(
                    b, position=position, axis_angles=_interpolated_angles(a, b, e, machine_profile, position)
                )
            if inside:
                crossed = True
                report.removed_length += (e - s) * float(np.linalg.norm(pb - pa))
                out.append(replace(end_point, extrusion=0.0))
                report.points_added += 0 if e >= 1.0 - _T_EPS else 1
                if e < 1.0 - _T_EPS:
                    # contrat step9 : une capsule exige deux points extrudes consecutifs
                    out.append(replace(end_point, extrusion=b.extrusion))
                    report.points_added += 1
            else:
                out.append(end_point)
                report.points_added += 0 if e >= 1.0 - _T_EPS else 1
        if crossed:
            report.crossings += 1

    report.applied = True
    return ToolpathMachine(points=out), report
