"""Etape 8 : lissage temporel des angles (limiteur de taux, filtre médian).

Découpé de `step8_machine_toolpath.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step8_machine_toolpath`.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import numpy as np
from core import _cpp
from core.mesh_types import MachineProfile, RotationAxisSpec


def _is_full_turn_axis(spec: RotationAxisSpec) -> bool:
    lo, hi = spec.angle_range_deg
    return (hi - lo) >= 360.0


def _unwrap_angle_sequence(angles_deg: np.ndarray, is_full_turn: bool) -> np.ndarray:
    """
    Convertit une sequence brute en sequence continue (deplie) via le delta
    au plus court chemin - UNIQUEMENT pour un axe a tour complet (359deg ->
    1deg doit se lire comme un pas de 2deg, pas 358deg). Une plage bornee
    strictement (ex: B in [0,180]) ne peut par construction pas boucler :
    depli desactive (no-op), le depliage y introduirait une discontinuite
    qui n'existe pas physiquement.
    """
    if len(angles_deg) == 0:
        return angles_deg
    if not is_full_turn:
        return angles_deg.copy()

    unwrapped = np.empty_like(angles_deg)
    unwrapped[0] = angles_deg[0]
    for i in range(1, len(angles_deg)):
        delta = ((angles_deg[i] - angles_deg[i - 1] + 180.0) % 360.0) - 180.0
        unwrapped[i] = unwrapped[i - 1] + delta
    return unwrapped


def _rewrap_angle_sequence(
    unwrapped_deg: np.ndarray, angle_range_deg: tuple[float, float], is_full_turn: bool
) -> np.ndarray:
    """Inverse de _unwrap_angle_sequence : replie dans la plage physique valide apres lissage."""
    if len(unwrapped_deg) == 0:
        return unwrapped_deg
    if is_full_turn:
        lo, _hi = angle_range_deg
        return lo + np.mod(unwrapped_deg - lo, 360.0)
    lo, hi = angle_range_deg
    return np.clip(unwrapped_deg, lo, hi)


def _clamp_rate(
    unwrapped_deg: np.ndarray,
    dt_per_point: np.ndarray,
    max_velocity_rad_s: float,
    max_accel_rad_s2: float,
) -> np.ndarray:
    """
    Limiteur de taux causal (avant uniquement) : applique sequentiellement
    les contraintes |delta_angle/delta_t| <= max_velocity puis
    |delta_vitesse/delta_t| <= max_accel entre points consecutifs.

    Prefere a une spline : les limites machine sont des contraintes DURES
    que ce filtre respecte par construction (clamp direct), alors qu'une
    spline necessiterait un ajustement iteratif sans garantie de
    convergence sous contrainte de derivees. `dt_per_point[i]` est le temps
    entre le point i-1 et i (issu du feedrate/distance, pas suppose
    uniforme) - le filtre reste donc generique sur un espacement de points
    quelconque, coherent avec l'ordre gele et non re-echantillonne de
    PrintSchedule.

    Si max_velocity_rad_s ou max_accel_rad_s2 <= 0 (non specifie), aucune
    limite n'est appliquee sur la grandeur correspondante - une machine
    pour laquelle ces limites ne sont pas encore connues ne doit pas voir
    sa trajectoire tronquee a zero mouvement.
    """
    n = len(unwrapped_deg)
    if n <= 1:
        return unwrapped_deg.copy()

    clamped = unwrapped_deg.copy()
    previous_velocity_deg_s = 0.0

    for i in range(1, n):
        dt = max(float(dt_per_point[i]), 1e-9)
        desired_delta_deg = clamped[i] - clamped[i - 1]
        desired_velocity_deg_s = desired_delta_deg / dt

        if max_velocity_rad_s > 0:
            max_velocity_deg_s = np.degrees(max_velocity_rad_s)
            desired_velocity_deg_s = float(
                np.clip(desired_velocity_deg_s, -max_velocity_deg_s, max_velocity_deg_s)
            )

        if max_accel_rad_s2 > 0:
            max_accel_deg_s2 = np.degrees(max_accel_rad_s2)
            max_delta_velocity = max_accel_deg_s2 * dt
            desired_velocity_deg_s = float(
                np.clip(
                    desired_velocity_deg_s,
                    previous_velocity_deg_s - max_delta_velocity,
                    previous_velocity_deg_s + max_delta_velocity,
                )
            )

        clamped[i] = clamped[i - 1] + desired_velocity_deg_s * dt
        previous_velocity_deg_s = desired_velocity_deg_s

    return clamped


@dataclass
class _RawPoint:
    position: np.ndarray
    desired_direction: np.ndarray
    layer_id: int
    segment_index: int
    color_id: Optional[int]
    feedrate: float
    extrusion: float
    is_travel: bool


@dataclass(slots=True, eq=False)
class _RawPoints:
    """Les points bruts de la trajectoire en TABLEAUX (une ligne par point),
    au lieu de ~300 000 objets `_RawPoint` (death @ 0,8). Memes valeurs, au
    bit pres ; l'acces point par point (`raw[i]`, iteration) reste possible
    pour les chemins de repli, qui recoivent alors un `_RawPoint`."""

    position: np.ndarray            # (N,3)
    desired_direction: np.ndarray   # (N,3)
    layer_id: np.ndarray            # (N,) int64
    segment_index: np.ndarray       # (N,) int64
    color_id: list                  # None autorise
    feedrate: np.ndarray            # (N,)
    extrusion: np.ndarray           # (N,)
    is_travel: np.ndarray           # (N,) bool

    def __len__(self) -> int:
        return len(self.position)

    def __getitem__(self, i: int) -> _RawPoint:
        return _RawPoint(
            position=self.position[i], desired_direction=self.desired_direction[i],
            layer_id=int(self.layer_id[i]), segment_index=int(self.segment_index[i]),
            color_id=self.color_id[i], feedrate=float(self.feedrate[i]),
            extrusion=float(self.extrusion[i]), is_travel=bool(self.is_travel[i]),
        )

    def __iter__(self):
        return (self[i] for i in range(len(self)))


def _point_arrays(raw_points) -> tuple[np.ndarray, np.ndarray]:
    """(positions, directions) en tableaux, depuis `_RawPoints` (sans copie)
    ou une liste de `_RawPoint` (tests, appelants historiques)."""
    if isinstance(raw_points, _RawPoints):
        return raw_points.position, raw_points.desired_direction
    return (np.array([p.position for p in raw_points], dtype=float),
            np.array([p.desired_direction for p in raw_points], dtype=float))


def _median_filter_free_axes(resolved_angles, raw_points, axes, machine_profile, reference_direction, window):
    """Filtre median des axes non couples, par polyligne (voir
    Step8Params.free_axis_median_window). Une polyligne s'arrete a chaque
    trajet a vide et a chaque changement de (couche, segment). Rend (angles,
    ecarts residuels recalculés en radians)."""
    from scipy.ndimage import median_filter

    angles = np.array(resolved_angles, dtype=float)
    coupled = set(getattr(machine_profile, "position_coupled_axes", None) or [])
    free = [i for i in range(len(axes)) if i not in coupled]
    native_median = (isinstance(raw_points, _RawPoints) and _cpp.enabled() and window % 2 == 1
                     and hasattr(_cpp.module(), "median_filter_free_axes"))
    if native_median:
        layer_ids = segments = travels = None
    elif isinstance(raw_points, _RawPoints):
        layer_ids, segments, travels = (raw_points.layer_id.tolist(), raw_points.segment_index.tolist(),
                                        raw_points.is_travel.tolist())
    else:
        layer_ids = [p.layer_id for p in raw_points]
        segments = [p.segment_index for p in raw_points]
        travels = [p.is_travel for p in raw_points]
    if native_median:
        _cpp.module().median_filter_free_axes(
            angles, np.asarray(free, dtype=np.int64), np.ascontiguousarray(raw_points.layer_id, dtype=np.int64),
            np.ascontiguousarray(raw_points.segment_index, dtype=np.int64),
            np.ascontiguousarray(raw_points.is_travel, dtype=np.uint8), int(window))
        keys = None
    else:
        keys = list(zip(layer_ids, segments, travels))
    start = 0
    for i in range(1, len(raw_points) + 1 if keys is not None else 0):
        if i == len(raw_points) or keys[i] != keys[start] or travels[i]:
            if i - start > 2 and not travels[start]:
                for ax in free:
                    angles[start:i, ax] = median_filter(angles[start:i, ax], size=window, mode="nearest")
            start = i
    units = [np.asarray(s.axis, dtype=float) / np.linalg.norm(s.axis) for s in axes]
    _positions, desired = _point_arrays(raw_points)
    # Ecart residuel recalcule PAR LOTS (2026-10-06 : la boucle point par point
    # coutait ~12 s sur death avec la normale par point). Meme formule de Rodrigues,
    # memes operations dans le meme ordre que `_rodrigues_rotate_single` ; seul le
    # cos/sin vectorise de numpy peut differer d'un ULP (sans effet sur un seuil en degres).
    n = len(raw_points)
    d = np.tile(np.asarray(reference_direction, dtype=float), (n, 1))
    for k in range(len(units) - 1, -1, -1):
        u = units[k]
        rad = np.radians(angles[:, k])
        cos_a, sin_a = np.cos(rad), np.sin(rad)
        dot = d[:, 0] * u[0] + d[:, 1] * u[1] + d[:, 2] * u[2]
        cx = u[1] * d[:, 2] - u[2] * d[:, 1]
        cy = u[2] * d[:, 0] - u[0] * d[:, 2]
        cz = u[0] * d[:, 1] - u[1] * d[:, 0]
        factor = dot * (1.0 - cos_a)
        d = np.column_stack([d[:, 0] * cos_a + cx * sin_a + u[0] * factor,
                             d[:, 1] * cos_a + cy * sin_a + u[1] * factor,
                             d[:, 2] * cos_a + cz * sin_a + u[2] * factor])
    t = np.asarray(desired, dtype=float).reshape(n, 3)
    t = t / np.maximum(np.linalg.norm(t, axis=1), 1e-12)[:, None]
    residual = np.arccos(np.clip((d * t).sum(axis=1), -1.0, 1.0)) if n else np.zeros(0)
    return [list(a) for a in angles], [float(r) for r in residual]


def _smooth_all_axes(
    resolved_angles: list[list[float]],
    raw_points: list["_RawPoint"],
    axes: list[RotationAxisSpec],
    machine_profile: MachineProfile,
) -> list[list[float]]:
    """
    Lisse la serie temporelle d'angles, axe par axe, sur l'ENSEMBLE des
    points (pas par passe : le lissage a besoin de la serie complete pour
    respecter les limites de vitesse/acceleration sur toute la trajectoire).
    Depli/repli (wraparound) pour les axes a tour complet, clamp de
    taux, puis reprojection - voir _unwrap_angle_sequence/_clamp_rate/
    _rewrap_angle_sequence.

    `dt_per_point[i]` derive du feedrate et de la distance parcourue entre
    les points i-1 et i, donnant un pas de temps physiquement significatif
    plutot qu'un taux d'echantillonnage suppose fixe.
    """
    n = len(raw_points)
    if n == 0 or not axes:
        return [list(a) for a in resolved_angles]

    dt_per_point = np.ones(n)
    positions, _desired = _point_arrays(raw_points)
    feedrates = (raw_points.feedrate.tolist() if isinstance(raw_points, _RawPoints)
                 else [p.feedrate for p in raw_points])
    for i in range(1, n):
        distance = float(np.linalg.norm(positions[i] - positions[i - 1]))
        feedrate = max(feedrates[i], 1e-6)
        dt_per_point[i] = distance / feedrate

    angles_matrix = np.array(resolved_angles, dtype=float)  # (n, len(axes))
    smoothed_matrix = angles_matrix.copy()

    coupled = set(getattr(machine_profile, "position_coupled_axes", None) or [])
    for axis_index, spec in enumerate(axes):
        if axis_index in coupled:
            # angle impose par la position : le lisser deplacerait le point
            # atteint (limite #59) ; c'est la vitesse d'avance qui devra
            # respecter la vitesse du plateau (a traiter a l'export G-code)
            continue
        is_full_turn = _is_full_turn_axis(spec)
        column = angles_matrix[:, axis_index]

        unwrapped = _unwrap_angle_sequence(column, is_full_turn)
        clamped = _clamp_rate(
            unwrapped,
            dt_per_point,
            machine_profile.max_angular_velocity,
            machine_profile.max_angular_acceleration,
        )
        rewrapped = _rewrap_angle_sequence(clamped, spec.angle_range_deg, is_full_turn)
        smoothed_matrix[:, axis_index] = rewrapped

    return [list(row) for row in smoothed_matrix]
