"""
pipeline/step8_machine_toolpath.py

Etape 8 du pipeline : resolution de la trajectoire machine reelle a partir
de PrintSchedule (etape 7) et des couches recolorees qu'il reference.

Pour chaque point de trajectoire, ce fichier :
  1. Releve la position 2D locale vers le 3D (`step6._lift_to_3d`).
  2. Resout la direction d'impression desiree (normale du repere local de
     la couche) en angles concrets pour CHAQUE axe de
     `MachineProfile.rotation_axes` - generiquement, pour N axes, jamais 2
     codes en dur (voir CLAUDE.md limite de conception #18).
  3. Detecte et gere la singularite au centre de rotation (distance nulle a
     l'axe -> angle herite du point precedent plutot que recalcule).
  4. Signale (sans jamais lever d'exception) toute orientation inatteignable
     dans les plages d'axes disponibles - donnee pure vers l'etape 10.
  5. Lisse la serie temporelle d'angles pour respecter
     max_angular_velocity/max_angular_acceleration.
  6. Insere des trajets a vide (extrusion=0) aux frontieres de polyligne.

Ce fichier ne connait RIEN d'un firmware particulier : sa sortie est
position cartesienne + angles generiques + extrusion/feedrate. La
traduction vers une commande machine precise (ordre X/Y/Z/B, noms d'axes)
est repoussee vers un futur module d'export (io/gcode_export.py), hors
perimetre ici.

Depend de list[Layer] (retour de step7, toolpath_2d re-decoupe par couleur),
PrintSchedule (step7, ordre d'execution gele), MachineProfile.

Sortie : (ToolpathMachine, list[Step8Diagnostic]).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from core.geometry_utils import (
    _rodrigues_rotate,
    composed_rotation_angles_and_directions,
    refine_axis_angles_local_search,
)
from core.mesh_types import (
    Layer,
    MachineProfile,
    PrintPass,
    PrintSchedule,
    RotationAxisSpec,
    Step8Diagnostic,
    ToolpathMachine,
    ToolpathPoint,
)
from pipeline.step6_toolpath2d import _lift_to_3d


class ToolpathResolutionError(Exception):
    """
    Levee uniquement pour une entree structurellement invalide (toolpath_2d
    manquant, schedule pointant vers un segment/couche inexistant) - JAMAIS
    pour une orientation inatteignable, qui est un cas attendu couvert par
    Step8Diagnostic, pas une erreur bloquante.
    """


@dataclass
class Step8Params:
    # Nombre d'echantillons par axe pour la recherche du plus proche voisin
    # (composed_rotation_angles_and_directions) - cout total = ~samples**N.
    ik_samples_per_axis: int = 10
    ik_refine_iterations: int = 8
    # Multiplicateur du pas initial d'affinement, PAS une constante fixe en
    # degres : le pas de depart = ik_refine_step_multiplier * (360 /
    # ik_samples_per_axis), a l'echelle de l'ecart entre deux points voisins
    # de la grille NN (jusqu'a ~36deg pour samples_per_axis=10). Un pas fixe
    # trop petit laisse l'affinement piege pres du candidat NN au lieu de
    # converger vers le vrai optimum (mesure : residu ~15-20deg avec un pas
    # fixe de 2deg ; <0.1deg avec ce pas derive - voir tests).
    ik_refine_step_multiplier: float = 0.5
    # Distance en-deca de laquelle la position courante est consideree SUR
    # l'axe de rotation (singularite : angle indetermine, herite du point
    # precedent plutot que recalcule).
    singularity_radius_epsilon: float = 1e-3
    # Feedrate des trajets a vide inseres aux frontieres de polyligne. 0 =>
    # reutilise le feedrate du point suivant (pas de vitesse magique par
    # defaut non documentee).
    travel_feedrate: float = 0.0
    default_feedrate: float = 1000.0
    # Au-dela de ce seuil, l'erreur residuelle d'un point est jugee assez
    # grande pour meriter un Step8Diagnostic (signal vers l'etape 10).
    max_unreachable_angular_error_deg: float = 5.0


# ---------------------------------------------------------------------------
# Direction desiree et resolution cinematique inverse generique
# ---------------------------------------------------------------------------

def _layer_normal(local_frame: tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    """
    Normale du repere local d'une couche = direction d'impression desiree en
    tout point de cette couche (meme lecture que step3/step4/step5 : le
    champ scalaire encode une direction d'impression par construction).
    """
    _origin, axis_u, axis_v = local_frame
    normal = np.cross(axis_u, axis_v)
    norm = np.linalg.norm(normal)
    if norm < 1e-12:
        return np.array([0.0, 0.0, 1.0])
    return normal / norm


def _resolve_axis_angles(
    desired_direction: np.ndarray,
    axes: list[RotationAxisSpec],
    reference_direction: np.ndarray,
    params: Step8Params,
) -> tuple[list[float], float]:
    """
    Resout `desired_direction` en un angle par axe de `axes`, generiquement
    pour N axes (0, 1, 2, ...), sans cinematique inverse analytique
    (impossible a ecrire une fois pour toute configuration d'axes
    arbitraire). Deux etapes :

      1. Recherche du plus proche voisin parmi les directions atteignables
         par composed_rotation_angles_and_directions (grille produit).
      2. Affinement local generique (refine_axis_angles_local_search),
         avec un pas initial derive de l'ecart de la grille NN - voir
         Step8Params.ik_refine_step_multiplier.

    Retourne (angles_deg, erreur_angulaire_residuelle_rad) - l'erreur est
    l'angle entre `desired_direction` et la direction reellement atteinte
    par les angles resolus (recalculee par cinematique directe), jamais
    supposee nulle.
    """
    if not axes:
        return [], 0.0

    directions, angles_grid = composed_rotation_angles_and_directions(
        axes, reference_direction, samples_per_axis=max(params.ik_samples_per_axis, 1)
    )

    desired = np.asarray(desired_direction, dtype=float)
    desired_norm = np.linalg.norm(desired)
    if desired_norm > 1e-12:
        desired = desired / desired_norm

    similarity = directions @ desired
    best_index = int(np.argmax(similarity))
    initial_angles = angles_grid[best_index]

    grid_spacing_deg = 360.0 / max(params.ik_samples_per_axis, 1)
    initial_step_deg = params.ik_refine_step_multiplier * grid_spacing_deg

    refined_angles, residual_error_rad = refine_axis_angles_local_search(
        axes,
        reference_direction,
        initial_angles,
        desired,
        iterations=params.ik_refine_iterations,
        initial_step_deg=initial_step_deg,
    )

    return [float(a) for a in refined_angles], float(residual_error_rad)


# ---------------------------------------------------------------------------
# Singularite au centre de rotation
# ---------------------------------------------------------------------------

def _axis_chain_world_lines(
    axes: list[RotationAxisSpec], angles_deg: list[float], platform_origin: np.ndarray
) -> list[np.ndarray]:
    """
    Direction monde (unitaire) de l'axe de rotation de chaque DOF, apres
    application des rotations des DOF qui le precedent dans la chaine, aux
    angles courants `angles_deg`. Toutes les rotations partagent la meme
    origine `platform_origin` (RotationAxisSpec n'a pas d'origine propre par
    axe dans ce modele) : seule la DIRECTION de chaque droite de rotation
    varie le long de la chaine, pas son point d'ancrage.

    Necessaire generiquement (pas seulement pour l'axe 0) car un axe interne
    (ex: la tete B) peut lui-meme devenir singulier selon la geometrie de la
    machine decrite - voir _detect_singular_axes.
    """
    lines: list[np.ndarray] = []
    accumulated_rotation_angles: list[float] = []
    accumulated_axes: list[np.ndarray] = []

    for axis_index, spec in enumerate(axes):
        axis = np.asarray(spec.axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        # applique les rotations des DOF PRECEDENTS (exterieurs) a cet axe
        world_axis = axis.reshape(1, 3)
        for prior_axis, prior_angle_deg in zip(
            reversed(accumulated_axes), reversed(accumulated_rotation_angles)
        ):
            world_axis = _rodrigues_rotate(
                world_axis, prior_axis, np.array([np.radians(prior_angle_deg)])
            )
        lines.append(world_axis[0])
        accumulated_axes.append(axis)
        accumulated_rotation_angles.append(angles_deg[axis_index])

    return lines


def _detect_singular_axes(
    position: np.ndarray,
    platform_origin: Optional[np.ndarray],
    axes: list[RotationAxisSpec],
    angles_deg: list[float],
    epsilon: float,
) -> list[int]:
    """
    Indices des axes dont la droite de rotation (ancree a `platform_origin`,
    orientee selon `_axis_chain_world_lines`) passe a moins de `epsilon` de
    `position` - generalise "R≈0 => C indetermine" a "distance a l'axe i≈0
    => angle i indetermine", pour n'importe quel axe de la chaine, pas
    seulement le premier.
    """
    if platform_origin is None or not axes:
        return []

    lines = _axis_chain_world_lines(axes, angles_deg, platform_origin)
    relative = np.asarray(position, dtype=float) - np.asarray(platform_origin, dtype=float)

    singular: list[int] = []
    for axis_index, line_direction in enumerate(lines):
        projection = np.dot(relative, line_direction) * line_direction
        perpendicular_distance = float(np.linalg.norm(relative - projection))
        if perpendicular_distance < epsilon:
            singular.append(axis_index)
    return singular


def _position_to_radius(
    position: np.ndarray, platform_origin: Optional[np.ndarray], axes: list[RotationAxisSpec]
) -> Optional[float]:
    """
    Distance perpendiculaire de `position` a la droite portee par l'axe de
    rotation LE PLUS EXTERNE (axes[0]), passant par `platform_origin` -
    concept de "rayon" bien defini des qu'il existe au moins un axe de
    rotation ancre a une origine. Pour une machine a 0 axe de rotation
    (orientation libre), aucun rayon n'a de sens : retourne None plutot
    qu'une valeur inventee.
    """
    if platform_origin is None or not axes:
        return None

    axis = np.asarray(axes[0].axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    relative = np.asarray(position, dtype=float) - np.asarray(platform_origin, dtype=float)
    projection = np.dot(relative, axis) * axis
    return float(np.linalg.norm(relative - projection))


def _inherit_angle_at_singularity(
    singular_axes: list[int], resolved_angles: list[float], previous_angles: Optional[list[float]]
) -> list[float]:
    """
    Remplace l'angle de chaque axe signale singulier par l'angle du POINT
    PRECEDENT plutot que par la valeur recalculee (potentiellement instable,
    type atan2(0,0)) - generalise l'heritage a autant d'axes simultanement
    singuliers que necessaire, pas seulement l'axe 0.
    """
    if not singular_axes or previous_angles is None:
        return resolved_angles

    angles = list(resolved_angles)
    for axis_index in singular_axes:
        if axis_index < len(previous_angles):
            angles[axis_index] = previous_angles[axis_index]
    return angles


# ---------------------------------------------------------------------------
# Lissage temporel (vitesse/acceleration angulaire)
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Parcours du schedule et insertion des trajets a vide
# ---------------------------------------------------------------------------

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


def _iter_schedule_points(
    schedule: PrintSchedule, layers_by_id: dict[int, Layer], params: Step8Params
) -> list[_RawPoint]:
    """
    Parcourt `schedule.passes` DANS L'ORDRE (jamais re-derive : c'est
    l'ordre d'execution gele par step7), releve chaque segment 2D vers le
    3D, et insere un point de trajet a vide (`extrusion=0.0`) a chaque
    frontiere de polyligne - y compris entre deux passes non contigues
    (changement de couleur ou saut de couche). Une frontiere est detectee
    des que la cle (layer_id, index) du segment courant differe de celle du
    segment precedent, quelle que soit son origine (meme couche, couche
    differente, passe differente).
    """
    raw_points: list[_RawPoint] = []
    previous_key: Optional[tuple[int, int]] = None

    def _segment_polyline(layer_id: int, segment_index: int) -> tuple[Layer, np.ndarray]:
        layer = layers_by_id.get(layer_id)
        if layer is None:
            raise ToolpathResolutionError(
                f"PrintSchedule reference la couche {layer_id}, absente de `layers`."
            )
        toolpath = layer.toolpath_2d
        if toolpath is None:
            raise ToolpathResolutionError(
                f"La couche {layer_id} n'a pas de toolpath_2d - "
                f"pipeline/step7_scheduling.py doit etre execute avant l'etape 8."
            )
        if segment_index >= len(toolpath):
            raise ToolpathResolutionError(
                f"PrintSchedule reference le segment {segment_index} de la couche "
                f"{layer_id}, hors limites ({len(toolpath)} segments disponibles)."
            )
        return layer, toolpath[segment_index]

    for pass_ in schedule.passes:
        for layer_id, segment_index in pass_.segments:
            layer, polyline_2d = _segment_polyline(layer_id, segment_index)
            if len(polyline_2d) == 0:
                continue

            points_3d = _lift_to_3d(polyline_2d, layer.local_frame)
            desired_direction = _layer_normal(layer.local_frame)
            key = (layer_id, segment_index)

            if previous_key is not None and key != previous_key:
                raw_points.append(
                    _RawPoint(
                        position=points_3d[0],
                        desired_direction=desired_direction,
                        layer_id=layer_id,
                        segment_index=segment_index,
                        color_id=pass_.color_id,
                        feedrate=(
                            params.travel_feedrate
                            if params.travel_feedrate > 0
                            else params.default_feedrate
                        ),
                        extrusion=0.0,
                        is_travel=True,
                    )
                )

            for point_3d in points_3d:
                raw_points.append(
                    _RawPoint(
                        position=point_3d,
                        desired_direction=desired_direction,
                        layer_id=layer_id,
                        segment_index=segment_index,
                        color_id=pass_.color_id,
                        feedrate=params.default_feedrate,
                        extrusion=1.0,
                        is_travel=False,
                    )
                )

            previous_key = key

    return raw_points


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def run(
    layers: list[Layer],
    schedule: PrintSchedule,
    machine_profile: MachineProfile,
    params: Step8Params,
) -> tuple[ToolpathMachine, list[Step8Diagnostic]]:
    """
    Resout la trajectoire machine reelle a partir du planning d'impression.

    Args:
        layers: sortie de pipeline/step7_scheduling.py (`colored_layers`) -
            PAS les couches d'entree de step7 : `PrintPass.segments` pointe
            vers ces couches-la (toolpath_2d re-decoupe par couleur).
        schedule: sortie de pipeline/step7_scheduling.py. `schedule.passes`
            est consomme dans l'ordre exact, jamais re-derive.
        machine_profile: cinematique de la machine (rotation_axes generique,
            longueur arbitraire) et origine du plateau.
        params: parametres de resolution IK, de singularite et de lissage.

    Returns:
        (ToolpathMachine, diagnostics) - diagnostics est une liste plate de
        Step8Diagnostic, une entree par point dont l'orientation desiree
        n'a pas pu etre atteinte a moins de
        `params.max_unreachable_angular_error_deg` pres. Jamais
        d'exception pour ce cas : c'est une donnee pure vers l'etape 10.

    Raises:
        ToolpathResolutionError: si le schedule reference une couche ou un
            segment absent/non rempli (donnee structurellement invalide,
            distinct d'une orientation inatteignable).
    """
    layers_by_id = {layer.id: layer for layer in layers}
    raw_points = _iter_schedule_points(schedule, layers_by_id, params)

    axes = machine_profile.rotation_axes or []
    reference_direction = np.array([0.0, 0.0, -1.0])
    platform_origin = machine_profile.platform_origin

    resolved_angles: list[list[float]] = []
    residual_errors: list[float] = []
    previous_angles: Optional[list[float]] = None

    for raw_point in raw_points:
        angles, residual_error_rad = _resolve_axis_angles(
            raw_point.desired_direction, axes, reference_direction, params
        )

        singular_axes = _detect_singular_axes(
            raw_point.position, platform_origin, axes, angles, params.singularity_radius_epsilon
        )
        angles = _inherit_angle_at_singularity(singular_axes, angles, previous_angles)

        resolved_angles.append(angles)
        residual_errors.append(residual_error_rad)
        previous_angles = angles

    smoothed_angles = _smooth_all_axes(resolved_angles, raw_points, axes, machine_profile)

    points: list[ToolpathPoint] = []
    diagnostics: list[Step8Diagnostic] = []

    for i, raw_point in enumerate(raw_points):
        points.append(
            ToolpathPoint(
                position=raw_point.position,
                axis_angles=smoothed_angles[i],
                extrusion=raw_point.extrusion,
                feedrate=raw_point.feedrate,
                color_id=raw_point.color_id,
            )
        )

        residual_error_deg = np.degrees(residual_errors[i])
        if residual_error_deg > params.max_unreachable_angular_error_deg:
            diagnostics.append(
                Step8Diagnostic(
                    point_index=i,
                    layer_id=raw_point.layer_id,
                    segment_index=raw_point.segment_index,
                    residual_angle_error_deg=float(residual_error_deg),
                    position=raw_point.position,
                )
            )

    return ToolpathMachine(points=points), diagnostics


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
    for i in range(1, n):
        distance = float(np.linalg.norm(raw_points[i].position - raw_points[i - 1].position))
        feedrate = max(raw_points[i].feedrate, 1e-6)
        dt_per_point[i] = distance / feedrate

    angles_matrix = np.array(resolved_angles, dtype=float)  # (n, len(axes))
    smoothed_matrix = angles_matrix.copy()

    for axis_index, spec in enumerate(axes):
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
