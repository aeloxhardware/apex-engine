"""Etape 8 : résolution des angles d'axes (grille, analytique), singularités, angle couplé.

Découpé de `step8_machine_toolpath.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step8_machine_toolpath`.
"""

from __future__ import annotations

from typing import Optional
import numpy as np
from core.geometry_utils import _rodrigues_rotate, composed_rotation_angles_and_directions, refine_axis_angles_local_search, resolve_axis_angles_analytic
from core.mesh_types import RotationAxisSpec
from element_ardor.step8_params import (  # noqa: F401
    Step8Params,
)


def _native_grid_args(params: "Step8Params") -> tuple[int, float, int]:
    if not params.native_grid_fallback:
        return (0, 0.0, 0)
    return (int(params.ik_samples_per_axis), float(params.ik_refine_step_multiplier), int(params.ik_refine_iterations))


def _resolve_axis_angles_via_grid_search(
    desired_direction: np.ndarray,
    axes: list[RotationAxisSpec],
    reference_direction: np.ndarray,
    params: Step8Params,
) -> tuple[list[float], float]:
    """
    Recherche par grille + affinement local - methode HISTORIQUE, conservee
    comme repli de `_resolve_axis_angles` (voir sa docstring) pour le cas
    rare ou la resolution analytique ne trouve geometriquement aucune
    solution (paire externe incapable d'atteindre l'angle cible, verifie
    reel et non hypothetique pendant la conception - CLAUDE.md/plan).

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


def _resolve_axis_angles(
    desired_direction: np.ndarray,
    axes: list[RotationAxisSpec],
    reference_direction: np.ndarray,
    params: Step8Params,
    previous_angles: Optional[list[float]] = None,
) -> tuple[list[float], float]:
    """
    Resout `desired_direction` en un angle par axe de `axes`, generiquement
    pour N axes (0, 1, 2, ...), par cinematique inverse ANALYTIQUE
    (`core.geometry_utils.resolve_axis_angles_analytic` - bissection + atan2,
    voir sa docstring pour la justification geometrique complete et
    CLAUDE.md "Piste chiffree pour step8" pour le diagnostic qui a motive
    ce fix : la grille historique avait une erreur mediane mesuree de
    9.26 deg contre un seuil de rejet de 5 deg, echouant structurellement
    son propre critere).

    Repli vers `_resolve_axis_angles_via_grid_search` (methode historique,
    conservee intacte) UNIQUEMENT quand la resolution analytique ne trouve
    geometriquement aucune solution pour la paire d'axes externe - cas
    confirme reel (pas hypothetique) pendant la conception sur des axes
    non orthogonaux, jamais rencontre sur la machine a 2 axes orthogonaux
    actuelle.

    `previous_angles` (angles resolus au point de trajectoire precedent,
    None au premier point) sert de valeur de repliement pour les axes
    internes au-dela de la paire externe (axes[2:]) et de repli en cas de
    singularite locale - voir docstring de resolve_axis_angles_analytic.

    Retourne (angles_deg, erreur_angulaire_residuelle_rad) - l'erreur est
    l'angle entre `desired_direction` et la direction reellement atteinte
    par les angles resolus (recalculee par cinematique directe), jamais
    supposee nulle.
    """
    if not axes:
        return [], 0.0

    def _fallback() -> tuple[list, float]:
        return _resolve_axis_angles_via_grid_search(
            desired_direction, axes, reference_direction, params
        )

    angles, residual_error_rad = resolve_axis_angles_analytic(
        axes, reference_direction, desired_direction, previous_angles, _fallback
    )

    return [float(a) for a in angles], float(residual_error_rad)


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


def coupled_axis_angle(
    position: np.ndarray,
    platform_origin: np.ndarray,
    axis: np.ndarray,
    radial_direction: np.ndarray,
    singularity_radius: float,
) -> Optional[float]:
    """
    Angle (degres, [0, 360)) de l'axe COUPLE A LA POSITION (limite #59) : la
    rotation autour de `axis` qui oriente `radial_direction` (+R a angles
    nuls) vers le point. La tete ne se deplacant que dans ce demi-plan, c'est
    le seul angle qui permet d'atteindre le point. `None` sur l'axe lui-meme
    (distance < singularity_radius) : l'angle y est indetermine (limite #19).
    """
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    rel = np.asarray(position, dtype=float) - np.asarray(platform_origin, dtype=float)
    rel = rel - axis * float(rel @ axis)
    if np.linalg.norm(rel) < singularity_radius:
        return None
    e1 = np.asarray(radial_direction, dtype=float)
    e1 = e1 - axis * float(e1 @ axis)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(axis, e1)
    return float(np.degrees(np.arctan2(rel @ e2, rel @ e1)) % 360.0)
