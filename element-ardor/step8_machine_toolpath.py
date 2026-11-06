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

from core import progress as _progress
from core import _cpp
from core.geometry_utils import (
    _rodrigues_rotate,
    composed_rotation_angles_and_directions,
    refine_axis_angles_local_search,
    resolve_axis_angles_analytic,
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

from pipeline.step8_params import (  # noqa: F401  (réexport : imports existants)
    Step8Params,
    ToolpathResolutionError,
)
from pipeline.step8_ik import (  # noqa: F401  (réexport : imports existants)
    _axis_chain_world_lines,
    _detect_singular_axes,
    _inherit_angle_at_singularity,
    _native_grid_args,
    _position_to_radius,
    _resolve_axis_angles,
    _resolve_axis_angles_via_grid_search,
    coupled_axis_angle,
)
from pipeline.step8_smoothing import (  # noqa: F401  (réexport : imports existants)
    _RawPoint,
    _RawPoints,
    _clamp_rate,
    _is_full_turn_axis,
    _median_filter_free_axes,
    _point_arrays,
    _rewrap_angle_sequence,
    _smooth_all_axes,
    _unwrap_angle_sequence,
)
from pipeline.step8_points import (  # noqa: F401  (réexport : imports existants)
    _apply_per_point_normals,
    _iter_schedule_points,
    _iter_schedule_points_native,
    _layer_normal,
    gradient_direction_field,
)


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def run(
    layers: list[Layer],
    schedule: PrintSchedule,
    machine_profile: MachineProfile,
    params: Step8Params,
    direction_field=None,
) -> tuple[ToolpathMachine, list[Step8Diagnostic]]:
    """
    Resout la trajectoire machine reelle a partir du planning d'impression.

    `direction_field` (optionnel, utilise si `params.per_point_normal`) :
    appelable `(M,3) positions -> (M,3)` directions unitaires +grad(phi),
    lignes NaN la ou le champ est inexploitable (repli sur la normale de
    couche). Construit par l'orchestrateur (`gradient_direction_field`).

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
    _progress.milestone(0.15)                # R-017 : points bruts assembles
    if params.per_point_normal and direction_field is not None and raw_points:
        _apply_per_point_normals(raw_points, layers_by_id, direction_field)
    _progress.milestone(0.3)                 # normales par point
    return resolve_raw_points(raw_points, machine_profile, params)


def resolve_raw_points(
    raw_points: "_RawPoints",
    machine_profile: MachineProfile,
    params: Step8Params,
) -> tuple[ToolpathMachine, list[Step8Diagnostic]]:
    """
    Points bruts (position + direction visee, extrusion, couche) -> trajectoire
    machine : IK (couplee a la position si la machine l'impose), filtre median
    des axes libres, lissage optionnel, diagnostics. Extrait de `run` le
    2026-10-07 sans changement de comportement, pour qu'un appelant qui produit
    ses points sans couches ni planning passe par exactement la meme resolution.
    """
    axes = machine_profile.rotation_axes or []
    # Direction a angles nuls = SENS D'IMPRESSION de la tete au repos (tete vers le bas, impression vers +Z) :
    # meme convention que step3/step7. Correction 2026-09-27 : la reference valait -Z (couche plate => B =
    # 179,9 deg, tete retournee, corps SOUS la buse).
    reference_direction = np.array([0.0, 0.0, 1.0])
    platform_origin = machine_profile.platform_origin

    if getattr(machine_profile, "position_coupled_axes", None) and axes:
        # limite #59 : premier axe impose par la position (machine R C Z B)
        resolved_angles, residual_errors = _resolve_all_axis_angles_coupled(
            raw_points, axes, reference_direction, machine_profile, params
        )
    else:
        resolved_angles, residual_errors = _resolve_all_axis_angles(
            raw_points, axes, reference_direction, platform_origin, params
        )

    _progress.milestone(0.85)                # cinematique inverse resolue (l'etape la plus longue, indivisible)
    if params.free_axis_median_window > 1 and axes:
        resolved_angles, residual_errors = _median_filter_free_axes(
            resolved_angles, raw_points, axes, machine_profile, reference_direction,
            params.free_axis_median_window,
        )

    if params.clamp_angle_rates:
        smoothed_angles = _smooth_all_axes(resolved_angles, raw_points, axes, machine_profile)
    else:
        # axes coordonnes par le firmware : angles exacts (voir clamp_angle_rates).
        # Listes neuves, propres a cet appel : aucune recopie necessaire.
        smoothed_angles = resolved_angles

    return _assemble_machine_toolpath(raw_points, smoothed_angles, residual_errors, params)


#: Issue du dernier appel a `_resolve_all_axis_angles` : `None` = chemin
#: Python (interrupteur eteint ou aucun point), sinon NOMBRE de points ayant
#: exige le repli grille cote Python (0 = natif de bout en bout).
#: Existe pour que les tests d'equivalence puissent VERIFIER que le natif a
#: bien produit le resultat - sans quoi un repli silencieux rendrait le test
#: vert en ne comparant que Python a lui-meme (trou reellement rencontre).
_LAST_NATIVE_FALLBACK_COUNT: Optional[int] = None


def _resolve_all_axis_angles_native(
    raw_points: list["_RawPoint"],
    axes: list[RotationAxisSpec],
    reference_direction: np.ndarray,
    platform_origin: Optional[np.ndarray],
    params: Step8Params,
) -> tuple[list[list[float]], list[float]]:
    """
    Chemin natif, avec REPRISE apres chaque repli.

    Le natif resout analytiquement et s'arrete au premier point que la
    resolution de la paire externe ne sait pas traiter (cible colineaire a
    l'axe externe : le residu du bracketing y est une tangence a zero, jamais
    un changement de signe - typiquement un depot vertical). Ce point-la est
    resolu par la recherche par grille Python, deja ecrite et testee, puis le
    natif reprend au suivant.

    **Reprendre au lieu d'abandonner est l'essentiel de ce mecanisme.** Une
    premiere version rejouait TOUTE la boucle en Python des le premier repli :
    le premier point d'une trajectoire a fond plat depose verticalement, donc
    UN point sur 297 994 suffisait a annuler l'integralite du gain (mesure :
    89,5 s contre 89,7 s, soit rien). La chaine `previous_angles` est passee
    explicitement d'un segment au suivant, donc l'equivalence est preservee.
    """
    global _LAST_NATIVE_FALLBACK_COUNT

    native = _cpp.module()
    positions, directions = _point_arrays(raw_points)
    axes_spec = [(np.asarray(spec.axis, dtype=float), *spec.angle_range_deg) for spec in axes]
    origin = None if platform_origin is None else np.asarray(platform_origin, dtype=float)

    resolved_angles: list[list[float]] = []
    residual_errors: list[float] = []
    previous_angles: Optional[list[float]] = None
    cursor = 0
    fallback_count = 0

    while cursor < len(raw_points):
        angles_matrix, residuals, fallback_index, grid_done = native.resolve_axis_angles_batch(
            directions[cursor:],
            positions[cursor:],
            axes_spec,
            np.asarray(reference_direction, dtype=float),
            origin,
            params.singularity_radius_epsilon,
            None if previous_angles is None else np.asarray(previous_angles, dtype=float),
            *_native_grid_args(params),
        )
        fallback_count += int(grid_done)

        # Sur repli, les lignes AVANT `fallback_index` sont valides ; les
        # suivantes n'ont jamais ete calculees et ne doivent pas etre lues.
        resolved_count = len(angles_matrix) if fallback_index < 0 else int(fallback_index)
        for k in range(resolved_count):
            resolved_angles.append(list(angles_matrix[k]))
            residual_errors.append(float(residuals[k]))
        if resolved_count:
            previous_angles = resolved_angles[-1]

        if fallback_index < 0:
            break

        fallback_count += 1
        raw_point = raw_points[cursor + resolved_count]
        angles, residual_error_rad = _resolve_axis_angles(
            raw_point.desired_direction, axes, reference_direction, params, previous_angles
        )
        singular_axes = _detect_singular_axes(
            raw_point.position, platform_origin, axes, angles, params.singularity_radius_epsilon
        )
        angles = _inherit_angle_at_singularity(singular_axes, angles, previous_angles)

        resolved_angles.append(angles)
        residual_errors.append(residual_error_rad)
        previous_angles = angles
        cursor += resolved_count + 1

    _LAST_NATIVE_FALLBACK_COUNT = fallback_count
    return resolved_angles, residual_errors


def _resolve_all_axis_angles_coupled(
    raw_points: list["_RawPoint"],
    axes: list[RotationAxisSpec],
    reference_direction: np.ndarray,
    machine_profile: MachineProfile,
    params: Step8Params,
) -> tuple[list[list[float]], list[float]]:
    """
    IK quand le premier axe est IMPOSE PAR LA POSITION (machine actuelle :
    R C Z B, pas d'axe Y ; limite #59). Pour chaque point : l'angle couple
    vient de la position (coupled_axis_angle), puis la direction cible est
    ramenee dans le repere de cet axe et les axes RESTANTS la resolvent (meme
    `_resolve_axis_angles` que le cas libre). L'erreur residuelle est celle de
    la chaine complete : sur la machine actuelle, B ne penche que dans le plan
    radial, donc une cible qui en sort garde une erreur egale a son angle au
    plan -- elle devient un Step8Diagnostic, jamais une exception.

    Sur l'axe (singularite), l'angle couple est herite du point precedent
    (0 au premier point), comme le cas libre (limite #19).

    Chemin Python uniquement, dans les deux positions de APEX_USE_CPP : le
    natif ne connait pas ce couplage (resultats identiques par construction).
    """
    coupled = list(machine_profile.position_coupled_axes)
    if coupled != [0]:
        raise ValueError(
            f"position_coupled_axes={coupled} : seul le premier axe peut etre couple a la position"
        )
    if machine_profile.radial_direction is None or machine_profile.platform_origin is None:
        raise ValueError("position_coupled_axes exige radial_direction et platform_origin")
    outer = axes[0]
    outer_axis = np.asarray(outer.axis, dtype=float)
    outer_axis = outer_axis / np.linalg.norm(outer_axis)
    inner_axes = axes[1:]

    if _cpp.enabled() and raw_points:
        return _resolve_all_axis_angles_coupled_native(
            raw_points, axes, outer_axis, reference_direction, machine_profile, params
        )

    resolved_angles: list[list[float]] = []
    residual_errors: list[float] = []
    previous: Optional[list[float]] = None
    for raw_point in raw_points:
        c = coupled_axis_angle(
            raw_point.position, machine_profile.platform_origin, outer_axis,
            machine_profile.radial_direction, params.singularity_radius_epsilon,
        )
        if c is None:
            c = previous[0] if previous is not None else 0.0
        target = np.asarray(raw_point.desired_direction, dtype=float)
        target = target / max(np.linalg.norm(target), 1e-12)
        local_target = _rodrigues_rotate(target[None, :], outer_axis, np.array([np.radians(-c)]))[0]
        inner, _res = _resolve_axis_angles(
            local_target, inner_axes, reference_direction, params,
            None if previous is None else previous[1:],
        )
        angles = [c] + list(inner)
        achieved = reference_direction
        for spec, angle_deg in zip(reversed(axes), reversed(angles)):
            ax = np.asarray(spec.axis, dtype=float)
            achieved = _rodrigues_rotate(achieved[None, :], ax / np.linalg.norm(ax),
                                         np.array([np.radians(angle_deg)]))[0]
        residual = float(np.arccos(np.clip(achieved @ target, -1.0, 1.0)))
        resolved_angles.append(angles)
        residual_errors.append(residual)
        previous = angles
    return resolved_angles, residual_errors


def _resolve_all_axis_angles_coupled_native(
    raw_points: list["_RawPoint"],
    axes: list[RotationAxisSpec],
    outer_axis: np.ndarray,
    reference_direction: np.ndarray,
    machine_profile: MachineProfile,
    params: Step8Params,
) -> tuple[list[list[float]], list[float]]:
    """
    Pas 17 du portage C++ (limite #56) : meme calcul que la boucle Python de
    `_resolve_all_axis_angles_coupled`, reorganise en trois temps : 1. angle couple C de tous les points
    (`_coupled_axis_angles`) ; 2. cibles ramenees dans le repere de C (Rodrigues par lots) ; 3. axes internes
    resolus par l'IK NATIVE du pas 2, SANS detection de singularite (la boucle Python n'en fait pas non plus
    pour les axes internes), avec la meme reprise au point fautif : la grille Python ne traite que ce point.

    Equivalence : TOLERANCE declaree, comme au pas 2 (angles a ~1e-9 deg) ; les produits scalaires vectorises
    peuvent differer d'un ULP de ceux de la boucle Python.
    """
    global _LAST_NATIVE_FALLBACK_COUNT
    native = _cpp.module()
    positions, targets = _point_arrays(raw_points)
    targets = targets / np.maximum(np.linalg.norm(targets, axis=1), 1e-12)[:, None]
    c = _coupled_axis_angles(positions, outer_axis, machine_profile, params)

    local_targets = _rodrigues_rotate(targets, outer_axis, np.radians(-c))
    inner_axes = axes[1:]
    inner_spec = [(np.asarray(spec.axis, dtype=float), *spec.angle_range_deg) for spec in inner_axes]

    # blocs de lignes (natif) et lignes isolees (repli), concatenes une fois : aucune conversion ligne par
    # ligne (0,26 s sur death @ 0,8)
    inner_blocks: list[np.ndarray] = []
    previous: Optional[list[float]] = None
    cursor = 0
    fallback_count = 0
    n = len(raw_points)
    width = len(inner_axes)
    while cursor < n:
        angles_matrix, _res, fallback_index, grid_done = native.resolve_axis_angles_batch(
            local_targets[cursor:], positions[cursor:], inner_spec,
            np.asarray(reference_direction, dtype=float), None,
            params.singularity_radius_epsilon,
            None if previous is None else np.asarray(previous, dtype=float),
            *_native_grid_args(params),
        )
        fallback_count += int(grid_done)
        count = len(angles_matrix) if fallback_index < 0 else int(fallback_index)
        if count:
            block = np.asarray(angles_matrix[:count], dtype=float).reshape(count, width)
            inner_blocks.append(block)
            previous = list(block[-1])
        if fallback_index < 0:
            break
        fallback_count += 1
        k = cursor + count
        inner, _r = _resolve_axis_angles(local_targets[k], inner_axes, reference_direction, params, previous)
        inner_blocks.append(np.asarray(list(inner), dtype=float).reshape(1, width))
        previous = list(inner_blocks[-1][-1])
        cursor = k + 1
    _LAST_NATIVE_FALLBACK_COUNT = fallback_count

    inner_matrix = np.concatenate(inner_blocks, axis=0) if inner_blocks else np.zeros((0, width))
    angles = np.column_stack([c, inner_matrix.reshape(n, width)])
    achieved = np.tile(np.asarray(reference_direction, dtype=float), (n, 1))
    for k in range(len(axes) - 1, -1, -1):
        ax = np.asarray(axes[k].axis, dtype=float)
        achieved = _rodrigues_rotate(achieved, ax / np.linalg.norm(ax), np.radians(angles[:, k]))
    residuals = np.arccos(np.clip(np.einsum("ij,ij->i", achieved, targets), -1.0, 1.0))
    return angles.tolist(), residuals.tolist()

def _resolve_all_axis_angles(
    raw_points: list["_RawPoint"],
    axes: list[RotationAxisSpec],
    reference_direction: np.ndarray,
    platform_origin: Optional[np.ndarray],
    params: Step8Params,
) -> tuple[list[list[float]], list[float]]:
    """
    Resout l'angle de chaque axe pour TOUS les points de la trajectoire.

    La boucle porte une dependance SEQUENTIELLE reelle : `previous_angles`
    sert de valeur de repliement pour les axes internes (`axes[2:]`), de
    repli en cas de singularite locale, et d'angle herite quand un axe est
    detecte singulier. Elle n'est donc PAS parallelisable telle quelle -
    contrairement a ce que laissait croire son allure de boucle par point.

    Chemin natif (portage C++, voir cpp/src/core/kinematics.cpp) quand
    l'interrupteur est actif ; sinon la boucle Python historique, conservee
    intacte juste en dessous.

    Le natif ne porte PAS la recherche par grille de repli
    (`_resolve_axis_angles_via_grid_search`) : quand la resolution analytique
    de la paire externe echoue, il rend la main en signalant l'index fautif,
    et TOUTE la boucle est rejouee en Python. Rejouer tout, et non le seul
    point concerne, est necessaire a l'equivalence : reprendre au milieu avec
    une chaine `previous_angles` differente divergerait silencieusement.
    Ce cas n'est jamais atteint sur la machine a 2 axes orthogonaux actuelle.
    """
    global _LAST_NATIVE_FALLBACK_COUNT
    _LAST_NATIVE_FALLBACK_COUNT = None

    if _cpp.enabled() and raw_points:
        return _resolve_all_axis_angles_native(
            raw_points, axes, reference_direction, platform_origin, params
        )

    resolved_angles: list[list[float]] = []
    residual_errors: list[float] = []
    previous_angles: Optional[list[float]] = None

    for raw_point in raw_points:
        angles, residual_error_rad = _resolve_axis_angles(
            raw_point.desired_direction, axes, reference_direction, params, previous_angles
        )

        singular_axes = _detect_singular_axes(
            raw_point.position, platform_origin, axes, angles, params.singularity_radius_epsilon
        )
        angles = _inherit_angle_at_singularity(singular_axes, angles, previous_angles)

        resolved_angles.append(angles)
        residual_errors.append(residual_error_rad)
        previous_angles = angles

    return resolved_angles, residual_errors


def _assemble_machine_toolpath(
    raw_points: "_RawPoints", smoothed_angles, residual_errors, params: Step8Params
) -> tuple[ToolpathMachine, list[Step8Diagnostic]]:
    """ToolpathPoint + diagnostics d'orientation inatteignable a partir des angles resolus."""
    diagnostics: list[Step8Diagnostic] = []
    positions = raw_points.position
    extrusion = raw_points.extrusion.tolist()
    feedrate = raw_points.feedrate.tolist()
    colors = raw_points.color_id
    points = [
        ToolpathPoint(
            position=positions[i],
            axis_angles=smoothed_angles[i],
            extrusion=extrusion[i],
            feedrate=feedrate[i],
            color_id=colors[i],
        )
        for i in range(len(raw_points))
    ]
    # conversion en degres VECTORISEE : 364 192 appels numpy evites sur death @ 0,8
    residual_deg = np.degrees(np.asarray(residual_errors, dtype=float))
    for i in np.nonzero(residual_deg > params.max_unreachable_angular_error_deg)[0]:
        diagnostics.append(
            Step8Diagnostic(
                point_index=int(i),
                layer_id=int(raw_points.layer_id[i]),
                segment_index=int(raw_points.segment_index[i]),
                residual_angle_error_deg=float(residual_deg[i]),
                position=positions[i],
            )
        )

    return ToolpathMachine(points=points, layer_ids=np.asarray(raw_points.layer_id, dtype=np.int64).copy()), diagnostics


def _coupled_axis_angles(positions, outer_axis, machine_profile: MachineProfile, params: Step8Params):
    """Angle couple C de TOUS les points d'un coup. Sa seule dependance sequentielle est l'heritage sur l'axe
    (singularite) : un remplissage vers l'avant, 0 au debut."""
    rel = positions - np.asarray(machine_profile.platform_origin, dtype=float)
    rel = rel - np.outer(rel @ outer_axis, outer_axis)
    on_axis = np.linalg.norm(rel, axis=1) < params.singularity_radius_epsilon
    e1 = np.asarray(machine_profile.radial_direction, dtype=float)
    e1 = e1 - outer_axis * float(e1 @ outer_axis)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(outer_axis, e1)
    c = np.degrees(np.arctan2(rel @ e2, rel @ e1)) % 360.0
    if on_axis.any():
        c = c.copy()
        last = 0.0
        for i in range(len(c)):
            if on_axis[i]:
                c[i] = last
            last = c[i]
    return c
