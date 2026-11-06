"""
pipeline/step10_adjustment.py

Etape 10 du pipeline : boucle d'ajustement iterative sur les collisions
rapportees par l'etape 9 et les orientations inatteignables rapportees par
l'etape 8.

Pour chaque probleme, tente dans l'ordre :
  1. REORIENTATION LOCALE : cherche un jeu d'angles voisin qui elimine la
     collision tout en restant sous le seuil de surplomb - via la MEME
     machinerie de recherche que step8 (composed_rotation_angles_and_directions
     + refine_axis_angles_local_search), jamais une nouvelle methode.
  2. RE-SEGMENTATION CIBLEE : emet une ResegmentationRequest vers l'appelant
     (hook injecte), avec la plage d'angles a exclure. Voir la docstring de
     ResegmentationRequest pour pourquoi l'etape 10 ne re-joue PAS
     elle-meme step3 -> step8.
  3. SUPPORT DE SECOURS : marque le point, avec log explicite.

Apres chaque ajustement, relance la verification de collision UNIQUEMENT
sur la portion modifiee (voir _collisions_at_points), pas la simulation
complete.

Sortie : (ToolpathMachine ajustee, AdjustmentReport).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable, Optional

import numpy as np
import trimesh

from core import _cpp
from core import progress as _progress
from core.geometry_utils import _rodrigues_rotate_single, composed_rotation_angles_and_directions
from core.mesh_types import (
    AdjustmentLogEntry,
    AdjustmentReport,
    CollisionEvent,
    CollisionReport,
    MachineProfile,
    ResegmentationRequest,
    RotationAxisSpec,
    Step8Diagnostic,
    ToolpathMachine,
    ToolpathPoint,
)
from pipeline.step9_collision_sim import (
    CollisionSimParams,
    _evaluate_point,
    _grow_printed_index_to,
    printed_margin_of,
    _filter_nozzle_contact_zone,
    _PrintedGeometryIndex,
    _surface_sample_points,
    run as collision_run,
)

from pipeline.step10_params import (  # noqa: F401  (réexport : imports existants)
    ResegmentationHook,
    Step10Params,
)
from pipeline.step10_reorient import (  # noqa: F401  (réexport : imports existants)
    _ReorientationContext,
    _build_reorientation_context,
    _collisions_at_points,
    _coupled_candidates,
    _direction_for_angles,
    _fast_checker,
    _fast_events,
    _try_local_reorientation,
    _violates_rate_limits,
    _with_angles,
    clamp_head_height,
    reorient_coupled_tiers,
)


# ---------------------------------------------------------------------------
# Niveau 2 - Demande de re-segmentation ciblee
# ---------------------------------------------------------------------------

def _excluded_ranges_for_point(
    point: ToolpathPoint, axes: list[RotationAxisSpec], half_width_deg: float
) -> list[Optional[tuple[float, float]]]:
    """
    Plage d'angles a exclure par axe : une fenetre centree sur l'angle
    courant, celui qui a mene a la collision. Formule ainsi (et pas comme
    "un angle unique interdit") parce que la segmentation raisonne en
    directions candidates continues, pas en valeurs discretes.
    """
    return [
        (angle - half_width_deg, angle + half_width_deg)
        for angle, _spec in zip(point.axis_angles, axes)
    ]


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def run(
    toolpath: ToolpathMachine,
    collision_report: CollisionReport,
    machine_profile: MachineProfile,
    sim_params: CollisionSimParams,
    params: Step10Params,
    step8_diagnostics: Optional[list[Step8Diagnostic]] = None,
    resegmentation_hook: Optional[ResegmentationHook] = None,
) -> tuple[ToolpathMachine, AdjustmentReport]:
    """
    Boucle d'ajustement sur les collisions (etape 9) et les orientations
    inatteignables (etape 8).

    Args:
        toolpath: sortie de l'etape 8 (jamais mutee - une copie ajustee est
            retournee).
        collision_report: sortie de l'etape 9 sur cette trajectoire.
        machine_profile: cinematique + `nozzle_model` deja charge.
        sim_params: parametres de l'etape 9, reutilises tels quels pour
            chaque re-verification (meme marge, meme echantillonnage).
        params: parametres propres a l'ajustement.
        step8_diagnostics: orientations inatteignables signalees par
            l'etape 8. Elles sautent le niveau 1 : la recherche IK locale
            a deja ete epuisee par l'etape 8 pour ces points, la relancer
            a l'identique ne peut rien donner de nouveau.
        resegmentation_hook: appele quand la reorientation locale echoue.
            Absent => le niveau 2 est saute et l'on passe directement au
            support de secours.

    Returns:
        (trajectoire ajustee, AdjustmentReport). Ne leve JAMAIS pour une
        collision non resolue : a l'epuisement de `max_adjustment_iterations`,
        les collisions restantes basculent en support de secours et
        `converged` passe a False - signal explicite, coherent avec
        `Step8Diagnostic` et `CollisionReport.chuck_checked`.
    """
    current = ToolpathMachine(points=list(toolpath.points))
    entries: list[AdjustmentLogEntry] = []
    fallback_marked: set[int] = set()
    axes = machine_profile.rotation_axes or []

    tracked_indices = sorted({event.point_index for event in collision_report.events})
    diagnostics = list(step8_diagnostics or [])
    threshold = params.orientation_support_min_residual_deg
    if threshold is None:
        threshold = float("inf")
    if diagnostics:
        # Ecart d'orientation sous le seuil : journalise, jamais un support
        # (voir Step10Params.orientation_support_min_residual_deg).
        logged: set[int] = set()
        for diagnostic in diagnostics:
            if diagnostic.residual_angle_error_deg <= threshold and diagnostic.point_index not in logged:
                logged.add(diagnostic.point_index)
                entries.append(
                    AdjustmentLogEntry(
                        point_index=diagnostic.point_index,
                        method="orientation_deviation",
                        detail=(
                            f"ecart d'orientation {diagnostic.residual_angle_error_deg:.1f} deg "
                            f"<= {threshold:g} deg : depot incline, sans support"
                        ),
                    )
                )
        diagnostics = [d for d in diagnostics if d.residual_angle_error_deg > threshold]
    diagnostic_indices = sorted({diagnostic.point_index for diagnostic in diagnostics})

    pending_events = list(collision_report.events)
    iterations_used = 0
    converged = False

    for iteration in range(1, max(params.max_adjustment_iterations, 0) + 1):
        iterations_used = iteration

        problem_indices = [
            index
            for index in sorted({event.point_index for event in pending_events} | set(diagnostic_indices))
            if index not in fallback_marked and index < len(current.points)
        ]
        if not problem_indices:
            converged = True
            break

        full_recheck_needed = False

        # Construit une seule fois par iteration externe, grossi ensuite de
        # facon monotone au fil de `problem_indices` (deja trie CROISSANT
        # ci-dessus) - voir _ReorientationContext.
        reorientation_context = _build_reorientation_context(current, machine_profile, sim_params)

        # Ensembles, pas listes : avec la normale par point, death compte ~115 000
        # diagnostics d'etape 8, et un `in` sur la liste triee rendait cette boucle
        # quadratique (73 s sur 81 s de step10, 2026-10-06). Meme resultat.
        diagnostic_set = set(diagnostic_indices)
        pending_event_indices = {event.point_index for event in pending_events}
        # R-017 : le nombre d'iterations n'est pas connu d'avance. La premiere occupe [0, 0,8] de la barre,
        # chaque suivante la moitie de ce qui reste : la fraction est croissante et ne repose que sur des
        # points deja traites.
        lo = 0.0 if iteration == 1 else 1.0 - 0.2 * 0.5 ** (iteration - 2)
        hi = 0.8 if iteration == 1 else 1.0 - 0.2 * 0.5 ** (iteration - 1)
        for done, point_index in enumerate(problem_indices):
            _progress.report(done, len(problem_indices), span=(lo, hi))
            is_diagnostic_only = point_index in diagnostic_set and point_index not in pending_event_indices

            # --- Niveau 1 : reorientation locale ------------------------------
            if not is_diagnostic_only:
                new_angles = _try_local_reorientation(
                    current, point_index, machine_profile, sim_params, params, reorientation_context
                )
                if new_angles is not None:
                    current = _with_angles(current, point_index, new_angles)
                    if reorientation_context.angles is not None:
                        reorientation_context.angles[point_index] = new_angles
                    entries.append(
                        AdjustmentLogEntry(
                            point_index=point_index,
                            method="local_reorientation",
                            detail=f"angles ajustes vers {[round(a, 3) for a in new_angles]}",
                        )
                    )
                    continue

            # --- Niveau 2 : re-segmentation ciblee ----------------------------
            if resegmentation_hook is not None:
                request = ResegmentationRequest(
                    point_index=point_index,
                    position=current.points[point_index].position,
                    excluded_axis_ranges_deg=_excluded_ranges_for_point(
                        current.points[point_index], axes, params.excluded_range_half_width_deg
                    ),
                    reason=(
                        "orientation inatteignable (etape 8)"
                        if is_diagnostic_only
                        else "collision non resolue par reorientation locale"
                    ),
                )
                new_toolpath = resegmentation_hook(request)
                if new_toolpath is not None:
                    current = new_toolpath
                    entries.append(
                        AdjustmentLogEntry(
                            point_index=point_index,
                            method="region_resegmentation",
                            detail=request.reason,
                        )
                    )
                    # Une re-segmentation invalide TOUTE la trajectoire (champ
                    # scalaire re-resolu globalement) : la re-verification
                    # partielle ne s'applique plus, il faut une passe complete.
                    full_recheck_needed = True
                    diagnostic_indices = []
                    break

            # --- Niveau 3 : support de secours --------------------------------
            fallback_marked.add(point_index)
            entries.append(
                AdjustmentLogEntry(
                    point_index=point_index,
                    method="fallback_support",
                    detail=(
                        "ni reorientation locale ni re-segmentation n'ont resolu "
                        "la collision - support classique requis"
                    ),
                )
            )

        if full_recheck_needed:
            pending_events = [
                event
                for event in collision_run(current, machine_profile, sim_params).events
                if event.point_index not in fallback_marked
            ]
            tracked_indices = sorted({event.point_index for event in pending_events})
        else:
            pending_events = [
                event
                for event in _collisions_at_points(
                    current, machine_profile, sim_params, tracked_indices
                )
                if event.point_index not in fallback_marked
            ]

        if not pending_events and not [
            index for index in diagnostic_indices if index not in fallback_marked
        ]:
            converged = True
            break

    # Epuisement des iterations avec des collisions encore ouvertes : bascule
    # forcee en support de secours plutot qu'une exception (voir docstring).
    unresolved = [event for event in pending_events if event.point_index not in fallback_marked]
    for event in unresolved:
        fallback_marked.add(event.point_index)
        entries.append(
            AdjustmentLogEntry(
                point_index=event.point_index,
                method="fallback_support",
                detail=(
                    f"non resolu apres {iterations_used} iteration(s) "
                    f"(max_adjustment_iterations atteint) - support classique requis"
                ),
            )
        )
    if unresolved:
        converged = False

    return current, AdjustmentReport(
        entries=entries,
        converged=converged,
        unresolved_collisions=unresolved,
        iterations_used=iterations_used,
    )
