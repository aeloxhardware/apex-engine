"""
pipeline/step10_loop.py - boucle d'ajustement de l'etape 10 (decoupee de `step10_adjustment.run`, R-030 point 1).

`step10_adjustment.run` garde le contrat public (signature, docstring) et delegue ici ; les primitives
(reorientation locale, re-verification, contexte) restent dans `step10_adjustment`. Ce module n'est importe
que par `run`, donc sans cycle a l'import.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from core import progress as _progress
from core.mesh_types import (
    AdjustmentLogEntry,
    AdjustmentReport,
    CollisionReport,
    MachineProfile,
    ResegmentationRequest,
    Step8Diagnostic,
    ToolpathMachine,
)
from pipeline.step9_collision_sim import CollisionSimParams, run as collision_run
from pipeline.step10_adjustment import (
    ResegmentationHook,
    Step10Params,
    _build_reorientation_context,
    _collisions_at_points,
    _excluded_ranges_for_point,
    _try_local_reorientation,
    _with_angles,
)


@dataclass
class _Ctx:
    """Configuration fixe de la boucle d'ajustement."""
    machine_profile: MachineProfile
    sim_params: CollisionSimParams
    params: Step10Params
    resegmentation_hook: Optional[ResegmentationHook]
    axes: list


@dataclass
class _LoopState:
    """Etat mutable de la boucle d'ajustement (trajectoire courante et journal)."""
    current: ToolpathMachine
    entries: list
    fallback_marked: set
    diagnostic_indices: list
    full_recheck_needed: bool = False


def _log_orientation_deviations(diagnostics: list, min_residual_deg: Optional[float], entries: list) -> list:
    """Ecart d'orientation sous le seuil : journalise, jamais un support (voir
    Step10Params.orientation_support_min_residual_deg). Rend les diagnostics restants."""
    threshold = float("inf") if min_residual_deg is None else min_residual_deg
    if not diagnostics:
        return diagnostics
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
    return [d for d in diagnostics if d.residual_angle_error_deg > threshold]


def _adjust_point(state: _LoopState, ctx: _Ctx, context, point_index: int, is_diagnostic_only: bool) -> bool:
    """Niveaux 1 a 3 pour un point. Rend True quand la boucle sur les points doit s'arreter
    (re-segmentation reussie : toute la trajectoire est invalidee)."""
    # --- Niveau 1 : reorientation locale ------------------------------
    if not is_diagnostic_only:
        new_angles = _try_local_reorientation(
            state.current, point_index, ctx.machine_profile, ctx.sim_params, ctx.params, context
        )
        if new_angles is not None:
            state.current = _with_angles(state.current, point_index, new_angles)
            if context.angles is not None:
                context.angles[point_index] = new_angles
            state.entries.append(
                AdjustmentLogEntry(
                    point_index=point_index,
                    method="local_reorientation",
                    detail=f"angles ajustes vers {[round(a, 3) for a in new_angles]}",
                )
            )
            return False

    # --- Niveau 2 : re-segmentation ciblee ----------------------------
    if ctx.resegmentation_hook is not None:
        request = ResegmentationRequest(
            point_index=point_index,
            position=state.current.points[point_index].position,
            excluded_axis_ranges_deg=_excluded_ranges_for_point(
                state.current.points[point_index], ctx.axes, ctx.params.excluded_range_half_width_deg
            ),
            reason=(
                "orientation inatteignable (etape 8)"
                if is_diagnostic_only
                else "collision non resolue par reorientation locale"
            ),
        )
        new_toolpath = ctx.resegmentation_hook(request)
        if new_toolpath is not None:
            state.current = new_toolpath
            state.entries.append(
                AdjustmentLogEntry(point_index=point_index, method="region_resegmentation", detail=request.reason)
            )
            # Une re-segmentation invalide TOUTE la trajectoire (champ scalaire re-resolu globalement) :
            # la re-verification partielle ne s'applique plus, il faut une passe complete.
            state.full_recheck_needed = True
            state.diagnostic_indices = []
            return True

    # --- Niveau 3 : support de secours --------------------------------
    state.fallback_marked.add(point_index)
    state.entries.append(
        AdjustmentLogEntry(
            point_index=point_index,
            method="fallback_support",
            detail=(
                "ni reorientation locale ni re-segmentation n'ont resolu "
                "la collision - support classique requis"
            ),
        )
    )
    return False


def _recheck(state: _LoopState, ctx: _Ctx, tracked_indices: list):
    """Re-verification apres une passe : complete apres une re-segmentation, sinon limitee aux points
    suivis. -> (pending_events, tracked_indices)."""
    if state.full_recheck_needed:
        pending = [
            event
            for event in collision_run(state.current, ctx.machine_profile, ctx.sim_params).events
            if event.point_index not in state.fallback_marked
        ]
        return pending, sorted({event.point_index for event in pending})
    pending = [
        event
        for event in _collisions_at_points(state.current, ctx.machine_profile, ctx.sim_params, tracked_indices)
        if event.point_index not in state.fallback_marked
    ]
    return pending, tracked_indices


def _force_fallback(state: _LoopState, pending_events: list, iterations_used: int) -> list:
    """Epuisement des iterations avec des collisions encore ouvertes : bascule forcee en support de secours
    plutot qu'une exception (voir la docstring de `run`). Rend les collisions non resolues."""
    unresolved = [event for event in pending_events if event.point_index not in state.fallback_marked]
    for event in unresolved:
        state.fallback_marked.add(event.point_index)
        state.entries.append(
            AdjustmentLogEntry(
                point_index=event.point_index,
                method="fallback_support",
                detail=(
                    f"non resolu apres {iterations_used} iteration(s) "
                    f"(max_adjustment_iterations atteint) - support classique requis"
                ),
            )
        )
    return unresolved



def run_loop(
    toolpath: ToolpathMachine,
    collision_report: CollisionReport,
    machine_profile: MachineProfile,
    sim_params: CollisionSimParams,
    params: Step10Params,
    step8_diagnostics: Optional[list[Step8Diagnostic]] = None,
    resegmentation_hook: Optional[ResegmentationHook] = None,
) -> tuple[ToolpathMachine, AdjustmentReport]:
    """Boucle d'ajustement (voir `step10_adjustment.run` pour le contrat)."""
    state = _LoopState(
        current=ToolpathMachine(points=list(toolpath.points)), entries=[], fallback_marked=set(), diagnostic_indices=[]
    )
    ctx = _Ctx(machine_profile, sim_params, params, resegmentation_hook, machine_profile.rotation_axes or [])

    tracked_indices = sorted({event.point_index for event in collision_report.events})
    diagnostics = _log_orientation_deviations(
        list(step8_diagnostics or []), params.orientation_support_min_residual_deg, state.entries
    )
    state.diagnostic_indices = sorted({diagnostic.point_index for diagnostic in diagnostics})

    pending_events = list(collision_report.events)
    iterations_used = 0
    converged = False

    for iteration in range(1, max(params.max_adjustment_iterations, 0) + 1):
        iterations_used = iteration

        problem_indices = [
            index
            for index in sorted({event.point_index for event in pending_events} | set(state.diagnostic_indices))
            if index not in state.fallback_marked and index < len(state.current.points)
        ]
        if not problem_indices:
            converged = True
            break

        state.full_recheck_needed = False

        # Construit une seule fois par iteration externe, grossi ensuite de facon monotone au fil de
        # `problem_indices` (deja trie CROISSANT ci-dessus) - voir _ReorientationContext.
        reorientation_context = _build_reorientation_context(state.current, machine_profile, sim_params)

        # Ensembles, pas listes : avec la normale par point, death compte ~115 000 diagnostics d'etape 8, et un
        # `in` sur la liste triee rendait cette boucle quadratique (73 s sur 81 s de step10, 2026-10-06).
        diagnostic_set = set(state.diagnostic_indices)
        pending_event_indices = {event.point_index for event in pending_events}
        # R-017 : le nombre d'iterations n'est pas connu d'avance. La premiere occupe [0, 0,8] de la barre,
        # chaque suivante la moitie de ce qui reste : la fraction est croissante et ne repose que sur des
        # points deja traites.
        lo = 0.0 if iteration == 1 else 1.0 - 0.2 * 0.5 ** (iteration - 2)
        hi = 0.8 if iteration == 1 else 1.0 - 0.2 * 0.5 ** (iteration - 1)
        for done, point_index in enumerate(problem_indices):
            _progress.report(done, len(problem_indices), span=(lo, hi))
            is_diagnostic_only = point_index in diagnostic_set and point_index not in pending_event_indices
            if _adjust_point(state, ctx, reorientation_context, point_index, is_diagnostic_only):
                break

        pending_events, tracked_indices = _recheck(state, ctx, tracked_indices)

        if not pending_events and not [
            index for index in state.diagnostic_indices if index not in state.fallback_marked
        ]:
            converged = True
            break

    unresolved = _force_fallback(state, pending_events, iterations_used)
    if unresolved:
        converged = False

    return state.current, AdjustmentReport(
        entries=state.entries,
        converged=converged,
        unresolved_collisions=unresolved,
        iterations_used=iterations_used,
    )
