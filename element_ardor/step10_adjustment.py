"""
element_ardor/step10_adjustment.py

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
from element_ardor.step9_collision_sim import (
    CollisionSimParams,
    _evaluate_point,
    _grow_printed_index_to,
    printed_margin_of,
    _filter_nozzle_contact_zone,
    _PrintedGeometryIndex,
    _surface_sample_points,
    run as collision_run,
)

from element_ardor.step10_params import (  # noqa: F401  (réexport : imports existants)
    ResegmentationHook,
    Step10Params,
)
from element_ardor.step10_reorient import (  # noqa: F401  (réexport : imports existants)
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

    `toolpath` (sortie de l'etape 8) n'est jamais mutee. `sim_params` est reutilise tel quel pour chaque
    re-verification. `step8_diagnostics` (orientations inatteignables) saute le niveau 1 : la recherche IK
    locale est deja epuisee pour ces points. `resegmentation_hook` absent => niveau 2 saute.

    Retourne (trajectoire ajustee, AdjustmentReport). Ne leve JAMAIS pour une collision non resolue : a
    l'epuisement de `max_adjustment_iterations` les collisions restantes basculent en support de secours et
    `converged` passe a False (signal explicite, comme `Step8Diagnostic`).
    """
    adjuster = _Adjuster(toolpath, collision_report, machine_profile, sim_params, params, resegmentation_hook)
    adjuster.set_diagnostics(step8_diagnostics)
    for iteration in range(1, max(params.max_adjustment_iterations, 0) + 1):
        if adjuster.run_iteration(iteration):
            break
    return adjuster.finish()


class _Adjuster:
    """Etat de la boucle d'ajustement de `run` : trajectoire courante, journal, points deja bascules en
    support de secours, indices a traiter et evenements de collision en attente."""

    def __init__(self, toolpath: ToolpathMachine, collision_report: CollisionReport,
                 machine_profile: MachineProfile, sim_params: CollisionSimParams, params: Step10Params,
                 resegmentation_hook: Optional[ResegmentationHook]):
        self.machine_profile, self.sim_params, self.params = machine_profile, sim_params, params
        self.hook = resegmentation_hook
        self.axes = machine_profile.rotation_axes or []
        self.current = ToolpathMachine(points=list(toolpath.points))
        self.entries: list[AdjustmentLogEntry] = []
        self.fallback_marked: set[int] = set()
        self.tracked_indices = sorted({event.point_index for event in collision_report.events})
        self.diagnostic_indices: list[int] = []
        self.pending_events = list(collision_report.events)
        self.iterations_used = 0
        self.converged = False

    def set_diagnostics(self, step8_diagnostics: Optional[list[Step8Diagnostic]]) -> None:
        """Ecart d'orientation sous le seuil : journalise, jamais un support
        (voir Step10Params.orientation_support_min_residual_deg)."""
        diagnostics = list(step8_diagnostics or [])
        threshold = self.params.orientation_support_min_residual_deg
        if threshold is None:
            threshold = float("inf")
        if diagnostics:
            logged: set[int] = set()
            for diagnostic in diagnostics:
                if diagnostic.residual_angle_error_deg <= threshold and diagnostic.point_index not in logged:
                    logged.add(diagnostic.point_index)
                    self.entries.append(
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
        self.diagnostic_indices = sorted({diagnostic.point_index for diagnostic in diagnostics})

    def handle_point(self, point_index: int, is_diagnostic_only: bool, context: _ReorientationContext) -> bool:
        """Niveaux 1 a 3 pour un point. Retourne True si une re-segmentation a remplace la trajectoire
        (la boucle des points s'arrete, une passe complete est necessaire)."""
        # --- Niveau 1 : reorientation locale ----------------------------------------------
        if not is_diagnostic_only:
            new_angles = _try_local_reorientation(
                self.current, point_index, self.machine_profile, self.sim_params, self.params, context
            )
            if new_angles is not None:
                self.current = _with_angles(self.current, point_index, new_angles)
                if context.angles is not None:
                    context.angles[point_index] = new_angles
                self.entries.append(
                    AdjustmentLogEntry(
                        point_index=point_index,
                        method="local_reorientation",
                        detail=f"angles ajustes vers {[round(a, 3) for a in new_angles]}",
                    )
                )
                return False

        # --- Niveau 2 : re-segmentation ciblee ---------------------------------------------
        if self.hook is not None:
            request = ResegmentationRequest(
                point_index=point_index,
                position=self.current.points[point_index].position,
                excluded_axis_ranges_deg=_excluded_ranges_for_point(
                    self.current.points[point_index], self.axes, self.params.excluded_range_half_width_deg
                ),
                reason=(
                    "orientation inatteignable (etape 8)"
                    if is_diagnostic_only
                    else "collision non resolue par reorientation locale"
                ),
            )
            new_toolpath = self.hook(request)
            if new_toolpath is not None:
                self.current = new_toolpath
                self.entries.append(
                    AdjustmentLogEntry(point_index=point_index, method="region_resegmentation",
                                       detail=request.reason)
                )
                # Une re-segmentation invalide TOUTE la trajectoire (champ scalaire re-resolu globalement) :
                # la re-verification partielle ne s'applique plus, il faut une passe complete.
                self.diagnostic_indices = []
                return True

        # --- Niveau 3 : support de secours -------------------------------------------------
        self.fallback_marked.add(point_index)
        self.entries.append(
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

    def run_iteration(self, iteration: int) -> bool:
        """Une iteration externe. Retourne True quand la boucle doit s'arreter (convergence)."""
        self.iterations_used = iteration
        problem_indices = [
            index
            for index in sorted({event.point_index for event in self.pending_events} | set(self.diagnostic_indices))
            if index not in self.fallback_marked and index < len(self.current.points)
        ]
        if not problem_indices:
            self.converged = True
            return True

        # Construit une seule fois par iteration externe, grossi ensuite de facon monotone au fil de
        # `problem_indices` (deja trie CROISSANT) - voir _ReorientationContext.
        context = _build_reorientation_context(self.current, self.machine_profile, self.sim_params)

        # Ensembles, pas listes : death compte ~115 000 diagnostics d'etape 8, et un `in` sur la liste triee
        # rendait cette boucle quadratique (73 s sur 81 s de step10, 2026-10-06). Meme resultat.
        diagnostic_set = set(self.diagnostic_indices)
        pending_event_indices = {event.point_index for event in self.pending_events}
        # R-017 : le nombre d'iterations n'est pas connu d'avance. La premiere occupe [0, 0,8] de la barre,
        # chaque suivante la moitie de ce qui reste : la fraction est croissante.
        lo = 0.0 if iteration == 1 else 1.0 - 0.2 * 0.5 ** (iteration - 2)
        hi = 0.8 if iteration == 1 else 1.0 - 0.2 * 0.5 ** (iteration - 1)
        full_recheck_needed = False
        for done, point_index in enumerate(problem_indices):
            _progress.report(done, len(problem_indices), span=(lo, hi))
            is_diagnostic_only = point_index in diagnostic_set and point_index not in pending_event_indices
            if self.handle_point(point_index, is_diagnostic_only, context):
                full_recheck_needed = True
                break
        self._recheck(full_recheck_needed)
        if not self.pending_events and not [
            index for index in self.diagnostic_indices if index not in self.fallback_marked
        ]:
            self.converged = True
            return True
        return False

    def _recheck(self, full_recheck_needed: bool) -> None:
        if full_recheck_needed:
            self.pending_events = [
                event
                for event in collision_run(self.current, self.machine_profile, self.sim_params).events
                if event.point_index not in self.fallback_marked
            ]
            self.tracked_indices = sorted({event.point_index for event in self.pending_events})
        else:
            self.pending_events = [
                event
                for event in _collisions_at_points(
                    self.current, self.machine_profile, self.sim_params, self.tracked_indices
                )
                if event.point_index not in self.fallback_marked
            ]

    def finish(self) -> tuple[ToolpathMachine, AdjustmentReport]:
        """Epuisement des iterations avec des collisions encore ouvertes : bascule forcee en support de
        secours plutot qu'une exception (voir docstring de `run`)."""
        unresolved = [event for event in self.pending_events if event.point_index not in self.fallback_marked]
        for event in unresolved:
            self.fallback_marked.add(event.point_index)
            self.entries.append(
                AdjustmentLogEntry(
                    point_index=event.point_index,
                    method="fallback_support",
                    detail=(
                        f"non resolu apres {self.iterations_used} iteration(s) "
                        f"(max_adjustment_iterations atteint) - support classique requis"
                    ),
                )
            )
        if unresolved:
            self.converged = False
        return self.current, AdjustmentReport(
            entries=self.entries,
            converged=self.converged,
            unresolved_collisions=unresolved,
            iterations_used=self.iterations_used,
        )

from element_ardor.step10_reorient import (  # noqa: E402,F401  (reexport : imports et patchs existants)
    _head_zlow_table,
    _clamp_apply_native,
    _TierSetup,
    _reorient_tiers_native,
    _reorient_tiers_python,
    _native_coupled_reorientation,
)
