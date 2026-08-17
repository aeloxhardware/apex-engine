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

from core.geometry_utils import _rodrigues_rotate, composed_rotation_angles_and_directions
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
    _PrintedGeometryIndex,
    _surface_sample_points,
    run as collision_run,
)

# Un hook recoit une demande de re-segmentation et retourne la nouvelle
# trajectoire complete (step3 -> step8 re-joues par l'appelant), ou None
# s'il ne peut pas satisfaire la demande.
ResegmentationHook = Callable[[ResegmentationRequest], Optional[ToolpathMachine]]


@dataclass
class Step10Params:
    max_adjustment_iterations: int = 3
    # Grille de recherche de la reorientation locale - meme parametrage que
    # la recherche NN de step8 (grille produit sur les plages d'axes).
    reorientation_samples_per_axis: int = 12
    # Seuil de surplomb : ecart angulaire maximal tolere entre la direction
    # de depot d'origine et celle du candidat de reorientation. Un candidat
    # qui elimine la collision mais depasse ce seuil reintroduirait un
    # surplomb - il est rejete, c'est exactement la condition "tout en
    # respectant encore le seuil de surplomb" de la specification.
    max_direction_deviation_deg: float = 45.0
    # Demi-largeur de la plage d'angles signalee comme problematique dans
    # une ResegmentationRequest (autour de l'angle courant de chaque axe).
    excluded_range_half_width_deg: float = 15.0
    # Nombre de voisins de part et d'autre du point modifie sur lesquels
    # valider les limites de vitesse/acceleration angulaire.
    rate_limit_neighbor_margin: int = 1


# ---------------------------------------------------------------------------
# Direction de depot pour un jeu d'angles donne
# ---------------------------------------------------------------------------

def _direction_for_angles(
    axes: list[RotationAxisSpec], angles_deg: list[float], reference_direction: np.ndarray
) -> np.ndarray:
    """
    Direction de depot resultant de `angles_deg`, avec la MEME convention de
    composition que composed_rotation_angles_and_directions / step8 /
    step9 (rotations appliquees en ordre INVERSE de la liste) - jamais une
    convention propre a l'etape 10.
    """
    direction = np.asarray(reference_direction, dtype=float).reshape(1, 3)
    for spec, angle_deg in zip(reversed(axes), reversed(list(angles_deg))):
        axis = np.asarray(spec.axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        direction = _rodrigues_rotate(direction, axis, np.array([np.radians(angle_deg)]))
    return direction[0]


# ---------------------------------------------------------------------------
# Re-verification partielle (question 1 du plan)
# ---------------------------------------------------------------------------

def _collisions_at_points(
    toolpath: ToolpathMachine,
    machine_profile: MachineProfile,
    sim_params: CollisionSimParams,
    point_indices: list[int],
) -> list[CollisionEvent]:
    """
    Re-teste la collision UNIQUEMENT aux indices demandes, en reutilisant
    les primitives de l'etape 9 (jamais une seconde implementation de la
    detection).

    Pourquoi c'est exact et pourquoi la "portion modifiee" se reduit aux
    seuls points reorientes : l'index de matiere imprimee de l'etape 9 est
    construit a partir de `points[i].position` UNIQUEMENT, jamais de
    `axis_angles` (voir _grow_printed_index_to). Une reorientation locale ne
    touche que les angles, donc la geometrie imprimee est strictement
    identique avant/apres, et aucun point en aval ne peut voir son etat de
    collision changer. Seul le point reoriente doit etre re-teste.

    L'index est reconstruit jusqu'au plus grand indice demande (croissance
    strictement monotone, meme invariant que l'etape 9 : jamais en avance
    sur le point teste, cordon en cours de depot exclu).
    """
    nozzle_model = machine_profile.nozzle_model
    if not isinstance(nozzle_model, trimesh.Trimesh) or len(nozzle_model.vertices) == 0:
        return []
    if not point_indices:
        return []

    # toolhead_sample_spacing (echelle de la TETE), jamais nozzle_diameter
    # (echelle du CORDON) - meme distinction que step9 (CLAUDE.md limite
    # #26). Utiliser nozzle_diameter ici reproduirait exactement le meme
    # blocage (17734 points de test par appel au lieu de ~1540, mesure
    # concretement sur aetheris-4.glb), simplement deplace de l'etape 9 a
    # l'etape 10 plutot qu'elimine.
    samples = _surface_sample_points(nozzle_model, sim_params.toolhead_sample_spacing)
    printed_index = _PrintedGeometryIndex(
        machine_profile.nozzle_diameter, sim_params.grid_cell_size_multiplier
    )

    events: list[CollisionEvent] = []
    last_grown_idx = 0
    for point_idx in sorted(set(point_indices)):
        if point_idx >= len(toolpath.points):
            continue
        last_grown_idx = _grow_printed_index_to(
            printed_index, toolpath.points, point_idx, last_grown_idx
        )
        events.extend(
            _evaluate_point(
                point_idx,
                samples,
                toolpath.points[point_idx],
                machine_profile,
                printed_index,
                sim_params.safety_margin,
            )
        )
    return events


# ---------------------------------------------------------------------------
# Niveau 1 - Reorientation locale
# ---------------------------------------------------------------------------

def _violates_rate_limits(
    points: list[ToolpathPoint],
    point_idx: int,
    candidate_angles: list[float],
    machine_profile: MachineProfile,
    neighbor_margin: int,
) -> bool:
    """
    Une reorientation ponctuelle peut creer un saut d'angle que le lissage
    temporel de l'etape 8 avait precisement elimine. On rejette donc tout
    candidat qui violerait `max_angular_velocity` vis-a-vis de ses voisins
    immediats - sans cette verification, l'etape 10 "resoudrait" une
    collision en produisant une trajectoire mecaniquement inexecutable.

    Limite non traitee (assumee) : seule la vitesse est verifiee, pas
    l'acceleration (qui demanderait de raisonner sur trois points
    consecutifs et de re-lisser toute la serie).
    """
    max_velocity = machine_profile.max_angular_velocity
    if max_velocity <= 0:
        return False  # limite non specifiee : aucune contrainte a faire respecter

    max_velocity_deg_s = np.degrees(max_velocity)

    for offset in range(1, max(neighbor_margin, 1) + 1):
        for neighbor_idx in (point_idx - offset, point_idx + offset):
            if neighbor_idx < 0 or neighbor_idx >= len(points):
                continue
            distance = float(
                np.linalg.norm(points[neighbor_idx].position - points[point_idx].position)
            )
            feedrate = max(points[neighbor_idx].feedrate, 1e-6)
            dt = max(distance / feedrate, 1e-9)
            for current, neighbor in zip(candidate_angles, points[neighbor_idx].axis_angles):
                if abs(current - neighbor) / dt > max_velocity_deg_s:
                    return True
    return False


def _with_angles(toolpath: ToolpathMachine, point_idx: int, angles: list[float]) -> ToolpathMachine:
    """Copie de la trajectoire avec les angles d'un seul point remplaces (jamais de mutation en place)."""
    points = list(toolpath.points)
    points[point_idx] = replace(points[point_idx], axis_angles=list(angles))
    return ToolpathMachine(points=points)


@dataclass
class _ReorientationContext:
    """
    Etat partage entre TOUS les points traites par le niveau 1 (reorientation
    locale) d'une meme iteration externe de `run()` - jamais reconstruit par
    point ni par candidat.

    Les echantillons de tete ne dependent que du modele (fixe pour toute la
    boucle). L'index de matiere imprimee ne depend que de `position` (jamais
    de `axis_angles`, voir `_collisions_at_points`) et grandit de facon
    strictement monotone au fil de `problem_indices` traites en ordre
    CROISSANT par `run()` - exactement le meme invariant que la croissance
    de l'etape 9 (`_grow_printed_index_to` dans `step9_collision_sim.run`).

    Sans cet etat partage, `_try_local_reorientation` reconstruisait son
    propre index depuis 0 a CHAQUE point traite (jusqu'a 3588 fois sur la
    trajectoire de test), chacun re-inserant toutes les capsules jusqu'a son
    propre `point_idx` - mesure concretement : le taux d'appel chutait de
    ~2.5 appels/s au debut a 0.39 appels/s au 200e appel (elapsed 512s),
    une degradation croissante avec la position dans la trajectoire,
    incompatible avec le cout lineaire vise.
    """

    toolhead_samples: Optional[np.ndarray]
    printed_index: Optional[_PrintedGeometryIndex]
    last_grown_idx: int = 0

    def grow_to(self, points: list[ToolpathPoint], point_idx: int) -> None:
        if self.printed_index is None:
            return
        self.last_grown_idx = _grow_printed_index_to(
            self.printed_index, points, point_idx, self.last_grown_idx
        )


def _build_reorientation_context(
    toolpath: ToolpathMachine, machine_profile: MachineProfile, sim_params: CollisionSimParams
) -> _ReorientationContext:
    """
    Construit l'etat partage une seule fois par iteration externe de `run()`.
    `last_grown_idx=0` : l'appelant doit parcourir les points en ordre
    CROISSANT et appeler `grow_to` avant chaque test, jamais en desordre -
    meme invariant "jamais en avance sur le point teste" que l'etape 9.
    """
    nozzle_model = machine_profile.nozzle_model
    has_toolhead = isinstance(nozzle_model, trimesh.Trimesh) and len(nozzle_model.vertices) > 0
    if not has_toolhead:
        return _ReorientationContext(toolhead_samples=None, printed_index=None)

    # toolhead_sample_spacing (echelle de la TETE), jamais nozzle_diameter
    # (echelle du CORDON) - meme distinction que step9 (CLAUDE.md limite #26).
    toolhead_samples = _surface_sample_points(nozzle_model, sim_params.toolhead_sample_spacing)
    printed_index = _PrintedGeometryIndex(
        machine_profile.nozzle_diameter, sim_params.grid_cell_size_multiplier
    )
    return _ReorientationContext(toolhead_samples=toolhead_samples, printed_index=printed_index)


def _try_local_reorientation(
    toolpath: ToolpathMachine,
    point_idx: int,
    machine_profile: MachineProfile,
    sim_params: CollisionSimParams,
    params: Step10Params,
    context: _ReorientationContext,
) -> Optional[list[float]]:
    """
    Cherche un jeu d'angles voisin qui elimine la collision au point
    `point_idx` sans depasser le seuil de surplomb.

    Reutilise la grille de recherche de l'etape 8
    (`composed_rotation_angles_and_directions`) plutot qu'une nouvelle
    exploration : les candidats sont exactement les orientations
    atteignables par la cinematique decrite, pour un nombre d'axes
    quelconque. Ils sont essayes par ecart angulaire CROISSANT vis-a-vis de
    la direction d'origine (le plus petit changement d'abord), et la
    recherche s'arrete des que l'ecart depasse `max_direction_deviation_deg`
    - au-dela, le candidat reintroduirait un surplomb.

    Retourne les angles retenus, ou None si aucun candidat ne convient.

    Performance : `context` (echantillons de tete + index de matiere
    imprimee) est fourni par l'appelant, construit UNE SEULE FOIS par
    iteration externe de `run()` et grossi de facon incrementale au fil des
    points traites (voir `_ReorientationContext`) - jamais reconstruit ici,
    ni par point ni par candidat.
    """
    axes = machine_profile.rotation_axes or []
    if not axes:
        return None  # aucune liberte de rotation : la reorientation n'a pas de sens

    point = toolpath.points[point_idx]
    reference_direction = np.array([0.0, 0.0, -1.0])
    original_direction = _direction_for_angles(axes, point.axis_angles, reference_direction)

    directions, angles_grid = composed_rotation_angles_and_directions(
        axes, reference_direction, samples_per_axis=max(params.reorientation_samples_per_axis, 2)
    )

    cosines = np.clip(directions @ original_direction, -1.0, 1.0)
    deviations_deg = np.degrees(np.arccos(cosines))

    has_toolhead = context.printed_index is not None
    if has_toolhead:
        context.grow_to(toolpath.points, point_idx)
    toolhead_samples = context.toolhead_samples
    printed_index = context.printed_index

    for candidate_idx in np.argsort(deviations_deg):
        if deviations_deg[candidate_idx] > params.max_direction_deviation_deg:
            break  # tries par ecart croissant : tout le reste depasse aussi le seuil

        candidate = [float(a) for a in angles_grid[candidate_idx]]
        if np.allclose(candidate, point.axis_angles, atol=1e-9):
            continue  # angles courants : c'est l'etat qui collisionne deja

        if _violates_rate_limits(
            toolpath.points, point_idx, candidate, machine_profile, params.rate_limit_neighbor_margin
        ):
            continue

        if not has_toolhead:
            return candidate  # aucun modele de tete : rien a re-verifier, meme garde que _collisions_at_points

        candidate_point = replace(point, axis_angles=list(candidate))
        events = _evaluate_point(
            point_idx, toolhead_samples, candidate_point, machine_profile, printed_index, sim_params.safety_margin
        )
        if not events:
            return candidate

    return None


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
    diagnostic_indices = sorted(
        {diagnostic.point_index for diagnostic in (step8_diagnostics or [])}
    )

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

        for point_index in problem_indices:
            is_diagnostic_only = point_index in diagnostic_indices and not any(
                event.point_index == point_index for event in pending_events
            )

            # --- Niveau 1 : reorientation locale ------------------------------
            if not is_diagnostic_only:
                new_angles = _try_local_reorientation(
                    current, point_index, machine_profile, sim_params, params, reorientation_context
                )
                if new_angles is not None:
                    current = _with_angles(current, point_index, new_angles)
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
