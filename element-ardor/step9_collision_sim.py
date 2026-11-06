"""
pipeline/step9_collision_sim.py

Etape 9 du pipeline : simulation de collision de la tete d'impression
(MachineProfile.nozzle_model) le long de ToolpathMachine (etape 8), contre
la piece deja imprimee (etat INCREMENTAL, jamais le mesh final complet),
le plateau a son angle courant, et le mandrin/mors (non represente
aujourd'hui - signal explicite plutot que test silencieusement absent).

Ce que ce fichier fait, dans l'ordre :
  1. Echantillonne ToolpathMachine.points selon `path_sampling_step`.
  2. Pour chaque point echantillonne, positionne une copie du modele de
     tete (rotation N-axes generique + translation - meme convention de
     composition que step3/step7/step8, jamais recodee).
  3. Fait grandir un index spatial de la matiere deja imprimee AU FUR ET A
     MESURE (jamais en avance sur le point teste), sous forme de capsules
     (segment + rayon), pas d'union booleene (cout prohibitif).
  4. Teste la tete positionnee contre cet index (distance point-capsule),
     et contre le plateau (ProximityQuery) - broad-phase AABB commun aux
     deux avant tout calcul de distance precis.
  5. Enregistre chaque collision : point concerne, type d'obstacle,
     profondeur de penetration.

Depend de ToolpathMachine (etape 8), MachineProfile (nozzle_model deja
charge en mm via io/mesh_loader.py, rotation_axes generique).

Sortie : CollisionReport.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import trimesh
import trimesh.proximity

from core import _cpp
from core import progress as _progress
from core.geometry_utils import _rodrigues_rotate, point_has_support_along_direction
from core.mesh_types import (
    AABB,
    CollisionEvent,
    CollisionReport,
    MachineProfile,
    RotationAxisSpec,
    ToolpathMachine,
    ToolpathPoint,
)

from pipeline.step9_params import (  # noqa: F401  (réexport : imports existants)
    CollisionSimParams,
    printed_margin_of,
)
from pipeline.step9_primitives import (  # noqa: F401  (réexport : imports existants)
    PLATFORM_CENTER_OFFSET,
    PLATFORM_HALF_HEIGHT,
    PLATFORM_RADIUS,
    _EMPTY_ENDPOINTS,
    _PrintedGeometryIndex,
    _PrintedSegment,
    _aabb_overlaps,
    _cell_coords_for_aabb,
    _check_platform_disc,
    _check_printed_part,
    _check_static_obstacle,
    _evaluate_point,
    _filter_nozzle_contact_zone,
    _grow_printed_index_to,
    _platform_mesh_at_angle,
    _point_to_segment_distance,
    _points_to_segment_distances,
    _positioned_toolhead_vertices,
    _positioned_toolhead_vertices_batch,
    _sample_point_indices,
    _segment_aabb,
    _surface_sample_points,
    _toolhead_aabb,
)
from pipeline.step9_fast import (  # noqa: F401  (réexport : imports existants)
    FastCollisionChecker,
    _FAST_CHECKER_CACHE,
    _NATIVE_MIN_TRANCHE,
    _NATIVE_TRANCHES,
    _densified_printed_material,
    _platform_disc_batch,
    _run_fast_kdtree,
    fast_batch_events,
    fast_checker_for,
)


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------


def run(
    toolpath: ToolpathMachine,
    machine_profile: MachineProfile,
    params: CollisionSimParams,
) -> CollisionReport:
    """
    Simule les collisions de la tete d'impression le long de la trajectoire
    machine, contre la piece deja imprimee (etat incremental), le plateau,
    et le mandrin/mors (jamais teste a ce stade - voir `chuck_checked`).

    Args:
        toolpath: sortie de pipeline/step8_machine_toolpath.py.
        machine_profile: cinematique + `nozzle_model` (Trimesh deja charge
            en mm, ex: io.mesh_loader.load_toolhead_mesh - PAS charge ici,
            step9 ne connait pas le systeme de fichiers).
        params: pas d'echantillonnage, marge de securite.

    Returns:
        CollisionReport. `nozzle_model_checked=False` si `machine_profile.nozzle_model`
        est absent/non exploitable (aucun test n'a alors eu lieu, jamais
        d'exception). `chuck_checked` reste toujours `False` : aucune
        geometrie de mandrin/mors n'est disponible dans le projet a ce
        stade - signal explicite, pas un skip silencieux.
    """
    nozzle_model = machine_profile.nozzle_model
    if not isinstance(nozzle_model, trimesh.Trimesh) or len(nozzle_model.vertices) == 0:
        return CollisionReport(events=[], chuck_checked=False, nozzle_model_checked=False)

    # Sommets + points le long des aretes, calcules UNE FOIS en repere local
    # (voir _surface_sample_points) : un obstacle fin peut traverser une face
    # sans passer pres d'aucun sommet, les sommets seuls manqueraient cette
    # collision. Le pas d'echantillonnage suit `toolhead_sample_spacing`
    # (echelle de la TETE, pas `nozzle_diameter` - echelle du CORDON, un
    # concept different : confondre les deux produit des dizaines de
    # milliers de points de test par point de trajectoire sur un modele de
    # tete reel, voir CollisionSimParams.toolhead_sample_spacing).
    toolhead_base_samples = _surface_sample_points(nozzle_model, params.toolhead_sample_spacing)
    toolhead_base_samples = _filter_nozzle_contact_zone(
        toolhead_base_samples, params.nozzle_contact_clearance
    )

    sample_indices = _sample_point_indices(toolpath.points, params.path_sampling_step)
    if params.fast_kdtree:
        events = _run_fast_kdtree(toolpath, machine_profile, params, toolhead_base_samples, sample_indices)
        return CollisionReport(events=events, chuck_checked=False, nozzle_model_checked=True)
    printed_index = _PrintedGeometryIndex(
        machine_profile.nozzle_diameter, params.grid_cell_size_multiplier
    )

    events: list[CollisionEvent] = []
    last_grown_idx = 0
    for done, point_idx in enumerate(sample_indices):
        _progress.report(done, len(sample_indices))   # R-017 : points deja testes
        last_grown_idx = _grow_printed_index_to(
            printed_index, toolpath.points, point_idx, last_grown_idx
        )
        events.extend(
            _evaluate_point(
                point_idx,
                toolhead_base_samples,
                toolpath.points[point_idx],
                machine_profile,
                printed_index,
                params.safety_margin,
                printed_margin_of(params),
            )
        )

    return CollisionReport(events=events, chuck_checked=False, nozzle_model_checked=True)
