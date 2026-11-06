"""
pipeline/step7_scheduling.py

Etape 7 du pipeline : application de la carte de couleur peinte par
l'utilisateur sur les trajectoires de step6, puis planification de l'ordre
d'impression global.

Depend de Layer (step5/step6), MeshData (step1, kdtree + bvh), ColorMap
(peinture utilisateur), RegionDependencyGraph (step3), ScalarField (step4).

Objectif produit : exploiter la liberte 4/5 axes pour imprimer une couleur
sur PLUSIEURS couches d'affilee (repassage), en laissant des TROUS
TEMPORAIRES la ou une autre couleur doit aller, combles plus tard - plutot
que de purger la buse a chaque changement de couleur a chaque couche.

Ce que ce fichier fait, dans l'ordre :
  1. Projette ColorMap.per_face_color sur chaque segment de trajectoire, via
     le pont kdtree surface/volume deja utilise par step4/step5/step6, mais
     adapte sur deux points (voir _edge_colors_for_polyline) : relevement
     2D -> 3D prealable, et exploitation de la DISTANCE retournee par le
     kdtree (ignoree par step6) pour distinguer matiere visible et matiere
     interne.
  2. Lisse la carte projetee selon ColorMap.merge_tolerance (filtre 1D en
     longueur d'arc le long de la trajectoire).
  3. Scinde les polylignes aux transitions de couleur, pour garantir
     l'invariant "une polyligne = une couleur" dont depend tout
     l'ordonnancement aval.
  4. Construit le graphe (region, color) -> dependances physiques.
  5. Calcule la deadline d'accessibilite de chaque trou candidat (a partir
     de quelle couche la buse ne peut plus l'atteindre).
  6. Arbitre trou temporaire vs purge classique (cout + faisabilite +
     max_open_holes), puis ordonnance en earliest-deadline-first entrelace
     avec l'ordre topologique de base.

Sortie : (list[Layer] recolorees/resegmentees, PrintSchedule).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Optional

import networkx as nx
import numpy as np

import trimesh
import trimesh.proximity

from core import progress as _progress
from core.geometry_utils import (
    composed_rotation_directions,
    fibonacci_sphere_directions,
)
from core.mesh_types import (
    ColorMap,
    Layer,
    MeshData,
    PrintPass,
    PrintSchedule,
    RegionDependencyGraph,
    RotationAxisSpec,
    ScalarField,
    TetMesh,
)
from core.spatial import KDTree
from pipeline.step6_toolpath2d import _lift_to_3d

from pipeline.step7_params import (  # noqa: F401  (réexport : imports existants)
    HoleCostModel,
    PurgeCostModel,
    SchedulingError,
    Step7Params,
)
from pipeline.step7_colors import (  # noqa: F401  (réexport : imports existants)
    _ColorGroup,
    _apply_colors_to_layers,
    _color_runs,
    _dominant_color,
    _edge_colors_for_polyline,
    _edge_lengths,
    _edge_midpoints_2d,
    _group_segments_by_color,
    _smooth_edge_colors,
    _split_polyline_by_color,
    build_color_dependency_graph,
)


    # NOTE: `merge_tolerance` n'est deliberement PAS duplique ici - il vit
    # sur ColorMap (contrat UI), seule source de verite. Deux sources qui
    # divergent en silence est exactement le piege deja rencontre avec
    # `platform_z` a l'etape 4.


# ---------------------------------------------------------------------------
# Question 2 - Deadline d'accessibilite d'un trou temporaire
# ---------------------------------------------------------------------------

def _escape_directions(params: Step7Params) -> np.ndarray:
    """
    Directions d'echappement candidates pour tester l'accessibilite d'un trou.

    Dispatch generique sur la cinematique decrite par `rotation_axes` (meme
    principe que step3._generate_candidate_directions) : ni None/vide
    (orientation libre, sphere complete) ni une liste d'axes quelconque
    n'est un cas special code en dur - composed_rotation_directions absorbe
    uniformement 0, 1 ou N axes. Tester une sphere complete quand la
    cinematique reelle est plus restreinte surestimerait l'accessibilite et
    planifierait des repassages impossibles (voir CLAUDE.md limite #18).
    """
    if not params.rotation_axes:
        return fibonacci_sphere_directions(params.escape_direction_count)

    up = np.array([0.0, 0.0, 1.0])
    return composed_rotation_directions(
        params.rotation_axes,
        reference_direction=up,
        samples_per_axis=_escape_samples_per_axis(
            params.escape_direction_count, len(params.rotation_axes)
        ),
    )


def _escape_samples_per_axis(escape_direction_count: int, axis_count: int) -> int:
    """Meme conversion count->samples_per_axis que step3._samples_per_axis (grille produit)."""
    if axis_count <= 0:
        return escape_direction_count
    return max(int(round(escape_direction_count ** (1.0 / axis_count))), 2)


def _blocking_layer_along_ray(
    origin: np.ndarray,
    direction: np.ndarray,
    hole_layer_index: int,
    context: DeadlineContext,
    samples: int,
) -> int:
    """
    Indice de la PREMIERE couche qui scelle le couloir d'echappement partant
    de `origin` dans la direction `direction`.

    Seule la matiere deposee APRES le trou peut le rendre inaccessible : les
    echantillons d'indice <= `hole_layer_index` sont soit deja presents au
    moment ou le trou est cree (ils ne scellent rien de nouveau), soit le
    trou lui-meme. Ne pas les exclure rendait le mecanisme structurellement
    inatteignable pour les trous de SURFACE - c'est-a-dire exactement le cas
    d'usage, la couleur etant une propriete de surface : le rayon sortant
    lateralement quitte le solide quasi immediatement, tous les echantillons
    retombaient sur la couche du trou, et le report etait toujours refuse.

    Le couloir est scelle des que N'IMPORTE QUELLE section transversale est
    remplie, donc par le MINIMUM des indices de couche futurs rencontres.
    """
    mesh_data, field_tree, field_values = context.mesh_data, context.field_tree, context.field_values
    sorted_isovalues, layer_count = context.sorted_isovalues, context.layer_count
    epsilon = 1e-6
    hits, _index_ray, _index_tri = mesh_data.bvh.intersects_location(
        origin.reshape(1, 3), direction.reshape(1, 3)
    )
    if len(hits) == 0:
        return layer_count  # aucune matiere devant : jamais scelle

    distances = np.linalg.norm(hits - origin, axis=1)
    forward = distances[distances > epsilon]
    if len(forward) == 0:
        return layer_count

    exit_distance = float(forward.min())
    sample_distances = np.linspace(epsilon, exit_distance, samples)
    sample_points = origin[None, :] + direction[None, :] * sample_distances[:, None]

    _, nearest_nodes = field_tree.query(sample_points)
    phi = field_values[np.asarray(nearest_nodes, dtype=np.int64)]
    layer_indices = np.searchsorted(sorted_isovalues, phi)

    future = layer_indices[layer_indices > hole_layer_index]
    if len(future) == 0:
        return layer_count  # rien de futur sur ce couloir : jamais scelle

    return int(future.min())


def compute_hole_deadline(
    position_3d: np.ndarray,
    hole_layer_index: int,
    context: DeadlineContext,
    directions: np.ndarray,
    params: Step7Params,
) -> int:
    """
    Derniere couche a partir de laquelle la buse peut encore atteindre le
    trou situe en `position_3d`.

    Un trou reste comblable tant qu'il existe au moins UN couloir
    d'echappement libre vers l'exterieur ; il meurt quand la matiere imprimee
    entre-temps les a tous scelles. La deadline est donc le MAXIMUM sur les
    directions (la meilleure route survit le plus longtemps), moins 1.

    Deux approximations assumees :
      - `phi` en un point 3D quelconque est lu au noeud du TetMesh le plus
        proche (kdtree), plutot que par localisation tetraedrique +
        interpolation barycentrique : machinerie neuve pour un gain nul, une
        deadline d'ordonnancement etant une heuristique par nature.
      - Le test RIGOUREUX d'accessibilite, avec le vrai modele 3D de la tete
        (buse + porte-buse + carriage) et l'etat incremental de la piece, est
        explicitement le travail de l'etape 9, sur laquelle l'etape 10 boucle.
        L'etape 7 produit un ordre plausible ; l'etape 9 le valide. Meme
        statut que le graphe de dependances de step3, lui aussi heuristique.
    """
    best_blocking = -1
    for direction in directions:
        blocking = _blocking_layer_along_ray(
            position_3d,
            direction,
            hole_layer_index,
            context,
            params.deadline_ray_samples,
        )
        best_blocking = max(best_blocking, blocking)

    deadline = best_blocking - 1
    return int(np.clip(deadline, hole_layer_index, context.layer_count - 1))


def _group_position_3d(layer: Layer, group: _ColorGroup) -> np.ndarray:
    """Position representative d'un groupe : centroide de ses segments, en 3D."""
    toolpath = layer.toolpath_2d or []
    points = [toolpath[i] for i in group.segment_indices if i < len(toolpath)]
    if not points:
        origin, _u, _v = layer.local_frame
        return np.asarray(origin, dtype=float)
    # centroide sur autant de colonnes que porte toolpath_2d (2 = u,v ;
    # 3 = u,v,offset hors-plan - cf. step6_toolpath2d._to_toolpath_3col) :
    # reshape(1, -1) reste correct dans les deux cas, contrairement a un
    # reshape(1, 2) fige qui plantait des que la 3e colonne a ete ajoutee.
    centroid = np.concatenate(points, axis=0).mean(axis=0)
    return _lift_to_3d(centroid.reshape(1, -1), layer.local_frame)[0]


# ---------------------------------------------------------------------------
# Decision 6 - Arbitrage trou temporaire vs purge classique
# ---------------------------------------------------------------------------

def _hole_cost(open_layer_span: int, params: Step7Params) -> float:
    model = params.hole_cost_model
    return model.reheat_penalty + model.adhesion_risk_weight * float(max(open_layer_span, 0))


def _purge_cost(params: Step7Params) -> float:
    model = params.purge_cost_model
    return model.time_per_purge + model.volume_per_purge


def _should_defer(
    open_layer_span: int, deadline: int, fill_layer_index: int, open_hole_count: int, params: Step7Params
) -> bool:
    """
    Un trou temporaire n'est retenu que s'il est a la fois REALISABLE et
    RENTABLE. Bascule en purge classique sinon :

      1. `deadline < fill_layer_index` : le trou serait deja inaccessible au
         moment du repassage - c'est la raison d'etre du calcul de deadline,
         sans lui le systeme planifierait des repassages impossibles.
      2. `max_open_holes` deja atteint.
      3. `hole_cost > purge_cost` : pas rentable.
    """
    if deadline < fill_layer_index:
        return False
    if open_hole_count >= params.max_open_holes:
        return False
    return _hole_cost(open_layer_span, params) <= _purge_cost(params)


# ---------------------------------------------------------------------------
# Decision 7 - Ordonnancement earliest-deadline-first entrelace
# ---------------------------------------------------------------------------

@dataclass
class _PendingHole:
    group: _ColorGroup
    deadline: int


def _make_pass(
    groups: list[_ColorGroup], color_id: Optional[int], is_hole_fill: bool
) -> Optional[PrintPass]:
    segments = [
        (group.layer_id, segment_index)
        for group in groups
        for segment_index in group.segment_indices
    ]
    if not segments:
        return None

    layer_ids = [layer_id for layer_id, _ in segments]
    regions = sorted({r for group in groups for r in group.region_ids})

    return PrintPass(
        segments=segments,
        color_id=color_id,
        is_hole_fill=is_hole_fill,
        layer_range=(min(layer_ids), max(layer_ids)),
        deferred_regions=regions if is_hole_fill else [],
    )


def _schedule_passes(
    groups_per_layer: list[list[_ColorGroup]],
    deadline_of: "callable",
    params: Step7Params,
) -> PrintSchedule:
    """
    Parcourt les couches dans l'ordre topologique de base (l'ordre des Layer de step5 est deja valide), en
    maintenant une couleur active aussi longtemps que possible. Les couleurs non actives d'une couche sont
    soit differees en trou temporaire, soit imprimees immediatement au prix d'une purge (Decision 6). Les
    trous en attente sont combles en earliest-deadline-first : des qu'un trou approche sa deadline, ou que
    `max_open_holes` est atteint, la passe de comblement est inseree.
    """
    scheduler = _PassScheduler(deadline_of, params)
    for layer_index, groups in enumerate(groups_per_layer):
        scheduler.schedule_layer(layer_index, groups)
    scheduler.flush_holes(len(groups_per_layer), force_all=True)
    return PrintSchedule(passes=scheduler.passes)


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def _native_schedule(layers, mesh_data, tet_mesh, color_map, scalar_field, params):
    """Pas 21c du portage (cpp/src/core/scheduling.cpp) : toute l'etape en
    natif. Tolerances aux seules frontieres (distance a la surface, sommes de
    longueurs : voir apex/scheduling.hpp). None si le noyau natif est inactif
    ou si les polylignes n'ont pas toutes le meme nombre de colonnes."""
    from core import _cpp
    if not _cpp.enabled():
        return None
    module = _cpp.module()
    if module is None or not hasattr(module, "schedule_layers"):
        return None
    packed = []
    for layer in layers:
        polys = [np.ascontiguousarray(p, dtype=float) for p in (layer.toolpath_2d or [])]
        columns = {p.shape[1] for p in polys if p.ndim == 2}
        if len(columns) > 1 or any(p.ndim != 2 for p in polys):
            return None
        origin, axis_u, axis_v = layer.local_frame
        frame = np.concatenate([np.asarray(origin, float), np.asarray(axis_u, float), np.asarray(axis_v, float)])
        packed.append((int(layer.id), frame, polys, columns.pop() if columns else 2,
                       [int(r) for r in layer.region_ids],
                       None if layer.isovalue is None else float(layer.isovalue)))
    face_count = mesh_data.face_count
    face_color = np.zeros(0, dtype=np.int64)
    painted = np.zeros(0, dtype=np.uint8)
    merge_tolerance = 0.0
    if color_map is not None and color_map.per_face_color:
        face_color = np.zeros(face_count, dtype=np.int64)
        painted = np.zeros(face_count, dtype=np.uint8)
        for face, color in color_map.per_face_color.items():
            if 0 <= int(face) < face_count and color is not None:
                face_color[int(face)] = int(color)
                painted[int(face)] = 1
        merge_tolerance = float(color_map.merge_tolerance)
    out_layers, passes = module.schedule_layers(
        packed,
        np.ascontiguousarray(mesh_data.vertices, dtype=float),
        np.ascontiguousarray(mesh_data.faces, dtype=np.int64),
        np.ascontiguousarray(mesh_data.face_normals, dtype=float),
        face_color, painted, merge_tolerance,
        np.ascontiguousarray(tet_mesh.tet_vertices, dtype=float),
        np.ascontiguousarray(scalar_field.values, dtype=float),
        np.ascontiguousarray(_escape_directions(params), dtype=float),
        float(params.surface_color_distance),
        float(params.hole_cost_model.reheat_penalty), float(params.hole_cost_model.adhesion_risk_weight),
        float(params.purge_cost_model.volume_per_purge), float(params.purge_cost_model.time_per_purge),
        int(params.max_open_holes), int(params.deadline_ray_samples),
    )
    colored_layers = [
        replace(layer, toolpath_2d=list(polys), segment_colors=list(colors), color_id=dominant)
        for layer, (polys, colors, dominant) in zip(layers, out_layers)
    ]
    schedule = PrintSchedule(passes=[
        PrintPass(segments=[tuple(seg) for seg in segs], color_id=color, is_hole_fill=bool(fill),
                  layer_range=tuple(rng), deferred_regions=list(deferred))
        for segs, color, fill, rng, deferred in passes
    ])
    return colored_layers, schedule


def run(
    layers: list[Layer],
    mesh_data: MeshData,
    tet_mesh: TetMesh,
    color_map: Optional[ColorMap],
    region_dependency_graph: RegionDependencyGraph,
    scalar_field: ScalarField,
    params: Step7Params,
) -> tuple[list[Layer], PrintSchedule]:
    """
    Applique la peinture utilisateur aux trajectoires et produit la sequence d'impression globale.

    Args:
        layers: sortie de step6 (toolpath_2d rempli).
        mesh_data: sortie de step1 (kdtree pour la projection de couleur, bvh pour les rayons d'accessibilite).
        tet_mesh: sortie de step1. Requis (et pas seulement `scalar_field`) : `ScalarField.values` est indexe
            par NOEUD du TetMesh, sans sa geometrie une valeur de champ ne se relie pas a une position 3D.
        color_map: peinture utilisateur, ou None (pipeline sans couleur : cas d'usage par defaut).
        region_dependency_graph: sortie de step3.
        scalar_field: sortie de step4, pour convertir une position en indice de couche (deadlines).
        params: parametres de projection, de cout et d'ordonnancement.

    Returns:
        (layers recolorees et re-segmentees, PrintSchedule). Les couches retournees sont de NOUVEAUX objets
        (`toolpath_2d` re-decoupe aux transitions de couleur, `segment_colors` renseigne) ; les `PrintPass`
        pointent vers CES couches-la.

    Raises:
        SchedulingError: si une couche a un `toolpath_2d` non rempli (step6 n'a pas ete execute).
    """
    for layer in layers:
        if layer.toolpath_2d is None:
            raise SchedulingError(
                f"La couche {layer.id} n'a pas de toolpath_2d - "
                f"pipeline/step6_toolpath2d.py doit etre execute avant l'etape 7."
            )

    native = _native_schedule(layers, mesh_data, tet_mesh, color_map, scalar_field, params)
    if native is not None:
        return native
    return _python_schedule(layers, mesh_data, tet_mesh, color_map, region_dependency_graph, scalar_field, params)


@dataclass(frozen=True)
class DeadlineContext:
    """Donnees partagees par le calcul des deadlines de trous (maillage, champ, couches)."""
    mesh_data: MeshData
    field_tree: KDTree
    field_values: np.ndarray
    sorted_isovalues: np.ndarray
    layer_count: int


class _PassScheduler:
    """Etat de `_schedule_passes` : passes emises, trous en attente et couleur active."""

    def __init__(self, deadline_of: "callable", params: Step7Params):
        self.deadline_of, self.params = deadline_of, params
        self.passes: list[PrintPass] = []
        self.pending: list[_PendingHole] = []
        self.active_color: Optional[int] = None

    def flush_holes(self, up_to_layer_index: int, force_all: bool = False) -> None:
        """Comble les trous dont la deadline est atteinte (EDF)."""
        self.pending.sort(key=lambda hole: hole.deadline)
        remaining: list[_PendingHole] = []
        for hole in self.pending:
            must_fill = force_all or hole.deadline <= up_to_layer_index
            if not must_fill:
                remaining.append(hole)
                continue
            fill_pass = _make_pass([hole.group], hole.group.color_id, is_hole_fill=True)
            if fill_pass is not None:
                self.passes.append(fill_pass)
                self.active_color = hole.group.color_id
        self.pending[:] = remaining

    def schedule_layer(self, layer_index: int, groups: list[_ColorGroup]) -> None:
        free_groups = [g for g in groups if g.color_id is None]
        colored_groups = [g for g in groups if g.color_id is not None]

        if not colored_groups:
            free_pass = _make_pass(free_groups, self.active_color, is_hole_fill=False)
            if free_pass is not None:
                self.passes.append(free_pass)
            self.flush_holes(layer_index)
            return

        # Continuer la couleur active si elle est presente dans cette couche, sinon prendre la couleur la
        # plus representee (moins de switches).
        present_colors = {g.color_id: g for g in colored_groups}
        if self.active_color in present_colors:
            chosen = present_colors[self.active_color]
        else:
            chosen = max(colored_groups, key=lambda g: len(g.segment_indices))

        # La matiere invisible part avec la passe courante : imprimable avec n'importe quelle couleur chargee.
        printed_pass = _make_pass([chosen, *free_groups], chosen.color_id, is_hole_fill=False)
        if printed_pass is not None:
            self.passes.append(printed_pass)
        self.active_color = chosen.color_id

        for group in colored_groups:
            if group is chosen:
                continue
            deadline = self.deadline_of(layer_index, group)
            if _should_defer(
                open_layer_span=deadline - layer_index,
                deadline=deadline,
                fill_layer_index=layer_index + 1,
                open_hole_count=len(self.pending),
                params=self.params,
            ):
                self.pending.append(_PendingHole(group=group, deadline=deadline))
            else:
                purge_pass = _make_pass([group], group.color_id, is_hole_fill=False)
                if purge_pass is not None:
                    self.passes.append(purge_pass)
                self.active_color = group.color_id

        self.flush_holes(layer_index)


def _python_schedule(layers, mesh_data, tet_mesh, color_map, region_dependency_graph, scalar_field, params):
    """Chemin Python : couleurs projetees puis passes ordonnancees."""
    colored_layers = _apply_colors_to_layers(layers, mesh_data, color_map, params)
    _progress.milestone(0.5)                 # R-017 : couleurs projetees ; reste l'ordonnancement
    groups_per_layer = _group_segments_by_color(colored_layers)

    # Graphe (region, color) : sortie intermediaire exploitable par l'UI et garde-fou de coherence.
    build_color_dependency_graph(groups_per_layer, region_dependency_graph)

    directions = _escape_directions(params)
    layer_count = len(colored_layers)

    # phi est defini par NOEUD du TetMesh : le kdtree de lecture du champ se construit donc sur `tet_vertices`,
    # dont `scalar_field.values` partage exactement l'indexation.
    field_tree = KDTree(tet_mesh.tet_vertices)
    field_values = np.asarray(scalar_field.values, dtype=float)

    isovalues = [layer.isovalue for layer in colored_layers]
    has_isovalues = layer_count > 0 and all(value is not None for value in isovalues)
    sorted_isovalues = (
        np.sort(np.asarray(isovalues, dtype=float)) if has_isovalues else np.zeros(0)
    )

    deadline_cache: dict[tuple[int, Optional[int]], int] = {}

    def deadline_of(layer_index: int, group: _ColorGroup) -> int:
        if not has_isovalues:
            # Sans isovaleur par couche, aucune conversion phi -> indice de couche n'est possible : on ne
            # contraint aucun trou plutot que d'inventer une deadline arbitraire.
            return layer_count - 1
        key = (layer_index, group.color_id)
        if key not in deadline_cache:
            position = _group_position_3d(colored_layers[layer_index], group)
            deadline_cache[key] = compute_hole_deadline(
                position,
                layer_index,
                DeadlineContext(mesh_data, field_tree, field_values, sorted_isovalues, layer_count),
                directions,
                params,
            )
        return deadline_cache[key]

    return colored_layers, _schedule_passes(groups_per_layer, deadline_of, params)
