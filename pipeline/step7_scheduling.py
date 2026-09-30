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


class SchedulingError(Exception):
    """Levee si le planning ne peut pas etre construit (donnees incoherentes)."""


@dataclass
class HoleCostModel:
    """Cout d'un trou temporaire (contrat UI : SchedulingParams.hole_cost_model)."""
    reheat_penalty: float = 1.0
    adhesion_risk_weight: float = 1.0


@dataclass
class PurgeCostModel:
    """Cout d'une purge classique (contrat UI : SchedulingParams.purge_cost_model)."""
    volume_per_purge: float = 1.0
    time_per_purge: float = 1.0


@dataclass
class Step7Params:
    # Distance max entre un point de trajectoire et la surface du mesh
    # au-dela de laquelle ce point est considere INTERNE (invisible), donc
    # sans couleur imposee. Typiquement de l'ordre d'une largeur de buse.
    surface_color_distance: float = 0.5
    hole_cost_model: HoleCostModel = field(default_factory=HoleCostModel)
    purge_cost_model: PurgeCostModel = field(default_factory=PurgeCostModel)
    max_open_holes: int = 4
    # Directions d'echappement testees pour la deadline d'accessibilite.
    escape_direction_count: int = 16
    # Cinematique de rotation de la machine : si None/vide, les directions
    # d'echappement sont prises sur la sphere complete (orientation libre).
    # Si fourni, decrit une cinematique arbitraire (N axes composes, chacun
    # avec sa propre plage) via composed_rotation_directions - PAS figee sur
    # une seule machine (voir CLAUDE.md limite de conception #18). Un seul
    # axe a plage complete degenere vers l'ancien comportement a 1 DOF.
    rotation_axes: Optional[list[RotationAxisSpec]] = None
    # Nombre d'echantillons le long d'un rayon d'echappement.
    deadline_ray_samples: int = 24

    # NOTE: `merge_tolerance` n'est deliberement PAS duplique ici - il vit
    # sur ColorMap (contrat UI), seule source de verite. Deux sources qui
    # divergent en silence est exactement le piege deja rencontre avec
    # `platform_z` a l'etape 4.


# ---------------------------------------------------------------------------
# Question 1 - Projection ColorMap (faces surface) -> segments de trajectoire
# ---------------------------------------------------------------------------

def _edge_midpoints_2d(polyline_2d: np.ndarray) -> np.ndarray:
    """Milieu de chaque arete de la polyligne (K,2) -> (K-1,2)."""
    return 0.5 * (polyline_2d[:-1] + polyline_2d[1:])


def _edge_lengths(polyline_2d: np.ndarray) -> np.ndarray:
    return np.linalg.norm(polyline_2d[1:] - polyline_2d[:-1], axis=1)


def _edge_colors_for_polyline(
    polyline_2d: np.ndarray,
    local_frame: tuple[np.ndarray, np.ndarray, np.ndarray],
    proximity: trimesh.proximity.ProximityQuery,
    color_map: ColorMap,
    params: Step7Params,
) -> list[Optional[int]]:
    """
    Couleur peinte de chaque ARETE de la polyligne (une couleur par arete,
    pas par point : une arete est l'unite reellement extrudee, et raisonner
    par arete evite de produire des sous-polylignes degenerees a 1 point lors
    du decoupage par couleur).

    Deux adaptations par rapport au pont surface/volume de step6 :

    1. `toolpath_2d` vit dans le repere local 2D de la couche, pas en 3D :
       il faut le relever d'abord (`step6_toolpath2d._lift_to_3d`).
    2. Contrairement a `contour_3d` (par construction SUR la surface du mesh
       source), la majorite de `toolpath_2d` est INTERNE au volume (infill,
       parois internes). Une couleur peinte etant une propriete de surface
       VISIBLE, il faut une vraie DISTANCE A LA SURFACE pour trancher.

    IMPORTANT - pourquoi PAS `mesh_data.kdtree` ici : ce kdtree est construit
    sur les CENTRES DE FACE (`step1_ingest`), donc la distance qu'il retourne
    est une distance a un centre de face, PAS a la surface. Mesure sur un
    cube (12 grands triangles) : des points situes exactement SUR la surface
    ressortent a une distance de 0.370 a 3.901 du centre de face le plus
    proche, alors que leur distance reelle a la surface est 0.000 - avec un
    seuil de 0.5, seuls 3 points sur 35 etaient correctement classes visibles.
    Step6 n'etait pas affecte car il n'utilise que l'INDICE ("quelle face est
    la plus proche"), jamais la distance. Voir CLAUDE.md limite #15.

    `ProximityQuery.on_surface` retourne les DEUX grandeurs correctes : la
    distance reelle a la surface et l'identifiant du triangle reellement le
    plus proche (l'indice du kdtree devient lui aussi douteux pres des aretes
    sur un maillage grossier).

    Une arete a plus de `surface_color_distance` de la surface recoit `None`
    : matiere invisible, imprimable avec n'importe quelle couleur deja
    chargee. Ce n'est pas un cas degrade mais le principal levier de
    reduction des purges - sans lui, l'infill forcerait des changements de
    couleur qui ne se voient pas.
    """
    if len(polyline_2d) < 2:
        return []

    midpoints_3d = _lift_to_3d(_edge_midpoints_2d(polyline_2d), local_frame)
    _closest, distances, triangle_ids = proximity.on_surface(midpoints_3d)
    distances = np.atleast_1d(np.asarray(distances, dtype=float))
    triangle_ids = np.atleast_1d(np.asarray(triangle_ids, dtype=np.int64))

    colors: list[Optional[int]] = []
    for distance, face_id in zip(distances, triangle_ids):
        if distance > params.surface_color_distance:
            colors.append(None)  # matiere interne : couleur libre
            continue
        colors.append(color_map.per_face_color.get(int(face_id)))
    return colors


# ---------------------------------------------------------------------------
# Decision 4 - Lissage de la carte de couleur (merge_tolerance)
# ---------------------------------------------------------------------------

def _smooth_edge_colors(
    edge_colors: list[Optional[int]], edge_lengths: np.ndarray, merge_tolerance: float
) -> list[Optional[int]]:
    """
    Filtre morphologique 1D LE LONG de la trajectoire (parametree en longueur
    d'arc) : tout run de couleur plus court que `merge_tolerance` est absorbe
    par son voisin le plus long.

    Lissage le long du chemin plutot que sur le maillage de surface : ce
    qu'on cherche a borner est le nombre de CHANGEMENTS DE COULEUR le long du
    trajet (chacun coutant une purge ou un trou), pas la finesse de la
    peinture sur le mesh. Interpretation directe et verifiable : "ne jamais
    changer de couleur pour moins de X mm d'extrusion".
    """
    if merge_tolerance <= 0 or len(edge_colors) == 0:
        return edge_colors

    colors = list(edge_colors)

    while True:
        runs = _color_runs(colors)
        if len(runs) <= 1:
            return colors

        run_lengths = [float(edge_lengths[start:stop].sum()) for start, stop, _c in runs]
        shortest = int(np.argmin(run_lengths))
        if run_lengths[shortest] >= merge_tolerance:
            return colors

        start, stop, _color = runs[shortest]
        previous_length = run_lengths[shortest - 1] if shortest > 0 else -1.0
        next_length = run_lengths[shortest + 1] if shortest + 1 < len(runs) else -1.0
        absorbing = runs[shortest - 1] if previous_length >= next_length else runs[shortest + 1]

        for i in range(start, stop):
            colors[i] = absorbing[2]


def _color_runs(colors: list[Optional[int]]) -> list[tuple[int, int, Optional[int]]]:
    """Decoupe une liste de couleurs en runs contigus : [(start, stop, color), ...]."""
    runs: list[tuple[int, int, Optional[int]]] = []
    start = 0
    for i in range(1, len(colors) + 1):
        if i == len(colors) or colors[i] != colors[start]:
            runs.append((start, i, colors[start]))
            start = i
    return runs


def _split_polyline_by_color(
    polyline_2d: np.ndarray, edge_colors: list[Optional[int]]
) -> list[tuple[np.ndarray, Optional[int]]]:
    """
    Scinde une polyligne en sous-polylignes de couleur constante, pour
    garantir l'invariant "une polyligne = une couleur" dont depend tout
    l'ordonnancement aval (une passe de couleur ne peut pointer que vers des
    segments entiers).

    Les sous-polylignes partagent leur point de jonction (le dernier point de
    l'une est le premier de la suivante), donc la trajectoire reste continue.
    """
    if not edge_colors:
        return [(polyline_2d, None)]

    return [
        (polyline_2d[start : stop + 1], color)
        for start, stop, color in _color_runs(edge_colors)
    ]


def _apply_colors_to_layers(
    layers: list[Layer],
    mesh_data: MeshData,
    color_map: Optional[ColorMap],
    params: Step7Params,
) -> list[Layer]:
    """
    Projette la couleur sur chaque segment, lisse, puis re-segmente chaque
    couche. Retourne de NOUVEAUX Layer (jamais de mutation en place, meme
    style pur que step6).
    """
    result: list[Layer] = []

    # Construit UNE SEULE FOIS la structure de proximite surfacique (et
    # seulement si une peinture existe reellement) : la reconstruire par
    # polyligne serait le cout dominant de l'etape.
    proximity = None
    if color_map is not None and color_map.per_face_color:
        surface = trimesh.Trimesh(
            vertices=mesh_data.vertices, faces=mesh_data.faces, process=False
        )
        proximity = trimesh.proximity.ProximityQuery(surface)

    for layer in layers:
        toolpath = layer.toolpath_2d or []

        if color_map is None or not color_map.per_face_color:
            # Pipeline sans peinture : cas d'usage par defaut, doit rester
            # pleinement fonctionnel (toute la matiere est "libre").
            result.append(
                replace(
                    layer,
                    toolpath_2d=list(toolpath),
                    segment_colors=[None] * len(toolpath),
                    color_id=None,
                )
            )
            continue

        new_toolpath: list[np.ndarray] = []
        new_colors: list[Optional[int]] = []

        for polyline in toolpath:
            edge_colors = _edge_colors_for_polyline(
                polyline, layer.local_frame, proximity, color_map, params
            )
            edge_colors = _smooth_edge_colors(
                edge_colors, _edge_lengths(polyline), color_map.merge_tolerance
            )
            for sub_polyline, color in _split_polyline_by_color(polyline, edge_colors):
                if len(sub_polyline) < 2:
                    continue
                new_toolpath.append(sub_polyline)
                new_colors.append(color)

        result.append(
            replace(
                layer,
                toolpath_2d=new_toolpath,
                segment_colors=new_colors,
                color_id=_dominant_color(new_colors),
            )
        )

    return result


def _dominant_color(colors: list[Optional[int]]) -> Optional[int]:
    """Couleur la plus representee (hors None) - purement indicatif pour l'UI."""
    counts: dict[int, int] = {}
    for color in colors:
        if color is not None:
            counts[color] = counts.get(color, 0) + 1
    if not counts:
        return None
    return max(counts, key=counts.get)


# ---------------------------------------------------------------------------
# Decision 5 - Graphe (region, color) -> dependances physiques
# ---------------------------------------------------------------------------

@dataclass
class _ColorGroup:
    """Les segments d'une meme couleur dans une meme couche."""
    layer_index: int              # index dans la liste ordonnee des couches
    layer_id: int
    color_id: Optional[int]
    segment_indices: list[int]
    region_ids: list[int]


def _group_segments_by_color(layers: list[Layer]) -> list[list[_ColorGroup]]:
    """Regroupe les segments de chaque couche par couleur, en preservant l'ordre."""
    groups_per_layer: list[list[_ColorGroup]] = []

    for layer_index, layer in enumerate(layers):
        colors = layer.segment_colors or []
        by_color: dict[Optional[int], list[int]] = {}
        for segment_index, color in enumerate(colors):
            by_color.setdefault(color, []).append(segment_index)

        groups_per_layer.append(
            [
                _ColorGroup(
                    layer_index=layer_index,
                    layer_id=layer.id,
                    color_id=color,
                    segment_indices=indices,
                    region_ids=list(layer.region_ids),
                )
                for color, indices in by_color.items()
            ]
        )

    return groups_per_layer


def build_color_dependency_graph(
    groups_per_layer: list[list[_ColorGroup]],
    region_dependency_graph: RegionDependencyGraph,
) -> nx.DiGraph:
    """
    Graphe oriente sur les noeuds `(region_id, color_id)` reellement presents.

    Deux familles d'aretes :
      - HERITEES : chaque arete `a -> b` de `RegionDependencyGraph` (step3) se
        releve en `(a, *) -> (b, *)` : la contrainte geometrique entre regions
        vaut quelles que soient les couleurs.
      - SUPPORT : une couleur presente en couche `k` a besoin de la matiere de
        la couche `k-1` sous elle. Modelise a la granularite (couche, couleur)
        et non segment-a-segment : le segment-a-segment exploserait en cout
        pour un gain d'ordonnancement nul a ce stade.

    Expose publiquement (pas de prefixe `_`) car c'est une sortie
    intermediaire utile a l'UI et aux etapes aval pour expliquer un ordre.
    """
    graph = nx.DiGraph()

    for groups in groups_per_layer:
        for group in groups:
            for region_id in group.region_ids:
                graph.add_node((region_id, group.color_id))

    for region_a, region_b in region_dependency_graph.dependencies:
        for node_a in [n for n in graph.nodes if n[0] == region_a]:
            for node_b in [n for n in graph.nodes if n[0] == region_b]:
                graph.add_edge(node_a, node_b)

    for layer_index in range(1, len(groups_per_layer)):
        below = groups_per_layer[layer_index - 1]
        above = groups_per_layer[layer_index]
        for group_above in above:
            for group_below in below:
                shared_regions = set(group_above.region_ids) & set(group_below.region_ids)
                for region_id in shared_regions:
                    node_below = (region_id, group_below.color_id)
                    node_above = (region_id, group_above.color_id)
                    if node_below == node_above or not graph.has_node(node_below):
                        continue
                    graph.add_edge(node_below, node_above)

    return graph


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
    mesh_data: MeshData,
    field_tree: KDTree,
    field_values: np.ndarray,
    sorted_isovalues: np.ndarray,
    layer_count: int,
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
    mesh_data: MeshData,
    field_tree: KDTree,
    field_values: np.ndarray,
    sorted_isovalues: np.ndarray,
    layer_count: int,
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
            mesh_data,
            field_tree,
            field_values,
            sorted_isovalues,
            layer_count,
            params.deadline_ray_samples,
        )
        best_blocking = max(best_blocking, blocking)

    deadline = best_blocking - 1
    return int(np.clip(deadline, hole_layer_index, layer_count - 1))


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
    Parcourt les couches dans l'ordre topologique de base (l'ordre des Layer
    de step5 est deja valide), en maintenant une couleur active aussi
    longtemps que possible. Les couleurs non actives d'une couche sont soit
    differees en trou temporaire, soit imprimees immediatement au prix d'une
    purge (Decision 6).

    Les trous en attente sont combles en earliest-deadline-first : des qu'un
    trou approche sa deadline, ou que `max_open_holes` est atteint, la passe
    de comblement est inseree.
    """
    passes: list[PrintPass] = []
    pending: list[_PendingHole] = []
    active_color: Optional[int] = None

    def flush_holes(up_to_layer_index: int, force_all: bool = False) -> None:
        """Comble les trous dont la deadline est atteinte (EDF)."""
        nonlocal active_color
        pending.sort(key=lambda hole: hole.deadline)
        remaining: list[_PendingHole] = []
        for hole in pending:
            must_fill = force_all or hole.deadline <= up_to_layer_index
            if not must_fill:
                remaining.append(hole)
                continue
            fill_pass = _make_pass([hole.group], hole.group.color_id, is_hole_fill=True)
            if fill_pass is not None:
                passes.append(fill_pass)
                active_color = hole.group.color_id
        pending[:] = remaining

    for layer_index, groups in enumerate(groups_per_layer):
        free_groups = [g for g in groups if g.color_id is None]
        colored_groups = [g for g in groups if g.color_id is not None]

        if not colored_groups:
            free_pass = _make_pass(free_groups, active_color, is_hole_fill=False)
            if free_pass is not None:
                passes.append(free_pass)
            flush_holes(layer_index)
            continue

        # Continuer la couleur active si elle est presente dans cette couche,
        # sinon prendre la couleur la plus representee (moins de switches).
        present_colors = {g.color_id: g for g in colored_groups}
        if active_color in present_colors:
            chosen = present_colors[active_color]
        else:
            chosen = max(colored_groups, key=lambda g: len(g.segment_indices))

        # La matiere invisible part avec la passe courante : elle est
        # imprimable avec n'importe quelle couleur deja chargee.
        printed_pass = _make_pass([chosen, *free_groups], chosen.color_id, is_hole_fill=False)
        if printed_pass is not None:
            passes.append(printed_pass)
        active_color = chosen.color_id

        for group in colored_groups:
            if group is chosen:
                continue
            deadline = deadline_of(layer_index, group)
            if _should_defer(
                open_layer_span=deadline - layer_index,
                deadline=deadline,
                fill_layer_index=layer_index + 1,
                open_hole_count=len(pending),
                params=params,
            ):
                pending.append(_PendingHole(group=group, deadline=deadline))
            else:
                purge_pass = _make_pass([group], group.color_id, is_hole_fill=False)
                if purge_pass is not None:
                    passes.append(purge_pass)
                active_color = group.color_id

        flush_holes(layer_index)

    flush_holes(len(groups_per_layer), force_all=True)

    return PrintSchedule(passes=passes)


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

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
    Applique la peinture utilisateur aux trajectoires et produit la sequence
    d'impression globale.

    Args:
        layers: sortie de pipeline/step6_toolpath2d.py (toolpath_2d rempli).
        mesh_data: sortie de step1 (kdtree pour la projection de couleur,
            bvh pour les rayons d'accessibilite).
        tet_mesh: sortie de step1. Requis (et pas seulement `scalar_field`)
            car `ScalarField.values` est indexe par NOEUD du TetMesh :
            sans la geometrie des noeuds, une valeur de champ ne peut pas
            etre reliee a une position 3D le long d'un rayon.
        color_map: peinture utilisateur, ou None (pipeline sans couleur -
            cas d'usage par defaut, doit rester fonctionnel).
        region_dependency_graph: sortie de step3.
        scalar_field: sortie de step4, pour convertir une position en indice
            de couche lors du calcul des deadlines.
        params: parametres de projection, de cout et d'ordonnancement.

    Returns:
        (layers recolorees et re-segmentees, PrintSchedule). Les couches
        retournees sont de NOUVEAUX objets : `toolpath_2d` y est re-decoupe
        aux transitions de couleur et `segment_colors` est renseigne. Les
        `PrintPass` pointent vers CES couches-la, pas vers celles d'entree.

    Raises:
        SchedulingError: si une couche a un `toolpath_2d` non rempli
            (step6 n'a pas ete execute).
    """
    for layer in layers:
        if layer.toolpath_2d is None:
            raise SchedulingError(
                f"La couche {layer.id} n'a pas de toolpath_2d - "
                f"pipeline/step6_toolpath2d.py doit etre execute avant l'etape 7."
            )

    colored_layers = _apply_colors_to_layers(layers, mesh_data, color_map, params)
    groups_per_layer = _group_segments_by_color(colored_layers)

    # Construit le graphe (region, color) : sortie intermediaire exploitable
    # par l'UI pour expliquer un ordre, et garde-fou de coherence.
    build_color_dependency_graph(groups_per_layer, region_dependency_graph)

    directions = _escape_directions(params)
    layer_count = len(colored_layers)

    # phi est defini par NOEUD du TetMesh : le kdtree de lecture du champ se
    # construit donc sur `tet_vertices`, dont `scalar_field.values` partage
    # exactement l'indexation.
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
            # Sans isovaleur par couche, aucune conversion phi -> indice de
            # couche n'est possible : on ne contraint alors aucun trou plutot
            # que d'inventer une deadline arbitraire.
            return layer_count - 1
        key = (layer_index, group.color_id)
        if key not in deadline_cache:
            position = _group_position_3d(colored_layers[layer_index], group)
            deadline_cache[key] = compute_hole_deadline(
                position,
                layer_index,
                mesh_data,
                field_tree,
                field_values,
                sorted_isovalues,
                layer_count,
                directions,
                params,
            )
        return deadline_cache[key]

    schedule = _schedule_passes(groups_per_layer, deadline_of, params)

    return colored_layers, schedule
