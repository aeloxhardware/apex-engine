"""Etape 7 : couleur par arête, lissage, découpage par couleur et groupes de segments.

Découpé de `step7_scheduling.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step7_scheduling`.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional
import networkx as nx
import numpy as np
import trimesh
import trimesh.proximity
from core.mesh_types import ColorMap, Layer, MeshData, RegionDependencyGraph
from element_ardor.step6_toolpath2d import _lift_to_3d
from element_ardor.step7_params import (  # noqa: F401
    Step7Params,
)


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
