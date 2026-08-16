"""
pipeline/step6_toolpath2d.py

Etape 6 du pipeline : remplissage 2D de chaque couche (parois + infill),
dans le repere local de la couche produit par step5.

Depend de Layer (step5, contour_3d + local_frame), MeshData (step1, kdtree
pour le pont surface/volume), OverhangMap (step2, is_bridge +
bridge_span_direction).

Ce que ce fichier fait, dans l'ordre, par Layer :
  1. Projette chaque boucle de contour_3d dans le plan local (local_frame).
  2. Classe les boucles projetees en contours "matiere" / trous par
     imbrication geometrique (point_in_polygon_2d), PAS par signe d'aire -
     voir CLAUDE.md limite #11 et la Decision 3bis de step5.
  3. Genere les parois par offset polygonal successif vers l'interieur
     (shapely.buffer), wall_count niveaux, en s'arretant proprement si
     l'offset fait disparaitre une zone trop fine avant d'atteindre
     wall_count.
  4. Genere l'infill a l'interieur du dernier offset.
  5. Pour les boucles associees a une zone is_bridge (step2, via le meme
     pont kdtree surface/volume que step4/step5) : ajoute des parois
     d'ancrage supplementaires sur le perimetre puis reoriente l'infill
     pour traverser perpendiculairement a bridge_span_direction.

Sortie : la meme liste de Layer, avec toolpath_2d rempli (repere local 2D,
pas 3D - step8 se chargera de le relever au repere machine).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum

import numpy as np
import shapely
import shapely.affinity
from shapely.geometry import LineString, Polygon

from core.geometry_utils import point_in_polygon_2d
from core.mesh_types import Layer, MeshData, OverhangMap


class ToolpathGenerationError(Exception):
    """Levee si une couche ne peut produire aucune trajectoire exploitable."""


class InfillPattern(Enum):
    LINES = "lines"
    GRID = "grid"


@dataclass
class Step6Params:
    nozzle_width: float = 0.4
    wall_count: int = 2
    infill_density: float = 0.2          # dans (0, 1]
    infill_pattern: InfillPattern = InfillPattern.LINES
    # seuil de PRESENCE (pas de majorite) : proportion minimale de points
    # d'une boucle de contour dont la face de surface la plus proche est
    # is_bridge, au-dela de laquelle la boucle entiere est traitee comme
    # trajectoire de pont. Un pont n'est qu'une petite portion locale du
    # perimetre total d'une couche (CLAUDE.md : zones is_bridge "qui
    # retombent dans cette couche", pas la couche entiere) - un seuil de
    # majorite ne peut structurellement jamais se declencher (mesure sur la
    # fixture pont : ratio reel entre 0.06 et 0.21 sur toutes les couches
    # concernees, jamais proche de 0.5). D'ou un defaut bas.
    bridge_contour_ratio_threshold: float = 0.05
    # parois supplementaires resserrees sur le perimetre d'un contour pont,
    # en plus de wall_count, avant de generer l'infill balaye - lecture
    # directe de CLAUDE.md "augmente la densite d'ancrage sur les bords"
    # (renfort de perimetre, pas de densite d'infill - cf. plan step6).
    bridge_anchor_extra_walls: int = 2
    # aire residuelle (relative a l'aire du contour d'origine) sous laquelle
    # un offset est considere comme ayant fait disparaitre la zone.
    degenerate_area_ratio: float = 1e-6


# ---------------------------------------------------------------------------
# Projection / relevement repere local <-> 3D
# ---------------------------------------------------------------------------

def _project_to_local_2d(
    points_3d: np.ndarray, local_frame: tuple[np.ndarray, np.ndarray, np.ndarray]
) -> np.ndarray:
    origin, axis_u, axis_v = local_frame
    relative = points_3d - origin
    u = relative @ axis_u
    v = relative @ axis_v
    return np.stack([u, v], axis=1)


def _lift_to_3d(
    points_2d: np.ndarray, local_frame: tuple[np.ndarray, np.ndarray, np.ndarray]
) -> np.ndarray:
    origin, axis_u, axis_v = local_frame
    return origin + points_2d[:, 0:1] * axis_u + points_2d[:, 1:2] * axis_v


# ---------------------------------------------------------------------------
# Decision 1 - Classification contour exterieur / trou par imbrication
# ---------------------------------------------------------------------------

def _signed_area_2d(polygon: np.ndarray) -> float:
    x = polygon[:, 0]
    y = polygon[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def _classify_contours(
    contours_2d: list[np.ndarray],
) -> list[tuple[np.ndarray, list[np.ndarray]]]:
    """
    Classe les boucles projetees en groupes (contour_matiere, [trous]) par
    imbrication geometrique (point_in_polygon_2d), PAS par signe d'aire
    (hypothese refutee empiriquement - cf. CLAUDE.md limite #11).

    Regle de parite : pour une boucle donnee, le nombre de boucles (parmi
    les plus grandes qu'elle) qui la contiennent determine son role -
    pair (0, 2, 4...) = matiere, impair (1, 3...) = trou. Generalise le cas
    a 2 boucles a un nombre arbitraire de boucles imbriquees (anneau dans
    anneau).

    Chaque trou est ensuite rattache au plus petit contour "matiere" qui le
    contient directement (son parent immediat), pas necessairement le plus
    grand contour de la couche - important des qu'il y a plusieurs ilots
    matiere disjoints.
    """
    if not contours_2d:
        return []

    # une seule boucle sans point de fermeture duplique pour les calculs
    open_loops = [c[:-1] if np.allclose(c[0], c[-1]) else c for c in contours_2d]
    areas = [abs(_signed_area_2d(loop)) for loop in open_loops]
    order = sorted(range(len(open_loops)), key=lambda i: areas[i], reverse=True)

    # containers[i] = liste des indices de boucles (parmi celles deja
    # traitees, donc plus grandes) qui contiennent la boucle i
    containers: dict[int, list[int]] = {i: [] for i in range(len(open_loops))}
    for rank, i in enumerate(order):
        probe_point = open_loops[i][0]
        for j in order[:rank]:
            if point_in_polygon_2d(probe_point, open_loops[j]):
                containers[i].append(j)

    is_hole = {i: (len(containers[i]) % 2 == 1) for i in range(len(open_loops))}

    # rattache chaque trou a son parent immediat = le plus petit
    # englobant (aire minimale parmi ses conteneurs)
    groups: dict[int, list[int]] = {i: [] for i in range(len(open_loops)) if not is_hole[i]}
    for i in range(len(open_loops)):
        if not is_hole[i]:
            continue
        parents = containers[i]
        if not parents:
            # trou sans parent matiere (topologiquement incoherent pour une
            # isosurface bien formee) : ignore plutot que de planter -
            # meme esprit defensif que step5 (_map_tets_to_regions).
            continue
        parent = min(parents, key=lambda p: areas[p])
        groups[parent].append(i)

    return [
        (open_loops[outer_idx], [open_loops[h] for h in holes])
        for outer_idx, holes in groups.items()
    ]


# ---------------------------------------------------------------------------
# Decision 3 - Offset polygonal progressif (parois), arret propre si degenere
# ---------------------------------------------------------------------------

def _to_shapely_polygon(outer: np.ndarray, holes: list[np.ndarray]) -> Polygon:
    return Polygon(outer, holes=[h for h in holes])


def _offset_walls(
    outer: np.ndarray, holes: list[np.ndarray], nozzle_width: float, wall_count: int, degenerate_area_ratio: float
) -> tuple[list[np.ndarray], Polygon | None]:
    """
    Genere les parois d'un contour matiere par offset progressif vers
    l'interieur (n = wall_count .. 1, puis le contour original n=0),
    ordonnees du centre vers l'exterieur dans la liste retournee.

    S'arrete proprement (sans erreur) des qu'un offset fait disparaitre la
    zone - une paroi plus fine que nozzle_width * wall_count ne peut pas
    recevoir wall_count parois pleines, cf. plan step6 Question 3. Retourne
    aussi le polygone du dernier offset valide (pour servir de frontiere a
    l'infill), ou None si meme n=1 a disparu.
    """
    base_polygon = _to_shapely_polygon(outer, holes)
    if base_polygon.is_empty or base_polygon.area <= 0:
        return [], None

    reference_area = base_polygon.area
    walls: list[tuple[int, np.ndarray]] = []
    innermost: Polygon | None = None

    for n in range(wall_count, 0, -1):
        offset = base_polygon.buffer(-nozzle_width * n, join_style="mitre")
        polygon = _largest_polygon(offset)
        if polygon is None or polygon.area <= degenerate_area_ratio * reference_area:
            continue
        walls.append((n, np.asarray(polygon.exterior.coords)))
        if innermost is None:
            innermost = polygon

    # parois du centre vers l'exterieur : n decroissant a deja ce sens
    ordered = [pts for _n, pts in walls]
    ordered.append(outer)  # n=0, contour original, en dernier (le plus exterieur)

    return ordered, innermost


def _largest_polygon(geometry) -> Polygon | None:
    """
    `buffer` peut retourner un Polygon, un MultiPolygon (zone separee en
    plusieurs morceaux), ou une geometrie vide. On ne garde que le plus
    grand morceau : les morceaux secondaires, plus petits, sont un signe
    de la meme degenerescence que Question 3 (zone trop fine) sur une
    sous-partie du contour - traites de la meme facon (silencieusement
    ignores), pas une erreur.
    """
    if geometry.is_empty:
        return None
    if geometry.geom_type == "Polygon":
        return geometry
    if geometry.geom_type == "MultiPolygon":
        parts = list(geometry.geoms)
        if not parts:
            return None
        return max(parts, key=lambda p: p.area)
    return None


# ---------------------------------------------------------------------------
# Decision 4 - Infill
# ---------------------------------------------------------------------------

def _generate_infill(
    boundary: Polygon, density: float, pattern: InfillPattern, sweep_angle_rad: float = 0.0
) -> list[np.ndarray]:
    """
    Genere l'infill a l'interieur de `boundary` par un balayage de lignes
    paralleles, clippees par le polygone (shapely.intersection). L'espacement
    entre lignes est derive de density (1.0 = lignes espacees d'une largeur
    de buse, densite plus faible = espacement plus grand) - premiere passe
    simple, suffisante pour l'objectif actuel (pas d'optimisation de trajet).

    sweep_angle_rad force l'orientation du balayage (utilise pour les ponts,
    cf. Decision 2 du plan : traverser perpendiculairement a la portee).
    """
    if boundary is None or boundary.is_empty or density <= 0:
        return []

    minx, miny, maxx, maxy = boundary.bounds
    diagonal = float(np.hypot(maxx - minx, maxy - miny))
    if diagonal <= 0:
        return []

    spacing = max(diagonal * 0.02, diagonal / max(density, 1e-6) * 0.05)
    # spacing decroit quand density augmente ; borne basse pour eviter un
    # nombre de lignes explosif sur des couches degenerees en aire.
    spacing = max(spacing, 1e-6)

    center = np.array([(minx + maxx) / 2.0, (miny + maxy) / 2.0])
    half_extent = diagonal  # sur-couvre largement pour rester valide apres rotation

    lines: list[np.ndarray] = []
    offsets = np.arange(-half_extent, half_extent, spacing)
    cos_a, sin_a = np.cos(sweep_angle_rad), np.sin(sweep_angle_rad)
    direction = np.array([cos_a, sin_a])
    perpendicular = np.array([-sin_a, cos_a])

    for offset in offsets:
        base = center + perpendicular * offset
        p0 = base - direction * half_extent
        p1 = base + direction * half_extent
        segment = LineString([p0, p1])
        clipped = segment.intersection(boundary)
        if clipped.is_empty:
            continue
        if clipped.geom_type == "LineString":
            coords = np.asarray(clipped.coords)
            if len(coords) >= 2:
                lines.append(coords)
        elif clipped.geom_type == "MultiLineString":
            for part in clipped.geoms:
                coords = np.asarray(part.coords)
                if len(coords) >= 2:
                    lines.append(coords)

    return lines


# ---------------------------------------------------------------------------
# Decision 2 - Pont surface (is_bridge) <-> boucle 2D
# ---------------------------------------------------------------------------

def _bridge_direction_for_loop(
    contour_3d_loop: np.ndarray,
    mesh_data: MeshData,
    overhang_map: OverhangMap,
    local_frame: tuple[np.ndarray, np.ndarray, np.ndarray],
    threshold: float,
) -> np.ndarray | None:
    """
    Determine si une boucle de contour 3D touche une zone is_bridge (step2),
    via le meme pont kdtree surface/volume deja utilise par step4 (ancrage
    plateau) et step5 (association region). Retourne la direction de portee
    du pont projetee dans le plan local de la couche (unitaire), ou None si
    la boucle ne touche aucune zone de pont.

    Seuil de PRESENCE, pas de majorite : un pont est une sous-region locale
    de la matiere d'une couche, jamais une majorite du perimetre total de sa
    boucle de contour (mesure directe sur la fixture pont : ratio entre 0.06
    et 0.21 sur toutes les couches concernees) - voir CORRECTION 1 du plan.
    """
    _, nearest_faces = mesh_data.kdtree.query(contour_3d_loop)
    nearest_faces = np.asarray(nearest_faces, dtype=np.int64)

    is_bridge_per_point = overhang_map.is_bridge[nearest_faces]
    ratio = float(is_bridge_per_point.mean()) if len(is_bridge_per_point) else 0.0
    if ratio < threshold:
        return None

    span_vectors = overhang_map.bridge_span_direction[nearest_faces[is_bridge_per_point]]
    valid = ~np.isnan(span_vectors).any(axis=1)
    if not np.any(valid):
        return None

    mean_span = span_vectors[valid].mean(axis=0)
    norm = np.linalg.norm(mean_span)
    if norm < 1e-12:
        return None
    mean_span = mean_span / norm

    _origin, axis_u, axis_v = local_frame
    span_2d = np.array([np.dot(mean_span, axis_u), np.dot(mean_span, axis_v)])
    span_norm = np.linalg.norm(span_2d)
    if span_norm < 1e-12:
        return None
    return span_2d / span_norm


def _reorient_bridge_toolpath(
    outer: np.ndarray,
    holes: list[np.ndarray],
    span_direction_2d: np.ndarray,
    params: Step6Params,
) -> list[np.ndarray]:
    """
    Trajectoire d'un contour pont : parois normales (wall_count), PUIS
    bridge_anchor_extra_walls parois supplementaires resserrees sur le
    perimetre (renfort d'ancrage sur les bords, cf. plan step6 Decision
    explicite : parois supplementaires plutot que densite d'infill
    augmentee), PUIS infill balaye perpendiculairement a span_direction_2d
    (la trajectoire doit TRAVERSER le pont, donc etre perpendiculaire a son
    axe de plus grande portee - balayer selon la portee reviendrait a
    longer le vide plutot qu'a le franchir).
    """
    total_wall_levels = params.wall_count + params.bridge_anchor_extra_walls
    walls, innermost = _offset_walls(
        outer, holes, params.nozzle_width, total_wall_levels, params.degenerate_area_ratio
    )
    if not walls:
        return [outer]

    sweep_angle = float(np.arctan2(span_direction_2d[1], span_direction_2d[0])) + np.pi / 2.0
    infill = _generate_infill(innermost, params.infill_density, params.infill_pattern, sweep_angle)

    return infill + walls


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

def run(
    layers: list[Layer],
    mesh_data: MeshData,
    overhang_map: OverhangMap,
    params: Step6Params,
) -> list[Layer]:
    """
    Genere la trajectoire 2D remplie de chaque couche (parois + infill,
    gestion des ponts), dans le repere local (2D) de chaque Layer.

    Args:
        layers: sortie de pipeline/step5_layers.py.
        mesh_data: sortie de pipeline/step1_ingest.py (kdtree, pont
            surface/volume pour la detection de pont).
        overhang_map: sortie de pipeline/step2_overhang.py (is_bridge,
            bridge_span_direction).
        params: parametres de parois/infill/ponts.

    Returns:
        Nouvelle liste de Layer (pas de mutation en place), avec
        toolpath_2d rempli - liste de polylignes (K,2) dans le repere local
        de la couche, ordonnees du centre vers l'exterieur.
    """
    result: list[Layer] = []

    for layer in layers:
        contours_2d = [_project_to_local_2d(c, layer.local_frame) for c in layer.contour_3d]
        groups = _classify_contours(contours_2d)

        toolpath: list[np.ndarray] = []

        for outer, holes in groups:
            bridge_direction = None
            for contour_3d_loop in layer.contour_3d:
                loop_2d = _project_to_local_2d(contour_3d_loop, layer.local_frame)
                if not _loop_matches(loop_2d, outer):
                    continue
                bridge_direction = _bridge_direction_for_loop(
                    contour_3d_loop, mesh_data, overhang_map, layer.local_frame,
                    params.bridge_contour_ratio_threshold,
                )
                break

            if bridge_direction is not None:
                toolpath.extend(
                    _reorient_bridge_toolpath(outer, holes, bridge_direction, params)
                )
                continue

            walls, innermost = _offset_walls(
                outer, holes, params.nozzle_width, params.wall_count, params.degenerate_area_ratio
            )
            if not walls:
                continue
            infill = _generate_infill(innermost, params.infill_density, params.infill_pattern)
            toolpath.extend(infill + walls)

        result.append(replace(layer, toolpath_2d=toolpath))

    return result


def _loop_matches(loop_2d: np.ndarray, outer_open: np.ndarray) -> bool:
    """
    Rapproche une boucle projetee (avec point de fermeture duplique) du
    contour "matiere" ouvert issu de _classify_contours (meme boucle,
    representation differente) - comparaison sur le premier point, unique
    par construction (les boucles d'une meme couche ne se croisent pas).
    """
    open_loop = loop_2d[:-1] if np.allclose(loop_2d[0], loop_2d[-1]) else loop_2d
    return open_loop.shape == outer_open.shape and np.allclose(open_loop, outer_open)
