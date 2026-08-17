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

from core.geometry_utils import _rodrigues_rotate
from core.mesh_types import (
    AABB,
    CollisionEvent,
    CollisionReport,
    MachineProfile,
    RotationAxisSpec,
    ToolpathMachine,
    ToolpathPoint,
)


@dataclass
class CollisionSimParams:
    path_sampling_step: float
    safety_margin: float
    # Non utilise par step9 lui-meme - porte pour le contrat CollisionSimParams
    # (CLAUDE.md), consomme par l'etape 10 qui boucle sur step9.
    max_adjustment_iterations: int = 3
    # Taille de cellule de la grille spatiale des capsules imprimees, en
    # multiple de nozzle_diameter (voir _PrintedGeometryIndex).
    grid_cell_size_multiplier: float = 4.0
    # Pas d'echantillonnage des ARETES du modele de tete (_surface_sample_points),
    # en mm - DELIBEREMENT distinct de nozzle_diameter (largeur de buse,
    # rayon des capsules imprimees : un concept different, la taille du
    # CORDON, pas la taille de la TETE). Les deux etaient confondus a
    # l'origine (limite #21) ; mesure concretement en executant main.py sur
    # le modele de tete reel (aetheris-4.glb, aretes de plusieurs dizaines
    # de mm) avec nozzle_diameter=0.4 par defaut : 17734 points de test par
    # point de trajectoire echantillonne, simulation totale de l'ordre de
    # la dizaine de minutes sur un modele trivial - inutilisable. Defaut
    # choisi a l'echelle typique d'une tete d'impression (quelques mm),
    # PAS a l'echelle d'un cordon d'extrusion.
    toolhead_sample_spacing: float = 5.0


# ---------------------------------------------------------------------------
# Echantillonnage de la trajectoire
# ---------------------------------------------------------------------------

def _sample_point_indices(points: list[ToolpathPoint], path_sampling_step: float) -> list[int]:
    """
    Indices de points a tester, espaces d'au moins `path_sampling_step` en
    distance cumulee le long de la trajectoire - tester TOUS les points
    d'une trajectoire de plusieurs milliers de points serait un cout inutile
    pour une detection de collision qui n'a pas besoin de la resolution
    complete de l'extrusion. Le premier et le dernier point sont toujours
    inclus.
    """
    if not points:
        return []
    if path_sampling_step <= 0:
        return list(range(len(points)))

    indices = [0]
    accumulated = 0.0
    for i in range(1, len(points)):
        accumulated += float(np.linalg.norm(points[i].position - points[i - 1].position))
        if accumulated >= path_sampling_step:
            indices.append(i)
            accumulated = 0.0

    if indices[-1] != len(points) - 1:
        indices.append(len(points) - 1)
    return indices


# ---------------------------------------------------------------------------
# Positionnement du modele de tete (question 2 du plan)
# ---------------------------------------------------------------------------

def _positioned_toolhead_vertices(
    base_vertices: np.ndarray,
    position: np.ndarray,
    axis_angles: list[float],
    rotation_axes: list[RotationAxisSpec],
) -> np.ndarray:
    """
    Positionne les sommets du modele de tete (repere local, origine =
    pointe de buse, voir CLAUDE.md) en espace machine, pour un nombre
    arbitraire d'axes de rotation - meme convention de composition que
    composed_rotation_angles_and_directions/step8._axis_chain_world_lines :
    rotations appliquees en ordre INVERSE de la liste (l'axe le plus
    interne, axis_angles[-1], en premier ; le plus externe, axis_angles[0],
    en dernier), jamais une nouvelle convention inventee ici.

    Translation par `position` APRES toutes les rotations : l'origine
    locale du modele EST la pointe de buse (repere verifie
    geometriquement), donc translater apres place exactement la pointe a
    `position`. Translater avant serait un pivot faux.

    N=0 axe : la boucle ne fait rien, translation seule - degenere
    correctement sans cas special.
    """
    vertices = base_vertices
    for spec, angle_deg in zip(reversed(rotation_axes), reversed(axis_angles)):
        axis = np.asarray(spec.axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        vertices = _rodrigues_rotate(vertices, axis, np.full(len(vertices), np.radians(angle_deg)))
    return vertices + np.asarray(position, dtype=float)


def _surface_sample_points(mesh: trimesh.Trimesh, sample_spacing: float) -> np.ndarray:
    """
    Points de test contre les obstacles : sommets du mesh DE TETE PLUS des
    points echantillonnes le long de chaque arete (pas ~`sample_spacing`).

    Les sommets seuls ne suffisent pas : un obstacle fin (ex: un mur
    imprime) peut traverser une FACE du mesh sans passer pres d'aucun
    sommet, ce qui manquerait une collision reelle (faux negatif mesure
    concretement : une boite de test traversee en son milieu par un mur
    n'etait detectee par aucun des 8 sommets, tous a plusieurs mm du mur).
    Echantillonner les aretes couvre ce cas sans introduire de calcul
    point-triangle (plus couteux, pas necessaire pour un mesh de tete de
    taille modeste comme aetheris-4.glb, ~377 sommets/~1000 aretes).
    """
    edges = mesh.edges_unique
    starts = mesh.vertices[edges[:, 0]]
    ends = mesh.vertices[edges[:, 1]]
    lengths = np.linalg.norm(ends - starts, axis=1)

    extra_points = [mesh.vertices]
    for start, end, length in zip(starts, ends, lengths):
        if length < 1e-9:
            continue
        steps = max(int(np.ceil(length / max(sample_spacing, 1e-6))), 1)
        t = np.linspace(0.0, 1.0, steps + 1)[1:-1]  # exclut les sommets, deja inclus
        if len(t) > 0:
            extra_points.append(start[None, :] + t[:, None] * (end - start)[None, :])

    return np.concatenate(extra_points, axis=0)


# ---------------------------------------------------------------------------
# Geometrie imprimee incrementale (question 1 du plan)
# ---------------------------------------------------------------------------

@dataclass
class _PrintedSegment:
    start: np.ndarray
    end: np.ndarray
    radius: float
    aabb: AABB


def _segment_aabb(start: np.ndarray, end: np.ndarray, radius: float) -> AABB:
    lo = np.minimum(start, end) - radius
    hi = np.maximum(start, end) + radius
    return AABB(min_corner=lo, max_corner=hi)


def _cell_coords_for_aabb(aabb: AABB, cell_size: float) -> list[tuple[int, int, int]]:
    """Toutes les cellules (grille uniforme) que l'AABB touche, pour l'indexation par ajout."""
    lo = np.floor(aabb.min_corner / cell_size).astype(int)
    hi = np.floor(aabb.max_corner / cell_size).astype(int)
    return [
        (x, y, z)
        for x in range(lo[0], hi[0] + 1)
        for y in range(lo[1], hi[1] + 1)
        for z in range(lo[2], hi[2] + 1)
    ]


class _PrintedGeometryIndex:
    """
    Index spatial de la matiere deja imprimee, sous forme de capsules
    (segment + rayon = nozzle_diameter/2), jamais une union booleenne
    (cout prohibitif sur plusieurs milliers de points - voir plan).

    Grille maison (dict cellule -> indices de capsules) plutot qu'un rtree :
    insertion O(1) par capsule, pas de restructuration de noeuds, suffisant
    pour une densite de trajectoire a peu pres uniforme. Sert aussi de
    broad-phase pour les tests de collision (question 3 du plan).
    """

    def __init__(self, nozzle_diameter: float, cell_size_multiplier: float = 4.0):
        self._radius = max(nozzle_diameter, 1e-6) * 0.5
        self._cell_size = max(nozzle_diameter, 1e-6) * cell_size_multiplier
        self._segments: list[_PrintedSegment] = []
        self._grid: dict[tuple[int, int, int], list[int]] = {}

    def append_segment(self, start: np.ndarray, end: np.ndarray) -> None:
        aabb = _segment_aabb(start, end, self._radius)
        segment = _PrintedSegment(start=start, end=end, radius=self._radius, aabb=aabb)
        index = len(self._segments)
        self._segments.append(segment)
        for cell in _cell_coords_for_aabb(aabb, self._cell_size):
            self._grid.setdefault(cell, []).append(index)

    def candidates_near(self, aabb: AABB, margin: float) -> list[_PrintedSegment]:
        dilated = AABB(min_corner=aabb.min_corner - margin, max_corner=aabb.max_corner + margin)
        seen: set[int] = set()
        for cell in _cell_coords_for_aabb(dilated, self._cell_size):
            for index in self._grid.get(cell, []):
                seen.add(index)
        return [self._segments[i] for i in seen]


def _grow_printed_index_to(
    index: _PrintedGeometryIndex, points: list[ToolpathPoint], up_to_point_idx: int, last_grown_idx: int
) -> int:
    """
    Ajoute au plus une capsule par paire consecutive de points en extrusion
    (`extrusion > 0`) strictement entre `last_grown_idx` et `up_to_point_idx`
    EXCLU - jamais la capsule qui se termine AU point actuellement teste.

    Cette exclusion n'est pas seulement "jamais en avance" (l'invariant
    central de l'etat incremental) : elle evite aussi un faux positif
    structurel. La buse, au point ou elle vient de deposer de la matiere,
    est par definition quasi confondue avec le cordon qu'elle est en train
    d'extruder - le tester contre CE cordon la declasserait a tort comme
    collision a chaque point en extrusion. Le cordon qui se termine au
    point teste devient une collision possible seulement au(x) point(s)
    suivant(s), une fois que la buse s'en est reellement eloignee.
    """
    grown = last_grown_idx
    for i in range(max(last_grown_idx, 1), up_to_point_idx):
        if points[i].extrusion > 0.0 and points[i - 1].extrusion > 0.0:
            index.append_segment(points[i - 1].position, points[i].position)
        grown = i
    return grown


# ---------------------------------------------------------------------------
# Detection de collision (question 3 du plan)
# ---------------------------------------------------------------------------

def _toolhead_aabb(vertices: np.ndarray) -> AABB:
    return AABB(min_corner=vertices.min(axis=0), max_corner=vertices.max(axis=0))


def _aabb_overlaps(a: AABB, b: AABB, margin: float) -> bool:
    return bool(
        np.all(a.min_corner - margin <= b.max_corner)
        and np.all(b.min_corner - margin <= a.max_corner)
    )


def _point_to_segment_distance(point: np.ndarray, start: np.ndarray, end: np.ndarray) -> float:
    """Distance fermee point-segment (projection scalaire clampee sur [0,1]), sans bibliotheque."""
    segment = end - start
    length_sq = float(np.dot(segment, segment))
    if length_sq < 1e-18:
        return float(np.linalg.norm(point - start))
    t = float(np.clip(np.dot(point - start, segment) / length_sq, 0.0, 1.0))
    closest = start + t * segment
    return float(np.linalg.norm(point - closest))


def _points_to_segment_distances(points: np.ndarray, start: np.ndarray, end: np.ndarray) -> np.ndarray:
    """
    Version vectorisee de _point_to_segment_distance : distance de CHAQUE
    ligne de `points` (V,3) au segment [start,end] en une seule passe
    numpy, au lieu d'une boucle Python par sommet.

    Necessaire en pratique, pas seulement en theorie : avec ~1500
    points de test par tete positionnee (voir toolhead_sample_spacing) et
    plusieurs dizaines de capsules candidates par point de trajectoire
    teste, la boucle Python pure (`for vertex in toolhead_vertices`)
    domine largement le temps total de l'etape - mesure concretement en
    profilant main.py sur un modele trivial (24 faces) avec le modele de
    tete reel : plusieurs minutes pour quelques milliers de points de
    trajectoire, alors que la meme operation vectorisee est de l'ordre de
    la seconde. Le calcul lui-meme est inchange (meme formule que
    _point_to_segment_distance, projection scalaire clampee), seulement
    applique a tout le tableau de points a la fois.
    """
    segment = end - start
    length_sq = float(np.dot(segment, segment))
    if length_sq < 1e-18:
        return np.linalg.norm(points - start, axis=1)
    t = np.clip((points - start) @ segment / length_sq, 0.0, 1.0)
    closest = start[None, :] + t[:, None] * segment[None, :]
    return np.linalg.norm(points - closest, axis=1)


def _check_printed_part(
    toolhead_vertices: np.ndarray, printed_index: _PrintedGeometryIndex, safety_margin: float
) -> Optional[float]:
    """
    Plus grande penetration (safety_margin - distance minimale) contre la
    matiere deja imprimee, ou None si aucune capsule n'est assez proche.
    Broad-phase AABB via l'index avant tout calcul de distance precis.
    """
    aabb = _toolhead_aabb(toolhead_vertices)
    candidates = printed_index.candidates_near(aabb, safety_margin)
    if not candidates:
        return None

    min_distance = np.inf
    for segment in candidates:
        if not _aabb_overlaps(aabb, segment.aabb, safety_margin):
            continue
        distances = _points_to_segment_distances(toolhead_vertices, segment.start, segment.end)
        distances -= segment.radius
        candidate_min = float(distances.min())
        if candidate_min < min_distance:
            min_distance = candidate_min

    if min_distance == np.inf or min_distance >= safety_margin:
        return None
    return safety_margin - min_distance


def _check_static_obstacle(
    toolhead_vertices: np.ndarray,
    obstacle_mesh: Optional[trimesh.Trimesh],
    safety_margin: float,
) -> Optional[float]:
    """
    Penetration contre un obstacle statique (plateau, mandrin) via
    ProximityQuery - jamais mesh_data.kdtree (limite #15 : distance a un
    centre de face, pas a la surface), meme pattern deja etabli par step7.
    """
    if obstacle_mesh is None or len(obstacle_mesh.vertices) == 0:
        return None

    obstacle_aabb = AABB(
        min_corner=np.asarray(obstacle_mesh.bounds[0]), max_corner=np.asarray(obstacle_mesh.bounds[1])
    )
    toolhead_aabb = _toolhead_aabb(toolhead_vertices)
    if not _aabb_overlaps(toolhead_aabb, obstacle_aabb, safety_margin):
        return None

    proximity = trimesh.proximity.ProximityQuery(obstacle_mesh)
    _closest, distances, _triangle_ids = proximity.on_surface(toolhead_vertices)
    min_distance = float(np.min(distances))

    if min_distance >= safety_margin:
        return None
    return safety_margin - min_distance


PLATFORM_RADIUS = 150.0      # mm, ordre de grandeur d'un plateau, faute de mieux fourni
PLATFORM_HALF_HEIGHT = 0.5   # mm, moitie de l'epaisseur du disque (cylinder(height=1.0))


def _check_platform_disc(
    toolhead_vertices: np.ndarray,
    platform_origin: Optional[np.ndarray],
    safety_margin: float,
) -> Optional[float]:
    """
    Penetration contre le plateau, calculee en FORME FERMEE (distance
    point-disque vectorisee) plutot que par `ProximityQuery.on_surface` sur
    un maillage de disque.

    Meme raisonnement que la limite #27 (boucle point-capsule vectorisee) :
    le plateau n'est pas une geometrie arbitraire recue de l'exterieur, c'est
    une primitive que CE fichier synthetise lui-meme
    (`_platform_mesh_at_angle`). Lui appliquer une requete de proximite
    maillage generique fait payer le prix d'un cas general dont on n'a pas
    besoin. Cout mesure sur le modele de tete reel (1540 points de test) :
    53.02 ms par appel via `on_surface`, contre une fraction de milliseconde
    ici - et cet appel a lieu une fois PAR CANDIDAT de reorientation a
    l'etape 10 (jusqu'a 144 par point en collision), ou il representait 69%
    du temps total mesure.

    Precision : la forme fermee decrit le disque EXACT, la ou le maillage
    etait un polygone inscrit a 32 cotes - la version analytique est donc
    legerement plus grande pres du bord, et plus fidele a un vrai plateau
    circulaire, jamais moins conservatrice.

    Le disque est suppose d'axe Z et centre sur `platform_origin` (exactement
    ce que construit `_platform_mesh_at_angle`, dont la geometrie ne depend
    d'aucun angle - voir sa docstring).
    """
    if platform_origin is None:
        return None

    local = toolhead_vertices - np.asarray(platform_origin, dtype=float)
    radial = np.hypot(local[:, 0], local[:, 1])
    # Distance a un cylindre plein : ecart radial hors du rayon, ecart axial
    # hors de la demi-epaisseur, combines - nulle a l'interieur (penetration).
    dr = np.maximum(radial - PLATFORM_RADIUS, 0.0)
    dz = np.maximum(np.abs(local[:, 2]) - PLATFORM_HALF_HEIGHT, 0.0)
    min_distance = float(np.hypot(dr, dz).min())

    if min_distance >= safety_margin:
        return None
    return safety_margin - min_distance


def _platform_mesh_at_angle(
    machine_profile: MachineProfile, axis_angles: list[float]
) -> Optional[trimesh.Trimesh]:
    """
    Plateau approxime par un disque plat positionne selon l'angle courant
    du premier axe de rotation (si applicable) - aucune geometrie de
    plateau reelle n'est fournie a ce stade (contrairement a nozzle_model),
    donc une primitive simple plutot qu'une hypothese de forme non verifiee.
    None si aucune information de plateau n'est disponible (pas d'origine,
    pas d'axe de rotation) : le test est alors simplement saute (pas
    d'evenement produit), sans faux positif ni exception.
    """
    if machine_profile.platform_origin is None or not machine_profile.rotation_axes:
        return None

    disc = trimesh.creation.cylinder(
        radius=PLATFORM_RADIUS, height=PLATFORM_HALF_HEIGHT * 2.0, sections=32
    )
    disc.apply_translation(np.asarray(machine_profile.platform_origin, dtype=float))
    return disc


# ---------------------------------------------------------------------------
# Orchestration par point
# ---------------------------------------------------------------------------

def _evaluate_point(
    point_idx: int,
    toolhead_base_samples: np.ndarray,
    point: ToolpathPoint,
    machine_profile: MachineProfile,
    printed_index: _PrintedGeometryIndex,
    safety_margin: float,
) -> list[CollisionEvent]:
    events: list[CollisionEvent] = []

    toolhead_samples = _positioned_toolhead_vertices(
        toolhead_base_samples, point.position, point.axis_angles, machine_profile.rotation_axes
    )

    printed_penetration = _check_printed_part(toolhead_samples, printed_index, safety_margin)
    if printed_penetration is not None:
        events.append(
            CollisionEvent(
                point_index=point_idx, obstacle_type="printed_part", penetration_depth=printed_penetration
            )
        )

    # Forme fermee plutot que ProximityQuery sur un maillage de disque : le
    # plateau est une primitive synthetisee ici meme, pas une geometrie
    # arbitraire (voir _check_platform_disc). `_check_static_obstacle` reste
    # la voie pour tout VRAI maillage d'obstacle (mandrin/mors le jour ou une
    # geometrie sera fournie).
    platform_penetration = None
    if machine_profile.platform_origin is not None and machine_profile.rotation_axes:
        platform_penetration = _check_platform_disc(
            toolhead_samples, machine_profile.platform_origin, safety_margin
        )
    if platform_penetration is not None:
        events.append(
            CollisionEvent(
                point_index=point_idx, obstacle_type="platform", penetration_depth=platform_penetration
            )
        )

    return events


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

    sample_indices = _sample_point_indices(toolpath.points, params.path_sampling_step)
    printed_index = _PrintedGeometryIndex(
        machine_profile.nozzle_diameter, params.grid_cell_size_multiplier
    )

    events: list[CollisionEvent] = []
    last_grown_idx = 0
    for point_idx in sample_indices:
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
            )
        )

    return CollisionReport(events=events, chuck_checked=False, nozzle_model_checked=True)
