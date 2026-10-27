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
    # Rayon (mm, repere local du modele de tete) autour de la pointe de buse
    # (origine locale) a l'interieur duquel les points d'echantillonnage sont
    # EXCLUS des tests de collision. La pointe de buse est par construction
    # toujours au contact du cordon qu'elle vient de deposer et proche du
    # plateau sur les premieres couches - sans cette zone de garde, ce seul
    # point declenche une collision a CHAQUE point de trajectoire teste,
    # sans exception (mesure sur l_shape_overhang.stl : exactement 2
    # evenements par point echantillonne, printed_part + platform - pas un
    # motif de collision reel, un artefact structurel du test). A l'echelle
    # de la buse/marge de securite (quelques mm), PAS de nozzle_diameter
    # (echelle du cordon, meme piege que toolhead_sample_spacing ci-dessus,
    # limites #21/#26). Defaut non nul : ce champ corrige un bug de faux
    # positif systematique, pas une fonctionnalite optionnelle a activer.
    nozzle_contact_clearance: float = 3.0
    # Controle RAPIDE (opt-in, porte de la recherche le 2026-09-30) : la matiere
    # imprimee est densifiee le long des cordons (pas `fast_densify_step`) dans
    # UN arbre k-d construit une fois, chaque point densifie portant l'indice
    # du point de trajectoire qui le depose ; un point teste ne voit que la
    # matiere deposee AVANT lui (meme invariant que _grow_printed_index_to).
    # Remplace la boucle de capsules, dont le debit s'effondre avec la matiere
    # accumulee (limite #27) : 30 min -> ~1 min sur le Benchy en recherche.
    # Densification : la distance a l'axe du cordon est SURESTIMEE d'au plus
    # un demi-pas (0,125 mm) -- a couvrir par safety_margin (1 mm).
    fast_kdtree: bool = False
    fast_densify_step: float = 0.25


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

    # longueurs de segments vectorisees (2026-10-01 : 400 000 appels a
    # np.linalg.norm prenaient 1 s sur death). `norm(axis=1)` peut differer
    # d'un ULP de l'appel unitaire (BLAS ddot) : un seuil franchi pile a un
    # ULP pres peut deplacer un echantillon, sans consequence pour un
    # echantillonnage tous les `path_sampling_step` mm.
    positions = np.array([p.position for p in points], dtype=float)
    lengths = np.linalg.norm(np.diff(positions, axis=0), axis=1).tolist()
    indices = [0]
    accumulated = 0.0
    for i in range(1, len(points)):
        accumulated += lengths[i - 1]
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


def _positioned_toolhead_vertices_batch(
    base_vertices: np.ndarray, positions: np.ndarray, axis_angles: np.ndarray,
    rotation_axes: list[RotationAxisSpec],
) -> np.ndarray:
    """(P,S,3) : `_positioned_toolhead_vertices` pour P poses d'un coup, par
    matrices de rotation composees dans le MEME ordre (axe le plus interne en
    premier). Ecart a la version pose par pose : quelques ULP (produit
    matriciel au lieu de Rodrigues), sans consequence pour un controle de
    distance a 0,125 mm pres. Cree pour step9 rapide, ou le positionnement
    pose par pose prenait 14,5 s sur 20 (death @ 0,8)."""
    positions = np.asarray(positions, dtype=float)
    count = len(positions)
    total = np.broadcast_to(np.eye(3), (count, 3, 3)).copy()
    angles = np.asarray(axis_angles, dtype=float).reshape(count, len(rotation_axes))
    for k, spec in enumerate(rotation_axes):
        axis = np.asarray(spec.axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        theta = np.radians(angles[:, k])
        c, s_ = np.cos(theta)[:, None, None], np.sin(theta)[:, None, None]
        cross = np.array([[0.0, -axis[2], axis[1]], [axis[2], 0.0, -axis[0]], [-axis[1], axis[0], 0.0]])
        total = total @ (c * np.eye(3) + s_ * cross + (1.0 - c) * np.outer(axis, axis))
    return np.matmul(base_vertices[None, :, :], total.transpose(0, 2, 1)) + positions[:, None, :]


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


def _filter_nozzle_contact_zone(samples: np.ndarray, clearance: float) -> np.ndarray:
    """
    Retire de `samples` (repere LOCAL du modele de tete) les points a moins
    de `clearance` de l'origine locale - l'origine EST la pointe de buse
    (convention verifiee, voir docstring du module), zone de contact
    attendue par construction avec la matiere/le plateau, pas un obstacle.

    `clearance <= 0` desactive le filtre (tableau inchange) - comportement
    historique explicitement accessible, pas seulement par defaut absent.

    Si TOUS les points tombent dans la zone de garde (tete plus petite que
    `clearance` - un cas synthetique/degenere, jamais rencontre sur le
    modele de tete reel dont le corps monte a ~77mm), le filtre est annule
    et le tableau d'origine est retourne intact : mieux vaut re-tester la
    zone de contact que de ne plus rien tester du tout (un tableau vide
    ferait planter `_toolhead_aabb` en aval sur `.min(axis=0)`, et
    desactiverait silencieusement toute detection de collision pour cette
    tete - pire que le faux positif qu'on corrige).

    Applique UNE FOIS sur l'echantillonnage de base (avant positionnement
    par point de trajectoire) : le mesh de tete et son echantillonnage sont
    invariants le long de la trajectoire, refaire ce filtre par point
    repeterait un cout inutile pour un resultat identique a chaque fois.
    """
    if clearance <= 0.0 or len(samples) == 0:
        return samples
    distances = np.linalg.norm(samples, axis=1)
    filtered = samples[distances >= clearance]
    if len(filtered) == 0:
        return samples
    return filtered


# ---------------------------------------------------------------------------
# Geometrie imprimee incrementale (question 1 du plan)
# ---------------------------------------------------------------------------

#: Tableau `(0,3)` rendu par `candidates_endpoints` quand aucun candidat n'est
#: trouve - la forme importe, les appelants indexent la 2e dimension.
_EMPTY_ENDPOINTS = np.zeros((0, 3), dtype=float)


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
        # Noyau natif (portage C++ en cours) : `None` = chemin Python
        # historique, strictement inchange. Voir core/_cpp.py pour
        # l'interrupteur et pourquoi le Python reste le defaut.
        self._native = None
        if _cpp.enabled():
            self._native = _cpp.module().PrintedGrid(
                cell_size=self._cell_size, radius=self._radius, mode=_cpp.grid_mode()
            )

    def append_segment(self, start: np.ndarray, end: np.ndarray) -> None:
        if self._native is not None:
            # AUCUN objet Python n'est cree ici. Le natif detient les
            # extremites ; les `_PrintedSegment` sont construits a la demande
            # dans `candidates_near`, pour les ~690 candidats reellement
            # rendus et non pour les 744 324 capsules inserees. Materialiser
            # a l'insertion coutait 1,79 s de `_segment_aabb` mesurees, pour
            # des objets dont la quasi-totalite n'est jamais lue.
            self._native.append(start, end)
            return
        aabb = _segment_aabb(start, end, self._radius)
        segment = _PrintedSegment(start=start, end=end, radius=self._radius, aabb=aabb)
        index = len(self._segments)
        self._segments.append(segment)
        for cell in _cell_coords_for_aabb(aabb, self._cell_size):
            self._grid.setdefault(cell, []).append(index)

    def append_segments(self, starts: np.ndarray, ends: np.ndarray) -> None:
        """
        Insertion en LOT de N capsules.

        Le chemin historique traversait la frontiere Python/C++ une fois par
        capsule — 744 324 appels sur death.stl pour une operation triviale.
        La granularite de l'appel coutait plus cher que le calcul, exactement
        comme pour le test de support (voir `has_support`).
        """
        if len(starts) == 0:
            return
        if self._native is not None:
            self._native.append_batch(
                np.asarray(starts, dtype=float), np.asarray(ends, dtype=float)
            )
            return
        for i in range(len(starts)):
            self.append_segment(starts[i], ends[i])

    def candidates_near(self, aabb: AABB, margin: float) -> list[_PrintedSegment]:
        if self._native is not None:
            starts, ends = self._native.candidates_endpoints(
                aabb.min_corner - margin, aabb.max_corner + margin
            )
            radius = self._radius
            return [
                _PrintedSegment(
                    start=starts[i], end=ends[i], radius=radius,
                    aabb=_segment_aabb(starts[i], ends[i], radius),
                )
                for i in range(len(starts))
            ]
        dilated = AABB(min_corner=aabb.min_corner - margin, max_corner=aabb.max_corner + margin)
        seen: set[int] = set()
        for cell in _cell_coords_for_aabb(dilated, self._cell_size):
            for index in self._grid.get(cell, []):
                seen.add(index)
        return [self._segments[i] for i in seen]

    def has_support(
        self,
        aabb: AABB,
        position: np.ndarray,
        directions: np.ndarray,
        horizontal_radius: float,
        vertical_height: float,
    ) -> bool:
        """
        `position` a-t-il un appui, dans la matiere indexee, le long d'AU MOINS
        UNE des `directions` (D,3) ? La recherche est bornee a `aabb`.

        FUSIONNE la requete spatiale et le test geometrique. C'est le point
        d'entree qui rapporte le plus du portage : le chemin historique
        materialisait deux tableaux `(C,3)` (~690 candidats en mediane) puis
        executait ~15 operations numpy dessus pour ne rendre qu'un BOOLEEN,
        sans jamais court-circuiter. Mesure sur death.stl : 11,97 s dans
        `point_has_support_along_direction` + ~10 s de dispatch numpy + 3,95 s
        d'allocation des tableaux, pour 273 146 appels.

        Le chemin natif lit les extremites directement dans la grille et
        s'arrete au PREMIER appui trouve. L'arret anticipe ne change pas le
        resultat : `np.any(...)` est vrai des qu'un candidat satisfait les
        conditions.

        Le chemin Python reproduit exactement le comportement historique.
        """
        if self._native is not None:
            return self._native.has_support(
                aabb.min_corner, aabb.max_corner, position, directions,
                horizontal_radius, vertical_height,
            )

        starts, ends = self.candidates_endpoints(aabb, margin=0.0)
        if len(starts) == 0:
            return False
        for direction in directions:
            if point_has_support_along_direction(
                position, direction, starts, ends, horizontal_radius, vertical_height
            ):
                return True
        return False

    def candidates_endpoints(
        self, aabb: AABB, margin: float
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Extremites des capsules candidates, en deux tableaux `(C,3)`.

        Chemin dedie aux trois consommateurs qui reconstruisaient exactement
        ces tableaux par comprehension sur les objets rendus par
        `candidates_near` (`step9b._has_support_below`,
        `step9b._has_support_in_any_direction`, `step5._layer_support_mask`).
        C'est la source des 41 s de `numpy.array` mesurees au profilage de
        death.stl : 817 653 requetes pour ~432 candidats en mediane, soit deux
        comprehensions Python et deux allocations a chaque appel.

        Le chemin natif les remplit directement ; le chemin Python reproduit
        le comportement historique a l'identique, pour que la sortie soit la
        meme quel que soit l'interrupteur.
        """
        if self._native is not None:
            return self._native.candidates_endpoints(
                aabb.min_corner - margin, aabb.max_corner + margin
            )
        candidates = self.candidates_near(aabb, margin)
        if not candidates:
            return _EMPTY_ENDPOINTS, _EMPTY_ENDPOINTS
        starts = np.array([segment.start for segment in candidates], dtype=float)
        ends = np.array([segment.end for segment in candidates], dtype=float)
        return starts, ends


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
# SURFACE du plateau a z = platform_origin.z (corrige le 2026-10-06) : les pieces posent a
# z = platform_z (step3) et l'export G-code compte le plateau en z = 0. Le disque
# etait jusque-la CENTRE sur l'origine, donc sa surface etait a +0,5 mm et step9 etait 0,5 mm plus
# strict que la marge annoncee (limite #70). Le centre du disque est desormais a origine - 0,5 mm.
PLATFORM_CENTER_OFFSET = -PLATFORM_HALF_HEIGHT


def _check_platform_disc(
    toolhead_vertices: np.ndarray,
    platform_origin: Optional[np.ndarray],
    safety_margin: float,
    platform_radius: Optional[float] = None,
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
    radius = PLATFORM_RADIUS if platform_radius is None else float(platform_radius)
    dr = np.maximum(radial - radius, 0.0)
    dz = np.maximum(np.abs(local[:, 2] - PLATFORM_CENTER_OFFSET) - PLATFORM_HALF_HEIGHT, 0.0)
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
    disc.apply_translation(np.asarray(machine_profile.platform_origin, dtype=float)
                           + np.array([0.0, 0.0, PLATFORM_CENTER_OFFSET]))
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
            toolhead_samples, machine_profile.platform_origin, safety_margin,
            getattr(machine_profile, "platform_radius", None),
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

def _densified_printed_material(points: list[ToolpathPoint], step: float) -> tuple[np.ndarray, np.ndarray]:
    """Points densifies le long de chaque cordon (meme regle que
    _grow_printed_index_to : segment i-1 -> i quand les deux points sont en
    extrusion) et indice `i` du point qui termine ce cordon."""
    pos = np.array([p.position for p in points], dtype=float)
    ext = np.array([p.extrusion for p in points], dtype=float)
    if len(pos) < 2:
        return np.zeros((0, 3)), np.zeros(0, dtype=int)
    ends = np.flatnonzero((ext[1:] > 0.0) & (ext[:-1] > 0.0)) + 1
    if len(ends) == 0:
        return np.zeros((0, 3)), np.zeros(0, dtype=int)
    a, b = pos[ends - 1], pos[ends]
    counts = np.maximum(1, np.ceil(np.linalg.norm(b - a, axis=1) / step).astype(int))
    seg = np.repeat(np.arange(len(ends)), counts + 1)
    # meme t que np.arange(c + 1) / c par segment, sans boucle Python
    offset = np.arange(len(seg)) - np.repeat(np.cumsum(counts + 1) - (counts + 1), counts + 1)
    t = offset / np.repeat(counts, counts + 1)
    dense = a[seg] + t[:, None] * (b[seg] - a[seg])
    return dense, ends[seg]


class FastCollisionChecker:
    """Controle rapide (voir CollisionSimParams.fast_kdtree), reutilisable
    point par point -- par `run` et par la reorientation de l'etape 10, qui
    re-teste d'autres angles au MEME point (la matiere ne depend que des
    positions, jamais des angles).

    Par point teste : les K plus proches points de matiere de chaque sommet de
    tete dans le rayon `reach`, filtres a la matiere deja deposee. Si les K
    sont tous dans le rayon (liste saturee) sans matiere passee, la requete
    exacte par boule est refaite pour ce sommet : la tete survole souvent de
    la matiere FUTURE, qui ne doit masquer aucune matiere passee.
    scipy directement (k > 1, que l'enveloppe core.spatial n'expose pas)."""

    K_NEAR = 16

    def __init__(self, toolpath: ToolpathMachine, machine_profile: MachineProfile, params: CollisionSimParams,
                 toolhead_base_samples: np.ndarray):
        from scipy.spatial import cKDTree

        self.machine, self.params, self.samples = machine_profile, params, toolhead_base_samples
        self.reach = params.safety_margin + machine_profile.nozzle_diameter / 2.0
        self.dense, self.owner = _densified_printed_material(toolpath.points, params.fast_densify_step)
        self._tree = None       # arbre scipy construit a la demande : inutile en mode natif (2026-10-06)
        self.owner_pad = np.r_[self.owner, np.iinfo(np.int64).max]   # indice "absent" de cKDTree
        # Pas 19 du portage C++ : meme minimum exact, en natif, avec elagage
        # des sous-arbres entierement futurs (CLAUDE.md limite #56).
        self.native = None
        if len(self.dense) and _cpp.enabled() and hasattr(_cpp.module(), "PastMaterialIndex"):
            self.native = _cpp.module().PastMaterialIndex(
                np.ascontiguousarray(self.dense, dtype=float), np.ascontiguousarray(self.owner, dtype=np.int64))

    @property
    def tree(self):
        if self._tree is None and len(self.dense):
            from scipy.spatial import cKDTree
            self._tree = cKDTree(self.dense)
        return self._tree

    def penetrations(self, point_idx: int, position: np.ndarray, axis_angles: list[float], verts=None):
        """(penetration plateau ou None, penetration piece ou None). `verts` :
        tete deja positionnee (calcul par lots de `_run_fast_kdtree`)."""
        if verts is None:
            verts = _positioned_toolhead_vertices(self.samples, position, axis_angles, self.machine.rotation_axes)
        part = None
        if self.native is not None:
            best = self.native.min_past_distance(np.ascontiguousarray(verts, dtype=float), int(point_idx),
                                                 float(self.reach))
            if best < self.reach:
                part = self.reach - best
        elif self.tree is not None:
            k = min(self.K_NEAR, len(self.dense))
            dist, idx = self.tree.query(verts, k=k, distance_upper_bound=self.reach)
            dist, idx = dist.reshape(len(verts), k), idx.reshape(len(verts), k)
            past = self.owner_pad[idx] < point_idx
            d_past = np.where(past, dist, np.inf).min(axis=1)
            saturated = np.isfinite(dist[:, -1]) & ~np.isfinite(d_past)
            for vi in np.flatnonzero(saturated):
                near = np.asarray(self.tree.query_ball_point(verts[vi], self.reach), dtype=int)
                near = near[self.owner[near] < point_idx]
                if len(near):
                    d_past[vi] = float(np.linalg.norm(self.dense[near] - verts[vi], axis=1).min())
            best = float(d_past.min())
            if best < self.reach:
                part = self.reach - best
        plate = None
        if self.machine.platform_origin is not None and self.machine.rotation_axes:
            plate = _check_platform_disc(verts, self.machine.platform_origin, self.params.safety_margin,
                                         getattr(self.machine, "platform_radius", None))
        return plate, part


def _platform_disc_batch(verts_block: np.ndarray, machine_profile: MachineProfile,
                         params: CollisionSimParams) -> list:
    """`_check_platform_disc` pour un bloc de poses (P,S,3) : memes operations
    element par element, minimum par pose - memes valeurs que l'appel unitaire."""
    if machine_profile.platform_origin is None or not machine_profile.rotation_axes:
        return [None] * len(verts_block)
    margin = params.safety_margin
    origin = np.asarray(machine_profile.platform_origin, dtype=float)
    # Raccourci EXACT : distance >= ecart axial ; une pose dont tous les points
    # sont a plus de `margin` au-dessus (ou au-dessous) de l'epaisseur du
    # plateau ne peut pas le toucher. Seules les autres paient le calcul
    # complet (2 s sur 2 sur death sans ce tri : 86 M d'hypot).
    dz_all = np.maximum(np.abs(verts_block[..., 2] - origin[2] - PLATFORM_CENTER_OFFSET) - PLATFORM_HALF_HEIGHT, 0.0)
    near = np.flatnonzero(dz_all.min(axis=1) < margin)
    result: list = [None] * len(verts_block)
    if len(near) == 0:
        return result
    local = verts_block[near] - origin
    radial = np.hypot(local[..., 0], local[..., 1])
    radius_value = getattr(machine_profile, "platform_radius", None)
    radius = PLATFORM_RADIUS if radius_value is None else float(radius_value)
    dr = np.maximum(radial - radius, 0.0)
    dz = np.maximum(np.abs(local[..., 2] - PLATFORM_CENTER_OFFSET) - PLATFORM_HALF_HEIGHT, 0.0)
    minimum = np.hypot(dr, dz).min(axis=1)
    for k, m in zip(near, minimum):
        result[k] = None if float(m) >= margin else margin - float(m)
    return result


_FAST_CHECKER_CACHE: dict = {}


def fast_checker_for(toolpath: ToolpathMachine, machine_profile: MachineProfile, params: CollisionSimParams,
                     toolhead_base_samples: np.ndarray) -> "FastCollisionChecker":
    """FastCollisionChecker REUTILISE tant que la matiere est la meme (2026-10-06). La matiere ne
    depend que des positions et des extrusions, jamais des angles : step9 puis les iterations de
    step10 (qui ne changent que des angles) partagent le meme controleur au lieu de le
    reconstruire (5 constructions, 7,3 s sur death). Cache a une entree, compare par valeur
    (positions, extrusions, echantillons de tete, parametres) et par identite du profil machine."""
    pos = np.array([p.position for p in toolpath.points], dtype=float)
    ext = np.array([p.extrusion for p in toolpath.points], dtype=float)
    samples = np.asarray(toolhead_base_samples, dtype=float)
    hit = _FAST_CHECKER_CACHE.get("entry")
    if (hit is not None and hit[0] is machine_profile and hit[1] == params
            and hit[2].shape == pos.shape and np.array_equal(hit[2], pos) and np.array_equal(hit[3], ext)
            and hit[4].shape == samples.shape and np.array_equal(hit[4], samples)):
        return hit[5]
    checker = FastCollisionChecker(toolpath, machine_profile, params, samples)
    _FAST_CHECKER_CACHE["entry"] = (machine_profile, params, pos, ext, samples, checker)
    return checker


def fast_batch_events(checker: "FastCollisionChecker", toolpath: ToolpathMachine,
                      indices: list[int]) -> Optional[list[CollisionEvent]]:
    """Evenements du controle rapide pour une liste de points, en UN appel natif parallele
    (pose de la tete, plateau, matiere passee). None si le noyau natif ne s'applique pas
    (l'appelant retombe sur la boucle Python). Meme ordre que la boucle : par point,
    piece puis plateau."""
    if checker.native is None or not _cpp.enabled() or not hasattr(_cpp.module(), "fast_penetrations_batch"):
        return None
    axes = checker.machine.rotation_axes or []
    if not axes or not indices:
        return None
    pts = toolpath.points
    positions = np.array([pts[i].position for i in indices], dtype=float)
    angles = np.array([pts[i].axis_angles for i in indices], dtype=float).reshape(len(indices), len(axes))
    units = np.array([np.asarray(a.axis, dtype=float) / np.linalg.norm(a.axis) for a in axes])
    origin = checker.machine.platform_origin
    radius = getattr(checker.machine, "platform_radius", None)
    plate, part = _cpp.module().fast_penetrations_batch(
        checker.native, np.ascontiguousarray(checker.samples, dtype=float), np.ascontiguousarray(units),
        np.ascontiguousarray(positions), np.ascontiguousarray(angles), np.asarray(indices, dtype=np.int64),
        float(checker.reach), None if origin is None else np.asarray(origin, dtype=float),
        float(PLATFORM_RADIUS if radius is None else radius), float(PLATFORM_HALF_HEIGHT),
        float(PLATFORM_CENTER_OFFSET), float(checker.params.safety_margin))
    events: list[CollisionEvent] = []
    for k, point_idx in enumerate(indices):
        if not np.isnan(part[k]):
            events.append(CollisionEvent(point_index=point_idx, obstacle_type="printed_part",
                                         penetration_depth=float(part[k])))
        if not np.isnan(plate[k]):
            events.append(CollisionEvent(point_index=point_idx, obstacle_type="platform",
                                         penetration_depth=float(plate[k])))
    return events


def _run_fast_kdtree(
    toolpath: ToolpathMachine,
    machine_profile: MachineProfile,
    params: CollisionSimParams,
    toolhead_base_samples: np.ndarray,
    sample_indices: list[int],
) -> list[CollisionEvent]:
    """Variante rapide de la boucle de `run` (voir FastCollisionChecker)."""
    checker = fast_checker_for(toolpath, machine_profile, params, toolhead_base_samples)
    native_events = fast_batch_events(checker, toolpath, list(sample_indices))
    if native_events is not None:
        return native_events
    events: list[CollisionEvent] = []
    axes = machine_profile.rotation_axes or []
    chunk = 1000
    for start in range(0, len(sample_indices), chunk):
        block = sample_indices[start:start + chunk]
        pts = [toolpath.points[i] for i in block]
        verts_block = _positioned_toolhead_vertices_batch(
            toolhead_base_samples, np.array([p.position for p in pts], dtype=float),
            np.array([p.axis_angles for p in pts], dtype=float), axes,
        ) if axes else None
        batch = None
        if verts_block is not None and checker.native is not None and hasattr(checker.native, "min_past_distance_batch"):
            # tout le bloc d'un coup, en parallele cote natif : meme calcul pose
            # par pose que `penetrations`, donc memes evenements au bit
            best = checker.native.min_past_distance_batch(
                np.ascontiguousarray(verts_block, dtype=float), np.asarray(block, dtype=np.int64), float(checker.reach))
            plates = _platform_disc_batch(verts_block, machine_profile, params)
            batch = (best, plates)
        for k, point_idx in enumerate(block):
            point = pts[k]
            verts = verts_block[k] if verts_block is not None else None
            if batch is not None:
                best_k = float(batch[0][k])
                part = checker.reach - best_k if best_k < checker.reach else None
                plate = batch[1][k]
            else:
                plate, part = checker.penetrations(point_idx, point.position, point.axis_angles, verts)
            if part is not None:
                events.append(CollisionEvent(point_index=point_idx, obstacle_type="printed_part",
                                             penetration_depth=part))
            if plate is not None:
                events.append(CollisionEvent(point_index=point_idx, obstacle_type="platform",
                                             penetration_depth=plate))
    return events


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
