"""
pipeline/step5_layers.py

Etape 5 du pipeline : extraction des couches courbes (isosurfaces) a partir
du ScalarField resolu par step4, sur le TetMesh de step1.

Depend de TetMesh + MeshData (step1), ScalarField (step4),
RegionDependencyGraph (step3).

Ce que ce fichier fait, dans l'ordre :
  1. Convertit `layer_height_base` (une distance physique) en un espacement
     de valeurs de champ, via la mediane ponderee par volume de ||grad(phi)||
     sur tout le TetMesh (une seule valeur representative - un champ non
     eikonal ne permet pas un espacement exactement constant partout, cf.
     Decision 2 du plan).
  2. Pour chaque isovaleur, extrait le patch triangule par marching
     tetrahedra (implementation maison - aucune bibliotheque du projet ne
     fournit cette operation sur un maillage tetraedrique arbitraire).
  3. En extrait les boucles de contour fermees via le trace de bord du patch
     (`trimesh.Trimesh.outline()`) : comme le TetMesh est conforme a la
     frontiere (tetgen appele avec le switch `p`), le bord du patch EST deja
     geometriquement l'intersection avec le mesh original - pas besoin d'une
     deuxieme operation geometrique explicite (cf. Decision 3 du plan).
  4. Associe chaque couche aux regions de step3 qu'elle traverse, et
     l'ordonne par rang topologique maximal (la contrainte la plus stricte
     parmi les regions traversees, cf. Decision 4) puis par isovaleur - ou,
     si `split_disjoint_clusters` est actif, par un ordre DFS de continuite
     spatiale a l'interieur de chaque rang (cf. Decision 8) : necessaire des
     qu'une meme region step3 (une nappe de surface connexe) couvre
     plusieurs composantes spatialement disjointes du solide (ex: pied et
     pointe de bras d'un T partageant la meme region "verticale" - sinon le
     tri par isovaleur seule les fait alterner a chaque couche, CLAUDE.md).

Sortie : liste ordonnee de Layer.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import networkx as nx
import numpy as np
import trimesh

from core.geometry_utils import (
    composed_rotation_directions,
    point_has_support_along_direction,
    point_in_polygon_2d,
)
from core.mesh_types import (
    AABB,
    Layer,
    MeshData,
    RegionDependencyGraph,
    ScalarField,
    TetMesh,
)
from core.spatial import KDTree


class LayerExtractionError(Exception):
    """Levee si le ScalarField ne peut pas etre converti en couches (gradient manquant, champ degenere)."""


@dataclass
class LayerParams:
    layer_height_base: float = 0.2
    # tolerance de fusion des sommets du patch triangule, avant extraction du
    # contour, relative a l'etendue du TetMesh - evite des micro-trous de
    # bord dus au bruit flottant plutot qu'a une vraie discontinuite
    # geometrique (cf. Decision 5 du plan).
    weld_tolerance_factor: float = 1e-6
    # Scinde une couche en plusieurs Layer distincts (meme isovalue,
    # region_ids et local_frame recalcules par sous-ensemble) quand ses
    # boucles de contour forment plusieurs composantes spatialement
    # disjointes ET non imbriquees (ex: jonction en T ou en etoile, ou une
    # isovaleur traverse simultanement le pied et l'extremite d'un bras
    # eloigne - CLAUDE.md). None (defaut) = DESACTIVE, comportement
    # historique inchange (un seul Layer par isovaleur, local_frame moyenne
    # sur toutes les boucles) - meme discipline que
    # direction_blend_radius=0.0 (step4), near_platform_resolution=None
    # (step1), max_region_merge_angle_deg=None (step3). Valeur = facteur
    # multiplicatif de layer_height_base utilise comme tolerance de
    # distance (cf. _cluster_contours) : deux boucles non imbriquees
    # distantes de moins de split_disjoint_clusters * layer_height_base
    # restent dans le meme cluster (evite une discontinuite artificielle
    # juste avant qu'elles ne fusionnent geometriquement a une isovaleur
    # voisine).
    #
    # Active AUSSI, automatiquement, le tri spatial par continuite (DFS sur
    # un graphe de proximite entre clusters d'isovaleurs successives,
    # respectant STRICTEMENT le rang topologique de region comme tri
    # primaire - jamais traverse, cf. Decision 8) : sans ce tri, deux
    # clusters partageant le meme rang ET les memes region_ids (typiquement
    # pied et pointe de bras d'un T, meme region topologique au sens step3
    # malgre une topologie SOLIDE separee) restent ordonnes uniquement par
    # isovaleur croissante, ce qui les fait alterner a chaque isovaleur
    # (pied->pointe->pied->pointe...) et produit un faisceau de trajets a
    # vide entre les deux a chaque couche - mesure directement sur
    # res/test_models/l-triple-axis-overhanged-shape.stl (CLAUDE.md). Pas de
    # parametre separe : ce tri n'a de sens/d'effet que si le clustering
    # produit deja plusieurs Layer par isovaleur, et degenere exactement
    # vers le tri historique (rank, isovalue) en leur absence - un flag
    # independant ajouterait une combinatoire sans benefice reel.
    split_disjoint_clusters: Optional[float] = None

    # Ordonne les couches par SUPPORT REEL plutot que par isovaleur croissante
    # (voir `_support_driven_order`). None (defaut) = DESACTIVE, tri historique
    # inchange - meme discipline que `split_disjoint_clusters` ci-dessus.
    #
    # Le tri historique `(rang de region, isovaleur)` traite le rang comme un
    # ordre TOTAL et, A L'INTERIEUR d'un rang, ne connait que l'isovaleur
    # croissante. Or l'isovaleur n'est PAS le support : une coque rattachee au
    # reste du modele UNIQUEMENT PAR LE HAUT voit ses isovaleurs basses passer
    # en premier, donc sa matiere deposee avant l'unique matiere qui la porte.
    # Ce mecanisme pose directement la bonne question, couche par couche : "si
    # j'imprime cette couche maintenant, ses points ont-ils de la matiere DEJA
    # PLACEE en dessous ?", avec exactement la geometrie de test qui sert deja
    # a COMPTER ces defauts en step9b (`point_has_support_along_direction`) -
    # pour DECIDER l'ordre au lieu de seulement CONSTATER les degats.
    #
    # Valeur = fraction minimale de points de contour devant etre supportes
    # pour qu'une couche soit ELIGIBLE (0.9 = 90%, valeur de reference
    # CLAUDE.md). Le support est la cle PRIMAIRE ; proximite a la couche
    # precedente puis altitude ne departagent que des couches DEJA jugees
    # imprimables - c'est ce qui protege les surplombs par construction (une
    # piece pendante est plus basse que son porteur, mais tant que le porteur
    # n'existe pas elle n'est pas supportee, donc jamais choisie).
    #
    # NON UNIVERSEL, et c'est pourquoi il reste opt-in : mesure CLAUDE.md,
    # -44% de defauts d'ordre sur une piece massive (death.stl) mais une
    # DEGRADATION sur une piece fine/ajouree (l-triple-axis, 113/141 couches
    # `support_unresolved`) - le test porte sur les CONTOURS seuls, l'infill
    # n'existant qu'a partir de step6, donc le support est SOUS-ESTIME des que
    # les contours de couches successives ne se recouvrent pas.
    support_driven_order: Optional[float] = None

    # --- ORDRE PAR DISTANCE GEODESIQUE AU PLATEAU -------------------------
    # Separe deux questions que le pipeline confond depuis l'origine, parce
    # qu'un seul champ scalaire y repond :
    #   - quelle FORME doit avoir une couche ? -> il faut un champ lisse et
    #     regulierement espace, c'est le role du FEM de step4 ;
    #   - dans quel ORDRE les imprimer ? -> il faut savoir par ou la matiere
    #     rejoint le plateau, ce qui est un probleme de PLUS COURT CHEMIN.
    # Rien n'oblige la meme quantite a bien repondre aux deux. Ici, la forme
    # reste celle des isosurfaces de phi (inchangee) et seul l'ORDRE est
    # re-derive de la distance geodesique au plateau A TRAVERS LA MATIERE
    # (Dijkstra multi-source sur les aretes du TetMesh, ~220ms sur death.stl).
    #
    # Pourquoi ce critere ne peut PAS reproduire l'echec de
    # `support_driven_order` : celui-ci posait une question BINAIRE et
    # gloutonne ("cette couche est-elle posee sur de la matiere deja
    # ancree ?"), qui detruit la continuite de proche en proche des surplombs
    # progressifs. La distance geodesique varie CONTINUMENT le long de la
    # matiere : une coque suspendue est ordonnee depuis son attache vers son
    # extremite, de facon monotone, sans saut possible par construction.
    #
    # False (defaut) = DESACTIVE, ordre historique strictement inchange.
    geodesic_order: bool = False

    # --- CROISSANCE DE FRONT (Dijkstra sur l'adjacence des couches) -------
    # Contrainte plus forte que `geodesic_order`, et c'est tout l'ecart :
    # celui-ci TRIE les couches par distance geodesique, ce qui n'interdit
    # nullement d'imprimer deux couches voisines en distance mais situees aux
    # antipodes du modele (mesure : x3,7 de degradation). Ici, une couche
    # n'est imprimable QUE si elle touche de la matiere DEJA imprimee.
    #
    # Parcours de Dijkstra sur le graphe d'adjacence des couches, amorce sur
    # la couche la plus proche du plateau : a chaque etape, on choisit parmi
    # les couches ATTEIGNABLES (adjacentes au deja-imprime) celle de plus
    # petite distance geodesique. Le front de matiere est connexe a chaque
    # instant, par construction.
    #
    # Consequence recherchee, impossible a exprimer avec l'ordre par isovaleur :
    # un pilier suspendu, relie au reste par un pont situe au-dessus de lui,
    # n'est atteignable qu'APRES ce pont - et se construit alors du haut vers
    # le bas. Mesure sur `real_overhang.stl` en ordre historique : ce pilier
    # demarre a 29% de l'impression, son point le plus bas a 6,87mm au-dessus
    # du plateau, sans rien dessous ni a cote.
    #
    # LIMITE ASSUMEE : une seule composante progresse a la fois. Une piece a
    # branches reellement simultanees (helice multi-pales) demanderait un front
    # multiple - non traite, et signale plutot que masque.
    #
    # False (defaut) = DESACTIVE, ordre historique strictement inchange.
    frontier_order: bool = False

    # Largeur de cordon, REQUISE des que `support_driven_order` est actif
    # (rayon des capsules de matiere + echelle des rayons de recherche). Aucun
    # defaut arbitraire : step5 ne connait pas `MachineProfile`, et deviner
    # cette valeur creerait une seconde source de verite divergeant en silence
    # de `MachineProfile.nozzle_diameter` - meme discipline que `platform_z`,
    # passe explicitement a step3/step9b.
    nozzle_width: Optional[float] = None

    # Rayons de recherche du support, en multiples de `nozzle_width` et de
    # `layer_height_base` - memes grandeurs et memes defauts que
    # `SupportCheckParams` (step9b), pour que la couche jugee imprimable ici et
    # le defaut compte la-bas partagent exactement le meme critere.
    support_horizontal_radius_factor: float = 2.5
    support_vertical_height_factor: float = 7.5

    # Tolerance "pose sur le plateau", en multiples de `nozzle_width` - meme
    # grandeur et meme defaut que `SupportCheckParams.platform_contact_tolerance_factor`
    # (step9b). Amorce l'ANCRAGE (voir `_support_driven_order`) sans avoir
    # besoin de `platform_z` : un point assez proche du Z minimal de toutes
    # les couches est structurellement au contact du plateau.
    platform_contact_tolerance_factor: float = 2.5

    # Cinematique de la machine (RotationAxisSpec), pour tester le support le
    # long des directions REELLEMENT atteignables, pas seulement la
    # verticale monde - meme fix que step9b (limite "Validite du test
    # vertical", voir sa docstring de module). Bug reel trouve en session sur
    # res/test_models/l_shape_overhang.stl (surplomb a 90 degres, cas simple) :
    # un test purement vertical fait grimper `support_unresolved` a 55.5% des
    # couches, alors que le bras incline REPOSE bien sur le pied, juste pas
    # verticalement au-dessus de lui - la meme cascade genait deja death.stl,
    # en pire (94%). None (defaut) = comportement historique inchange
    # (verticale seule) : step5 ne connait pas `MachineProfile`, l'appelant
    # (main.py) doit passer explicitement `MachineProfile.rotation_axes` s'il
    # veut ce fix - meme discipline que `nozzle_width`, jamais une seconde
    # source de verite devinee.
    support_rotation_axes: Optional[list] = None

    # Nombre de directions candidates testees quand `support_rotation_axes`
    # est actif - meme grandeur et meme defaut que
    # `SegmentationParams.direction_candidate_count` (step3) et
    # `SupportCheckParams.direction_candidate_count` (step9b), jamais une
    # nouvelle densite inventee.
    support_direction_candidate_count: int = 26


# ---------------------------------------------------------------------------
# Decision 1 - Marching tetrahedra (3 cas geometriques)
# ---------------------------------------------------------------------------

def _interpolate_edge(
    p_from: np.ndarray, v_from: float, p_to: np.ndarray, v_to: float, isovalue: float
) -> np.ndarray:
    t = (isovalue - v_from) / (v_to - v_from)
    return p_from + t * (p_to - p_from)


def _orient_triangle(
    p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, reference_direction: np.ndarray
) -> np.ndarray:
    """
    Oriente le triangle (p0,p1,p2) pour que sa normale ait un produit
    scalaire positif avec `reference_direction`. Utilise pour toutes les
    variantes de marching tetrahedra (cf. Decision 3bis du plan) : garantit
    une orientation coherente sur tout le patch extrait, sans dependre de
    l'ordre local (arbitraire) des sommets du tetraedre source.
    """
    normal = np.cross(p1 - p0, p2 - p0)
    if np.dot(normal, reference_direction) < 0.0:
        return np.array([p0, p2, p1])
    return np.array([p0, p1, p2])


def _march_tetrahedron(
    tet_points: np.ndarray, tet_values: np.ndarray, isovalue: float
) -> list[np.ndarray]:
    """
    Extrait les triangles de l'isosurface `isovalue` pour UN tetraedre.

    tet_points: (4,3), tet_values: (4,). Retourne une liste de triangles
    (chacun (3,3)), toujours orientes normale pointant du cote "en dessous"
    vers le cote "au dessus" de l'isovaleur.

    La direction de reference pour l'orientation (`p[hi] - p[lo]`, ou hi/lo
    sont les sommets de valeur max/min du tetraedre) est GARANTIE d'avoir un
    produit scalaire positif avec le vrai gradient local : le champ est P1
    (lineaire par tetraedre), donc phi(hi) - phi(lo) = grad(phi) . (p[hi] -
    p[lo]) exactement, et phi(hi) > phi(lo) par construction. C'est ce qui
    rend l'orientation independante de l'ordre local (arbitraire, cf.
    step4_field_solver._compute_shape_gradients qui prend deja abs(volume)
    pour la meme raison) des sommets du tetraedre.

    Convention de degenerescence : un sommet dont la valeur egale exactement
    l'isovaleur est classe "au-dessus" (>=), pour eliminer toute ambiguite
    de classement sans biaiser les valeurs elles-memes.
    """
    above = tet_values >= isovalue
    n_above = int(above.sum())

    if n_above == 0 or n_above == 4:
        return []

    hi_idx = int(np.argmax(tet_values))
    lo_idx = int(np.argmin(tet_values))
    if tet_values[hi_idx] == tet_values[lo_idx]:
        # tetraedre a champ constant (degenere pour cette isovaleur) : pas
        # de direction de reference exploitable, on l'ignore proprement.
        return []
    reference_direction = tet_points[hi_idx] - tet_points[lo_idx]

    if n_above in (1, 3):
        lone_is_above = n_above == 1
        lone_idx = int(np.where(above == lone_is_above)[0][0])
        rest_idx = [i for i in range(4) if i != lone_idx]

        pts = []
        for r in rest_idx:
            if lone_is_above:
                p = _interpolate_edge(
                    tet_points[r], tet_values[r], tet_points[lone_idx], tet_values[lone_idx], isovalue
                )
            else:
                p = _interpolate_edge(
                    tet_points[lone_idx], tet_values[lone_idx], tet_points[r], tet_values[r], isovalue
                )
            pts.append(p)

        return [_orient_triangle(pts[0], pts[1], pts[2], reference_direction)]

    # n_above == 2 : quadrilatere -> 2 triangles. Les 4 aretes qui
    # traversent l'isovaleur relient chacun des 2 sommets "au-dessus" a
    # chacun des 2 sommets "en-dessous" (les aretes above-above et
    # below-below ne traversent pas).
    a1, a2 = [int(i) for i in np.where(above)[0]]
    b1, b2 = [int(i) for i in np.where(~above)[0]]

    p_a1b1 = _interpolate_edge(tet_points[b1], tet_values[b1], tet_points[a1], tet_values[a1], isovalue)
    p_a1b2 = _interpolate_edge(tet_points[b2], tet_values[b2], tet_points[a1], tet_values[a1], isovalue)
    p_a2b1 = _interpolate_edge(tet_points[b1], tet_values[b1], tet_points[a2], tet_values[a2], isovalue)
    p_a2b2 = _interpolate_edge(tet_points[b2], tet_values[b2], tet_points[a2], tet_values[a2], isovalue)

    return [
        _orient_triangle(p_a1b1, p_a1b2, p_a2b2, reference_direction),
        _orient_triangle(p_a1b1, p_a2b2, p_a2b1, reference_direction),
    ]


# ---------------------------------------------------------------------------
# Decision 1bis - Version VECTORISEE de _march_tetrahedron, appliquee a TOUS
# les tetraedres candidats d'une isovaleur a la fois (zero boucle Python sur
# les tetraedres/triangles). Remplace le goulot mesure : _march_tetrahedron
# etait appelee ~321000 fois sur death.stl (89% du temps de step5), dont
# 421545 appels a np.cross sur des vecteurs 3x3 individuels via
# _orient_triangle - meme anti-pattern de dispatch numpy generique deja
# corrige pour l'IK (CLAUDE.md limite #37). Logique mathematique
# STRICTEMENT identique a _march_tetrahedron/_interpolate_edge/
# _orient_triangle (conservees, jamais supprimees - _march_tetrahedron sert
# d'oracle de test), seule l'application passe du scalaire au batch.
# ---------------------------------------------------------------------------

# exclude_table[lone_idx] -> les 3 autres indices locaux, ordre croissant -
# remplace `[i for i in range(4) if i != lone_idx]` sans boucle Python.
_CAP_EXCLUDE_TABLE = np.array(
    [[j for j in range(4) if j != i] for i in range(4)], dtype=np.int64
)

_QUAD_BIT_WEIGHTS = np.array([1, 2, 4, 8], dtype=np.int64)


def _build_quad_tables() -> tuple[np.ndarray, np.ndarray]:
    """
    table_a[code]/table_b[code] -> les 2 indices locaux "au-dessus"/"en-dessous"
    (ordre croissant), pour le motif booleen `above` encode en entier 4 bits
    (code = somme des poids 2**k des bits a True). Les 10 codes dont le
    popcount != 2 (n_above != 2) sont laisses a -1 : ne sont JAMAIS indexes
    en pratique (code n'est calcule que sur des lignes deja filtrees par
    quad_mask, donc n_above==2 garanti) - le -1 n'est qu'un garde-fou
    defensif qui romprait bruyamment un futur usage incorrect.
    """
    table_a = np.full((16, 2), -1, dtype=np.int64)
    table_b = np.full((16, 2), -1, dtype=np.int64)
    for pattern in range(16):
        bits = [(pattern >> k) & 1 for k in range(4)]
        if sum(bits) != 2:
            continue
        table_a[pattern] = [k for k in range(4) if bits[k] == 1]
        table_b[pattern] = [k for k in range(4) if bits[k] == 0]
    return table_a, table_b


_QUAD_TABLE_A, _QUAD_TABLE_B = _build_quad_tables()


def _interpolate_edge_batch(
    p_from: np.ndarray, v_from: np.ndarray, p_to: np.ndarray, v_to: np.ndarray, isovalue: float
) -> np.ndarray:
    """Forme batchee de _interpolate_edge - meme formule, aucun changement
    de logique. p_from/p_to (M,3), v_from/v_to (M,). Retourne (M,3)."""
    t = (isovalue - v_from) / (v_to - v_from)
    return p_from + t[:, None] * (p_to - p_from)


def _orient_triangle_batch(
    p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, reference_direction: np.ndarray
) -> np.ndarray:
    """Forme batchee de _orient_triangle - meme formule, aucun changement
    de logique. p0/p1/p2/reference_direction (M,3). Retourne (M,3,3) :
    ordre des sommets par ligne = (p0,p1,p2) ou (p0,p2,p1)."""
    normal = np.cross(p1 - p0, p2 - p0)
    dot = np.einsum("ij,ij->i", normal, reference_direction)
    flip = dot < 0.0

    tri = np.stack([p0, p1, p2], axis=1)
    tri_flipped = np.stack([p0, p2, p1], axis=1)
    return np.where(flip[:, None, None], tri_flipped, tri)


def _march_tetrahedra_batch(
    tet_points_batch: np.ndarray, tet_values_batch: np.ndarray, isovalue: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Equivalent vectorise d'appeler _march_tetrahedron une fois par ligne de
    (tet_points_batch, tet_values_batch) - voir sa docstring pour la
    justification geometrique complete (champ P1, reference_direction
    garantie coherente avec le gradient reel), inchangee ici.

    Args:
        tet_points_batch: (M,4,3).
        tet_values_batch: (M,4).
        isovalue: valeur de champ de l'isosurface.

    Returns:
        (vertices, local_tet_idx) : vertices (3K,3), triplets consecutifs
        de sommets de triangle (K = nb total de triangles extraits) ;
        local_tet_idx (K,) indice de LIGNE dans le batch d'entree (PAS
        l'indice tetraedre global - l'appelant doit le remapper via
        candidates[local_tet_idx], voir _extract_isosurface_triangles).
        Une ligne qui ne produit aucun triangle (n_above in (0,4), ou champ
        constant degenere) ne contribue AUCUNE ligne en sortie - jamais un
        placeholder/NaN.
    """
    m = tet_points_batch.shape[0]
    if m == 0:
        return np.zeros((0, 3)), np.zeros(0, dtype=np.int64)

    above = tet_values_batch >= isovalue  # (M,4) - meme convention >= que le scalaire
    n_above = above.sum(axis=1)

    hi_idx = np.argmax(tet_values_batch, axis=1)
    lo_idx = np.argmin(tet_values_batch, axis=1)
    rows = np.arange(m)
    hi_val = tet_values_batch[rows, hi_idx]
    lo_val = tet_values_batch[rows, lo_idx]
    hi_pt = tet_points_batch[rows, hi_idx]
    lo_pt = tet_points_batch[rows, lo_idx]
    reference_direction = hi_pt - lo_pt

    # meme paire de branches "return []" que le scalaire, fusionnees en un
    # masque de validite : n_above in (0,4) OU champ constant (hi==lo).
    non_degenerate = hi_val != lo_val
    active = non_degenerate & (n_above != 0) & (n_above != 4)

    cap_mask = active & ((n_above == 1) | (n_above == 3))
    quad_mask = active & (n_above == 2)

    vertices_parts: list[np.ndarray] = []
    tet_idx_parts: list[np.ndarray] = []

    # ---------------- Cas "cap" : n_above in (1,3) -> 1 triangle/tet ------
    if np.any(cap_mask):
        idx = np.nonzero(cap_mask)[0]
        above_c = above[idx]
        n_above_c = n_above[idx]
        lone_is_above = n_above_c == 1

        match = above_c == lone_is_above[:, None]
        lone_idx = np.argmax(match, axis=1)
        rest_idx = _CAP_EXCLUDE_TABLE[lone_idx]  # (Mc,3)

        pts_c = tet_points_batch[idx]
        vals_c = tet_values_batch[idx]
        mc = len(idx)
        rows_c = np.arange(mc)

        lone_pt = pts_c[rows_c, lone_idx]
        lone_val = vals_c[rows_c, lone_idx]

        rest_pt = np.take_along_axis(pts_c, rest_idx[:, :, None], axis=1)  # (Mc,3,3)
        rest_val = np.take_along_axis(vals_c, rest_idx, axis=1)  # (Mc,3)

        # Aplati la dimension "3 aretes par tet" dans la dimension batch
        # pour un seul appel a _interpolate_edge_batch. Ordre C tet-majeur/
        # arete-mineur : rest_pt.reshape(3*Mc,3) donne
        # [tet0_e0,tet0_e1,tet0_e2, tet1_e0,...] et np.repeat(lone_pt,3,axis=0)
        # donne [tet0,tet0,tet0, tet1,...] - meme ordre, alignement correct
        # (verifie empiriquement pendant la conception, 7129 tetraedres
        # synthetiques compares a l'oracle scalaire, 0 ecart).
        rest_pt_flat = rest_pt.reshape(3 * mc, 3)
        rest_val_flat = rest_val.reshape(3 * mc)
        lone_pt_rep = np.repeat(lone_pt, 3, axis=0)
        lone_val_rep = np.repeat(lone_val, 3)
        lone_is_above_rep = np.repeat(lone_is_above, 3)

        # Direction d'interpolation identique au scalaire :
        #   lone_is_above=True  -> from=rest, to=lone
        #   lone_is_above=False -> from=lone, to=rest
        from_pt = np.where(lone_is_above_rep[:, None], rest_pt_flat, lone_pt_rep)
        from_val = np.where(lone_is_above_rep, rest_val_flat, lone_val_rep)
        to_pt = np.where(lone_is_above_rep[:, None], lone_pt_rep, rest_pt_flat)
        to_val = np.where(lone_is_above_rep, lone_val_rep, rest_val_flat)

        edge_pts = _interpolate_edge_batch(from_pt, from_val, to_pt, to_val, isovalue)
        edge_pts = edge_pts.reshape(mc, 3, 3)

        p0, p1, p2 = edge_pts[:, 0], edge_pts[:, 1], edge_pts[:, 2]
        ref_c = reference_direction[idx]
        tri_c = _orient_triangle_batch(p0, p1, p2, ref_c)

        vertices_parts.append(tri_c.reshape(3 * mc, 3))
        tet_idx_parts.append(idx)

    # ---------------- Cas "quad" : n_above == 2 -> 2 triangles/tet --------
    if np.any(quad_mask):
        idx = np.nonzero(quad_mask)[0]
        above_q = above[idx]
        code = (above_q.astype(np.int64) * _QUAD_BIT_WEIGHTS).sum(axis=1)

        ab = _QUAD_TABLE_A[code]  # (Mq,2) : a1,a2 croissant
        bb = _QUAD_TABLE_B[code]  # (Mq,2) : b1,b2 croissant

        pts_q = tet_points_batch[idx]
        vals_q = tet_values_batch[idx]
        mq = len(idx)
        rows_q = np.arange(mq)

        a1_pt, a2_pt = pts_q[rows_q, ab[:, 0]], pts_q[rows_q, ab[:, 1]]
        b1_pt, b2_pt = pts_q[rows_q, bb[:, 0]], pts_q[rows_q, bb[:, 1]]
        a1_val, a2_val = vals_q[rows_q, ab[:, 0]], vals_q[rows_q, ab[:, 1]]
        b1_val, b2_val = vals_q[rows_q, bb[:, 0]], vals_q[rows_q, bb[:, 1]]

        # 4 aretes par tet (chaque sommet "au-dessus" vers chaque sommet
        # "en-dessous"), toujours interpolees from=below, to=above - meme
        # convention que le scalaire.
        from_pt_stack = np.stack([b1_pt, b2_pt, b1_pt, b2_pt], axis=1)
        to_pt_stack = np.stack([a1_pt, a1_pt, a2_pt, a2_pt], axis=1)
        from_val_stack = np.stack([b1_val, b2_val, b1_val, b2_val], axis=1)
        to_val_stack = np.stack([a1_val, a1_val, a2_val, a2_val], axis=1)

        edge_pts = _interpolate_edge_batch(
            from_pt_stack.reshape(4 * mq, 3), from_val_stack.reshape(4 * mq),
            to_pt_stack.reshape(4 * mq, 3), to_val_stack.reshape(4 * mq),
            isovalue,
        ).reshape(mq, 4, 3)

        p_a1b1, p_a1b2, p_a2b1, p_a2b2 = (
            edge_pts[:, 0], edge_pts[:, 1], edge_pts[:, 2], edge_pts[:, 3]
        )
        ref_q = reference_direction[idx]

        tri1 = _orient_triangle_batch(p_a1b1, p_a1b2, p_a2b2, ref_q)
        tri2 = _orient_triangle_batch(p_a1b1, p_a2b2, p_a2b1, ref_q)

        # Interleave tri1[i], tri2[i] consecutifs pour le meme tet - parite
        # avec l'ordre d'ajout du scalaire (tri1 avant tri2 par tet).
        both = np.stack([tri1, tri2], axis=1)  # (Mq,2,3,3)
        vertices_parts.append(both.reshape(2 * mq * 3, 3))
        tet_idx_parts.append(np.repeat(idx, 2))

    if not vertices_parts:
        return np.zeros((0, 3)), np.zeros(0, dtype=np.int64)

    vertices = np.concatenate(vertices_parts, axis=0)
    local_tet_idx = np.concatenate(tet_idx_parts, axis=0)
    return vertices, local_tet_idx


def _tet_volumes(tet_mesh: TetMesh) -> np.ndarray:
    """Volume de chaque tetraedre (meme formule que step4_field_solver, non partagee - trop petit pour justifier un utilitaire commun)."""
    verts = tet_mesh.tet_vertices[tet_mesh.tetrahedra]  # (T,4,3)
    edge_matrices = np.transpose(verts[:, 1:, :] - verts[:, 0:1, :], (0, 2, 1))  # (T,3,3)
    return np.abs(np.linalg.det(edge_matrices)) / 6.0


def _extract_isosurface_triangles(
    tet_mesh: TetMesh,
    values: np.ndarray,
    isovalue: float,
    tet_min: np.ndarray,
    tet_max: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Extrait le patch triangule complet de l'isosurface `isovalue` sur tout
    le TetMesh. Filtre d'abord vectoriellement les tetraedres candidats
    (`tet_min <= isovalue <= tet_max`), pour ne classifier triangle par
    triangle que les tetraedres realistiquement traverses.

    Retourne (vertices (3N,3), faces (N,3), source_tet_ids (N,)) - un
    triplet de sommets et l'indice du tetraedre source par triangle extrait
    (necessaire pour l'association couche->regions de la decision 4).

    Classification/extraction VECTORISEE sur tous les tetraedres candidats
    a la fois (`_march_tetrahedra_batch`) - remplace l'ancienne boucle
    Python par tetraedre (goulot mesure : 89% du temps de step5 sur un
    modele dense, cf. CLAUDE.md). `_march_tetrahedron` reste la reference
    scalaire (oracle de test), non appelee ici.
    """
    candidates = np.where((tet_min <= isovalue) & (isovalue <= tet_max))[0]

    if len(candidates) == 0:
        return np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int64), np.zeros(0, dtype=np.int64)

    tets = tet_mesh.tetrahedra[candidates]
    tet_points_batch = tet_mesh.tet_vertices[tets]
    tet_values_batch = values[tets]

    vertices, local_tet_idx = _march_tetrahedra_batch(tet_points_batch, tet_values_batch, isovalue)

    if len(vertices) == 0:
        return np.zeros((0, 3)), np.zeros((0, 3), dtype=np.int64), np.zeros(0, dtype=np.int64)

    faces = np.arange(len(vertices), dtype=np.int64).reshape(-1, 3)
    source_tets = candidates[local_tet_idx]
    return vertices, faces, source_tets


def _weld_vertices(
    vertices: np.ndarray, faces: np.ndarray, tolerance: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Fusionne les sommets proches (a `tolerance` pres) avant extraction du
    contour. Deux tetraedres adjacents interpolent independamment le meme
    point d'arete partagee ; ce devrait etre bit-identique en theorie (memes
    entrees, meme formule), mais fusionner par tolerance absorbe tout ecart
    residuel sans hypothese fragile sur l'arithmetique flottante.
    """
    if tolerance <= 0 or len(vertices) == 0:
        return vertices, faces

    rounded = np.round(vertices / tolerance).astype(np.int64)
    unique_rounded, unique_indices, inverse = np.unique(
        rounded, axis=0, return_index=True, return_inverse=True
    )
    welded_vertices = vertices[unique_indices]
    remapped_faces = inverse[faces].reshape(faces.shape)
    return welded_vertices, remapped_faces


# ---------------------------------------------------------------------------
# Decision 3 - Contours fermes via trace de bord (trimesh.outline())
# ---------------------------------------------------------------------------

def _extract_closed_contours(vertices: np.ndarray, faces: np.ndarray) -> list[np.ndarray]:
    """
    Extrait les boucles de contour fermees du patch triangule (vertices,
    faces), via le trace de bord de trimesh. `process=False` : on ne laisse
    pas trimesh reordonner/fusionner les sommets a notre place (deja fait
    explicitement par `_weld_vertices`), pour garder le controle complet de
    la correspondance triangle -> tetraedre source en amont.

    Gere nativement les patchs a plusieurs boucles disjointes (ex: paroi
    exterieure + paroi interieure d'une cavite a la meme hauteur de champ) -
    chaque boucle devient un element de la liste retournee.
    """
    # Triangles DEGENERES (deux sommets confondus par `_weld_vertices`) : quand
    # un noeud est a moins de la tolerance de soudure du plan, les points de
    # coupe de toutes ses aretes fusionnent. Sans aire, ils ne portent aucun
    # bord, mais leurs aretes repetees (multiplicite 6 a 22 mesuree) cassent le
    # trace de bord : boucles ouvertes, donc jetees. Mesure sur un Benchy Bird
    # (noeuds a 2,6e-5 et 9,2e-5 mm du plan) : couche entiere perdue (192 mm2)
    # et 23 % d'une autre.
    if len(faces):
        faces = faces[(faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2])
                      & (faces[:, 0] != faces[:, 2])]
    if len(faces) == 0:
        return []

    try:
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
        outline = mesh.outline()
    except Exception:
        # patch degenere (triangles nuls, topologie non-variete locale) :
        # aucun contour exploitable plutot qu'un crash - cf. Decision 5.
        return []

    contours = []
    for entity in outline.entities:
        if not getattr(entity, "closed", True):
            continue
        points = np.asarray(entity.discrete(outline.vertices))
        if len(points) < 4:  # minimum : 3 sommets distincts + point de fermeture
            continue
        contours.append(points)

    return contours


# ---------------------------------------------------------------------------
# Decision 2 - Espacement des isovaleurs (mediane ponderee par volume)
# ---------------------------------------------------------------------------

def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    values_sorted = values[order]
    cum_weights = np.cumsum(weights[order])
    total = cum_weights[-1]
    if total <= 0.0:
        return float(np.median(values))
    idx = int(np.searchsorted(cum_weights, 0.5 * total))
    idx = min(idx, len(values_sorted) - 1)
    return float(values_sorted[idx])


def _representative_gradient_norm(scalar_field: ScalarField, tet_mesh: TetMesh) -> float:
    """
    Mediane ponderee par volume de ||grad(phi)|| sur tout le TetMesh - le
    scalaire "representatif" utilise pour convertir layer_height_base en
    espacement de champ (cf. Decision 2 du plan).

    Mediane plutot que moyenne : l'ancrage plateau de step4
    (`platform_anchor_weight`, poids eleve sur un sous-ensemble restreint de
    noeuds frontiere) tend a produire un gradient localement plus fort sur
    la fine couche de tetraedres adjacents a ces noeuds ancres. Une moyenne
    brute se laisse fausser par cette poignee de valeurs extremes ; la
    mediane ponderee reste representative du comportement typique du champ.
    """
    if scalar_field.gradient is None:
        raise LayerExtractionError(
            "ScalarField.gradient est requis (fourni par pipeline/step4_field_solver.py) "
            "pour convertir layer_height_base en espacement de champ."
        )

    norms = np.linalg.norm(scalar_field.gradient, axis=1)
    volumes = _tet_volumes(tet_mesh)
    return _weighted_median(norms, volumes)


def _generate_isovalues(scalar_field: ScalarField, tet_mesh: TetMesh, params: LayerParams) -> np.ndarray:
    """
    Genere la sequence d'isovaleurs espacees de `delta = layer_height_base *
    representatif(||grad(phi)||)`, entre les bornes REELLES du champ resolu
    (min/max mesures, pas supposees etre exactement [0, max] - l'ancrage
    plateau de step4 est une penalite molle, pas une elimination dure).
    """
    gradient_norm = _representative_gradient_norm(scalar_field, tet_mesh)
    if gradient_norm <= 1e-12:
        raise LayerExtractionError(
            "La norme representative du gradient du champ est quasi nulle - le "
            "champ scalaire semble degenere (constant), impossible d'en deriver "
            "un espacement d'isovaleurs."
        )

    delta = params.layer_height_base * gradient_norm
    v_min = float(scalar_field.values.min())
    v_max = float(scalar_field.values.max())

    if v_max - v_min < delta:
        # modele plus fin qu'une seule couche : une couche unique au centre
        return np.array([0.5 * (v_min + v_max)])

    first = v_min + 0.5 * delta
    count = int(np.floor((v_max - first) / delta)) + 1
    return first + delta * np.arange(count)


# ---------------------------------------------------------------------------
# Decision 4 - Association couche -> regions, rang topologique
# ---------------------------------------------------------------------------

def _map_tets_to_regions(
    tet_mesh: TetMesh, mesh_data: MeshData, region_dependency_graph: RegionDependencyGraph
) -> np.ndarray:
    """
    Pour chaque tetraedre, retrouve la region de sa face de surface la plus
    proche (via mesh_data.kdtree sur le centroide du tetraedre) - meme
    technique que step4_field_solver.build_static_system pour le meme
    pont surface/volume, pas reinventee ici.

    Retourne (T,) region_id, ou -1 si aucune region ne couvre la face
    trouvee (ne devrait pas arriver : step3 garantit que chaque face a
    exactement une region - garde-fou defensif, pas un cas attendu).
    """
    tet_centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra].mean(axis=1)
    _, nearest_face = mesh_data.kdtree.query(tet_centers)
    nearest_face = np.asarray(nearest_face, dtype=np.int64)

    face_to_region = np.full(mesh_data.face_count, -1, dtype=np.int64)
    for region in region_dependency_graph.regions.values():
        face_to_region[region.face_ids] = region.id

    return face_to_region[nearest_face]


def _compute_region_ranks(region_dependency_graph: RegionDependencyGraph) -> dict[int, int]:
    """
    Rang topologique de chaque region : plus long chemin depuis une region
    racine (aucune arete entrante) dans le DAG de dependances. Une region
    racine a rang 0. Accumulation en ordre topologique (meme outil -
    networkx - deja utilise par step3_segmentation pour la detection de
    cycle, nouvel usage ici).
    """
    graph = nx.DiGraph()
    graph.add_nodes_from(region_dependency_graph.regions.keys())
    graph.add_edges_from(region_dependency_graph.dependencies)

    ranks: dict[int, int] = {}
    for region_id in nx.topological_sort(graph):
        predecessors = list(graph.predecessors(region_id))
        ranks[region_id] = 0 if not predecessors else 1 + max(ranks[p] for p in predecessors)
    return ranks


# ---------------------------------------------------------------------------
# Decision 6 - Repere local (origine, axe_u, axe_v)
# ---------------------------------------------------------------------------

def _compute_local_frame(
    contours: list[np.ndarray], tet_gradients: np.ndarray, tet_volumes: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Repere local approximatif de la couche : origine = centroide des points
    de contour, normale = gradient moyen (pondere par volume) des
    tetraedres contribuant a la couche, axe_u/axe_v = base orthonormee du
    plan perpendiculaire (meme construction que
    core.geometry_utils.rotation_axis_circle_directions, reutilisee ici).
    """
    all_points = np.concatenate(contours, axis=0)
    origin = all_points.mean(axis=0)

    weighted_gradient = (tet_gradients * tet_volumes[:, None]).sum(axis=0)
    norm = np.linalg.norm(weighted_gradient)
    normal = np.array([0.0, 0.0, 1.0]) if norm < 1e-12 else weighted_gradient / norm

    reference = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(reference, normal)) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])

    axis_u = reference - np.dot(reference, normal) * normal
    axis_u = axis_u / np.linalg.norm(axis_u)
    axis_v = np.cross(normal, axis_u)

    return origin, axis_u, axis_v


# ---------------------------------------------------------------------------
# Decision 7 - Clustering spatial des boucles disjointes non imbriquees
# (CLAUDE.md : jonction en T/etoile, une isovaleur peut traverser
# simultanement deux parties du solide tres eloignees l'une de l'autre)
# ---------------------------------------------------------------------------

def _min_pairwise_distance(loop_a: np.ndarray, loop_b: np.ndarray) -> float:
    """Distance minimale point-a-point entre deux boucles 3D (enveloppes brutes, pas de projection)."""
    diffs = loop_a[:, None, :] - loop_b[None, :, :]
    return float(np.linalg.norm(diffs, axis=-1).min())


def _loops_connected(
    proj_a: np.ndarray,
    proj_b: np.ndarray,
    loop_a_3d: np.ndarray,
    loop_b_3d: np.ndarray,
    gap_threshold: float,
) -> bool:
    """
    Deux boucles sont "connectees" si elles sont soit IMBRIQUEES (une
    contient l'autre - ex: paroi exterieure + cavite interne, ne doit JAMAIS
    etre scinde, cf. CLAUDE.md), soit a une distance de bord-a-bord
    inferieure a `gap_threshold` (evite une discontinuite artificielle juste
    avant une fusion geometrique reelle a une isovaleur voisine).

    `proj_a`/`proj_b` : boucles projetees en 2D dans un repere de reference
    commun (utilise uniquement pour le test d'imbrication - la topologie
    relative de deux boucles franchement disjointes ne peut pas s'imbriquer
    accidentellement par un mauvais choix de plan, cf. CLAUDE.md).
    `loop_a_3d`/`loop_b_3d` : boucles brutes 3D (test de distance, pas de
    projection necessaire).

    Factorisee depuis `_cluster_contours` (clustering intra-isovaleur) pour
    reutilisation par `_build_continuity_graph` (continuite inter-isovaleur,
    Decision 8) - meme critere, une seule source de verite.
    """
    nested = (
        point_in_polygon_2d(proj_b[0], proj_a)
        or point_in_polygon_2d(proj_a[0], proj_b)
    )
    if nested:
        return True
    return _min_pairwise_distance(loop_a_3d, loop_b_3d) <= gap_threshold


def _cluster_contours(
    contours: list[np.ndarray],
    provisional_frame: tuple[np.ndarray, np.ndarray, np.ndarray],
    gap_threshold: float,
) -> list[list[int]]:
    """
    Partitionne `contours` (boucles fermees d'une meme couche) en clusters
    spatiaux, via `_loops_connected` (imbrication ou distance de bord-a-bord).

    `provisional_frame` sert UNIQUEMENT de reference de projection pour le
    test d'imbrication (calcule sur TOUTES les boucles de la couche, avant
    clustering - le meme plan "approximatif" qui motive ce fix en aval,
    mais les erreurs d'angle de ce plan moyen n'affectent pas la topologie
    relative de deux boucles tant que la projection ne degenere pas).

    Retourne une liste de clusters, chaque cluster etant une liste d'indices
    dans `contours` (triee, ordre stable).
    """
    n = len(contours)
    graph = nx.Graph()
    graph.add_nodes_from(range(n))

    origin, axis_u, axis_v = provisional_frame
    projected = []
    for c in contours:
        relative = c[:-1] - origin  # sans le point de fermeture duplique
        projected.append(np.stack([relative @ axis_u, relative @ axis_v], axis=1))

    for i in range(n):
        for j in range(i + 1, n):
            if _loops_connected(projected[i], projected[j], contours[i], contours[j], gap_threshold):
                graph.add_edge(i, j)

    return [sorted(component) for component in nx.connected_components(graph)]


def _tets_for_cluster(
    unique_tets: np.ndarray,
    tet_centers: np.ndarray,
    cluster_points: np.ndarray,
    tolerance: float,
) -> np.ndarray:
    """
    Sous-ensemble de `unique_tets` dont le centroide est proche (a
    `tolerance` pres) d'au moins un point du cluster - approximation de la
    correspondance triangle->tetraedre-source, PERDUE par
    `trimesh.outline()` (verifie : les entites de `outline.entities` ne
    portent aucune reference aux indices de `faces` d'entree). Un tetraedre
    ambigu (proche de plusieurs clusters) peut legitimement contribuer a
    plusieurs d'entre eux sans consequence : region_ids ne sert qu'a un
    rang topologique MAX pour le tri final (_sort_key), pas a une
    partition exacte.
    """
    if len(unique_tets) == 0:
        return unique_tets
    dists = np.linalg.norm(
        tet_centers[:, None, :] - cluster_points[None, :, :], axis=-1
    ).min(axis=1)
    return unique_tets[dists <= tolerance]


# ---------------------------------------------------------------------------
# Decision 8 - Continuite spatiale inter-isovaleur (ordre DFS des clusters)
#
# Le clustering (Decision 7) resout le cas ou une couche a une isovaleur
# DONNEE a plusieurs boucles disjointes. Mais deux clusters issus d'une
# meme region topologique step3 (une region est une nappe de surface
# connexe, PAS necessairement une seule composante spatiale du solide - ex:
# le pied et la pointe d'un bras d'un modele en T partagent souvent la
# region "verticale" unique, cf. CLAUDE.md) ont le meme rang ET les memes
# region_ids : le tri _sort_key=(rank, isovalue) les fait alterner a chaque
# isovaleur (pied->pointe->pied->pointe...), produisant un trajet a vide a
# chaque bascule. Ce mecanisme construit un graphe de continuite spatiale
# entre clusters d'isovaleurs successives, puis un parcours DFS qui termine
# une branche avant de sauter a une autre - generalise a un nombre
# arbitraire de branches, aux recombinaisons (deux branches qui
# fusionnent), sans jamais violer le rang topologique de region (partition
# stricte avant le DFS, pas une ponderation).
# ---------------------------------------------------------------------------

def _build_continuity_graph(
    raw_layers: list[tuple[float, list[np.ndarray], list[int], tuple]],
    gap_threshold: float,
    progress: Optional[np.ndarray] = None,
) -> nx.DiGraph:
    """
    Graphe oriente : un noeud par entree de `raw_layers` (indice dans cette
    liste), une arete i->j si le cluster i (isovaleur v_i) et le cluster j
    (isovaleur v_j) sont connectes (`_loops_connected`) ET v_i/v_j sont des
    isovaleurs distinctes ADJACENTES dans la sequence reellement presente
    dans `raw_layers` (pas la grille theorique de `_generate_isovalues` -
    une isovaleur sans contour exploitable est simplement absente et ne doit
    pas casser la continuite de ses voisines survivantes).

    Le seuil de connexion est `gap_threshold`, le MEME que celui utilise par
    `_cluster_contours` pour le clustering intra-isovaleur : deux isovaleurs
    successives sont physiquement proches en 3D (c'est la definition meme de
    "couches successives"), un cluster et son prolongement direct a
    l'isovaleur suivante doit etre au moins aussi proche que deux boucles de
    la meme couche qui ne fusionnent pas encore - pas une nouvelle constante
    arbitraire.

    La projection necessaire au test d'imbrication de `_loops_connected`
    reutilise le `local_frame` DEJA calcule par cluster dans `raw_layers`
    (contrairement a `_cluster_contours`, qui doit encore utiliser un plan
    provisoire commun avant que les clusters n'existent) : chaque cluster a
    deja son propre plan local coherent depuis le clustering. Projection
    faite dans le repere du cluster i (choix arbitraire mais coherent pour
    le test de distance ; le test d'imbrication teste deja les deux sens).

    `progress` (optionnel, une valeur par entree de `raw_layers`) remplace
    l'isovaleur pour ORIENTER les aretes, sans jamais toucher au critere
    d'ADJACENCE, qui reste l'isovaleur : "quelles couches se touchent" est une
    question de geometrie (deux isosurfaces successives sont voisines), tandis
    que "laquelle vient avant l'autre" est une question de progression
    d'impression. Confondre les deux est precisement ce qui fait echouer un
    tri global par distance geodesique seule (mesure : x4 de degradation,
    l'ordre entrelacant des couches spatialement eloignees qui n'ont en commun
    que leur distance au plateau).

    Avec `progress`, une arete peut donc etre INVERSEE par rapport a
    l'isovaleur : c'est exactement le cas d'une coque suspendue dont le champ
    croit de l'extremite vers l'attache, ou l'ordre correct est l'inverse de
    l'isovaleur. Une paire de progression strictement egale ne produit aucune
    arete (aucun des deux ne precede l'autre).

    Retourne un graphe acyclique par construction : les aretes vont toujours
    d'une valeur strictement inferieure vers une strictement superieure, que
    cette valeur soit l'isovaleur ou `progress`.
    """
    graph = nx.DiGraph()
    graph.add_nodes_from(range(len(raw_layers)))

    groups: dict[float, list[int]] = {}
    for idx, (isovalue, *_rest) in enumerate(raw_layers):
        groups.setdefault(round(isovalue, 6), []).append(idx)

    sorted_isovalues = sorted(groups.keys())

    for k in range(len(sorted_isovalues) - 1):
        indices_a = groups[sorted_isovalues[k]]
        indices_b = groups[sorted_isovalues[k + 1]]

        for i in indices_a:
            _, contours_i, _region_ids_i, frame_i = raw_layers[i]
            origin_i, axis_u_i, axis_v_i = frame_i
            points_i = np.concatenate(contours_i, axis=0)
            proj_i = np.stack(
                [
                    (points_i - origin_i) @ axis_u_i,
                    (points_i - origin_i) @ axis_v_i,
                ],
                axis=1,
            )

            for j in indices_b:
                _, contours_j, _region_ids_j, _frame_j = raw_layers[j]
                points_j = np.concatenate(contours_j, axis=0)
                proj_j = np.stack(
                    [
                        (points_j - origin_i) @ axis_u_i,
                        (points_j - origin_i) @ axis_v_i,
                    ],
                    axis=1,
                )
                if not _loops_connected(
                    proj_i, proj_j, points_i, points_j, gap_threshold
                ):
                    continue

                weight = _min_pairwise_distance(points_i, points_j)
                if progress is None:
                    graph.add_edge(i, j, weight=weight)
                elif progress[i] < progress[j]:
                    graph.add_edge(i, j, weight=weight)
                elif progress[j] < progress[i]:
                    # La progression d'impression s'oppose a l'isovaleur ici :
                    # l'arete est inversee, pas ignoree.
                    graph.add_edge(j, i, weight=weight)

    return graph


def _subtree_sizes(graph: nx.DiGraph, nodes: set[int]) -> dict[int, int]:
    """
    Taille du sous-arbre atteignable depuis chaque noeud (nombre de
    descendants, soi-meme inclus), restreint a `nodes`. Utilise
    `nx.descendants` (pas une somme naive des tailles d'enfants) car le
    graphe de continuite peut avoir des recombinaisons (in-degree >= 2) : un
    sous-arbre partage entre deux parents ne doit pas etre compte en double.
    Cout negligeable a l'echelle reelle de ces graphes (quelques dizaines de
    noeuds par couche).
    """
    subgraph = graph.subgraph(nodes)
    return {node: len(nx.descendants(subgraph, node)) + 1 for node in nodes}


def _dfs_order_within_rank(
    graph: nx.DiGraph,
    nodes: set[int],
    raw_layers: list[tuple[float, list[np.ndarray], list[int], tuple]],
    progress: Optional[np.ndarray] = None,
) -> list[int]:
    """
    Ordre de visite DFS des `nodes` (sous-ensemble de raw_layers appartenant
    a un meme rang topologique), respectant la continuite spatiale : une
    branche est entierement terminee avant de sauter a une autre.

    A un noeud de branchement, visite la branche de plus PETITE taille de
    sous-arbre en premier (`_subtree_sizes`) : termine les branches courtes
    d'abord, jamais de retour en arriere sur une grande branche restee en
    suspens - robuste a un nombre arbitraire d'enfants, contrairement a un
    critere de proximite au point de sortie du parent qui degenere en
    quasi-egalite au moment meme du split (cf. CLAUDE.md).

    Entre composantes faiblement connexes independantes d'un meme rang,
    visite par isovaleur minimale croissante (generalisation directe du tri
    historique par isovaleur, applique a l'unite "composante" plutot qu'a un
    cluster individuel).

    `progress` (optionnel) remplace l'isovaleur comme mesure de PROGRESSION
    d'impression - ordre des composantes et choix des racines. La continuite
    spatiale (structure du graphe, DFS, priorite aux petites branches) est
    strictement inchangee : ce parametre ne dit pas "ou aller", seulement "par
    quel bout commencer".

    Tri explicite partout (tie-break par indice de noeud) : deux appels
    identiques produisent toujours exactement le meme ordre.
    """
    def key_of(node: int) -> float:
        return float(progress[node]) if progress is not None else raw_layers[node][0]

    subgraph = graph.subgraph(nodes)
    components = sorted(
        nx.weakly_connected_components(subgraph),
        key=lambda comp: (min(key_of(n) for n in comp), min(comp)),
    )

    sizes = _subtree_sizes(graph, nodes)
    visited: set[int] = set()
    order: list[int] = []

    def visit(node: int) -> None:
        if node in visited:
            return
        visited.add(node)
        order.append(node)
        children = sorted(
            (c for c in subgraph.successors(node) if c not in visited),
            key=lambda c: (sizes[c], c),
        )
        for child in children:
            visit(child)

    for component in components:
        roots = sorted(
            (n for n in component if subgraph.in_degree(n) == 0),
            key=lambda n: (key_of(n), n),
        )
        if not roots:
            # garde-fou defensif (un DAG fini a toujours une source, ne
            # devrait pas arriver) : repli sur le noeud de plus petite
            # isovaleur (ou progression) de la composante.
            roots = [min(component, key=lambda n: (key_of(n), n))]
        for root in roots:
            visit(root)
        # filet de securite : un noeud non atteint depuis les racines
        # detectees est quand meme ajoute, en fin de composante, par
        # isovaleur croissante - ne perd jamais un Layer silencieusement.
        for node in sorted(component, key=lambda n: (raw_layers[n][0], n)):
            if node not in visited:
                visit(node)

    return order


def _spatial_dfs_order(
    raw_layers: list[tuple[float, list[np.ndarray], list[int], tuple]],
    region_ranks: dict[int, int],
    gap_threshold: float,
    progress: Optional[np.ndarray] = None,
) -> list[int]:
    """
    Ordre final des indices de `raw_layers` : rang topologique de region
    croissant (contrainte physique dure, JAMAIS traversee - partition
    stricte des noeuds par rang avant tout parcours, pas une ponderation),
    puis continuite spatiale DFS a l'interieur de chaque rang
    (`_dfs_order_within_rank`). Retourne une permutation complete des
    indices de `raw_layers`.
    """
    graph = _build_continuity_graph(raw_layers, gap_threshold, progress)

    by_rank: dict[int, set[int]] = {}
    for idx, (_isovalue, _contours, region_ids, _frame) in enumerate(raw_layers):
        rank = max((region_ranks.get(r, 0) for r in region_ids), default=0)
        by_rank.setdefault(rank, set()).add(idx)

    order: list[int] = []
    for rank in sorted(by_rank):
        order.extend(
            _dfs_order_within_rank(graph, by_rank[rank], raw_layers, progress)
        )
    return order


# ---------------------------------------------------------------------------
# Ordre pilote par le SUPPORT REEL (LayerParams.support_driven_order)
# ---------------------------------------------------------------------------

def _layer_support_points(contours: list[np.ndarray]) -> np.ndarray:
    """Tous les points de contour d'une couche, en un seul (P,3)."""
    usable = [np.asarray(c, dtype=float) for c in contours if len(c) > 0]
    if not usable:
        return np.empty((0, 3), dtype=float)
    return np.vstack(usable)


def _place_layer_in_index(
    index, contours: list[np.ndarray], anchored_mask: Optional[np.ndarray] = None
) -> int:
    """
    Ajoute la matiere d'une couche a l'index de capsules, un segment par
    arete de contour - c'est l'approximation "contours seuls" assumee du
    mecanisme (l'infill n'existe qu'a partir de step6, voir
    `LayerParams.support_driven_order`).

    `anchored_mask` (booleen, aligne sur la concatenation des contours -
    meme ordre et meme filtrage des contours vides que
    `_layer_support_points`) restreint l'ajout aux aretes dont les DEUX
    extremites sont ancrees. Une couche n'est ainsi jamais ancree "en bloc" :
    seule la part de sa matiere qui repose reellement sur de l'ancre le
    devient a son tour.

    C'est ce qui permet a une couche partiellement posee sur le plateau de
    servir d'appui pour la suivante, sans pretendre que sa part en
    porte-a-faux tient. Sans ce filtrage par ARETE, il fallait choisir entre
    deux comportements egalement faux : ancrer toute la couche (le faux
    support du defaut n°1 de la limite #51, qui laisse une coque suspendue
    s'auto-ancrer de proche en proche) ou n'ancrer que les couches
    integralement supportees (cascade bloquee des la premiere couche
    partiellement posee - mesure : 100% des couches signalees non resolues
    sur death.stl alors que la geometrie produite etait, elle, meilleure).

    Retourne le nombre d'aretes reellement ajoutees : 0 signifie que cette
    couche n'a rien ancre, et donc qu'aucune couche voisine n'a besoin d'etre
    reevaluee.
    """
    placed = 0
    offset = 0
    for contour in contours:
        points = np.asarray(contour, dtype=float)
        if len(points) == 0:
            continue  # meme filtrage que _layer_support_points : l'offset ne bouge pas
        for i in range(1, len(points)):
            if anchored_mask is not None and not (
                anchored_mask[offset + i - 1] and anchored_mask[offset + i]
            ):
                continue
            index.append_segment(points[i - 1], points[i])
            placed += 1
        offset += len(points)
    return placed


def _supported_mask(
    index,
    points: np.ndarray,
    up: np.ndarray,
    horizontal_radius: float,
    vertical_height: float,
    candidate_directions: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Masque booleen (P,) : pour chaque point, existe-t-il de la matiere deja
    placee qui le supporte le long de `-up`, OU le long d'au moins une des
    `candidate_directions` (D,3) si fournies ?

    `candidate_directions` reproduit le fix de step9b (docstring de module,
    "Validite du test vertical") : un point sans support VERTICAL peut quand
    meme reposer reellement sur de la matiere existante si la machine peut
    deposer selon une direction inclinee atteignable - bug reel trouve en
    session sur un simple surplomb a 90 degres (l_shape_overhang.stl) : le
    champ scalaire courbe fait "demarrer" un bras incline a une hauteur
    superieure a celle ou son porteur vertical s'arrete (~2.5mm d'ecart
    mesure), un test purement vertical ne voit alors JAMAIS le support qui
    existe bel et bien en diagonale, juste en dessous dans la direction
    inclinee reelle du bras - et comme aucune couche du bras ne s'ancre
    jamais, TOUT le reste bascule en cascade (55.5% des couches signalees a
    tort, mesure). None (defaut) = comportement historique, verticale seule.

    Une requete spatiale par point (AABB isotrope couvrant toutes les
    directions testees), puis EXACTEMENT la meme geometrie de test que
    step9b (`point_has_support_along_direction`, vectorisee sur les
    candidats pour chaque direction) - jamais une seconde implementation du
    critere de support. Court-circuite des la premiere direction trouvee.
    """
    supported = np.zeros(len(points), dtype=bool)
    directions = (
        candidate_directions if candidate_directions is not None else up.reshape(1, 3)
    )
    for i, position in enumerate(points):
        aabb = _search_aabb_below(position, up, horizontal_radius, vertical_height)
        if index.has_support(
            aabb, position, directions, horizontal_radius, vertical_height
        ):
            supported[i] = True
    return supported


def _search_aabb_below(
    position: np.ndarray, up: np.ndarray, horizontal_radius: float, vertical_height: float
) -> AABB:
    """
    Boite de recherche englobant la zone "sous" `position` le long de `-up`.

    Etendue isotrope de `max(horizontal_radius, vertical_height)` plutot que
    la boite verticale serree de step9b : `up` est ici une direction
    quelconque, et une boite alignee sur les axes ne peut pas la suivre. Un
    sur-ensemble strict des vrais candidats - jamais un faux negatif, seulement
    moins de filtrage grossier ; la geometrie precise reste faite par
    `point_has_support_along_direction`.
    """
    extent = max(horizontal_radius, vertical_height)
    return AABB(min_corner=position - extent, max_corner=position + extent)


def _support_driven_order(
    raw_layers: list[tuple],
    params: LayerParams,
    up_vector: np.ndarray,
) -> tuple[list[int], set[int]]:
    """
    Ordonne les couches par SUPPORT REEL ET CONNECTE AU PLATEAU, plutot que
    par isovaleur croissante.

    Glouton : a chaque tour, une couche est ELIGIBLE si au moins
    `support_driven_order` (fraction) de ses points de contour ont de la
    matiere ANCREE juste en dessous d'eux. Parmi les eligibles, on prend la
    plus proche de la couche precedente (anti-trajets a vide), puis la plus
    basse. Aucune eligible -> on place la plus basse restante et on la marque
    `support_unresolved`, jamais une exception ni un blocage (meme discipline
    que `CollisionReport.chuck_checked`).

    ANCRAGE, PAS SEULEMENT SUPPORT LOCAL - piege reel rencontre et corrige en
    cours de session : une premiere version testait juste "y a-t-il de la
    matiere deja placee en dessous", sans se soucier de savoir si CETTE
    matiere est elle-meme reliee au plateau. Sur une coque en surplomb
    imprimee du bas vers le haut, chaque couche se retrouve alors bien
    au-dessus de la precedente et passe le test a CHAQUE etape - meme si
    l'ensemble de la coque ne touche le reste du modele qu'en un point tout en
    haut, ce qui degenere silencieusement vers le tri par hauteur (mesure
    directement sur res/test_models/death.stl : ordre quasi strictement
    croissant en Z malgre le mecanisme actif). La geometrie DEJA PLACEE ne
    suffit donc pas : seule la geometrie ANCREE (au contact du plateau, ou
    supportee EXCLUSIVEMENT par de la matiere deja ancree) alimente l'index de
    test. Une couche choisie faute de mieux (`unresolved`) est bien imprimee,
    mais N'ENTRE PAS dans l'index ancre : elle ne peut pas servir de faux
    support a une couche suivante, exactement comme le veut la
    discipline "surplomb reel jamais masque".

    Amorcage de l'ancrage : une couche est ancree D'EMBLEE si un de ses points
    est a moins de `platform_contact_tolerance_factor * nozzle_width` du Z
    minimal de TOUTES les couches - la meme tolerance de contact plateau que
    `SupportCheckParams` (step9b), sans avoir besoin de cabler `platform_z`
    ici (source de verite = la geometrie des couches elles-memes).

    Le support est la cle PRIMAIRE ; proximite et altitude ne departagent que
    des couches deja jugees imprimables. C'est ce qui protege les surplombs
    PAR CONSTRUCTION : une piece pendante est plus basse que son porteur, mais
    tant que le porteur n'existe pas ANCRE elle n'est pas supportee, donc
    jamais choisie - la regle "ce qui pend s'imprime apres" n'est pas encodee
    a part, elle DECOULE du test.

    Performance : une reevaluation naive de toutes les couches a chaque tour
    coute O(couches^2 x points). On exploite la MONOTONIE (la matiere ancree
    ne disparait jamais, un point supporte le reste) : chaque couche ne
    conserve que ses points ENCORE non supportes, et seules les couches dont
    l'AABB touche celle qui vient d'etre ANCREE sont reevaluees.

    Returns:
        (permutation des indices de `raw_layers`, indices places sans support).
    """
    from pipeline.step3_segmentation import _samples_per_axis
    from pipeline.step9_collision_sim import _PrintedGeometryIndex

    if params.nozzle_width is None:
        raise LayerExtractionError(
            "LayerParams.support_driven_order est actif mais nozzle_width est None. "
            "Cette largeur de cordon n'a pas de defaut arbitraire : passez "
            "explicitement MachineProfile.nozzle_diameter, pour ne pas creer une "
            "seconde source de verite divergeant en silence."
        )

    nozzle_width = max(float(params.nozzle_width), 1e-6)
    horizontal_radius = params.support_horizontal_radius_factor * nozzle_width
    vertical_height = params.support_vertical_height_factor * max(
        float(params.layer_height_base), 1e-6
    )
    platform_tolerance = params.platform_contact_tolerance_factor * nozzle_width
    min_fraction = float(params.support_driven_order)

    up = np.asarray(up_vector, dtype=float)
    norm = np.linalg.norm(up)
    up = up / norm if norm > 1e-12 else np.array([0.0, 0.0, 1.0])

    # Directions candidates pour le test de support (voir docstring de
    # `_supported_mask`) - memes fonction et defaut que step3/step9b, jamais
    # une nouvelle densite inventee. None/vide -> verticale seule, repli
    # exact sur le comportement historique.
    candidate_directions: Optional[np.ndarray] = None
    if params.support_rotation_axes:
        samples_per_axis = _samples_per_axis(
            params.support_direction_candidate_count, len(params.support_rotation_axes)
        )
        candidate_directions = np.vstack(
            [up.reshape(1, 3), composed_rotation_directions(params.support_rotation_axes, up, samples_per_axis)]
        )

    all_points = [_layer_support_points(entry[1]) for entry in raw_layers]
    # Masque des points DEJA ANCRES, par couche (monotonie : la matiere
    # ancree ne disparait jamais, un point ancre le reste). Un MASQUE, et non
    # un simple compteur avec une liste de points restants : il faut savoir
    # QUELS points sont ancres pour n'ajouter a l'index que les aretes
    # correspondantes (voir `_place_layer_in_index`).
    anchored_masks = [np.zeros(len(points), dtype=bool) for points in all_points]
    total_counts = [len(points) for points in all_points]

    heights = [
        float(points[:, 2].mean()) if len(points) else 0.0 for points in all_points
    ]
    centroids = [
        points.mean(axis=0) if len(points) else np.zeros(3) for points in all_points
    ]
    aabbs = [
        (points.min(axis=0), points.max(axis=0)) if len(points) else (np.zeros(3), np.zeros(3))
        for points in all_points
    ]

    # Amorcage : les points au contact direct du plateau comptent d'emblee
    # comme "supportes par ancrage", sans attendre qu'une capsule ne les
    # couvre - exactement l'exemption "pose sur le plateau" de step9b.
    global_min_z = min((float(points[:, 2].min()) for points in all_points if len(points)), default=0.0)
    for idx, points in enumerate(all_points):
        if len(points) == 0:
            continue
        anchored_masks[idx] |= points[:, 2] <= global_min_z + platform_tolerance

    # Index de la matiere ANCREE UNIQUEMENT - jamais toute la matiere deja
    # imprimee (voir docstring : c'est precisement la distinction qui corrige
    # le defaut de support purement local). Sert au calcul du statut
    # `unresolved`, jamais au departage.
    anchored_index = _PrintedGeometryIndex(nozzle_width)
    remaining = set(range(len(raw_layers)))
    # Toutes les couches doivent etre evaluees au premier tour ; ensuite,
    # seules celles voisines de la derniere couche ANCREE.
    stale = set(remaining)
    order: list[int] = []
    unresolved: set[int] = set()
    previous_centroid: Optional[np.ndarray] = None

    def _fraction(idx: int) -> float:
        if total_counts[idx] == 0:
            return 1.0  # couche sans point exploitable : jamais bloquante
        return float(anchored_masks[idx].sum()) / total_counts[idx]

    while remaining:
        for idx in stale & remaining:
            pending = ~anchored_masks[idx]
            if not pending.any():
                continue
            # Seuls les points ENCORE non ancres sont retestes (monotonie) ;
            # le resultat est reinjecte a leur position d'origine, pour que le
            # masque reste aligne sur la concatenation des contours.
            pending_indices = np.where(pending)[0]
            mask = _supported_mask(
                anchored_index, all_points[idx][pending_indices], up,
                horizontal_radius, vertical_height, candidate_directions,
            )
            anchored_masks[idx][pending_indices[mask]] = True
        stale = set()

        eligible = [idx for idx in remaining if _fraction(idx) >= min_fraction]

        if eligible:
            if previous_centroid is None:
                chosen = min(eligible, key=lambda i: (heights[i], i))
            else:
                chosen = min(
                    eligible,
                    key=lambda i: (
                        float(np.linalg.norm(centroids[i] - previous_centroid)),
                        heights[i],
                        i,
                    ),
                )
        else:
            # Aucune couche n'est ANCREE (statut inchange, reste honnete) -
            # mais parmi les couches suspendues, on privilegie quand meme
            # celle la plus proche de la MATIERE DEJA DEPOSEE (ancree ou
            # non), plutot que la seule altitude : reduit les sauts a vide
            # entre couches suspendues qui se recouvrent geometriquement
            # (ex: un dome qui se referme progressivement sur lui-meme sans
            # jamais retoucher le plateau) sans jamais faire croire qu'une
            # couche suspendue est supportee - `unresolved` reste base
            # UNIQUEMENT sur `anchored_index`, inchange.
            if previous_centroid is None:
                chosen = min(remaining, key=lambda i: (heights[i], i))
            else:
                chosen = min(
                    remaining,
                    key=lambda i: (
                        float(np.linalg.norm(centroids[i] - previous_centroid)),
                        heights[i],
                        i,
                    ),
                )
            unresolved.add(chosen)

        order.append(chosen)
        remaining.discard(chosen)
        previous_centroid = centroids[chosen]

        # La couche est imprimee quoi qu'il arrive, mais seule la part de sa
        # matiere REELLEMENT ANCREE entre dans l'index - arete par arete, les
        # deux extremites devant etre ancrees. Une couche suspendue n'ancre
        # donc rien (son masque est vide) et ne peut jamais servir de faux
        # support ; une couche partiellement posee ancre sa seule part posee,
        # ce qui laisse la cascade progresser. Voir `_place_layer_in_index`.
        placed_edges = _place_layer_in_index(
            anchored_index, raw_layers[chosen][1], anchored_masks[chosen]
        )

        if placed_edges > 0:
            # Seules les couches dont l'AABB touche celle qui vient d'ancrer
            # de la matiere peuvent avoir gagne du support - les autres
            # gardent leur compte.
            placed_min, placed_max = aabbs[chosen]
            reach = max(horizontal_radius, vertical_height)
            for idx in remaining:
                other_min, other_max = aabbs[idx]
                if np.all(other_min - reach <= placed_max) and np.all(
                    other_max + reach >= placed_min
                ):
                    stale.add(idx)

    return order, unresolved


# ---------------------------------------------------------------------------
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

# Multiple de `layer_height_base` servant de seuil de connexion spatiale quand
# `geodesic_order` est actif SANS `split_disjoint_clusters`. Reprend la valeur
# de ce dernier mesuree comme fonctionnelle sur le modele en T (CLAUDE.md,
# limite #39) plutot que d'introduire une constante independante : c'est le
# meme seuil, pour le meme graphe de continuite.
GEODESIC_ORDER_GAP_FACTOR = 3.0


def _frontier_order(
    raw_layers: list[tuple],
    progress: np.ndarray,
    gap_threshold: float,
) -> tuple[list[int], set[int]]:
    """
    Ordre d'impression par CROISSANCE DE FRONT : une couche n'est imprimable
    que si elle touche de la matiere deja imprimee.

    Retourne `(ordre, indices deposes a contre-sens)`. Le second ensemble
    contient les couches atteintes depuis un voisin d'isovaleur SUPERIEURE :
    on y descend le long du champ, donc le depot s'y fait selon `-grad(phi)`
    et la tete doit pointer vers le haut. C'est le parcours - et lui seul -
    qui connait cette information : il sait par quel voisin chaque couche a
    ete atteinte, ce qu'aucun tri a posteriori ne peut retrouver.

    Dijkstra sur le graphe d'ADJACENCE des couches (non oriente - deux couches
    voisines peuvent etre imprimees dans un sens ou dans l'autre, c'est
    justement ce qui permet de descendre le long d'un pilier suspendu une fois
    son pont construit). A chaque etape, parmi les couches ATTEIGNABLES depuis
    le deja-imprime, on prend celle de plus petite `progress`.

    Difference essentielle avec un tri par `progress` : le tri autorise deux
    couches consecutives situees aux antipodes du modele du moment qu'elles
    sont a distance geodesique voisine (mesure : x3,7 de degradation). Ici la
    connexite du front est une contrainte DURE, pas une preference.

    L'amorce est la couche de plus petite `progress` (la plus proche du
    plateau a travers la matiere). Si le front s'epuise alors que des couches
    restent - piece en plusieurs composantes deconnectees, ou adjacence non
    detectee sous `gap_threshold` - une nouvelle amorce est prise parmi les
    restantes, toujours a plus petite `progress`. Le saut est alors REEL et
    inevitable (rien ne relie ces couches), jamais un choix de l'algorithme.
    """
    import heapq

    graph = _build_continuity_graph(raw_layers, gap_threshold).to_undirected()

    remaining = set(range(len(raw_layers)))
    order: list[int] = []
    reversed_layers: set[int] = set()
    heap: list[tuple[float, int, int]] = []   # (progression, noeud, parent)

    def seed() -> None:
        start = min(remaining, key=lambda n: (progress[n], n))
        heapq.heappush(heap, (float(progress[start]), start, -1))

    seed()
    while remaining:
        if not heap:
            seed()
        _key, node, parent = heapq.heappop(heap)
        if node not in remaining:
            continue
        remaining.discard(node)
        order.append(node)
        # Atteinte depuis un voisin d'isovaleur SUPERIEURE : on descend le
        # champ, donc on depose sous la matiere deja posee. Une amorce
        # (parent = -1) n'a aucun voisin anterieur, donc aucun sens impose.
        if parent >= 0 and raw_layers[parent][0] > raw_layers[node][0]:
            reversed_layers.add(node)
        for neighbour in graph.neighbors(node):
            if neighbour in remaining:
                heapq.heappush(heap, (float(progress[neighbour]), neighbour, node))

    return order, reversed_layers


def _geodesic_layer_progress(
    raw_layers: list[tuple],
    tet_mesh: TetMesh,
    mesh_data: MeshData,
    region_dependency_graph: RegionDependencyGraph,
) -> np.ndarray:
    """
    PROGRESSION d'impression de chaque couche, derivee de la distance
    GEODESIQUE au plateau mesuree a travers la matiere - "par ou cette couche
    rejoint-elle le plateau", et non "a quelle hauteur est-elle".

    C'est une mesure de progression, PAS un ordre a elle seule : trier
    globalement sur cette valeur a ete mesure comme degradant d'un facteur 4
    (l'isosurface d'une distance geodesique est une coquille entiere, donc un
    tri global entrelace des couches spatialement eloignees qui n'ont en
    commun que leur distance au plateau). Elle est destinee a etre passee en
    `progress` au parcours de continuite spatiale, qui l'utilise pour ORIENTER
    et AMORCER, jamais pour sauter.

    Chaque couche est representee par la MEDIANE de la distance geodesique de
    ses points de contour, pas par leur minimum ni leur moyenne : le minimum
    ferait passer une couche entiere devant les autres pour un seul point
    proche de l'attache, et la moyenne serait tiree par les valeurs extremes
    d'une couche qui s'etale. La mediane represente ou se trouve l'essentiel
    de la matiere de la couche.

    Les primitives sont celles de step4 (`_geodesic_distance_to_platform`,
    `_identify_platform_nodes`), jamais une seconde implementation du meme
    calcul : une divergence silencieuse entre deux definitions du "plateau"
    est exactement le piege deja rencontre avec `platform_z`.

    Repli : si la geodesique est inexploitable (aucun noeud plateau, ou une
    couche dont aucun point n'atteint un noeud fini), les couches concernees
    recoivent `+inf` et sont donc repoussees a la fin, dans leur ordre
    d'origine - jamais une exception, et jamais un ordre invente en silence.
    """
    from pipeline.step4_field_solver import (
        _geodesic_distance_to_platform,
        _identify_platform_nodes,
    )

    platform_nodes = _identify_platform_nodes(
        tet_mesh, mesh_data, region_dependency_graph
    )
    geodesic = _geodesic_distance_to_platform(tet_mesh, platform_nodes)

    node_tree = KDTree(tet_mesh.tet_vertices)

    keys: list[float] = []
    for _isovalue, contours, _region_ids, _frame in raw_layers:
        points = [np.asarray(c, dtype=float) for c in contours if len(c) > 0]
        if not points:
            keys.append(float("inf"))
            continue
        _, nearest = node_tree.query(np.vstack(points))
        values = geodesic[np.asarray(nearest, dtype=np.int64)]
        finite = values[np.isfinite(values)]
        keys.append(float(np.median(finite)) if len(finite) else float("inf"))

    return np.asarray(keys, dtype=float)


def run(
    tet_mesh: TetMesh,
    mesh_data: MeshData,
    scalar_field: ScalarField,
    region_dependency_graph: RegionDependencyGraph,
    params: LayerParams,
) -> list[Layer]:
    """
    Extrait les couches courbes du ScalarField et les ordonne.

    Args:
        tet_mesh: sortie de pipeline/step1_ingest.py.
        mesh_data: sortie de pipeline/step1_ingest.py (kdtree pour le pont
            surface/volume).
        scalar_field: sortie de pipeline/step4_field_solver.py (gradient
            requis, pas seulement les valeurs).
        region_dependency_graph: sortie de pipeline/step3_segmentation.py.
        params: hauteur de couche cible et tolerances.

    Returns:
        Liste ordonnee de Layer (id reassigne selon l'ordre d'impression
        final : rang topologique des regions traversees, puis isovaleur).

    Raises:
        LayerExtractionError: si le champ n'a pas de gradient ou est
            degenere (norme de gradient quasi nulle partout).
    """
    isovalues = _generate_isovalues(scalar_field, tet_mesh, params)

    tet_values_per_tet = scalar_field.values[tet_mesh.tetrahedra]  # (T,4)
    tet_min = tet_values_per_tet.min(axis=1)
    tet_max = tet_values_per_tet.max(axis=1)

    tet_regions = _map_tets_to_regions(tet_mesh, mesh_data, region_dependency_graph)
    region_ranks = _compute_region_ranks(region_dependency_graph)
    volumes = _tet_volumes(tet_mesh)

    scale = float(np.ptp(tet_mesh.tet_vertices, axis=0).max()) or 1.0
    weld_tolerance = params.weld_tolerance_factor * scale

    # Seuil de connexion spatiale entre boucles, invariant de boucle (il ne
    # depend que des params) - hisse ici parce que `geodesic_order` en a besoin
    # meme sans `split_disjoint_clusters`, les deux mecanismes partageant le
    # graphe de continuite. Repli sur GEODESIC_ORDER_GAP_FACTOR : deux couches
    # successives sont distantes d'environ `layer_height_base`, donc quelques
    # hauteurs de couche separent ce qui est physiquement le prolongement
    # direct d'un cluster de ce qui ne l'est pas.
    gap_factor = (
        params.split_disjoint_clusters
        if params.split_disjoint_clusters is not None
        else GEODESIC_ORDER_GAP_FACTOR
    )
    gap_threshold = gap_factor * params.layer_height_base

    # taille caracteristique d'un tetraedre (mediane des longueurs d'arete
    # sur tout le TetMesh) - utilisee comme tolerance de proximite pour
    # _tets_for_cluster (Decision 7). PAS une fraction de `scale` (l'etendue
    # du modele entier) : mesure directement sur un cas reel (modele
    # tripode, CLAUDE.md) qu'une fraction de `scale` (~0.32 sur ce modele)
    # tombe EN DESSOUS de la distance reelle centroide-tetraedre -> point de
    # cluster le plus proche (mesuree a 0.328 sur ce meme modele), filtrant
    # ALORS TOUS les tetraedres et forcant un repli sur unique_tets entier
    # (region_ids identiques pour tous les clusters, bug reel mesure avant
    # ce fix). La taille de tetraedre typique (~2.6 unites sur ce modele)
    # est la bonne echelle : un tetraedre contribuant a un triangle de bord
    # d'un cluster a necessairement un centroide proche de ce cluster a
    # l'echelle d'UN tetraedre, pas de l'echelle globale du modele.
    if params.split_disjoint_clusters is not None:
        tet_verts_per_tet = tet_mesh.tet_vertices[tet_mesh.tetrahedra]
        edge_lengths = np.linalg.norm(
            tet_verts_per_tet[:, 1:, :] - tet_verts_per_tet[:, 0:1, :], axis=-1
        )
        cluster_tet_tolerance = float(np.median(edge_lengths)) * 2.0 if edge_lengths.size else scale
    else:
        cluster_tet_tolerance = 0.0

    raw_layers: list[tuple[float, list[np.ndarray], list[int], tuple[np.ndarray, np.ndarray, np.ndarray]]] = []

    for isovalue in isovalues:
        vertices, faces, source_tets = _extract_isosurface_triangles(
            tet_mesh, scalar_field.values, float(isovalue), tet_min, tet_max
        )
        if len(faces) == 0:
            continue

        welded_vertices, welded_faces = _weld_vertices(vertices, faces, weld_tolerance)
        contours = _extract_closed_contours(welded_vertices, welded_faces)
        if not contours:
            continue

        unique_tets = np.unique(source_tets)

        if params.split_disjoint_clusters is None or len(contours) < 2:
            # chemin historique inchange (mot pour mot) - garantit la
            # non-regression par construction, pas par coincidence.
            region_ids = sorted({int(r) for r in tet_regions[unique_tets] if r >= 0})
            local_frame = _compute_local_frame(
                contours, scalar_field.gradient[unique_tets], volumes[unique_tets]
            )
            raw_layers.append((float(isovalue), contours, region_ids, local_frame))
            continue

        provisional_frame = _compute_local_frame(
            contours, scalar_field.gradient[unique_tets], volumes[unique_tets]
        )
        clusters = _cluster_contours(contours, provisional_frame, gap_threshold)

        tet_centers = tet_mesh.tet_vertices[tet_mesh.tetrahedra[unique_tets]].mean(axis=1)

        for cluster_indices in clusters:
            cluster_contours = [contours[i] for i in cluster_indices]
            cluster_points = np.concatenate(cluster_contours, axis=0)
            cluster_tets = _tets_for_cluster(
                unique_tets, tet_centers, cluster_points, tolerance=cluster_tet_tolerance
            )
            if len(cluster_tets) == 0:
                # aucun tetraedre proche detecte (tolerance trop stricte sur
                # un maillage tres grossier) : repli sur unique_tets entier
                # plutot qu'un local_frame/region_ids vides.
                cluster_tets = unique_tets

            region_ids = sorted({int(r) for r in tet_regions[cluster_tets] if r >= 0})
            local_frame = _compute_local_frame(
                cluster_contours, scalar_field.gradient[cluster_tets], volumes[cluster_tets]
            )
            raw_layers.append((float(isovalue), cluster_contours, region_ids, local_frame))

    unresolved: set[int] = set()
    reversed_final: set[int] = set()

    if params.support_driven_order is not None:
        # Ordre pilote par le SUPPORT REEL - prioritaire sur les deux tris
        # ci-dessous, qui ne connaissent que l'isovaleur/la continuite
        # spatiale et jamais ce qui porte reellement la matiere. `up` = +Z,
        # convention deja actee dans tout le projet (le plateau est a
        # z=platform_z, l'impression progresse en Z croissant).
        support_order, unresolved = _support_driven_order(
            raw_layers, params, np.array([0.0, 0.0, 1.0])
        )
        raw_layers = [raw_layers[i] for i in support_order]
        # `unresolved` indexait `raw_layers` AVANT permutation : le reindexer
        # sur les positions finales, sinon le diagnostic designe des couches
        # arbitraires.
        unresolved = {new_i for new_i, old_i in enumerate(support_order) if old_i in unresolved}
    elif params.frontier_order:
        # Croissance de front : le deja-imprime doit toujours etre connexe.
        # Voir `LayerParams.frontier_order` - contrainte DURE d'adjacence, la
        # geodesique ne servant qu'a departager les couches atteignables.
        progress = _geodesic_layer_progress(
            raw_layers, tet_mesh, mesh_data, region_dependency_graph
        )
        order, reversed_before = _frontier_order(raw_layers, progress, gap_threshold)
        raw_layers = [raw_layers[i] for i in order]
        # Les indices du parcours designent les couches AVANT reordonnancement :
        # les retraduire dans la numerotation finale, sans quoi le drapeau
        # atterrirait sur une autre couche.
        reversed_final = {
            new_i for new_i, old_i in enumerate(order) if old_i in reversed_before
        }
    elif params.geodesic_order:
        # La PROGRESSION vient du plus court chemin jusqu'au plateau a travers
        # la matiere ; la CONTINUITE SPATIALE reste celle du parcours DFS
        # existant. Les deux sont necessaires : la progression seule entrelace
        # des couches eloignees (mesure x4 de degradation), la continuite
        # seule suit l'isovaleur meme quand celle-ci part de l'extremite d'une
        # coque suspendue au lieu de son attache.
        progress = _geodesic_layer_progress(
            raw_layers, tet_mesh, mesh_data, region_dependency_graph
        )
        geodesic_order = _spatial_dfs_order(
            raw_layers, region_ranks, gap_threshold, progress
        )
        raw_layers = [raw_layers[i] for i in geodesic_order]
    elif params.split_disjoint_clusters is not None:
        # tri spatial DFS (Decision 8) : degenere exactement vers le tri
        # historique (rank, isovalue) quand le clustering ne produit jamais
        # plus d'un cluster par isovaleur (graphe = chaine simple sans
        # branchement) - pas de parametre separe, voir docstring de
        # split_disjoint_clusters.
        spatial_order = _spatial_dfs_order(raw_layers, region_ranks, gap_threshold)
        raw_layers = [raw_layers[i] for i in spatial_order]
    else:
        def _sort_key(entry: tuple) -> tuple[int, float]:
            isovalue, _contours, region_ids, _frame = entry
            rank = max((region_ranks.get(r, 0) for r in region_ids), default=0)
            return (rank, isovalue)

        raw_layers.sort(key=_sort_key)

    return [
        Layer(
            id=i,
            contour_3d=contours,
            region_ids=region_ids,
            local_frame=local_frame,
            isovalue=isovalue,
            support_unresolved=i in unresolved,
            deposition_reversed=i in reversed_final,
        )
        for i, (isovalue, contours, region_ids, local_frame) in enumerate(raw_layers)
    ]
