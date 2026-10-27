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
from core import _cpp
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

    # --- TRI PAR ISOVALEUR SEULE (2026-10-22) ------------------------------
    # Le tri historique `(rang de region, isovaleur)` suppose que l'isovaleur
    # ne dit pas a elle seule l'ordre de croissance, et s'en remet au rang de
    # step3 pour le dire. C'est l'inverse quand le champ EST une distance au
    # plateau a travers la matiere (`use_heat_geodesic_field`) : phi est alors
    # la progression elle-meme, et le rang de region la casse. Mesure sur
    # death.stl (champ de la chaleur) : les couches de rang 2 (phi 23,6 a 39,6)
    # passaient avant celles de rang 3 a phi 8,3, donc des pieces imprimees
    # avant la matiere qui les relie au reste (correlation ordre/isovaleur 0,78).
    #
    # True (defaut) : le champ est la distance a travers la matiere. False =
    # tri (rang de region, isovaleur).
    isovalue_order: bool = True


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
    if _cpp.enabled() and hasattr(_cpp.module(), "tet_volumes"):
        return _cpp.module().tet_volumes(np.ascontiguousarray(tet_mesh.tet_vertices, dtype=float),
                                         np.ascontiguousarray(tet_mesh.tetrahedra, dtype=np.int64))
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

# Nombre de couches ou le trace natif a rendu la main a trimesh (bord
# non-variete) - mesure du taux de repli, verifie par les tests d'equivalence.
_NATIVE_OUTLINE_FALLBACKS = 0


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
    # trace de bord : boucles ouvertes, donc jetees. Mesure sur un Benchy
    # (noeuds a 2,6e-5 et 9,2e-5 mm du plan) : couche entiere perdue (192 mm2)
    # et 23 % d'une autre.
    if len(faces):
        faces = faces[(faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2])
                      & (faces[:, 0] != faces[:, 2])]
    if len(faces) == 0:
        return []

    # Pas 14 du portage C++ : trace de bord natif, identique a outline() dans le
    # cas des boucles simples ; None = bord non-variete -> trimesh ci-dessous.
    global _NATIVE_OUTLINE_FALLBACKS
    if _cpp.enabled() and hasattr(_cpp.module(), "boundary_loops"):
        loops = _cpp.module().boundary_loops(np.ascontiguousarray(faces, dtype=np.int64), len(vertices))
        if loops is not None:
            return [vertices[loop] for loop in loops if len(loop) >= 4]
        _NATIVE_OUTLINE_FALLBACKS += 1

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
# Point d'entree principal de l'etape
# ---------------------------------------------------------------------------

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
    tet_regions_c = np.ascontiguousarray(tet_regions, dtype=np.int64)
    gradient_c = np.ascontiguousarray(scalar_field.gradient, dtype=float)
    volumes_c = np.ascontiguousarray(volumes, dtype=float)

    scale = float(np.ptp(tet_mesh.tet_vertices, axis=0).max()) or 1.0
    weld_tolerance = params.weld_tolerance_factor * scale

    # Seuil de connexion spatiale entre boucles (clustering et tri par
    # continuite) : quelques hauteurs de couche separent le prolongement direct
    # d'un cluster de ce qui ne l'est pas. Valeur mesuree sur le modele en T.
    gap_threshold = (params.split_disjoint_clusters if params.split_disjoint_clusters is not None else 3.0) \
        * params.layer_height_base

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

    # Pas 13 du portage C++ (CLAUDE.md limite #56) : marching tetrahedra +
    # soudure en natif, construit UNE fois par champ. Equivalence au bit
    # (tests/test_cpp_equivalence.py) ; chemin Python inchange sinon.
    native_extractor = None
    if _cpp.enabled() and hasattr(_cpp.module(), "IsosurfaceExtractor"):
        native_extractor = _cpp.module().IsosurfaceExtractor(
            np.ascontiguousarray(tet_mesh.tet_vertices, dtype=float),
            np.ascontiguousarray(tet_mesh.tetrahedra, dtype=np.int64),
            np.ascontiguousarray(scalar_field.values, dtype=float),
        )

    for isovalue in isovalues:
        if (native_extractor is not None and (params.split_disjoint_clusters is None)
                and hasattr(native_extractor, "extract_layer")):
            # Chemin par defaut en un appel natif : isosurface, contours, regions et
            # repere local. Statut 2 (bord non simple) : on retombe sur trimesh ci-dessous.
            layer = native_extractor.extract_layer(float(isovalue), weld_tolerance, gradient_c, volumes_c,
                                                   tet_regions_c)
            if isinstance(layer, tuple):
                contours, region_ids, local_frame = layer
                raw_layers.append((float(isovalue), contours, list(region_ids), local_frame))
                continue
            if layer == 0:
                continue
        if native_extractor is not None:
            welded_vertices, welded_faces, source_tets = native_extractor.extract(float(isovalue), weld_tolerance)
            if len(welded_faces) == 0:
                continue
        else:
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

    if params.split_disjoint_clusters is not None:
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
            if params.isovalue_order:
                return (0, isovalue)
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
        )
        for i, (isovalue, contours, region_ids, local_frame) in enumerate(raw_layers)
    ]
