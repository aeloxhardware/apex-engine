"""Etape 5 : marching tetrahedra (isosurface par lots) et soudure des sommets.

Découpé de `step5_layers.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step5_layers`.
"""

from __future__ import annotations

import numpy as np
from core.mesh_types import TetMesh
from core import _cpp


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
