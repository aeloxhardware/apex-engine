"""
Recouvrement entre tetraedres NON VOISINS -- le repli "global" d'un maillage
deforme.

Idee de l'utilisateur : qu'un tetraedre chevauche ses voisins directs est
NORMAL pendant une deformation Bird (courber une plaque comprime sa face
interieure, il y a "de l'espace en moins"). Un vrai repli, c'est une partie
de la piece qui vient recouvrir une partie ELOIGNEE, avec laquelle elle n'a
aucun contact dans le maillage d'origine.

Mesure :
  - voisins directs (rang 1) = partagent au moins un sommet -> ignores ;
  - rang 2 = un sommet de l'un est un sommet, ou le voisin d'arete d'un
    sommet, de l'autre -> compte separe, pour verifier que le resultat ne
    depend pas trop de cette frontiere ;
  - test exact par axe separateur (SAT) entre deux tetraedres : 4 + 4
    normales de faces, 6 x 6 produits d'aretes = 44 axes ;
  - profondeur de penetration = plus petit recouvrement des projections sur
    les axes (normalises). Distingue un effleurement d'un recouvrement.

Controle : sur un maillage NON deforme, le compte doit etre exactement 0
(dans un maillage conforme, deux tetraedres sans sommet commun sont
disjoints).
"""
import numpy as np
from scipy.spatial import cKDTree

_FACES = np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]])
_EDGES = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])


def _candidate_pairs(points, cells):
    """Paires (i<j) dont les spheres englobantes se recoupent, puis filtre
    AABB. Sur-ensemble conservateur des paires en intersection."""
    tp = points[cells]                                   # (T,4,3)
    centers = tp.mean(axis=1)
    radii = np.linalg.norm(tp - centers[:, None, :], axis=2).max(axis=1)
    tree = cKDTree(centers)
    # Chaque tetraedre cherche a 2 x SON rayon : si r_i >= r_j, la paire est
    # trouvee depuis i (r_i + r_j <= 2 r_i). Sur-ensemble exact, sans le
    # rayon global du plus gros tetraedre (qui rendait tout quadratique).
    lists = tree.query_ball_point(centers, 2.0 * radii)
    rows = np.repeat(np.arange(len(cells)), [len(l) for l in lists])
    cols = np.concatenate([np.asarray(l, dtype=np.int64) for l in lists])
    keep = cols != rows
    rows, cols = np.minimum(rows[keep], cols[keep]), np.maximum(rows[keep], cols[keep])
    pairs = np.unique(rows * len(cells) + cols)
    rows, cols = pairs // len(cells), pairs % len(cells)
    keep = np.linalg.norm(centers[rows] - centers[cols], axis=1) <= radii[rows] + radii[cols]
    rows, cols = rows[keep], cols[keep]
    lo, hi = tp.min(axis=1), tp.max(axis=1)
    keep = np.all((lo[rows] <= hi[cols]) & (lo[cols] <= hi[rows]), axis=1)
    return rows[keep], cols[keep]


def _neighbour_rank(cells, n_points, rows, cols):
    """Rang de voisinage de chaque paire : 1 (sommet commun), 2 (sommet de
    l'un a une arete d'un sommet de l'autre), 3 (au-dela)."""
    # Rang 1 : une des 16 combinaisons de sommets est egale.
    ci, cj = cells[rows], cells[cols]                       # (K,4)
    vi = np.repeat(ci, 4, axis=1)                           # (K,16)
    vj = np.tile(cj, (1, 4))
    shared = (vi == vj).any(axis=1)
    # Rang 2 : une des 16 combinaisons est une arete du maillage. Aretes
    # encodees en entiers tries, test d'appartenance par recherche binaire.
    e = cells[:, _EDGES].reshape(-1, 2)
    keys = np.unique(np.minimum(e[:, 0], e[:, 1]) * n_points + np.maximum(e[:, 0], e[:, 1]))
    q = np.minimum(vi, vj) * n_points + np.maximum(vi, vj)
    pos = np.clip(np.searchsorted(keys, q), 0, len(keys) - 1)
    edge = (keys[pos] == q).any(axis=1)
    return np.where(shared, 1, np.where(edge, 2, 3)).astype(np.int8)


def _sat_penetration(P, Q, tol):
    """Profondeur de penetration par SAT, (K,) -- <= tol si separes."""
    def face_normals(X):
        a, b, c = X[:, _FACES[:, 0]], X[:, _FACES[:, 1]], X[:, _FACES[:, 2]]
        return np.cross(b - a, c - a)                   # (K,4,3)

    eP = P[:, _EDGES[:, 1]] - P[:, _EDGES[:, 0]]        # (K,6,3)
    eQ = Q[:, _EDGES[:, 1]] - Q[:, _EDGES[:, 0]]
    cross = np.cross(eP[:, :, None, :], eQ[:, None, :, :]).reshape(len(P), 36, 3)
    axes = np.concatenate([face_normals(P), face_normals(Q), cross], axis=1)  # (K,44,3)
    norms = np.linalg.norm(axes, axis=2)
    ref = np.maximum(norms.max(axis=1, keepdims=True), 1e-300)
    valid = norms > 1e-9 * ref
    axes = axes / np.where(valid, norms, 1.0)[:, :, None]

    pP = np.einsum("kai,kvi->kav", axes, P)
    pQ = np.einsum("kai,kvi->kav", axes, Q)
    overlap = np.minimum(pP.max(2), pQ.max(2)) - np.maximum(pP.min(2), pQ.min(2))
    overlap = np.where(valid, overlap, np.inf)
    return overlap.min(axis=1)


def non_neighbour_overlaps(points, cells, rel_tol=1e-6):
    """Retourne un dict de mesures. `rel_tol` : tolerance relative a la
    longueur d'arete mediane, en dessous de laquelle une penetration est
    un simple contact numerique."""
    points = np.asarray(points, dtype=np.float64)
    cells = np.asarray(cells, dtype=np.int64)
    tp = points[cells]
    edge_len = np.median(np.linalg.norm(tp[:, _EDGES[:, 1]] - tp[:, _EDGES[:, 0]], axis=2))
    tol = rel_tol * edge_len

    rows, cols = _candidate_pairs(points, cells)
    rank = _neighbour_rank(cells, len(points), rows, cols)
    far = rank >= 2
    rows, cols, rank = rows[far], cols[far], rank[far]

    depth = np.empty(len(rows))
    batch = 100_000
    for s in range(0, len(rows), batch):
        depth[s:s + batch] = _sat_penetration(tp[rows[s:s + batch]], tp[cols[s:s + batch]], tol)
    hit = depth > tol

    out = {"median_edge": float(edge_len)}
    for label, mask in (("rang>=2", hit), ("rang>=3", hit & (rank >= 3))):
        tets = np.unique(np.r_[rows[mask], cols[mask]])
        d = depth[mask]
        out[label] = {
            "pairs": int(mask.sum()),
            "tets": int(len(tets)),
            "depth_median": float(np.median(d)) if len(d) else 0.0,
            "depth_max": float(d.max()) if len(d) else 0.0,
            "tet_ids": tets,
        }
    return out
