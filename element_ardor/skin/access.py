"""
Accessibilite de P2 en temps raisonnable (2026-11-19) : grille de phi et faces representatives.

Avant : chaque point de la tete (quelques centaines par pose, trois coques pour la marge) etait localise dans les
tetraedres par un arbre k-d et des coordonnees barycentriques, pour CHAQUE face de la surface (372 000 sur un Benchy) :
32 s pour 2 000 faces, une heure et demie pour le Benchy entier sous Windows (pas de `fork`, un seul coeur).

Ici :
  * `PhiGrid` : phi aux centres d'une grille reguliere (0,5 mm), +inf hors de la piece, calcule UNE fois en testant
    chaque tetraedre contre les centres de sa boite englobante ; un point de la tete devient une lecture de tableau ;
  * `representatives` : les faces candidates sont regroupees par cellule (1 mm) et par orientation (normales a moins
    de ~15°) ; l'accessibilite est evaluee sur un representant par groupe et etendue aux autres. L'ecart tete /
    normale d'une face du groupe est majore par celui du representant plus l'angle entre leurs normales.

Derive de S4 Slicer (Joshua Bird, GPL-3.0, https://github.com/jyjblrd/S4_Slicer).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_PAIR_CHUNK = 2_000_000


@dataclass(frozen=True)
class PhiGrid:
    """phi de la piece echantillonne aux centres d'une grille : `values[i, j, k]` au point
    `origin + (i + 0,5, j + 0,5, k + 0,5) * cell`, +inf hors de la piece."""
    origin: np.ndarray
    cell: float
    values: np.ndarray

    @classmethod
    def build(cls, tet_vertices: np.ndarray, tetrahedra: np.ndarray, phi_nodes: np.ndarray, cell: float) -> "PhiGrid":
        lo = tet_vertices.min(axis=0) - cell
        shape = np.ceil((tet_vertices.max(axis=0) + cell - lo) / cell).astype(int)
        values = np.full(int(np.prod(shape)), np.inf)
        tets = tet_vertices[tetrahedra]
        a = tets[:, 0]
        m = np.stack([tets[:, 1] - a, tets[:, 2] - a, tets[:, 3] - a], axis=2)
        det = np.linalg.det(m)
        ok = np.abs(det) > 1e-12
        inv = np.zeros_like(m)
        inv[ok] = np.linalg.inv(m[ok])
        i0 = np.ceil((tets.min(axis=1) - lo) / cell - 0.5).astype(int)
        i1 = np.floor((tets.max(axis=1) - lo) / cell - 0.5).astype(int)
        span = np.maximum(i1 - i0 + 1, 0)
        span[~ok] = 0
        count = span.prod(axis=1)
        phi_t = phi_nodes[tetrahedra]
        order = np.flatnonzero(count > 0)
        start = 0
        while start < len(order):
            stop = start + max(1, int(np.searchsorted(np.cumsum(count[order[start:]]), _PAIR_CHUNK)))
            ids = order[start:stop]
            start = stop
            n = count[ids]
            t = np.repeat(ids, n)
            local = np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)
            sy, sz = span[t, 1], span[t, 2]
            ijk = i0[t] + np.stack([local // (sy * sz), (local // sz) % sy, local % sz], axis=1)
            centres = lo + (ijk + 0.5) * cell
            w = np.einsum("nij,nj->ni", inv[t], centres - a[t])
            bary = np.concatenate([1.0 - w.sum(axis=1, keepdims=True), w], axis=1)
            inside = bary.min(axis=1) >= -1e-9
            flat = np.ravel_multi_index(ijk[inside].T, shape)
            np.minimum.at(values, flat, np.einsum("nv,nv->n", bary[inside], phi_t[t[inside]]))
        return cls(origin=lo, cell=float(cell), values=values.reshape(shape))

    def lookup(self, points: np.ndarray) -> np.ndarray:
        """phi au centre de la cellule de chaque point (+inf hors de la piece ou de la grille)."""
        ijk = np.floor((points - self.origin) / self.cell).astype(np.int64)
        shape = np.array(self.values.shape)
        inside = np.all((ijk >= 0) & (ijk < shape), axis=-1)
        out = np.full(points.shape[:-1], np.inf)
        out[inside] = self.values[tuple(ijk[inside].T)]
        return out


def representatives(centers: np.ndarray, normals: np.ndarray, ids: np.ndarray, cell: float,
                    normal_bins: int = 4) -> tuple[np.ndarray, np.ndarray]:
    """Groupes de faces `ids` (meme cellule de `cell` mm, normales quantifiees a 1/`normal_bins`) : rend
    `(rep, group)` — `rep[g]` est la face de `ids` la plus proche du centre de son groupe, `group[i]` le groupe de
    `ids[i]`."""
    if len(ids) == 0:
        return np.zeros(0, dtype=int), np.zeros(0, dtype=int)
    key = np.concatenate([np.floor(centers[ids] / cell), np.round(normals[ids] * normal_bins)], axis=1).astype(np.int64)
    _u, group = np.unique(key, axis=0, return_inverse=True)
    group = group.ravel()
    n_groups = int(group.max()) + 1
    mean = np.zeros((n_groups, 3))
    np.add.at(mean, group, centers[ids])
    mean /= np.bincount(group, minlength=n_groups)[:, None]
    d = np.linalg.norm(centers[ids] - mean[group], axis=1)
    order = np.lexsort((d, group))
    first = np.ones(len(order), bool)
    first[1:] = group[order][1:] != group[order][:-1]
    rep = np.zeros(n_groups, dtype=int)
    rep[group[order][first]] = ids[order][first]
    return rep, group


def last_accessible(grid: PhiGrid, head: np.ndarray, rotations: np.ndarray, centers: np.ndarray,
                    isovalues: np.ndarray, margin: float, batch: int = 256) -> tuple[np.ndarray, np.ndarray]:
    """(plateau, dernier rang de couche accessible) pour chaque pose : la tete `head` (repere de la buse) tournee par
    `rotations[f]` et posee en `centers[f]` entre-t-elle dans le plateau (z < `margin`) ? Quel est le dernier rang
    de couche imprime avant que la matiere ({phi <= isovaleur}) n'atteigne un de ses points ?"""
    plate = np.zeros(len(centers), bool)
    last = np.full(len(centers), len(isovalues) - 1)
    for s in range(0, len(centers), batch):
        pts = np.einsum("fij,pj->fpi", rotations[s:s + batch], head) + centers[s:s + batch, None, :]
        plate[s:s + batch] = (pts[:, :, 2] < margin).any(axis=1)
        low = grid.lookup(pts).min(axis=1)
        ins = np.isfinite(low)
        part = last[s:s + batch]
        part[ins] = np.searchsorted(isovalues, low[ins] - 1e-9, side="left") - 1
        last[s:s + batch] = part
    return plate, last
