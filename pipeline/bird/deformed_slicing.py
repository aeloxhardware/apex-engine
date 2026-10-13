"""
Tranchage Bird DANS L'ESPACE DEFORME, puis retour de la trajectoire dans la
piece reelle (porte de research/S4_Slicer/bird_slice_deformed.py et
bird_to_apex.py).

Pourquoi : dans l'espace deforme, chaque couche est EXACTEMENT un plan
z = constante, donc les etapes 5 et 6 d'Apex y travaillent sans
approximation (l'astuce de Bird, qui tranche des plans avec Cura, sans Cura).
Retour : un point est localise dans son tetraedre deforme (coordonnees
barycentriques) et replace au meme endroit relatif dans le tetraedre
d'origine (exact : deformation lineaire par tetraedre), sauf dans les plis ou
deux tetraedres deformes se recouvrent ; une polyligne choisit alors la suite
de tetraedres la plus continue (Viterbi).

Ponts desactives (sans sens physique dans l'espace deforme).

Derive de S4 Slicer (Joshua Bird, GPL-3.0, https://github.com/jyjblrd/S4_Slicer).
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
import trimesh
from scipy.spatial import cKDTree

from core.mesh_types import OverhangMap, Region, RegionDependencyGraph, ScalarField, TetMesh
from pipeline.step1_ingest import _build_mesh_data
from pipeline.step5_layers import LayerParams, run as layers_run
from pipeline.step6_toolpath2d import Step6Params, _lift_to_3d, run as toolpath2d_run

_OPPOSITE = [((1, 2, 3), 0), ((0, 2, 3), 1), ((0, 1, 3), 2), ((0, 1, 2), 3)]


@dataclass
class DeformedSlicingParams:
    # Longueur maximale d'un segment avant son retour (bien plus court qu'une arete de tetraedre).
    max_segment_mm: float = 0.5
    # Un point a plus de cut_factor x largeur de buse HORS du volume deforme n'est pas depose
    # (fausse zone pleine d'un pli). 2,0 depuis le 2026-10-02 (CLAUDE.md limite #68).
    cut_factor: float = 2.0
    candidates: int = 24                 # tetraedres candidats par point (k plus proches centres)
    refine_rounds: int = 6
    layers: LayerParams = field(default_factory=LayerParams)
    toolpath2d: Step6Params = field(default_factory=Step6Params)


def barycentric(points, tets_xyz):
    """Coordonnees barycentriques de chaque point dans son tetraedre, (K,4)."""
    a = tets_xyz[:, 0]
    m = np.stack([tets_xyz[:, 1] - a, tets_xyz[:, 2] - a, tets_xyz[:, 3] - a], axis=2)
    lam = np.linalg.solve(m, (points - a)[:, :, None])[:, :, 0]
    return np.column_stack([1.0 - lam.sum(axis=1), lam])


def tet_gradients(nodes, tetrahedra, values):
    """Gradient exact d'un champ lineaire par tetraedre (P1), (T,3)."""
    p = nodes[tetrahedra]
    m = np.stack([p[:, 1] - p[:, 0], p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]], axis=1)
    dv = values[tetrahedra[:, 1:]] - values[tetrahedra[:, :1]]
    return np.linalg.solve(m, dv[:, :, None])[:, :, 0]


def boundary_faces(cells, reference_points):
    """Faces de tetraedres qui n'apparaissent qu'une fois, orientees vers
    l'exterieur sur la piece d'ORIGINE (toujours valide, contrairement a
    l'espace deforme ou quelques tetraedres se retournent)."""
    tri = np.concatenate([cells[:, list(f)] for f, _ in _OPPOSITE])
    opp = np.concatenate([cells[:, o] for _, o in _OPPOSITE])
    keys = np.sort(tri, axis=1)
    _, first, counts = np.unique(keys, axis=0, return_index=True, return_counts=True)
    keep = first[counts == 1]
    tri, opp = tri[keep], opp[keep]
    p = reference_points
    normal = np.cross(p[tri[:, 1]] - p[tri[:, 0]], p[tri[:, 2]] - p[tri[:, 0]])
    inward = np.einsum("ij,ij->i", normal, p[opp] - p[tri[:, 0]]) > 0
    tri[inward] = tri[inward][:, [0, 2, 1]]
    return tri


class DeformedToOriginal:
    """Ramene des points de l'espace deforme vers la piece d'origine."""

    def __init__(self, cells, deformed, original, k=24):
        self.cells, self.deformed, self.original, self.k = cells, deformed, original, k
        self.tets = deformed[cells]
        self.tree = cKDTree(self.tets.mean(axis=1))

    def __call__(self, points):
        _, cand = self.tree.query(points, k=self.k)
        best_score = np.full(len(points), -np.inf)
        best_bary = np.zeros((len(points), 4))
        best_tet = np.zeros(len(points), dtype=np.int64)
        containing = np.zeros(len(points), dtype=np.int64)
        for j in range(self.k):
            t = cand[:, j]
            with np.errstate(all="ignore"):
                bary = barycentric(points, self.tets[t])
            score = np.nan_to_num(bary.min(axis=1), nan=-np.inf)
            containing += score > 1e-6
            better = score > best_score
            best_score[better], best_bary[better], best_tet[better] = score[better], bary[better], t[better]
        mapped = np.einsum("ij,ijk->ik", best_bary, self.original[self.cells[best_tet]])
        return mapped, best_tet, best_score, containing

    def _candidates(self, points, outside_margin=0.1):
        """Pour chaque point : ses k tetraedres candidats, coordonnees barycentriques,
        points ramenes, score, validite. Calcul point par point (independant du lot) :
        requete k-NN sur tous les coeurs (`workers=-1`) et barycentriques par paquets
        de points dans des fils (numpy libere le GIL) -- resultat identique au bit
        (2026-10-13, x4 mesure sur 200 000 points)."""
        _, cand = self.tree.query(points, k=self.k, workers=-1)
        native = _barycentric_native()
        if native is not None:
            if getattr(self, "_original_tets", None) is None:
                self._original_tets = np.ascontiguousarray(self.original[self.cells], dtype=float)
                self._deformed_tets = np.ascontiguousarray(self.tets, dtype=float)
            score, mapped = native(np.ascontiguousarray(points, dtype=float), self._deformed_tets,
                                   self._original_tets, np.ascontiguousarray(cand, dtype=np.int64))
        else:
            score, mapped = _parallel_barycentric(points, self.tets, cand, self.original, self.cells)
        valid = score > -1e-6
        none = ~valid.any(axis=1)
        if none.any():
            best = score[none].max(axis=1, keepdims=True)
            valid[none] = score[none] >= best - outside_margin
        containing = (score > 1e-6).sum(axis=1)
        return cand, score, mapped, valid, containing

    def _viterbi(self, cand, score, mapped, valid, containing):
        n = len(mapped)
        native = _viterbi_native()
        if native is not None:
            choice = native(np.ascontiguousarray(mapped), np.ascontiguousarray(valid))
            rows = np.arange(n)
            return mapped[rows, choice], cand[rows, choice], score[rows, choice], containing
        cost = np.where(valid[0], 0.0, np.inf)
        back = np.zeros((n, self.k), dtype=np.int64)
        for i in range(1, n):
            step = np.linalg.norm(mapped[i - 1][:, None, :] - mapped[i][None, :, :], axis=2)
            total = cost[:, None] + step
            back[i] = total.argmin(axis=0)
            cost = np.where(valid[i], total.min(axis=0), np.inf)
        choice = np.zeros(n, dtype=np.int64)
        choice[-1] = int(cost.argmin())
        for i in range(n - 1, 0, -1):
            choice[i - 1] = back[i, choice[i]]
        rows = np.arange(n)
        return mapped[rows, choice], cand[rows, choice], score[rows, choice], containing

    def map_polyline(self, points, outside_margin=0.1):
        """Comme __call__, mais le tetraedre de chaque point est choisi pour que
        le chemin ramene soit le plus CONTINU (Viterbi exact) : dans un pli, un
        choix point par point sautait d'une epaisseur a l'autre (2 a 15 mm)."""
        return self._viterbi(*self._candidates(points, outside_margin))

    def map_polylines(self, polylines, outside_margin=0.1):
        """map_polyline sur plusieurs polylignes : candidats calcules en UN lot (le
        calcul de chaque point est independant), Viterbi polyligne par polyligne.
        Identique au bit a des appels separes."""
        polylines = [np.asarray(p, dtype=float) for p in polylines]
        if not polylines:
            return []
        cand, score, mapped, valid, containing = self._candidates(np.vstack(polylines), outside_margin)
        out, start = [], 0
        for p in polylines:
            sl = slice(start, start + len(p))
            out.append(self._viterbi(cand[sl], score[sl], mapped[sl], valid[sl], containing[sl]))
            start += len(p)
        return out

_BARY_CHUNK = 20_000
_BATCH_POINTS = 200_000


def _parallel_barycentric(points, tets, cand, original, cells):
    """(score (N,k), points ramenes (N,k,3)) : pour chaque candidat j, barycentriques
    barycentric(points, tets[cand[:, j]]), score = plus petite coordonnee (NaN -> -inf),
    point ramene = combinaison des sommets d'ORIGINE. Par paquets de points repartis sur
    des fils ; chaque ligne est calculee exactement comme en serie (meme resolution 3x3
    par LAPACK, memes operations elementaires) : identique au bit."""
    n, k = cand.shape
    score = np.empty((n, k))
    mapped = np.empty((n, k, 3))

    def block(sl):
        bary = np.empty((sl.stop - sl.start, k, 4))
        with np.errstate(all="ignore"):
            for j in range(k):
                bary[:, j] = barycentric(points[sl], tets[cand[sl, j]])
        score[sl] = np.nan_to_num(bary.min(axis=2), nan=-np.inf)
        mapped[sl] = np.einsum("njv,njvk->njk", bary, original[cells[cand[sl]]])

    chunks = [slice(s, min(n, s + _BARY_CHUNK)) for s in range(0, n, _BARY_CHUNK)]
    if len(chunks) <= 1:
        for sl in chunks:
            block(sl)
        return score, mapped
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=min(len(chunks), os.cpu_count() or 1)) as pool:
        list(pool.map(block, chunks))
    return score, mapped


def _barycentric_native():
    """apex_core.tet_barycentric_candidates si `APEX_USE_CPP=1`. TOLERANCE DECLAREE (pas le
    bit) : resolution 3x3 par LU transcrite au lieu de LAPACK, ecart mesure sur le Benchy
    passe 4 (200 000 points x 24 candidats) <= 4,5e-9 sur le score et 1,2e-10 mm sur le
    point ramene, aucune decision changee aux seuils +-1e-6 (CLAUDE.md #77).
    `APEX_NATIVE_BARY=0` garde LAPACK meme avec `APEX_USE_CPP=1` (mesure dos a dos ; lu
    aussi par les sous-processus de P2)."""
    from core import _cpp
    if os.environ.get("APEX_NATIVE_BARY", "1") == "0":
        return None
    module = _cpp.module() if _cpp.enabled() else None
    return getattr(module, "tet_barycentric_candidates", None) if module is not None else None


def _viterbi_native():
    """apex_core.viterbi_choice si `APEX_USE_CPP=1` (identique au bit a la boucle numpy)."""
    from core import _cpp
    module = _cpp.module() if _cpp.enabled() else None
    return getattr(module, "viterbi_choice", None) if module is not None else None


def densify(points, step):
    """Subdivise chaque segment en ceil(longueur / step) morceaux egaux. Version
    vectorisee de la boucle d'origine, memes operations : la longueur est calculee
    segment par segment comme np.linalg.norm d'un vecteur (produit scalaire), puis
    a + (b - a) * (k / n)."""
    points = np.asarray(points, dtype=float)
    if len(points) < 2:
        return np.vstack([points[:1]])
    a, b = points[:-1], points[1:]
    lengths = np.array([np.linalg.norm(v) for v in (b - a)])
    n = np.maximum(1, np.ceil(lengths / step).astype(np.int64))
    seg = np.repeat(np.arange(len(a)), n)
    k = np.arange(len(seg)) - np.repeat(np.cumsum(n) - n, n) + 1
    t = (k / n[seg])[:, None]
    return np.vstack([points[:1], a[seg] + (b - a)[seg] * t])


@dataclass
class DeformedSlicingResult:
    layers: list[dict]          # ordre d'impression : contour_3d, isovalue, toolpath_2d, toolpath_3d, toolpath_score
    stats: dict


def run(cells: np.ndarray, original: np.ndarray, deformed: np.ndarray, offset: np.ndarray,
        params: DeformedSlicingParams = DeformedSlicingParams(),
        carve: Optional[Callable[[list], list]] = None,
        first_layer: Optional[Callable[[list], list]] = None) -> DeformedSlicingResult:
    """`original`/`deformed` : points (N,3) dans le repere Bird ; `offset` :
    repere Bird -> repere de sortie. `carve` (creusement sous le lissage) et
    `first_layer` (motif de premiere couche) : points d'entree du lissage non
    planaire, appliques avant / apres l'etape 6 comme dans la recherche."""
    w = params.toolpath2d.nozzle_width
    # La recherche oriente les faces sur l'ORIGINE puis construit la surface sur le deforme.
    faces = boundary_faces(cells, original)
    surface = trimesh.Trimesh(deformed, faces, process=False)
    surface.remove_unreferenced_vertices()
    mesh_data = _build_mesh_data(surface)
    tet_mesh = TetMesh(tet_vertices=deformed, tetrahedra=cells)
    field_ = ScalarField(values=deformed[:, 2].copy(), gradient=np.tile([0.0, 0.0, 1.0], (len(cells), 1)))
    graph = RegionDependencyGraph(
        regions={0: Region(id=0, face_ids=np.arange(mesh_data.face_count), target_direction=np.array([0.0, 0.0, 1.0]))},
        dependencies=[])
    m = mesh_data.face_count
    no_bridges = OverhangMap(risk_score=np.zeros(m), is_overhang=np.zeros(m, bool), is_bridge=np.zeros(m, bool),
                             bridge_span_direction=np.full((m, 3), np.nan))

    layers = layers_run(tet_mesh, mesh_data, field_, graph, params.layers)
    if carve is not None:
        layers = carve(layers)
    layers = toolpath2d_run(layers, mesh_data, no_bridges, params.toolpath2d)
    if first_layer is not None:
        layers = first_layer(layers)

    back = DeformedToOriginal(cells, deformed, original, k=params.candidates)
    step = params.max_segment_mm

    def map_batched(polylines):
        """map_polylines par lots bornes (~_BATCH_POINTS points) : memoire maitrisee."""
        mapped_, group, group_size = [], [], 0
        for poly in polylines:
            group.append(poly)
            group_size += len(poly)
            if group_size >= _BATCH_POINTS:
                mapped_ += back.map_polylines(group)
                group, group_size = [], 0
        return mapped_ + back.map_polylines(group)

    def refine_all(paths, results):
        """Raffinement de chaque polyligne (milieu insere sur tout segment ramene trop long,
        jusqu'a refine_rounds tours) ; les polylignes a raffiner d'un MEME tour sont
        ramenees en un lot. Chaque polyligne evolue exactement comme seule (calcul par
        point independant, map_polylines identique au bit a des appels separes) ; avant
        le 2026-10-13, un appel par polyligne et par tour (4 383 appels sur le dino)."""
        paths, results = list(paths), list(results)
        active = list(range(len(paths)))
        for _ in range(params.refine_rounds):
            todo = []
            for i in active:
                lengths = np.linalg.norm(np.diff(results[i][0], axis=0), axis=1)
                long = np.flatnonzero(lengths > 2 * step)
                if len(long):
                    pts = paths[i]
                    paths[i] = np.insert(pts, long + 1, (pts[long] + pts[long + 1]) / 2, axis=0)
                    todo.append(i)
            if not todo:
                break
            for i, r in zip(todo, map_batched([paths[i] for i in todo])):
                results[i] = r
            active = todo
        return [(p, *r) for p, r in zip(paths, results)]

    def outside_distance(points, tets, score):
        dist = np.zeros(len(points))
        m_ = score < -1e-6
        if m_.any():
            corners = back.tets[tets[m_]]
            with np.errstate(all="ignore"):
                b = np.clip(barycentric(points[m_], corners), 0.0, None)
            b /= b.sum(axis=1, keepdims=True)
            dist[m_] = np.linalg.norm(points[m_] - np.einsum("nv,nvk->nk", b, corners), axis=1)
        return dist

    out_layers, n_pts, n_out, n_amb, n_cut, n_jump = [], 0, 0, 0, 0, 0
    # Premier retour de TOUTES les couches (trajectoires puis contours) en UN lot : le
    # calcul des candidats est independant par point (resultat identique au bit a des
    # lots par couche), et un grand lot occupe tous les coeurs (2026-10-13).
    all_dense = [[densify(_lift_to_3d(np.asarray(t), L.local_frame), step) for t in (L.toolpath_2d or [])]
                 for L in layers]
    flat = [p for paths in all_dense for p in paths] + [c for L in layers for c in L.contour_3d]
    flat_mapped = map_batched(flat)
    cursor = 0
    first_by_layer = []
    for paths in all_dense:
        first_by_layer.append(flat_mapped[cursor:cursor + len(paths)])
        cursor += len(paths)
    contours_by_layer = []
    for L in layers:
        contours_by_layer.append(flat_mapped[cursor:cursor + len(L.contour_3d)])
        cursor += len(L.contour_3d)
    flat_paths = [p for paths in all_dense for p in paths]
    refined_all = iter(refine_all(flat_paths, [r for firsts in first_by_layer for r in firsts]))
    for L, dense_paths, contour_maps in zip(layers, all_dense, contours_by_layer):
        tp3, tps = [], []
        for _ in dense_paths:
            refined, mapped, tet, score, containing = next(refined_all)
            score = outside_distance(refined, tet, score)
            n_pts += len(mapped)
            n_out += int((score > 0).sum())
            n_amb += int((containing >= 2).sum())
            keep = score <= w * params.cut_factor
            n_cut += int((~keep).sum())
            edges = np.flatnonzero(np.diff(np.concatenate([[0], keep.astype(np.int8), [0]])))
            for start, stop in zip(edges[::2], edges[1::2]):
                run_, run_score = mapped[start:stop], score[start:stop]
                jumps = np.flatnonzero(np.linalg.norm(np.diff(run_, axis=0), axis=1) > 2 * step) + 1
                n_jump += len(jumps)
                for piece, piece_score in zip(np.split(run_, jumps), np.split(run_score, jumps)):
                    if len(piece) >= 2:
                        tp3.append(piece + offset)
                        tps.append(piece_score)
        contours = [r[0] + offset for r in contour_maps]
        out_layers.append({"contour_3d": contours, "isovalue": L.isovalue,
                           "toolpath_2d": [np.asarray(t) for t in (L.toolpath_2d or [])],
                           "toolpath_3d": tp3, "toolpath_score": tps})
    stats = {"layers": len(out_layers), "points": n_pts, "outside": n_out, "ambiguous": n_amb, "cut": n_cut,
             "jumps": n_jump}
    return DeformedSlicingResult(layers=out_layers, stats=stats)
