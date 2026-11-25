"""
Etape 5e : contours des couches COURBES pris sur la surface du STL d'origine (2026-11-25, limite #92).

Le contour d'une couche est la courbe ou son isosurface phi = c coupe la surface de la piece. Ardor le prenait sur le
maillage tetraedrique, qui n'epouse le STL qu'a l'echelle de ses mailles (limite #86) : bords crenele des couches
inclinees, parois en dents de scie la ou une couche courbe rencontre une paroi. L'etape 5d le corrige pour les couches
PLANES en recoupant le STL par un plan ; ici, la meme idee pour toutes les autres : phi est evalue aux sommets du STL
(fonction lineaire du tetraedre qui contient le sommet, extrapolee depuis le plus proche sinon), puis la ligne de niveau
phi = c est extraite triangle par triangle sur la surface. Pour une couche plane (phi = z), c'est la section plane.

- Le STL est d'abord subdivise (aretes <= `subdivide`) : sur une grande face plane, phi n'est pas lineaire la ou les
  couches se plient, et la ligne de niveau interpolee entre trois sommets eloignes s'ecarte de la vraie (1,2 mm mesure
  sur le toit du Benchy en couches conformes).
- Chaque boucle de l'isosurface est remplacee par la boucle FERMEE de la ligne de niveau la plus proche (p99 des
  distances, pire des deux sens), acceptee sous `match_tol`, sous `loose_match_tol` si la suivante est au moins
  `ambiguity_ratio` fois plus loin, ou sous `far_match_tol` si elle l'est `far_ratio` fois.
- Ne comptent PAS comme rivales : une boucle deja prise (sous `match_tol`) par une autre boucle de la meme couche, et une
  boucle qui court surtout (> `rival_tangent_fraction`) sur des faces ou la couche est presque tangente a la surface
  (gradient de phi le long de la face < `tangent_gradient`) : la ligne de niveau y est mal conditionnee (dessus d'une
  couche conforme, parallele au toit).
- Une boucle sans correspondante garde le contour de l'isosurface. Le repere local et l'isovaleur ne changent pas.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import trimesh
from scipy.spatial import cKDTree

from core.mesh_types import Layer, ScalarField, TetMesh


@dataclass
class SurfaceLevelParams:
    subdivide: float = 1.0             # longueur d'arete maximale du STL avant evaluation de phi (mm) ; 0 = telle quelle
    match_tol: float = 0.5             # p99 des distances, pire des deux sens (mm)
    loose_match_tol: float = 1.0       # accepte si la suivante est au moins ambiguity_ratio fois plus loin (mm)
    ambiguity_ratio: float = 2.5
    far_match_tol: float = 3.0         # accepte si la suivante est au moins far_ratio fois plus loin (mm)
    far_ratio: float = 5.0
    tangent_gradient: float = 0.17     # |grad phi| le long d'une face sous lequel la couche lui est presque tangente
    rival_tangent_fraction: float = 0.3
    sample_step: float = 0.05          # pas d'echantillonnage des boucles pour la comparaison (mm)
    simplify_tol: float = 0.01         # Douglas-Peucker 3D des boucles retenues (mm)
    # Segments redecoupes a cette longueur apres simplification (mm) ; 0 = aucun. Sur une couche COURBE, seuls les
    # sommets sont poses sur la couche : un long segment droit en coupe la corde. Mesure (death, longueur deposee sans
    # appui ; isosurface : 2 097 mm) : simplifie sans borne 3 956 mm, 0,2 mm 2 014 mm. Avec les couches conformes : 1 mm
    # 2 076 mm, 0,5 mm 1 949 mm, 0,3 mm 1 862 mm ; G-code du Benchy +19 %, +47 %, +88 % de lignes. 0,3 : aucun des
    # quatre modeles en recul (l_shape 338 -> 337 mm, l-triple 672 -> 672 mm) ; 0,5 : l_shape +2,3 %, l-triple +2,9 %.
    max_segment: float = 0.3
    # ... seulement si la boucle s'ecarte du plan de sa couche de plus de ceci (mm) : une boucle plane (meme inclinee) a un
    # decalage hors plan constant, que les parois interieures interpolent exactement sans points de plus.
    curved_tol: float = 0.05
    close_tol: float = 1e-3            # une ligne dont les extremites sont plus loin que ceci est ouverte (mm)
    neighbours: int = 8                # tetraedres essayes par sommet du STL


@dataclass
class SurfaceLevelReport:
    loops: int = 0                     # boucles examinees (non recoupees par l'etape 5d)
    replaced: int = 0


def phi_at(points: np.ndarray, tet_mesh: TetMesh, field: ScalarField, neighbours: int = 8,
           chunk: int = 150_000) -> np.ndarray:
    """phi aux `points` : fonction lineaire du tetraedre (parmi les `neighbours` de centre le plus proche) qui les
    contient le mieux (plus grande coordonnee barycentrique minimale), extrapolee si aucun ne les contient."""
    tets = tet_mesh.tet_vertices[tet_mesh.tetrahedra]
    tree = cKDTree(tets.mean(axis=1))
    out = np.empty(len(points))
    for start in range(0, len(points), chunk):
        p = points[start:start + chunk]
        _d, idx = tree.query(p, k=neighbours, workers=-1)
        idx = idx.reshape(len(p), -1)
        a = tets[idx, 0]
        m = np.stack([tets[idx, 1] - a, tets[idx, 2] - a, tets[idx, 3] - a], axis=-1)
        rhs = p[:, None, :] - a
        try:
            lam = np.linalg.solve(m, rhs[..., None])[..., 0]
        except np.linalg.LinAlgError:       # tetraedre degenere parmi les voisins
            lam = np.einsum("nkij,nkj->nki", np.linalg.pinv(m), rhs)
        bary = np.concatenate([1.0 - lam.sum(-1, keepdims=True), lam], axis=-1)
        # un tetraedre plat (volume ~ 0) donne des coordonnees sans valeur : jamais retenu
        size = np.linalg.norm(m, axis=-2).max(-1)
        flat = np.abs(np.linalg.det(m)) <= 1e-9 * np.maximum(size, 1e-12) ** 3
        best = np.argmax(np.where(flat, -np.inf, np.nan_to_num(bary.min(-1), nan=-np.inf)), axis=1)
        rows = np.arange(len(p))
        values = field.values[tet_mesh.tetrahedra[idx[rows, best]]]
        out[start:start + chunk] = (bary[rows, best] * values).sum(-1)
    return out


def surface_gradient(vertices: np.ndarray, faces: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Norme du gradient de phi (P1) le long de chaque face."""
    a, b, c = vertices[faces[:, 0]], vertices[faces[:, 1]], vertices[faces[:, 2]]
    e1, e2 = b - a, c - a
    n = np.cross(e1, e2)
    area2 = np.maximum(np.einsum("ij,ij->i", n, n), 1e-30)
    d1, d2 = phi[faces[:, 1]] - phi[faces[:, 0]], phi[faces[:, 2]] - phi[faces[:, 0]]
    g = (np.cross(n, e2) * d1[:, None] - np.cross(n, e1) * d2[:, None]) / area2[:, None]
    return np.linalg.norm(g, axis=1)


def level_segments(vertices: np.ndarray, faces: np.ndarray, phi: np.ndarray, c: float):
    """(segments (S,2,3), faces coupees (S,)) de la ligne de niveau phi = c, interpolee sur les aretes."""
    pf = phi[faces] - c
    above = pf >= 0.0                      # un sommet pile sur c compte au-dessus : 0 ou 2 aretes coupees par face
    hit = np.flatnonzero(above.any(1) & ~above.all(1))
    f, p, a = faces[hit], pf[hit], above[hit]
    pts = np.zeros((len(hit), 2, 3))
    count = np.zeros(len(hit), dtype=int)
    for i, j in ((0, 1), (1, 2), (2, 0)):
        s = np.flatnonzero(a[:, i] != a[:, j])
        t = p[s, i] / (p[s, i] - p[s, j])
        pts[s, count[s]] = vertices[f[s, i]] + (vertices[f[s, j]] - vertices[f[s, i]]) * t[:, None]
        count[s] += 1
    ok = count == 2
    return pts[ok], hit[ok]


def closed_loops(segments: np.ndarray, close_tol: float) -> list[np.ndarray]:
    if len(segments) == 0:
        return []
    try:
        path = trimesh.load_path(segments)
    except Exception:                       # ligne degeneree : aucune boucle sure
        return []
    loops = [np.asarray(q, dtype=float) for q in path.discrete]
    return [q for q in loops if len(q) >= 4 and np.linalg.norm(q[0] - q[-1]) <= close_tol]


def simplify(loop: np.ndarray, tol: float) -> np.ndarray:
    """Douglas-Peucker 3D (iteratif) d'une boucle fermee (premier point = dernier)."""
    n = len(loop)
    if tol <= 0 or n < 5:
        return loop
    keep = np.zeros(n, dtype=bool)
    far = int(np.argmax(np.linalg.norm(loop - loop[0], axis=1)))
    keep[[0, far, n - 1]] = True
    stack = [(0, far), (far, n - 1)]
    while stack:
        i, j = stack.pop()
        if j - i < 2:
            continue
        ab = loop[j] - loop[i]
        q = loop[i + 1:j] - loop[i]
        t = np.clip(q @ ab / max(float(ab @ ab), 1e-24), 0.0, 1.0)
        d = np.linalg.norm(q - t[:, None] * ab, axis=1)
        m = int(np.argmax(d))
        if d[m] > tol:
            keep[i + 1 + m] = True
            stack += [(i, i + 1 + m), (i + 1 + m, j)]
    return loop[keep]


def _densify(loop: np.ndarray, step: float) -> np.ndarray:
    seg = np.diff(loop, axis=0)
    n = np.maximum(1, np.ceil(np.linalg.norm(seg, axis=1) / step).astype(int))
    i = np.repeat(np.arange(len(seg)), n)
    t = (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)) / np.repeat(n, n)
    return np.vstack([loop[:-1][i] + seg[i] * t[:, None], loop[-1:]])


def _finish(loop: np.ndarray, frame, params: SurfaceLevelParams) -> np.ndarray:
    """Boucle retenue : simplifiee, puis redecoupee a `max_segment` si elle est courbe dans le repere de sa couche."""
    loop = simplify(loop, params.simplify_tol)
    if params.max_segment <= 0:
        return loop
    if frame is not None and params.curved_tol > 0:
        origin, u, v = (np.asarray(f, dtype=float) for f in frame)
        n = np.cross(u, v)
        offset = (loop - origin) @ (n / max(np.linalg.norm(n), 1e-12))
        if np.ptp(offset) < params.curved_tol:
            return loop
    return _densify(loop, params.max_segment)


@dataclass
class _Candidates:
    loops: list
    dense: list
    trees: list
    tangent: list                    # part de chaque boucle sur des faces presque tangentes a la couche


def _candidates(segments, hit, tangent_face, params: SurfaceLevelParams) -> _Candidates:
    loops = closed_loops(segments, params.close_tol)
    dense = [_densify(q, params.sample_step) for q in loops]
    tangent = []
    if loops:
        mid = cKDTree(segments.mean(axis=1))
        tangent = [float(tangent_face[hit][mid.query(q)[1]].mean()) for q in loops]
    return _Candidates(loops, dense, [cKDTree(q) for q in dense], tangent)


def _scores(loop: np.ndarray, cand: _Candidates, params: SurfaceLevelParams) -> list[tuple[float, int]]:
    own = _densify(np.asarray(loop, dtype=float), params.sample_step)
    own_tree = cKDTree(own)
    out = []
    for i, (pts, tree) in enumerate(zip(cand.dense, cand.trees)):
        score = float(np.percentile(tree.query(own)[0], 99))
        if score <= params.far_match_tol:
            score = max(score, float(np.percentile(own_tree.query(pts)[0], 99)))
        out.append((score, i))
    return sorted(out)


def _accepted(scores: list, j: int, taken: dict, cand: _Candidates, params: SurfaceLevelParams) -> int | None:
    """Indice de la boucle retenue pour la boucle `j` de la couche, ou None."""
    if not scores:
        return None
    best, i = scores[0]
    rivals = [s for s, r in scores[1:]
              if cand.tangent[r] <= params.rival_tangent_fraction and taken.get(r, j) == j]
    second = rivals[0] if rivals else np.inf
    if (best <= params.match_tol
            or (best <= params.loose_match_tol and second >= params.ambiguity_ratio * best)
            or (best <= params.far_match_tol and second >= params.far_ratio * best)):
        return i
    return None


def run(layers: list[Layer], mesh_data, tet_mesh: TetMesh, scalar_field: ScalarField, skip: set,
        params: SurfaceLevelParams = SurfaceLevelParams()) -> tuple[list[Layer], SurfaceLevelReport]:
    """Remplace chaque boucle (k, j) hors de `skip` (boucles deja recoupees par l'etape 5d) par sa ligne de niveau."""
    report = SurfaceLevelReport()
    todo = {k: [j for j in range(len(L.contour_3d or [])) if (k, j) not in skip]
            for k, L in enumerate(layers) if L.contour_3d and L.isovalue is not None}
    todo = {k: js for k, js in todo.items() if js}
    if not todo:
        return layers, report
    vertices, faces = np.asarray(mesh_data.vertices, dtype=float), np.asarray(mesh_data.faces)
    if params.subdivide > 0:
        vertices, faces = trimesh.remesh.subdivide_to_size(vertices, faces, params.subdivide)
    phi = phi_at(vertices, tet_mesh, scalar_field, params.neighbours)
    tangent_face = surface_gradient(vertices, faces, phi) < params.tangent_gradient
    out = list(layers)
    for k, js in todo.items():
        segments, hit = level_segments(vertices, faces, phi, layers[k].isovalue)
        cand = _candidates(segments, hit, tangent_face, params)
        contours = list(layers[k].contour_3d)
        scores = {j: _scores(contours[j], cand, params) for j in js}
        taken = {s[0][1]: j for j, s in scores.items() if s and s[0][0] <= params.match_tol}
        report.loops += len(js)
        for j in js:
            i = _accepted(scores[j], j, taken, cand, params)
            if i is not None:
                contours[j] = _finish(cand.loops[i], layers[k].local_frame, params)
                report.replaced += 1
        out[k] = replace(layers[k], contour_3d=contours)
    return out, report
