"""
Etape 5c — fidelite de la coque : projection des contours sur la surface d'origine (R-030 point 8, piste D).

Les contours d'une couche (`Layer.contour_3d`) sont l'intersection de l'isosurface avec le BORD du maillage
volumique. Avec un mailleur dont le bord n'est qu'une approximation de la surface (TetraNet : reseau regulier ramene
a la peau), ces contours s'ecartent du STL (mesure du 2026-11-07 sur death et turkish : p99 0,23 a 0,30 mm, maximum
0,55 mm, 7 a 8 % des points au-dela de 0,1 mm ; TetGen, qui conserve la surface : p99 0,007 a 0,03 mm).

Ici : chaque sommet de contour est ramene sur le point le plus proche de la surface d'origine, puis, tant qu'un
milieu de segment s'ecarte de la surface de plus de `chord_tolerance`, le point projete de ce milieu est insere
(au plus `max_rounds` tours). Les points ne sont ajoutes que la ou la corde s'ecarte (courbure, arete vive), pas
partout : le nombre de points de trajectoire croit peu.

Opt-in : `PipelineParams.skin_projection = None` (defaut) ne fait rien.
Requete de distance : noyau natif (`apex_core.MeshDistance`) quand il est actif, sinon trimesh ; meme resultat
(point le plus proche exact), verifie par `tests/test_step5c_skin_projection.py`.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np
import trimesh

from core import _cpp
from core.mesh_types import Layer, MeshData


@dataclass
class SkinProjectionParams:
    chord_tolerance: float = 0.02     # mm : un milieu de segment a plus de cela de la surface recoit un point de plus
    max_rounds: int = 3               # tours de subdivision adaptative
    max_shift: float = 1.5            # mm : un point plus loin de la surface n'est PAS deplace (securite)
    crease_angle_deg: float = 8.0     # deux faces dont les normales different de plus : le point d'arete est insere
    crease_tolerance: float = 1e-4    # mm : le point d'arete doit etre a moins de cela de la surface
    batch: int = 200_000              # points par requete (memoire)
    # Deplacement DANS LE PLAN de la couche (2026-11-19) : le point va jusqu'au plan tangent de la surface au point le
    # plus proche, en restant a sa hauteur de couche. Le point le plus proche lui-meme sortait le contour de sa couche
    # sur toute paroi en pente (couches planes d'un Benchy : +-0,2 mm sur le chanfrein du fond, premiere couche de
    # z = 0 a 0,3 mm). Paroi verticale : meme resultat. Surface presque parallele a la couche : point laisse en place.
    in_layer: bool = True


@dataclass
class SkinProjectionReport:
    layers: int = 0
    points_in: int = 0
    points_out: int = 0
    moved: int = 0                    # sommets effectivement deplaces
    refused: int = 0                  # sommets restes en place car a plus de `max_shift` de la surface
    max_shift_applied: float = 0.0
    rounds_used: int = 0


class _Surface:
    """Point le plus proche sur la surface d'origine : natif si actif, trimesh sinon. Peut aussi rendre le triangle."""

    def __init__(self, vertices: np.ndarray, faces: np.ndarray, face_normals: np.ndarray, batch: int):
        self.batch = batch
        self.vertices = np.asarray(vertices, dtype=np.float64)
        self.faces = np.asarray(faces, dtype=np.int64)
        self.face_normals = np.asarray(face_normals, dtype=np.float64)
        self._native = None
        self._native_faces = False
        self._trimesh = None
        if _cpp.enabled() and hasattr(_cpp.module(), "MeshDistance"):
            self._native = _cpp.module().MeshDistance(np.ascontiguousarray(self.vertices),
                                                      np.ascontiguousarray(self.faces))
            try:     # module construit avant l'ajout de `return_faces` : repli trimesh pour les triangles seulement
                self._native.query(np.zeros((1, 3)), 1, 0, True)
                self._native_faces = True
            except TypeError:
                self._native_faces = False

    def _trimesh_query(self):
        if self._trimesh is None:
            self._trimesh = trimesh.proximity.ProximityQuery(
                trimesh.Trimesh(self.vertices, self.faces, process=False))
        return self._trimesh

    def closest(self, points: np.ndarray, want_faces: bool = False
                ) -> tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
        """(points les plus proches, distances, triangles ou None), par lots."""
        out_points, out_dist, out_faces = [], [], []
        use_native = self._native is not None and (self._native_faces or not want_faces)
        for start in range(0, len(points), self.batch):
            chunk = np.ascontiguousarray(points[start:start + self.batch], dtype=np.float64)
            if use_native:
                result = self._native.query(chunk, 1, 0, True) if want_faces else self._native.query(chunk)
                dist, closest = result[0], result[1]
                faces = result[2] if want_faces else None
            else:
                closest, dist, faces = self._trimesh_query().on_surface(chunk)
            out_points.append(np.asarray(closest, dtype=np.float64))
            out_dist.append(np.asarray(dist, dtype=np.float64))
            if want_faces:
                out_faces.append(np.asarray(faces, dtype=np.int64))
        if not out_points:
            return np.zeros((0, 3)), np.zeros(0), (np.zeros(0, dtype=np.int64) if want_faces else None)
        return (np.concatenate(out_points), np.concatenate(out_dist),
                np.concatenate(out_faces) if want_faces else None)

    def plane_offsets(self, face_ids: np.ndarray) -> np.ndarray:
        """Decalage `n . x = h` du plan de chaque triangle."""
        first = self.vertices[self.faces[face_ids, 0]]
        return np.einsum("ij,ij->i", self.face_normals[face_ids], first)


def _crease_points(surface: _Surface, a: np.ndarray, b: np.ndarray, mid: np.ndarray,
                   params: SkinProjectionParams, layer_normals: Optional[np.ndarray] = None
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Pour les segments [a, b] dont les extremites sont sur deux faces d'orientations differentes : point de la
    droite d'intersection des deux plans le plus proche du milieu. `ok` : ce point est bien sur la surface (arete
    reelle, pas le prolongement d'un plan)."""
    count = len(mid)
    points = np.zeros((count, 3))
    ok = np.zeros(count, dtype=bool)
    if count == 0:
        return points, ok
    fa = surface.closest(a, want_faces=True)[2]
    fb = surface.closest(b, want_faces=True)[2]
    na, nb = surface.face_normals[fa], surface.face_normals[fb]
    direction = np.cross(na, nb)
    sq = np.einsum("ij,ij->i", direction, direction)
    creased = sq > np.sin(np.radians(params.crease_angle_deg)) ** 2
    if not creased.any():
        return points, ok
    ha, hb = surface.plane_offsets(fa), surface.plane_offsets(fb)
    safe = np.where(creased, sq, 1.0)[:, None]
    base = (ha[:, None] * np.cross(nb, direction) + hb[:, None] * np.cross(direction, na)) / safe
    unit = direction / np.sqrt(safe)
    candidate = base + unit * np.einsum("ij,ij->i", mid - base, unit)[:, None]
    if layer_normals is not None:
        # point de l'arete dans le plan de la couche (par le milieu), si l'arete n'est pas parallele a la couche
        along = np.einsum("ij,ij->i", unit, layer_normals)
        cut = np.abs(along) > 0.2
        s_cut = np.einsum("ij,ij->i", mid - base, layer_normals) / np.where(cut, along, 1.0)
        candidate = np.where(cut[:, None], base + unit * s_cut[:, None], candidate)
    _c, dist, _f = surface.closest(candidate)
    near = np.linalg.norm(candidate - mid, axis=1) <= params.max_shift
    not_endpoint = (np.linalg.norm(candidate - a, axis=1) > 1e-6) & (np.linalg.norm(candidate - b, axis=1) > 1e-6)
    ok = creased & (dist <= params.crease_tolerance) & near & not_endpoint
    if layer_normals is not None:
        # une arete parallele a la couche (bord du fond, d'un dessus) ne la coupe pas : son point sortirait de la couche
        ok &= np.abs(np.einsum("ij,ij->i", candidate - mid, layer_normals)) <= params.crease_tolerance
    return candidate, ok


def _in_layer(points: np.ndarray, closest: np.ndarray, normals: np.ndarray, layer_normals: np.ndarray
              ) -> tuple[np.ndarray, np.ndarray]:
    """Points deplaces dans le plan de leur couche (normale `layer_normals`) jusqu'au plan tangent de la surface
    (normale `normals`) au point le plus proche, et longueur du deplacement. Surface a moins de ~11 deg de la couche
    (dessus, dessous) : point laisse en place, deplacement nul."""
    n = normals / np.maximum(np.linalg.norm(normals, axis=1), 1e-12)[:, None]
    g = layer_normals
    m = n - np.einsum("ij,ij->i", n, g)[:, None] * g
    denom = np.einsum("ij,ij->i", m, n)
    ok = denom > 0.04
    t = np.where(ok, np.einsum("ij,ij->i", closest - points, n) / np.where(ok, denom, 1.0), 0.0)
    moved = points + t[:, None] * m
    return moved, np.abs(t) * np.linalg.norm(m, axis=1)


def _layer_normals(layers: list[Layer], owners, starts: np.ndarray, count: int) -> np.ndarray:
    """Normale du repere local de la couche de chaque point (plan moyen de la couche)."""
    per_loop = np.array([np.cross(np.asarray(layers[o[0]].local_frame[1], float), np.asarray(layers[o[0]].local_frame[2], float))
                         for o in owners]).reshape(-1, 3)
    per_loop /= np.maximum(np.linalg.norm(per_loop, axis=1), 1e-12)[:, None]
    owner = np.searchsorted(starts, np.arange(count), side="right") - 1
    return per_loop[owner]


def _flatten(layers: list[Layer]) -> tuple[np.ndarray, np.ndarray, list[tuple[int, int, bool]]]:
    """Tous les sommets de contour dans un tableau, les polylignes dans un tableau de debuts.
    Les boucles fermees (premier point repete a la fin) perdent leur doublon, restitue a la fin."""
    chunks, owners = [], []
    for layer_index, layer in enumerate(layers):
        for loop_index, loop in enumerate(layer.contour_3d):
            loop = np.asarray(loop, dtype=np.float64)
            closed = len(loop) > 2 and np.array_equal(loop[0], loop[-1])
            chunks.append(loop[:-1] if closed else loop)
            owners.append((layer_index, loop_index, closed))
    if not chunks:
        return np.zeros((0, 3)), np.zeros(0, dtype=np.int64), owners
    starts = np.cumsum([0] + [len(c) for c in chunks[:-1]]).astype(np.int64)
    return np.concatenate(chunks), starts, owners


def _project_vertices(points: np.ndarray, surface: _Surface, params: SkinProjectionParams,
                      report: SkinProjectionReport, layer_normals: Optional[np.ndarray] = None) -> np.ndarray:
    closest, dist, faces = surface.closest(points, want_faces=layer_normals is not None)
    if layer_normals is not None:
        closest, dist = _in_layer(points, closest, surface.face_normals[faces], layer_normals)
    allowed = dist <= params.max_shift
    report.refused = int((~allowed).sum())
    report.moved = int((allowed & (dist > 1e-9)).sum())
    if allowed.any():
        report.max_shift_applied = float(dist[allowed].max())
    return np.where(allowed[:, None], closest, points)


def _refine_once(points: np.ndarray, starts: np.ndarray, closed: np.ndarray, surface: _Surface,
                 params: SkinProjectionParams, normals_of=None) -> tuple[np.ndarray, np.ndarray, int]:
    """Un tour : insere le point projete du milieu des segments qui s'ecartent de la surface."""
    n = len(points)
    ends = np.append(starts[1:], n)                       # fin (exclue) de chaque polyligne
    index = np.arange(n)
    owner = np.searchsorted(starts, index, side="right") - 1
    nxt = index + 1
    is_last = nxt >= ends[owner]
    nxt = np.where(is_last, starts[owner], nxt)           # le dernier point d'une boucle fermee rejoint le premier
    has_segment = ~is_last | closed[owner]
    seg = index[has_segment]
    mid = 0.5 * (points[seg] + points[nxt[seg]])
    closest, dist, faces = surface.closest(mid, want_faces=normals_of is not None)
    g = normals_of(seg) if normals_of is not None else None
    gap = dist
    if g is not None:
        closest, moved = _in_layer(mid, closest, surface.face_normals[faces], g)
        gap = np.where(moved > 0.0, moved, 0.0)
    insert = (gap > params.chord_tolerance) & (gap <= params.max_shift)
    if not insert.any():
        return points, starts, 0
    where = seg[insert]
    inserted = closest[insert]
    # un segment qui relie deux faces d'orientations differentes coupe une arete : on y place le point d'arete
    crease, ok = _crease_points(surface, points[where], points[nxt[where]], mid[insert], params,
                                g[insert] if g is not None else None)
    inserted = np.where(ok[:, None], crease, inserted)
    new_points = np.insert(points, where + 1, inserted, axis=0)
    added_before = np.zeros(n + 1, dtype=np.int64)
    np.add.at(added_before, where + 1, 1)
    shift = np.cumsum(added_before)[:-1]                  # points inseres avant chaque indice d'origine
    return new_points, starts + shift[starts], int(insert.sum())


def _rebuild(layers: list[Layer], points: np.ndarray, starts: np.ndarray, owners) -> list[Layer]:
    ends = np.append(starts[1:], len(points))
    contours: dict[int, list[Optional[np.ndarray]]] = {}
    for k, (layer_index, loop_index, closed) in enumerate(owners):
        loop = points[starts[k]:ends[k]]
        contours.setdefault(layer_index, [None] * len(layers[layer_index].contour_3d))[loop_index] = (
            np.vstack([loop, loop[:1]]) if closed else loop.copy())
    return [replace(layer, contour_3d=contours.get(i, layer.contour_3d)) for i, layer in enumerate(layers)]


def run(layers: list[Layer], mesh_data: MeshData, params: SkinProjectionParams
        ) -> tuple[list[Layer], SkinProjectionReport]:
    """Retourne les couches aux contours ramenes sur la surface d'origine (`mesh_data`) et un rapport.
    Fonction pure : les couches d'entree ne sont pas modifiees."""
    report = SkinProjectionReport(layers=len(layers))
    points, starts, owners = _flatten(layers)
    report.points_in = len(points)
    if len(points) == 0:
        return layers, report
    closed = np.array([o[2] for o in owners], dtype=bool)
    surface = _Surface(mesh_data.vertices, mesh_data.faces, mesh_data.face_normals, params.batch)
    normals = _layer_normals(layers, owners, starts, len(points)) if params.in_layer else None
    points = _project_vertices(points, surface, params, report, normals)
    for round_index in range(params.max_rounds):
        normals_of = None
        if params.in_layer:
            point_normals = _layer_normals(layers, owners, starts, len(points))
            normals_of = point_normals.__getitem__
        points, starts, inserted = _refine_once(points, starts, closed, surface, params, normals_of)
        if inserted == 0:
            break
        report.rounds_used = round_index + 1
    report.points_out = len(points)
    return _rebuild(layers, points, starts, owners), report
