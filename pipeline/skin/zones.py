"""
P2 sur la voie chaleur : quelles faces de la surface d'origine peuvent recevoir une couche de lissage ?

Meme principe que `pipeline/bird/skin/zones.py` (branche `bird_old`) : la matiere imprimee au moment t est
{phi <= phi(t)} (croissante, car les couches sont les isosurfaces de phi, ordonnees par isovaleur) ; une face est
retenue si la vraie tete, posee sur elle la buse vers -n, n'entre ni dans le plateau ni dans la piece avant que la
face soit portee. Ici phi est celui de la methode de la chaleur (`ScalarField.values`) et les isovaleurs sont celles
des couches de l'etape 5.

Statuts : 0 retenue, 1 inatteignable (tete a plus de `head_tolerance_deg` de la normale), 2 collision plateau,
3 collision objet, 4 contact plateau.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import trimesh
from scipy.spatial import cKDTree

from core.mesh_types import MachineProfile, ScalarField, TetMesh
from pipeline.step8_machine_toolpath import Step8Params, _RawPoint, _resolve_all_axis_angles_coupled
from pipeline.step9_collision_sim import CollisionSimParams, _filter_nozzle_contact_zone, _surface_sample_points

NOZZLE_WIDTH = 0.4
STATUS_RETAINED, STATUS_UNREACHABLE, STATUS_PLATE, STATUS_OBJECT, STATUS_PLATE_CONTACT = 0, 1, 2, 3, 4
STATUS_NAMES = {0: "retenue", 1: "inatteignable", 2: "collision plateau", 3: "collision objet", 4: "contact plateau"}
PLATE_CONTACT_TOL = 1e-3
FACE_BATCH = 64


@dataclass
class SkinZoneParams:
    max_edge: float = 5.0 * NOZZLE_WIDTH     # subdivision de la surface analysee (mm)
    head_tolerance_deg: float = 45.0         # ecart maximal tete / normale (decision provisoire, comme Bird)
    head_sample_spacing: float = 5.0         # comme l'etape 9
    face_batch: int = FACE_BATCH


@dataclass
class SkinZones:
    vertices: np.ndarray        # (V,3) surface re-decoupee
    faces: np.ndarray           # (F,3)
    normals: np.ndarray         # (F,3) vers l'exterieur
    areas: np.ndarray           # (F,)
    status: np.ndarray          # (F,) int8
    moment: np.ndarray          # (F,) rang de la couche porteuse
    deviation_deg: np.ndarray   # (F,) ecart tete / normale (pire de centre + sommets) ; 0 si non evalue
    evaluated: np.ndarray       # (F,) bool : faces dont le statut a ete calcule

    @property
    def retained(self) -> np.ndarray:
        return np.flatnonzero(self.status == STATUS_RETAINED)


def boundary_faces(tetrahedra: np.ndarray) -> np.ndarray:
    """Faces triangulaires n'appartenant qu'a un tetraedre."""
    quads = np.concatenate([tetrahedra[:, [0, 1, 2]], tetrahedra[:, [0, 1, 3]], tetrahedra[:, [0, 2, 3]],
                            tetrahedra[:, [1, 2, 3]]])
    key = np.sort(quads, axis=1)
    _u, first, counts = np.unique(key, axis=0, return_index=True, return_counts=True)
    return quads[first[counts == 1]]


def barycentric(points: np.ndarray, tets: np.ndarray) -> np.ndarray:
    a = tets[:, 0]
    m = np.stack([tets[:, 1] - a, tets[:, 2] - a, tets[:, 3] - a], axis=2)
    with np.errstate(all="ignore"):
        w = np.linalg.solve(m, (points - a)[..., None])[..., 0]
    return np.concatenate([1 - w.sum(axis=1, keepdims=True), w], axis=1)


class PrintedSolid:
    """La piece avec, en chaque point, sa valeur phi : la matiere imprimee a l'instant t est {phi <= phi(t)}."""

    def __init__(self, tet_vertices, tetrahedra, phi_nodes, k: int = 16):
        self.tets = tet_vertices[tetrahedra]
        self.phi = phi_nodes[tetrahedra]
        self.tree = cKDTree(self.tets.mean(axis=1))
        self.k = k
        self.reach = float(np.linalg.norm(self.tets - self.tets.mean(axis=1, keepdims=True), axis=2).max())
        self.lo, self.hi = tet_vertices.min(axis=0), tet_vertices.max(axis=0)
        self.tolerance = -1e-9

    def phi_at(self, points: np.ndarray) -> np.ndarray:
        """phi aux points (NaN hors de la piece)."""
        out = np.full(len(points), np.nan)
        idx_box = np.flatnonzero(np.all((points >= self.lo) & (points <= self.hi), axis=1))
        if not len(idx_box):
            return out
        dist, cand = self.tree.query(points[idx_box], k=self.k)
        pending = dist[:, 0] <= self.reach
        for j in range(self.k):
            idx = np.flatnonzero(pending)
            if not len(idx):
                break
            t = cand[idx, j]
            bary = barycentric(points[idx_box[idx]], self.tets[t])
            inside = bary.min(axis=1) >= self.tolerance
            hit = idx[inside]
            out[idx_box[hit]] = np.einsum("nv,nv->n", bary[inside], self.phi[t[inside]])
            pending[hit] = False
        return out


def rotation_matrices(angles_deg: np.ndarray, axes) -> np.ndarray:
    """(F,3,3), meme composition que `step9._positioned_toolhead_vertices` (axe le plus interne d'abord)."""
    total = np.tile(np.eye(3), (len(angles_deg), 1, 1))
    for i, spec in reversed(list(enumerate(axes))):
        k = np.asarray(spec.axis, float)
        k = k / np.linalg.norm(k)
        kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        a = np.radians(angles_deg[:, i])[:, None, None]
        total = (np.eye(3) + np.sin(a) * kx + (1 - np.cos(a)) * (kx @ kx)) @ total
    return total


def surface_mesh(tet_mesh: TetMesh, solid: PrintedSolid, max_edge: float):
    """Surface du volume re-decoupee a `max_edge`, normales sortantes : (vertices, faces, normals, areas, centers)."""
    verts = tet_mesh.tet_vertices
    surf_v, faces = trimesh.remesh.subdivide_to_size(verts, boundary_faces(tet_mesh.tetrahedra), max_edge)
    tri = surf_v[faces]
    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    normals = cross / np.maximum(2 * areas, 1e-12)[:, None]
    centers = tri.mean(axis=1)
    # le maillage de frontiere des tetraedres n'a pas de sens garanti : une sonde juste au-dessus du centre
    # doit tomber hors de la piece
    inside = ~np.isnan(solid.phi_at(centers + 1e-3 * normals))
    normals[inside] *= -1
    return surf_v, faces, normals, areas, centers


def compute(
    tet_mesh: TetMesh,
    scalar_field: ScalarField,
    isovalues: np.ndarray,
    machine: MachineProfile,
    sim: CollisionSimParams,
    params: Optional[SkinZoneParams] = None,
    face_sample: Optional[int] = None,
    seed: int = 0,
) -> SkinZones:
    """Statut de chaque face de la surface. `face_sample` : n'evaluer qu'un echantillon (graine fixe) de faces
    hors plateau ; les autres gardent `evaluated=False` et le statut 1 (jamais retenues)."""
    params = params or SkinZoneParams()
    verts = tet_mesh.tet_vertices
    isovalues = np.sort(np.unique(np.asarray(isovalues, dtype=float)))
    solid = PrintedSolid(verts, tet_mesh.tetrahedra, scalar_field.values)
    surf_v, faces, normals, areas, centers = surface_mesh(tet_mesh, solid, params.max_edge)

    solid.tolerance = -1e-6                  # sommets de surface : sur la frontiere
    surf_phi = solid.phi_at(surf_v)
    solid.tolerance = -1e-9
    lost = np.isnan(surf_phi)
    if lost.any():
        surf_phi[lost] = scalar_field.values[cKDTree(verts).query(surf_v[lost])[1]]
    tri = surf_v[faces]
    status = np.full(len(faces), STATUS_UNREACHABLE, dtype=np.int8)
    free = tri[:, :, 2].max(axis=1) >= verts[:, 2].min() + PLATE_CONTACT_TOL
    status[~free] = STATUS_PLATE_CONTACT
    moment = np.clip(np.searchsorted(isovalues, surf_phi[faces].max(axis=1) - 1e-9, side="left"), 0,
                     len(isovalues) - 1)

    candidates = np.flatnonzero(free)
    if face_sample is not None and face_sample < len(candidates):
        candidates = np.sort(np.random.default_rng(seed).choice(candidates, size=face_sample, replace=False))
    evaluated = np.zeros(len(faces), bool)
    evaluated[candidates] = True
    status[candidates] = STATUS_RETAINED

    axes = machine.rotation_axes
    step8 = Step8Params()
    ref = np.array([0.0, 0.0, 1.0])
    angles = np.zeros((len(faces), len(axes)))
    deviation = np.zeros(len(faces))
    raw = [_RawPoint(position=centers[i], desired_direction=normals[i], layer_id=0, segment_index=0, color_id=None,
                     feedrate=1.0, extrusion=0.0, is_travel=False) for i in candidates]
    a, err = _resolve_all_axis_angles_coupled(raw, axes, ref, machine, step8)
    angles[candidates], deviation[candidates] = np.asarray(a), np.asarray(err)
    corners = tri[candidates].reshape(-1, 3)
    raw_c = [_RawPoint(position=q, desired_direction=normals[candidates[j // 3]], layer_id=0, segment_index=0,
                       color_id=None, feedrate=1.0, extrusion=0.0, is_travel=False) for j, q in enumerate(corners)]
    _a, err_c = _resolve_all_axis_angles_coupled(raw_c, axes, ref, machine, step8)
    deviation[candidates] = np.maximum(deviation[candidates], np.asarray(err_c).reshape(-1, 3).max(axis=1))
    deviation = np.degrees(deviation)
    status[candidates[deviation[candidates] > params.head_tolerance_deg]] = STATUS_UNREACHABLE

    head = _surface_sample_points(machine.nozzle_model, params.head_sample_spacing)
    head = _filter_nozzle_contact_zone(head, sim.nozzle_contact_clearance)
    _c, _d, tid = trimesh.proximity.closest_point(machine.nozzle_model, head)
    hn = machine.nozzle_model.face_normals[tid]
    head = np.vstack([head, head + sim.safety_margin * hn, head - sim.safety_margin * hn])   # marge de securite
    reach = candidates[status[candidates] == STATUS_RETAINED]
    plate = np.zeros(len(reach), bool)
    last = np.full(len(reach), len(isovalues) - 1)
    for s in range(0, len(reach), params.face_batch):
        ids = reach[s:s + params.face_batch]
        rot = rotation_matrices(angles[ids], axes)
        pts = np.einsum("fij,pj->fpi", rot, head) + centers[ids][:, None, :]
        plate[s:s + len(ids)] = (pts[:, :, 2] < sim.safety_margin).any(axis=1)
        phi = solid.phi_at(pts.reshape(-1, 3)).reshape(len(ids), -1)
        low = np.min(np.where(np.isnan(phi), np.inf, phi), axis=1)
        ins = np.isfinite(low)
        last[s:s + len(ids)][ins] = np.searchsorted(isovalues, low[ins] - 1e-9, side="left") - 1
    status[reach[plate]] = STATUS_PLATE
    blocked = np.zeros(len(faces), bool)
    blocked[reach] = last < moment[reach]
    status[(status == STATUS_RETAINED) & blocked] = STATUS_OBJECT
    return SkinZones(vertices=surf_v, faces=faces, normals=normals, areas=areas, status=status, moment=moment,
                     deviation_deg=deviation, evaluated=evaluated)
