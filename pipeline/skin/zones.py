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
    angles: Optional[np.ndarray] = None       # (F, axes) angles de tete (buse vers -n)
    last_access: Optional[np.ndarray] = None  # (F,) dernier moment accessible ; -1 = bloquee des le debut
    phi: Optional[np.ndarray] = None          # (V,) phi aux sommets de la surface
    isovalues: Optional[np.ndarray] = None

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


def surface_mesh(tet_mesh: TetMesh, solid: PrintedSolid, max_edge: float, surface=None):
    """Surface du volume re-decoupee a `max_edge`, normales sortantes : (vertices, faces, normals, areas, centers).
    `surface` : (vertices, faces) de la surface d'origine (`MeshData`) ; sinon la frontiere des tetraedres, bien plus
    fine pres du plateau (raffinement local) donc beaucoup plus de faces a evaluer."""
    if surface is not None:
        surf_v, faces = trimesh.remesh.subdivide_to_size(np.asarray(surface[0], dtype=float),
                                                        np.asarray(surface[1], dtype=np.int64), max_edge)
    else:
        surf_v, faces = trimesh.remesh.subdivide_to_size(tet_mesh.tet_vertices, boundary_faces(tet_mesh.tetrahedra),
                                                        max_edge)
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


_ACCESS: dict = {}


def _access_chunk(ids: np.ndarray):
    a = _ACCESS
    plate = np.zeros(len(ids), bool)
    last = np.full(len(ids), len(a["isovalues"]) - 1)
    for s in range(0, len(ids), a["batch"]):
        part = ids[s:s + a["batch"]]
        rot = rotation_matrices(a["angles"][part], a["axes"])
        pts = np.einsum("fij,pj->fpi", rot, a["head"]) + a["centers"][part][:, None, :]
        plate[s:s + len(part)] = (pts[:, :, 2] < a["margin"]).any(axis=1)
        phi = a["solid"].phi_at(pts.reshape(-1, 3)).reshape(len(part), -1)
        low = np.min(np.where(np.isnan(phi), np.inf, phi), axis=1)
        ins = np.isfinite(low)
        last[s:s + len(part)][ins] = np.searchsorted(a["isovalues"], low[ins] - 1e-9, side="left") - 1
    return plate, last


def _last_accessible(head, centers, angles, axes, solid, isovalues, face_ids, margin, batch, workers):
    """(plateau, dernier moment accessible) par face de `face_ids` : la tete posee sur la face, buse vers -n, ne doit
    ni entrer dans le plateau ni, tant que la face n'est pas portee, dans la piece. Faces independantes : reparties
    sur des processus (fork ; la memoire est partagee, rien n'est copie)."""
    import os
    _ACCESS.update(head=head, centers=centers, angles=angles, axes=axes, solid=solid, isovalues=isovalues,
                   margin=margin, batch=batch)
    workers = workers if workers is not None else max(1, min(8, (os.cpu_count() or 2)))
    if workers <= 1 or len(face_ids) < 4 * batch or not hasattr(os, "fork"):
        return _access_chunk(face_ids)
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing
    chunks = np.array_split(face_ids, workers * 4)
    with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("fork")) as pool:
        parts = list(pool.map(_access_chunk, chunks))
    return np.concatenate([p for p, _l in parts]), np.concatenate([l for _p, l in parts])


def compute(
    tet_mesh: TetMesh,
    scalar_field: ScalarField,
    isovalues: np.ndarray,
    machine: MachineProfile,
    sim: CollisionSimParams,
    params: Optional[SkinZoneParams] = None,
    face_sample: Optional[int] = None,
    seed: int = 0,
    surface=None,
    workers: Optional[int] = None,
) -> SkinZones:
    """Statut de chaque face de la surface. `face_sample` : n'evaluer qu'un echantillon (graine fixe) de faces
    hors plateau ; les autres gardent `evaluated=False` et le statut 1 (jamais retenues)."""
    params = params or SkinZoneParams()
    verts = tet_mesh.tet_vertices
    isovalues = np.sort(np.unique(np.asarray(isovalues, dtype=float)))
    solid = PrintedSolid(verts, tet_mesh.tetrahedra, scalar_field.values)
    surf_v, faces, normals, areas, centers = surface_mesh(tet_mesh, solid, params.max_edge, surface)

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
    last_all = np.full(len(faces), -1)
    plate, last = _last_accessible(head, centers, angles, axes, solid, isovalues, reach, sim.safety_margin,
                                   params.face_batch, workers)
    last_all[reach] = last
    status[reach[plate]] = STATUS_PLATE
    blocked = np.zeros(len(faces), bool)
    blocked[reach] = last < moment[reach]
    status[(status == STATUS_RETAINED) & blocked] = STATUS_OBJECT
    return SkinZones(vertices=surf_v, faces=faces, normals=normals, areas=areas, status=status, moment=moment,
                     deviation_deg=deviation, evaluated=evaluated, angles=angles, last_access=last_all,
                     phi=surf_phi, isovalues=isovalues)


# --- partition en morceaux, bandes de moments et zones (fin de P2 de la voie Bird, inchangee) -------------------
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402
from pipeline.skin import surface_paths as sp  # noqa: E402

EDGE_FACE_ANGLE_DEG = 150.0               # arete vive : angle ENTRE LES FACES en dessous
MIN_ZONE_AREA = 4.0 * NOZZLE_WIDTH ** 2   # plus petit : pas de couche de surface
MICRO_INRADIUS = 1.5 * NOZZLE_WIDTH       # micro-zone : plus grande distance au bord en dessous
STATUS_SMALL, STATUS_MICRO = 5, 6
STATUS_NAMES[5] = "zone trop petite"
STATUS_NAMES[6] = "micro-zone non lissee"


def layer_rank(values, isovalues):
    """Rang de la couche qui porte la hauteur deformee `values`."""
    return np.clip(np.searchsorted(isovalues, values - 1e-9, side="left"), 0, len(isovalues) - 1)


def face_adjacency(faces):
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(faces)), 3)
    order = np.lexsort((edges[:, 1], edges[:, 0]))
    e, o = edges[order], owner[order]
    same = np.all(e[1:] == e[:-1], axis=1)
    a, b = o[:-1][same], o[1:][same]
    n = len(faces)
    return coo_matrix((np.ones(len(a) * 2), (np.r_[a, b], np.r_[b, a])), shape=(n, n)).tocsr()


def minimal_stabbing(lo, hi):
    """Plus petit ensemble de moments touchant tous les intervalles [lo, hi]
    (glouton exact : extremite droite croissante)."""
    points = []
    for i in np.argsort(hi, kind="stable"):
        if not points or points[-1] < lo[i]:
            points.append(int(hi[i]))
    return np.array(points)


def split_along_level(vertices, faces, phi, face_mask, level):
    """Coupe les faces de `face_mask` traversees par phi = level (phi lineaire
    sur chaque face : coupe exacte). Sommets crees partages entre faces
    voisines (cle = arete). Rend (vertices, faces, phi, parent) ; parent =
    index de la face coupee d'origine pour chaque nouvelle face."""
    s = phi[faces] - level
    pos = s >= 0
    cross = face_mask & (pos.any(axis=1) & ~pos.all(axis=1))
    parent = np.arange(len(faces))
    if not cross.any():
        return vertices, faces, phi, parent
    ci = np.flatnonzero(cross)
    f, sp = faces[ci], pos[ci]
    # sommet isole en premier (celui dont le signe differe des deux autres), orientation conservee
    alone = np.where(sp.sum(axis=1) == 1, np.argmax(sp, axis=1), np.argmin(sp, axis=1))
    roll = np.stack([alone, (alone + 1) % 3, (alone + 2) % 3], axis=1)
    f = np.take_along_axis(f, roll, axis=1)
    a, b, c = f[:, 0], f[:, 1], f[:, 2]
    new_v, new_phi, key_index = [], [], {}

    def edge_point(p, q):
        out = np.empty(len(p), int)
        for j, (u, w) in enumerate(zip(p, q)):
            key = (min(u, w), max(u, w))
            if key not in key_index:
                su, sw = phi[u] - level, phi[w] - level
                t = su / (su - sw)
                if t <= 1e-9:
                    key_index[key] = u
                elif t >= 1 - 1e-9:
                    key_index[key] = w
                else:
                    key_index[key] = len(vertices) + len(new_v)
                    new_v.append(vertices[u] + t * (vertices[w] - vertices[u]))
                    new_phi.append(level)
            out[j] = key_index[key]
        return out

    pab, pac = edge_point(a, b), edge_point(a, c)
    tri = np.concatenate([np.stack([a, pab, pac], 1), np.stack([pab, b, c], 1), np.stack([pab, c, pac], 1)])
    par = np.tile(ci, 3)
    good = (tri[:, 0] != tri[:, 1]) & (tri[:, 1] != tri[:, 2]) & (tri[:, 0] != tri[:, 2])
    keep = np.flatnonzero(~cross)
    vertices = np.vstack([vertices, np.array(new_v).reshape(-1, 3)])
    phi = np.concatenate([phi, np.array(new_phi)])
    return vertices, np.vstack([faces[keep], tri[good]]), phi, np.concatenate([keep, par[good]])




def partition(z: SkinZones, edge_face_angle: float = EDGE_FACE_ANGLE_DEG, facet_smoothing=None,
              micro_inradius: float = MICRO_INRADIUS, trim_points=None, trim_radius: float = 2.0) -> dict:
    """Morceaux (aretes vives), bandes de moments en nombre minimal, frontieres lisses le long des isovaleurs, zones,
    micro-zones : renvoie le dictionnaire que la voie Bird ecrivait dans `_skin_pass<k>.npz` (cles `vertices`,
    `faces`, `status`, `moment`, `zone`, `zone_moment`, `normals`, `areas`, `isovalues`, `nozzle_width`, ...).
    `trim_points` : points de collision cumules de P5 (rognage local), faces a moins de `trim_radius` retirees."""
    isovalues = z.isovalues
    surf_vertices, faces, surf_phi = z.vertices, z.faces, z.phi
    normals, angles, status, moment, last = z.normals, z.angles, z.status.copy(), z.moment, z.last_access.copy()
    n_faces = len(faces)
    tri = surf_vertices[faces]
    # faces non evaluees (echantillon) : jamais retenues
    status[~z.evaluated & (status == STATUS_RETAINED)] = STATUS_UNREACHABLE
    areas = z.areas
    centers = tri.mean(axis=1)
    if trim_points is not None:
        trim_pts = np.asarray(trim_points, dtype=float).reshape(-1, 3)
        if len(trim_pts):
            tree = cKDTree(trim_pts)
            near = np.minimum(tree.query(tri.reshape(-1, 3))[0].reshape(-1, 3).min(axis=1), tree.query(centers)[0])
            status[(status == 0) & (near <= trim_radius)] = STATUS_OBJECT
    return _partition_core(surf_vertices, faces, surf_phi, normals, angles, status, moment, last, areas, isovalues,
                           edge_face_angle, facet_smoothing, micro_inradius)


def _partition_core(surf_vertices, faces, surf_phi, normals, angles, status, moment, last, areas, isovalues,
                    edge_face_angle, facet_smoothing, micro_inradius) -> dict:
    n_faces = len(faces)
    kept = np.flatnonzero(status == 0)
    adj = face_adjacency(faces[kept]).tocoo()
    cos_n = np.einsum("ij,ij->i", normals[kept][adj.row], normals[kept][adj.col])
    deviation = np.degrees(np.arccos(np.clip(cos_n, -1, 1)))      # angle entre normales
    cut_above = 180.0 - edge_face_angle
    if facet_smoothing is not None:
        cut_above = max(cut_above, facet_smoothing)
    link = deviation <= cut_above
    graph = coo_matrix((np.ones(link.sum()), (adj.row[link], adj.col[link])), shape=(len(kept),) * 2)
    _n, piece_of_kept = connected_components(graph, directed=False)
    piece = np.full(n_faces, -1)
    piece[kept] = piece_of_kept
    n_edge_pieces = piece_of_kept.max() + 1 if len(kept) else 0

    band_moments = {}                       # morceau -> moments choisis (croissants)
    levels_needed = {}                      # isovaleur -> faces (d'origine) a couper
    for p in range(n_edge_pieces):
        ids = np.flatnonzero(piece == p)
        points = minimal_stabbing(moment[ids], last[ids])
        band_moments[p] = points
        for m in points[:-1]:
            levels_needed.setdefault(float(isovalues[m]), []).append(ids)

    # coupe exacte le long des lignes de niveau (frontieres lisses)
    v, f, phi = surf_vertices, faces, surf_phi.copy()
    origin = np.arange(n_faces)             # face d'origine de chaque face courante
    for level in sorted(levels_needed):
        mask_orig = np.zeros(n_faces, bool)
        mask_orig[np.concatenate(levels_needed[level])] = True
        v, f, phi, par = split_along_level(v, f, phi, mask_orig[origin], level)
        origin = origin[par]

    # attributs herites de la face d'origine ; moment porteur recalcule
    status2, piece2, last2 = status[origin], piece[origin], last[origin]
    moment2 = layer_rank(phi[f].max(axis=1), isovalues)
    band = np.full(len(f), -1)
    for p in range(n_edge_pieces):
        sel = np.flatnonzero(piece2 == p)
        pts = band_moments[p]
        # plus petit moment choisi >= moment porteur ET <= dernier moment accessible
        j = np.searchsorted(pts, moment2[sel], side="left")
        ok = j < len(pts)
        ok[ok] = pts[j[ok]] <= last2[sel][ok]
        # repli (rare) : premier moment choisi valide pour la face
        for q in np.flatnonzero(~ok):
            valid = np.flatnonzero((pts >= moment2[sel[q]]) & (pts <= last2[sel[q]]))
            j[q] = valid[0] if len(valid) else -1
        band[sel] = np.where(j >= 0, j, -1)
    # zones = composantes connexes de (morceau, bande)
    tri2 = v[f]
    areas2 = 0.5 * np.linalg.norm(np.cross(tri2[:, 1] - tri2[:, 0], tri2[:, 2] - tri2[:, 0]), axis=1)
    zone = np.full(len(f), -1)
    zone_moment = []
    kept2 = np.flatnonzero((status2 == 0) & (band >= 0))
    adj2 = face_adjacency(f[kept2]).tocoo()
    same = (piece2[kept2][adj2.row] == piece2[kept2][adj2.col]) & (band[kept2][adj2.row] == band[kept2][adj2.col])
    g2 = coo_matrix((np.ones(same.sum()), (adj2.row[same], adj2.col[same])), shape=(len(kept2),) * 2)
    _n, comp = connected_components(g2, directed=False)
    for c in range(comp.max() + 1 if len(kept2) else 0):
        ids = kept2[comp == c]
        if areas2[ids].sum() < MIN_ZONE_AREA:
            status2[ids] = 5
            continue
        zone[ids] = len(zone_moment)
        zone_moment.append(int(moment2[ids].max()))          # au plus tot
    zone_moment = np.array(zone_moment, dtype=int)
    # micro-zones : rattachees a une voisine si la fusion reste accessible,
    # sinon NON lissees (decision utilisateur du 2026-09-28) -- statut 6, hors zone,
    # donc ni creusees par P3 ni tracees par P4
    zone, merged_moment, zone_merges, micro_left = sp.merge_micro_zones(
        v, f, areas2, zone, zone_moment, last2, micro_inradius)
    for z in micro_left:
        status2[zone == z] = 6
        zone[zone == z] = -1
    live = sorted(set(zone[zone >= 0].tolist()))
    remap = np.full(len(zone_moment), -1)
    remap[live] = np.arange(len(live))
    zone = np.where(zone >= 0, remap[np.maximum(zone, 0)], -1)
    zone_moment = np.array([merged_moment[z] for z in live], dtype=int)
    zone_merges = np.array([(remap[q] if remap[q] >= 0 else -1, zz, m) for zz, q, _a, _b, m, _l in zone_merges],
                           dtype=int).reshape(-1, 3)
    normals2, angles2 = normals[origin], angles[origin]
    return dict(vertices=v, faces=f, status=status2, moment=moment2, zone=zone,
                zone_moment=zone_moment, angles=angles2, normals=normals2, areas=areas2, isovalues=isovalues,
                nozzle_width=NOZZLE_WIDTH, last_access=last2, origin=origin, edge_piece=piece2,
                surface_vertices=surf_vertices, surface_faces=faces, surface_status=status,
                edge_face_angle=edge_face_angle, zone_merges=zone_merges, micro_inradius=micro_inradius)
