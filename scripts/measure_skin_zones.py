"""
Jalon 1 du lissage non-planaire sur la voie CHALEUR : quelle part de la surface peut recevoir une couche de
surface (P2 de la voie Bird, `pipeline/bird/skin/zones.py` de la branche `bird_old`) ?

Pour chaque face de la surface (re-decoupee a 5 largeurs de buse comme P2) :
  - moment porteur  = rang de la couche qui porte son sommet le plus haut (phi de la methode de la chaleur) ;
  - orientation     = IK de l'etape 8 (C impose par la position), buse vers -n, ecart tete / normale juge au centre
                      et aux trois sommets ;
  - dernier moment accessible = la matiere imprimee est {phi <= phi(t)} (croissante, comme avec Bird) : la vraie
                      tete (marge de securite comprise) posee au centre de la face ne doit pas entrer dans la piece
                      avant que la face soit portee ; plateau = plan relevee de la marge.
  Statuts comme P2 : retenue / inatteignable / collision plateau / collision objet / contact plateau.

Echantillon de faces (--faces=N, graine fixe) : le calcul complet coute des minutes en Python.
Usage : APEX_USE_CPP=0 python scripts/measure_skin_zones.py <modele.stl> [--faces=300] [--tolerance=45]
"""
from __future__ import annotations

import dataclasses
import sys
import time

import numpy as np
import trimesh
from scipy.spatial import cKDTree

ROOT = __import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
sys.path.insert(0, ROOT)

import importlib  # noqa: E402

M = importlib.import_module("main")
from pipeline import step1_ingest as s1, step2_overhang as s2, step3_segmentation as s3  # noqa: E402
from pipeline import step4_field_solver as s4, step5_layers as s5  # noqa: E402
from pipeline.step8_machine_toolpath import Step8Params, _RawPoint, _resolve_all_axis_angles_coupled  # noqa: E402
from pipeline.step9_collision_sim import _filter_nozzle_contact_zone, _surface_sample_points  # noqa: E402

NOZZLE_WIDTH = 0.4
MAX_EDGE = 5.0 * NOZZLE_WIDTH
PLATE_CONTACT_TOL = 1e-3
HEAD_SAMPLE_SPACING = 5.0
FACE_BATCH = 64


def barycentric(points, tets):
    a = tets[:, 0]
    m = np.stack([tets[:, 1] - a, tets[:, 2] - a, tets[:, 3] - a], axis=2)
    with np.errstate(all="ignore"):
        w = np.linalg.solve(m, (points - a)[..., None])[..., 0]
    return np.concatenate([1 - w.sum(axis=1, keepdims=True), w], axis=1)


class PrintedSolid:
    def __init__(self, tet_vertices, tetrahedra, phi_nodes, k=16):
        self.tets = tet_vertices[tetrahedra]
        self.phi = phi_nodes[tetrahedra]
        self.tree = cKDTree(self.tets.mean(axis=1))
        self.k = k
        self.reach = float(np.linalg.norm(self.tets - self.tets.mean(axis=1, keepdims=True), axis=2).max())
        self.lo, self.hi = tet_vertices.min(axis=0), tet_vertices.max(axis=0)
        self.tolerance = -1e-9

    def phi_at(self, points):
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


def rotation_matrices(angles_deg, axes):
    total = np.tile(np.eye(3), (len(angles_deg), 1, 1))
    for i, spec in reversed(list(enumerate(axes))):
        k = np.asarray(spec.axis, float)
        k = k / np.linalg.norm(k)
        kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        a = np.radians(angles_deg[:, i])[:, None, None]
        total = (np.eye(3) + np.sin(a) * kx + (1 - np.cos(a)) * (kx @ kx)) @ total
    return total


def boundary_faces(tetrahedra):
    quads = np.concatenate([tetrahedra[:, [0, 1, 2]], tetrahedra[:, [0, 1, 3]], tetrahedra[:, [0, 2, 3]],
                            tetrahedra[:, [1, 2, 3]]])
    key = np.sort(quads, axis=1)
    _u, first, counts = np.unique(key, axis=0, return_index=True, return_counts=True)
    return quads[first[counts == 1]]


def main():
    stl = [a for a in sys.argv[1:] if not a.startswith("--")][0]
    opts = dict(a[2:].split("=") for a in sys.argv[1:] if a.startswith("--") and "=" in a)
    n_sample, tolerance = int(opts.get("faces", 300)), float(opts.get("tolerance", 45))
    p = M.PipelineParams()
    t0 = time.time()
    md, tm, rep = s1.run(stl, resolution=p.resolution, near_platform_resolution=p.near_platform_resolution,
                         near_platform_margin=p.near_platform_margin, place_on_platform=p.place_on_platform)
    om = s2.run(md, p.overhang)
    seg = dataclasses.replace(p.segmentation, platform_z=rep.platform_z_used) if p.segmentation.platform_z is None else p.segmentation
    graph = s3.run(md, om, seg)
    field = s4.run(tm, md, om, graph, M._effective_field_solver_params(p))
    layers = s5.run(tm, md, field, graph, p.layers)
    isovalues = np.sort(np.unique([l.isovalue for l in layers]))
    print(f"{stl}: {len(tm.tetrahedra)} tetra, {len(layers)} couches, phi max {field.values.max():.1f} ({time.time()-t0:.0f} s)", flush=True)

    machine, sim = p.machine_profile, p.collision
    verts = tm.tet_vertices
    surf_v, faces = trimesh.remesh.subdivide_to_size(verts, boundary_faces(tm.tetrahedra), MAX_EDGE)
    solid = PrintedSolid(verts, tm.tetrahedra, field.values)
    solid.tolerance = -1e-6
    surf_phi = solid.phi_at(surf_v)
    solid.tolerance = -1e-9
    lost = np.isnan(surf_phi)
    if lost.any():
        surf_phi[lost] = field.values[cKDTree(verts).query(surf_v[lost])[1]]
    tri = surf_v[faces]
    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    normals = cross / np.maximum(2 * areas, 1e-12)[:, None]
    # sens des normales : vers l'exterieur (le maillage de frontiere des tetraedres n'a pas d'orientation garantie)
    probe = tri.mean(axis=1) + 1e-3 * normals
    inside = ~np.isnan(solid.phi_at(probe))
    normals[inside] *= -1
    centers = tri.mean(axis=1)
    status = np.zeros(len(faces), dtype=np.int8)
    status[tri[:, :, 2].max(axis=1) < verts[:, 2].min() + PLATE_CONTACT_TOL] = 4
    moment = np.clip(np.searchsorted(isovalues, surf_phi[faces].max(axis=1) - 1e-9, side="left"), 0, len(isovalues) - 1)

    rng = np.random.default_rng(0)
    free = np.flatnonzero(status == 0)
    todo = np.sort(rng.choice(free, size=min(n_sample, len(free)), replace=False))
    print(f"{len(faces)} faces ({areas.sum():.0f} mm2), dont {len(free)} hors plateau ; echantillon {len(todo)} faces", flush=True)

    axes = machine.rotation_axes
    raw = [_RawPoint(position=centers[i], desired_direction=normals[i], layer_id=0, segment_index=0, color_id=None,
                     feedrate=1.0, extrusion=0.0, is_travel=False) for i in todo]
    a, err = _resolve_all_axis_angles_coupled(raw, axes, np.array([0.0, 0.0, 1.0]), machine, Step8Params())
    angles = np.zeros((len(faces), len(axes)))
    residual = np.zeros(len(faces))
    angles[todo], residual[todo] = np.asarray(a), np.asarray(err)
    corners = tri[todo].reshape(-1, 3)
    raw_c = [_RawPoint(position=q, desired_direction=normals[todo[j // 3]], layer_id=0, segment_index=0, color_id=None,
                       feedrate=1.0, extrusion=0.0, is_travel=False) for j, q in enumerate(corners)]
    _a, err_c = _resolve_all_axis_angles_coupled(raw_c, axes, np.array([0.0, 0.0, 1.0]), machine, Step8Params())
    residual[todo] = np.maximum(residual[todo], np.asarray(err_c).reshape(-1, 3).max(axis=1))
    dev = np.degrees(residual)
    status[todo[dev[todo] > tolerance]] = 1

    head = _surface_sample_points(machine.nozzle_model, HEAD_SAMPLE_SPACING)
    head = _filter_nozzle_contact_zone(head, sim.nozzle_contact_clearance)
    _c, _d, tid = trimesh.proximity.closest_point(machine.nozzle_model, head)
    hn = machine.nozzle_model.face_normals[tid]
    head = np.vstack([head, head + sim.safety_margin * hn, head - sim.safety_margin * hn])
    print(f"tete : {len(head)} points ; calcul de l'accessibilite...", flush=True)
    reach = todo[status[todo] == 0]
    plate = np.zeros(len(reach), bool)
    last = np.full(len(reach), len(isovalues) - 1)
    t1 = time.time()
    for s in range(0, len(reach), FACE_BATCH):
        ids = reach[s:s + FACE_BATCH]
        rot = rotation_matrices(angles[ids], axes)
        pts = np.einsum("fij,pj->fpi", rot, head) + centers[ids][:, None, :]
        plate[s:s + len(ids)] = (pts[:, :, 2] < sim.safety_margin).any(axis=1)
        phi = solid.phi_at(pts.reshape(-1, 3)).reshape(len(ids), -1)
        low = np.min(np.where(np.isnan(phi), np.inf, phi), axis=1)
        ins = np.isfinite(low)
        last[s:s + len(ids)][ins] = np.searchsorted(isovalues, low[ins] - 1e-9, side="left") - 1
    print(f"accessibilite : {time.time()-t1:.0f} s pour {len(reach)} faces", flush=True)
    status[reach[plate]] = 2
    obj = np.zeros(len(faces), bool)
    obj[reach] = last < moment[reach]
    status[(status == 0) & obj] = 3

    names = {0: "retenue", 1: "inatteignable (tete a plus de %g deg)" % tolerance, 2: "collision plateau",
             3: "collision objet (bloquee avant d'etre portee)", 4: "contact plateau"}
    sel = todo
    area = areas[sel].sum()
    print(f"\nSur {len(sel)} faces tirees ({area:.0f} mm2) :")
    for code in (0, 1, 2, 3):
        m = status[sel] == code
        print(f"  {names[code]:52s} {m.sum():5d} faces  {100*areas[sel][m].sum()/area:5.1f} % de l'aire")
    print(f"ecart tete/normale (pire de centre+sommets) : mediane {np.median(dev[sel]):.1f} deg, "
          f"p90 {np.percentile(dev[sel],90):.1f} deg ; moments : couches {len(isovalues)}")
    for tol in (5, 15, 30, 45, 60):
        print(f"  ecart <= {tol:2d} deg : {100*areas[sel][dev[sel] <= tol].sum()/area:5.1f} % de l'aire")


if __name__ == "__main__":
    main()
