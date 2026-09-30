"""Back to v4's blanket removal (best triangulation success rate: 21/23),
but handle the 2 non-planar failing loops properly: for a highly non-planar
loop, fall back to trimesh.repair.fill_holes ON JUST THAT SUBSET (extract
as isolated submesh, fill, reinsert) rather than a flat 2D triangulation
which assumes near-planarity.
"""
import sys, re
import numpy as np
import trimesh
import networkx as nx
from collections import Counter
from shapely.geometry import Polygon
from trimesh.creation import triangulate_polygon
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1


def _local_plane_frame(points):
    centroid = points.mean(axis=0)
    centered = points - centroid
    _, _, vt = np.linalg.svd(centered)
    normal = vt[-1]
    arbitrary = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(arbitrary, normal)) > 0.9:
        arbitrary = np.array([0.0, 1.0, 0.0])
    axis_u = np.cross(normal, arbitrary)
    axis_u /= np.linalg.norm(axis_u)
    axis_v = np.cross(normal, axis_u)
    return centroid, axis_u, axis_v


def _split_loop_at_pinch_points(loop_edges):
    loops = []
    remaining = set(tuple(sorted(e)) for e in loop_edges)
    while remaining:
        e0 = next(iter(remaining))
        remaining.discard(e0)
        loop = [e0[0], e0[1]]
        cur = e0[1]
        while True:
            candidates = [e for e in remaining if cur in e]
            if not candidates:
                break
            nxt_edge = candidates[0]
            nxt = nxt_edge[0] if nxt_edge[1] == cur else nxt_edge[1]
            remaining.discard(nxt_edge)
            if nxt == loop[0]:
                break
            loop.append(nxt)
            cur = nxt
        loops.append(loop)
    return loops


MODEL = "res/test_models/turkish.stl"
raw = trimesh.load(MODEL)
cleaned, report = s1._clean_and_repair(raw)
print(f"cleaned: faces={len(cleaned.faces)}, watertight={cleaned.is_watertight}")

log_path = r"C:/Users/prazl/AppData/Local/Temp/claude/r--apex-alg-apex-beta/d91361a4-152b-4fdf-8401-377effdfa068/scratchpad/bad_triangles.log"
text = open(log_path).read()
flagged = set()
for m in re.finditer(r'\[([\d,]+)\]', text):
    for tok in m.group(1).split(','):
        if tok.strip().isdigit():
            flagged.add(int(tok.strip()))

face_mask = np.isin(cleaned.faces, list(flagged)).any(axis=1)
removed_verts = set(cleaned.faces[face_mask].ravel().tolist())
face_mask2 = np.isin(cleaned.faces, list(removed_verts)).any(axis=1)
keep_mask = ~face_mask2
new_mesh = cleaned.copy()
new_mesh.update_faces(keep_mask)
new_mesh.remove_unreferenced_vertices()
print(f"after removal: faces={len(new_mesh.faces)}")

edges = new_mesh.edges_sorted
edge_counts = Counter(map(tuple, edges))
boundary_edges = [e for e, c in edge_counts.items() if c == 1]
g = nx.Graph()
g.add_edges_from(boundary_edges)
components = list(nx.connected_components(g))

new_verts = list(new_mesh.vertices)
new_faces = list(new_mesh.faces)
total_loops = 0
failed = []
for comp in components:
    comp_edges = [e for e in boundary_edges if e[0] in comp and e[1] in comp]
    loops = _split_loop_at_pinch_points(comp_edges)
    for loop in loops:
        if len(loop) < 3:
            continue
        total_loops += 1
        loop_pts = new_mesh.vertices[loop]
        centroid, axis_u, axis_v = _local_plane_frame(loop_pts)
        pts2d = np.column_stack([(loop_pts - centroid) @ axis_u, (loop_pts - centroid) @ axis_v])
        poly = Polygon(pts2d)
        if poly.is_valid and poly.is_simple:
            v2b, f2b = triangulate_polygon(poly, force_vertices=True)
            if len(v2b) == len(pts2d):
                for tri in f2b:
                    new_faces.append([loop[tri[0]], loop[tri[1]], loop[tri[2]]])
                continue
        # FALLBACK: non-planar / self-intersecting-when-flattened loop.
        # Use a fan triangulation from the loop's 3D centroid instead of a
        # flat 2D polygon fill - works for any simple closed loop regardless
        # of planarity (each triangle uses two consecutive loop vertices +
        # a synthetic center vertex), at the cost of a slightly non-optimal
        # triangulation (thin/irregular triangles allowed, since this is a
        # tiny defect patch, not user-visible print geometry quality).
        failed.append(len(loop))
        center_idx = len(new_verts)
        new_verts.append(loop_pts.mean(axis=0).tolist())
        for i in range(len(loop)):
            a = loop[i]
            b = loop[(i + 1) % len(loop)]
            new_faces.append([a, b, center_idx])

print(f"total loops: {total_loops}, fan-fallback used for: {failed}")

final_mesh = trimesh.Trimesh(vertices=np.array(new_verts), faces=np.array(new_faces), process=False)
trimesh.repair.fix_normals(final_mesh)
trimesh.repair.fix_winding(final_mesh)
print(f"final mesh: faces={len(final_mesh.faces)}, watertight={final_mesh.is_watertight}")
print(f"z range: {final_mesh.vertices[:,2].min():.2f} .. {final_mesh.vertices[:,2].max():.2f}")

import pickle
with open(r"C:/Users/prazl/AppData/Local/Temp/claude/r--apex-alg-apex-beta/d91361a4-152b-4fdf-8401-377effdfa068/scratchpad/final_mesh.pkl", "wb") as fh:
    pickle.dump((np.array(final_mesh.vertices), np.array(final_mesh.faces)), fh)

edges3 = final_mesh.edges_sorted
edge_counts3 = Counter(map(tuple, edges3))
boundary_edges3 = [e for e, c in edge_counts3.items() if c == 1]
print(f"remaining open edges: {len(boundary_edges3)}")
g3 = nx.Graph()
g3.add_edges_from(boundary_edges3)
comps3 = list(nx.connected_components(g3))
print(f"remaining boundary components: {len(comps3)}, sizes: {sorted([len(c) for c in comps3], reverse=True)}")
broken = trimesh.repair.broken_faces(final_mesh)
print(f"broken_faces: {len(broken)}")
