"""Root cause of huge non-planar loops found: the 1-ring vertex-sharing
expansion jumps across the mesh via a single shared corner vertex touched
by many unrelated faces (e.g. bbox corner [100,100,30.26]). Fix: bound
expansion by EUCLIDEAN distance from each flagged vertex, not topological
vertex-sharing, so removal never leaves its local neighborhood - even if
that means NOT removing every face nominally "touching" a flagged vertex
(a corner vertex shared by 50 faces should only pull in the few faces
actually near the defect, not all 50).
"""
import sys
import re
import numpy as np
import trimesh
import networkx as nx
from collections import Counter
from scipy.spatial import cKDTree
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
flagged = sorted(set(
    int(tok) for m in re.finditer(r'\[([\d,]+)\]', text)
    for tok in m.group(1).split(',') if tok.strip().isdigit()
))
print(f"flagged vertex indices: {len(flagged)}")

flagged_pts = cleaned.vertices[flagged]

# Distance-bounded expansion: a face is removed only if its CENTROID is
# within EXPAND_RADIUS of ANY flagged vertex. This can never jump across
# the model via a single far-shared vertex like the topological 1-ring did.
EXPAND_RADIUS = 2.0  # mm - just enough to clear the immediate defect patch
face_centroids = cleaned.vertices[cleaned.faces].mean(axis=1)
tree = cKDTree(flagged_pts)
dists, _ = tree.query(face_centroids, k=1)
faces_to_remove = dists < EXPAND_RADIUS
print(f"faces marked for removal (distance-bounded, radius={EXPAND_RADIUS}mm): {faces_to_remove.sum()}")

keep_mask = ~faces_to_remove
new_mesh = cleaned.copy()
new_mesh.update_faces(keep_mask)
new_mesh.remove_unreferenced_vertices()
print(f"after removal: faces={len(new_mesh.faces)}")
print(f"z range: {new_mesh.vertices[:,2].min():.2f} .. {new_mesh.vertices[:,2].max():.2f}")

edges = new_mesh.edges_sorted
edge_counts = Counter(map(tuple, edges))
boundary_edges = [e for e, c in edge_counts.items() if c == 1]
g = nx.Graph()
g.add_edges_from(boundary_edges)
components = list(nx.connected_components(g))
print(f"boundary components: {len(components)}, sizes: {sorted([len(c) for c in components], reverse=True)}")

new_verts = list(new_mesh.vertices)
new_faces = list(new_mesh.faces)
total_loops = 0
failed = []
out_of_plane_stats = []
for comp in components:
    comp_edges = [e for e in boundary_edges if e[0] in comp and e[1] in comp]
    loops = _split_loop_at_pinch_points(comp_edges)
    for loop in loops:
        if len(loop) < 3:
            continue
        total_loops += 1
        loop_pts = new_mesh.vertices[loop]
        centroid, axis_u, axis_v = _local_plane_frame(loop_pts)
        normal = np.cross(axis_u, axis_v)
        out_of_plane = np.abs((loop_pts - centroid) @ normal).max()
        out_of_plane_stats.append(out_of_plane)
        pts2d = np.column_stack([(loop_pts - centroid) @ axis_u, (loop_pts - centroid) @ axis_v])
        poly = Polygon(pts2d)
        if poly.is_valid and poly.is_simple:
            v2b, f2b = triangulate_polygon(poly, force_vertices=True)
            if len(v2b) == len(pts2d):
                for tri in f2b:
                    new_faces.append([loop[tri[0]], loop[tri[1]], loop[tri[2]]])
                continue
        failed.append((len(loop), out_of_plane))

print(f"total loops: {total_loops}, out-of-plane max per loop: {sorted(out_of_plane_stats, reverse=True)[:10]}")
print(f"failed (size, out_of_plane): {failed}")

final_mesh = trimesh.Trimesh(vertices=np.array(new_verts), faces=np.array(new_faces), process=False)
trimesh.repair.fix_normals(final_mesh)
trimesh.repair.fix_winding(final_mesh)
print(f"final mesh: faces={len(final_mesh.faces)}, watertight={final_mesh.is_watertight}")
print(f"z range: {final_mesh.vertices[:,2].min():.2f} .. {final_mesh.vertices[:,2].max():.2f}")
broken = trimesh.repair.broken_faces(final_mesh)
print(f"broken_faces: {len(broken)}")
