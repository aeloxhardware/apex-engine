"""Fix the loop-splitting bug: the previous edge-walker didn't properly
detach at pinch points (degree>2 vertices), producing self-touching loops.
New approach: at a pinch point, greedily close the CURRENT loop as soon as
we return to ANY already-visited vertex in this walk (not just the start),
splitting there - guarantees every emitted loop is a simple cycle.
"""
import sys, re
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


def _split_loop_at_pinch_points_v2(loop_edges):
    """Decompose a (possibly branching, at degree>2 pinch vertices) set of
    boundary edges into SIMPLE closed loops. Each emitted loop visits every
    vertex at most once."""
    remaining = set(tuple(sorted(e)) for e in loop_edges)
    adjacency = {}
    for a, b in remaining:
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)

    loops = []
    while remaining:
        start_edge = next(iter(remaining))
        remaining.discard(start_edge)
        path = [start_edge[0], start_edge[1]]
        visited_in_path = {start_edge[0]: 0, start_edge[1]: 1}
        cur = start_edge[1]
        closed = False
        while True:
            nbrs = [
                v for v in adjacency.get(cur, ())
                if tuple(sorted((cur, v))) in remaining
            ]
            if not nbrs:
                break
            nxt = nbrs[0]
            edge = tuple(sorted((cur, nxt)))
            remaining.discard(edge)
            if nxt in visited_in_path:
                # closes a sub-loop: emit path[visited_in_path[nxt]:] as a
                # simple loop, then continue walking from nxt with the
                # remainder of path truncated (pinch point handled).
                start_idx = visited_in_path[nxt]
                loops.append(path[start_idx:])
                # truncate path back to the pinch point and keep going
                path = path[:start_idx + 1]
                visited_in_path = {v: i for i, v in enumerate(path)}
                cur = nxt
                continue
            path.append(nxt)
            visited_in_path[nxt] = len(path) - 1
            cur = nxt
        if len(path) >= 3:
            loops.append(path)
    return [loop for loop in loops if len(loop) >= 3]


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
flagged_pts = cleaned.vertices[flagged]

EXPAND_RADIUS = 2.0
face_centroids = cleaned.vertices[cleaned.faces].mean(axis=1)
tree = cKDTree(flagged_pts)
dists, _ = tree.query(face_centroids, k=1)
faces_to_remove = dists < EXPAND_RADIUS
print(f"faces removed: {faces_to_remove.sum()}")

keep_mask = ~faces_to_remove
new_mesh = cleaned.copy()
new_mesh.update_faces(keep_mask)
new_mesh.remove_unreferenced_vertices()

edges = new_mesh.edges_sorted
edge_counts = Counter(map(tuple, edges))
boundary_edges = [e for e, c in edge_counts.items() if c == 1]
g = nx.Graph()
g.add_edges_from(boundary_edges)
components = list(nx.connected_components(g))
print(f"boundary components: {len(components)}")

new_verts = list(new_mesh.vertices)
new_faces = list(new_mesh.faces)
total_loops = 0
failed = []
for comp in components:
    comp_edges = [e for e in boundary_edges if e[0] in comp and e[1] in comp]
    loops = _split_loop_at_pinch_points_v2(comp_edges)
    for loop in loops:
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
        failed.append(len(loop))

print(f"total loops: {total_loops}, failed: {failed}")

final_mesh = trimesh.Trimesh(vertices=np.array(new_verts), faces=np.array(new_faces), process=False)
trimesh.repair.fix_normals(final_mesh)
trimesh.repair.fix_winding(final_mesh)
print(f"final mesh: faces={len(final_mesh.faces)}, watertight={final_mesh.is_watertight}")
print(f"z range: {final_mesh.vertices[:,2].min():.2f} .. {final_mesh.vertices[:,2].max():.2f}")
broken = trimesh.repair.broken_faces(final_mesh)
print(f"broken_faces: {len(broken)}")

import pickle
with open(r"C:/Users/prazl/AppData/Local/Temp/claude/r--apex-alg-apex-beta/d91361a4-152b-4fdf-8401-377effdfa068/scratchpad/final_mesh_v2.pkl", "wb") as fh:
    pickle.dump((np.array(final_mesh.vertices), np.array(final_mesh.faces)), fh)
