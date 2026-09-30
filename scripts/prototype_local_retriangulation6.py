"""Redo local repair PER SPATIAL CLUSTER (not one global blob) - process
each of the ~15 flagged neighborhoods independently with a tightly bounded
face-removal radius, so unrelated clusters never get bridged together into
one giant non-planar hole."""
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
flagged = set()
for m in re.finditer(r'\[([\d,]+)\]', text):
    for tok in m.group(1).split(','):
        if tok.strip().isdigit():
            flagged.add(int(tok.strip()))
flagged = sorted(flagged)
print(f"flagged vertex indices: {len(flagged)}")

# spatial clustering of flagged vertices (radius 3mm, per earlier diagnostic)
pts = cleaned.vertices[flagged]
tree = cKDTree(pts)
pairs = tree.query_pairs(r=3.0)
gclust = nx.Graph()
gclust.add_nodes_from(range(len(flagged)))
gclust.add_edges_from(pairs)
spatial_components = list(nx.connected_components(gclust))
clusters = [[flagged[i] for i in comp] for comp in spatial_components]
print(f"spatial clusters: {len(clusters)}")

# work on a mutable copy
work_mesh = cleaned.copy()
faces_to_remove = np.zeros(len(work_mesh.faces), dtype=bool)
patch_faces = []  # list of (v0,v1,v2) new triangles to add at the end

REMOVE_RADIUS = 4.0  # mm around each flagged vertex - local neighborhood only

for cluster_verts in clusters:
    cluster_pts = cleaned.vertices[cluster_verts]
    center = cluster_pts.mean(axis=0)
    # faces whose centroid falls within REMOVE_RADIUS of this cluster's center
    face_centroids = cleaned.vertices[cleaned.faces].mean(axis=1)
    dist = np.linalg.norm(face_centroids - center, axis=1)
    local_mask = dist < REMOVE_RADIUS
    faces_to_remove |= local_mask

print(f"total faces marked for removal (local radius {REMOVE_RADIUS}mm per cluster): {faces_to_remove.sum()}")

keep_mask = ~faces_to_remove
new_mesh = cleaned.copy()
new_mesh.update_faces(keep_mask)
new_mesh.remove_unreferenced_vertices()
print(f"after removal: faces={len(new_mesh.faces)}, watertight={new_mesh.is_watertight}")
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
total_new_faces = 0
failed_loops = []
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
        if not poly.is_valid or not poly.is_simple:
            failed_loops.append(len(loop))
            continue
        v2b, f2b = triangulate_polygon(poly, force_vertices=True)
        if len(v2b) != len(pts2d):
            failed_loops.append(len(loop))
            continue
        for tri in f2b:
            new_faces.append([loop[tri[0]], loop[tri[1]], loop[tri[2]]])
            total_new_faces += 1

print(f"total simple loops: {total_loops}, triangulated ok: {total_loops - len(failed_loops)}, failed: {failed_loops}")

final_mesh = trimesh.Trimesh(vertices=np.array(new_verts), faces=np.array(new_faces), process=False)
trimesh.repair.fix_normals(final_mesh)
trimesh.repair.fix_winding(final_mesh)
print(f"final mesh: faces={len(final_mesh.faces)}, watertight={final_mesh.is_watertight}")
print(f"z range: {final_mesh.vertices[:,2].min():.2f} .. {final_mesh.vertices[:,2].max():.2f}")
