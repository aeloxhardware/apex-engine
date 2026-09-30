"""Full prototype: remove flagged micro-patches, split each boundary
component at pinch points into simple loops, triangulate each in its own
local best-fit plane, stitch back in. Verify final mesh is watertight,
z-range intact, and TetGen accepts it."""
import sys
import re
import numpy as np
import trimesh
import networkx as nx
from collections import Counter
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1
from core.geometry_utils import point_in_polygon_2d  # noqa: F401 (not used yet, checking availability)


def _local_plane_frame(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    centroid = points.mean(axis=0)
    centered = points - centroid
    # best-fit plane normal via SVD (smallest singular vector)
    _, _, vt = np.linalg.svd(centered)
    normal = vt[-1]
    # build an orthonormal basis (axis_u, axis_v) spanning the plane
    arbitrary = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(arbitrary, normal)) > 0.9:
        arbitrary = np.array([0.0, 1.0, 0.0])
    axis_u = np.cross(normal, arbitrary)
    axis_u /= np.linalg.norm(axis_u)
    axis_v = np.cross(normal, axis_u)
    return centroid, axis_u, axis_v


def _split_loop_at_pinch_points(loop_edges: list[tuple[int, int]]) -> list[list[int]]:
    """Given boundary edges of one connected component (possibly with
    degree-4 pinch vertices), split into simple closed loops (each vertex
    visited once) by walking edges and detaching at any degree>2 vertex."""
    g = nx.MultiGraph()
    g.add_edges_from(loop_edges)
    loops = []
    remaining = set(tuple(sorted(e)) for e in loop_edges)
    while remaining:
        e0 = next(iter(remaining))
        remaining.discard(e0)
        loop = [e0[0], e0[1]]
        cur = e0[1]
        prev = e0[0]
        while True:
            candidates = [
                e for e in remaining
                if cur in e
            ]
            if not candidates:
                break
            # pick the edge that isn't going back the way we came, prefer any
            nxt_edge = candidates[0]
            nxt = nxt_edge[0] if nxt_edge[1] == cur else nxt_edge[1]
            remaining.discard(nxt_edge)
            if nxt == loop[0]:
                break
            loop.append(nxt)
            prev, cur = cur, nxt
        loops.append(loop)
    return loops


MODEL = "res/test_models/turkish.stl"
raw = trimesh.load(MODEL)
cleaned, report = s1._clean_and_repair(raw)
print(f"cleaned: faces={len(cleaned.faces)}, watertight={cleaned.is_watertight}")

log_path = r"C:/Users/prazl/AppData/Local/Temp/claude/r--apex-alg-apex-beta/d91361a4-152b-4fdf-8401-377effdfa068/scratchpad/bad_triangles.log"
text = open(log_path).read()
indices = set()
for m in re.finditer(r'\[([\d,]+)\]', text):
    for tok in m.group(1).split(','):
        if tok.strip().isdigit():
            indices.add(int(tok.strip()))

face_mask = np.isin(cleaned.faces, list(indices)).any(axis=1)
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
print(f"boundary components: {len(components)}")

new_verts = list(new_mesh.vertices)
new_faces = list(new_mesh.faces)

total_loops = 0
total_new_faces = 0
for comp in components:
    comp_edges = [e for e in boundary_edges if e[0] in comp and e[1] in comp]
    loops = _split_loop_at_pinch_points(comp_edges)
    for loop in loops:
        if len(loop) < 3:
            continue
        total_loops += 1
        loop_pts = new_mesh.vertices[loop]
        centroid, axis_u, axis_v = _local_plane_frame(loop_pts)
        pts2d = np.column_stack([
            (loop_pts - centroid) @ axis_u,
            (loop_pts - centroid) @ axis_v,
        ])
        try:
            from trimesh.creation import triangulate_polygon
            from shapely.geometry import Polygon
            poly = Polygon(pts2d)
            if not poly.is_valid:
                poly = poly.buffer(0)
            v2, f2 = triangulate_polygon(poly, force_vertices=False)
            # map triangulated (possibly Steiner-added) verts back: v2 are
            # the polygon's own coords possibly reordered/deduped - need
            # correspondence back to original loop vertex indices. Use
            # force_vertices=True to guarantee v2 == input coords, same
            # order, no Steiner points added.
        except Exception as e:
            print(f"  loop size {len(loop)} triangulation FAILED: {e}")
            continue

        v2b, f2b = triangulate_polygon(Polygon(pts2d), force_vertices=True)
        if len(v2b) != len(pts2d):
            print(f"  loop size {len(loop)}: force_vertices mismatch ({len(v2b)} vs {len(pts2d)})")
            continue
        for tri in f2b:
            new_faces.append([loop[tri[0]], loop[tri[1]], loop[tri[2]]])
            total_new_faces += 1

print(f"total simple loops found: {total_loops}, total new faces added: {total_new_faces}")

final_mesh = trimesh.Trimesh(vertices=np.array(new_verts), faces=np.array(new_faces), process=False)
print(f"final mesh: faces={len(final_mesh.faces)}, watertight={final_mesh.is_watertight}")
print(f"z range: {final_mesh.vertices[:,2].min():.2f} .. {final_mesh.vertices[:,2].max():.2f}")

trimesh.repair.fix_normals(final_mesh)
trimesh.repair.fix_winding(final_mesh)
print(f"after fix_normals/winding: watertight={final_mesh.is_watertight}")

edges3 = final_mesh.edges_sorted
edge_counts3 = Counter(map(tuple, edges3))
boundary_edges3 = [e for e, c in edge_counts3.items() if c == 1]
print(f"remaining open edges: {len(boundary_edges3)}")
g3 = nx.Graph()
g3.add_edges_from(boundary_edges3)
comps3 = list(nx.connected_components(g3))
print(f"remaining boundary components: {len(comps3)}, sizes: {sorted([len(c) for c in comps3], reverse=True)}")
