"""Inspect individual boundary loops left after local removal - are they
simple closed polygons (fixable by ear-clipping / Delaunay-in-plane), or
genuinely broken topology (multiple loops per hole, non-planar, etc)?"""
import sys
import re
import numpy as np
import trimesh
import networkx as nx
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1

MODEL = "res/test_models/turkish.stl"
raw = trimesh.load(MODEL)
cleaned, report = s1._clean_and_repair(raw)

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

edges = new_mesh.edges_sorted
from collections import Counter
edge_counts = Counter(map(tuple, edges))
boundary_edges = [e for e, c in edge_counts.items() if c == 1]

# build a graph of boundary edges, find connected loops, check simplicity
g = nx.Graph()
g.add_edges_from(boundary_edges)
components = list(nx.connected_components(g))
print(f"boundary components: {len(components)}")

degree_hist = Counter(dict(g.degree()).values())
print(f"degree histogram (should be all 2 for simple loops): {degree_hist}")

sizes = sorted([len(c) for c in components], reverse=True)
print(f"component sizes: {sizes[:20]}")

# check if any vertex has degree != 2 (non-simple loop / branching)
non_simple = [n for n, d in g.degree() if d != 2]
print(f"vertices with boundary degree != 2: {len(non_simple)}")
if non_simple:
    for n in non_simple[:10]:
        print(f"  vertex {n}: degree {g.degree(n)}, pos={new_mesh.vertices[n]}")
