"""Check ALL vertex indices mentioned in TetGen's self-intersection warnings
to see if the double-surface-patch pattern (near-duplicate coplanar-ish
micro triangles) repeats everywhere, or if it's isolated to one spot."""
import sys
import re
import numpy as np
import trimesh
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1

MODEL = "res/test_models/turkish.stl"
raw = trimesh.load(MODEL)
cleaned, report = s1._clean_and_repair(raw)

log_path = r"C:/Users/prazl/AppData/Local/Temp/claude/r--apex-alg-apex-beta/d91361a4-152b-4fdf-8401-377effdfa068/scratchpad/bad_triangles.log"
text = open(log_path).read()

# extract all vertex indices mentioned anywhere in the warnings
indices = set()
for m in re.finditer(r'\[([\d,]+)\]', text):
    for tok in m.group(1).split(','):
        if tok.strip().isdigit():
            indices.add(int(tok.strip()))

indices = sorted(indices)
print(f"total distinct vertex indices flagged: {len(indices)}")

pts = cleaned.vertices[indices]

# cluster them spatially (simple: connected components within 3mm)
from scipy.spatial import cKDTree
tree = cKDTree(pts)
pairs = tree.query_pairs(r=3.0)
import networkx as nx
g = nx.Graph()
g.add_nodes_from(range(len(indices)))
g.add_edges_from(pairs)
components = list(nx.connected_components(g))
print(f"spatial clusters (radius 3mm): {len(components)}")
for comp in sorted(components, key=len, reverse=True)[:15]:
    comp_pts = pts[list(comp)]
    centroid = comp_pts.mean(axis=0)
    print(f"  cluster size={len(comp)}, centroid={centroid}, z_range=[{comp_pts[:,2].min():.2f},{comp_pts[:,2].max():.2f}]")
