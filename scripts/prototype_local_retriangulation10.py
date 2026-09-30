"""Diagnose WHY the two large boundary loops (54, 83 verts) are so
non-planar: are they really one physical hole, or does the blanket 1-ring
expansion connect multiple SEPARATE clusters via a shared face/vertex,
merging them into one giant boundary that spans much of the model?"""
import sys
import re
import numpy as np
import trimesh
import networkx as nx
from collections import Counter, defaultdict
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1

MODEL = "res/test_models/turkish.stl"
raw = trimesh.load(MODEL)
cleaned, report = s1._clean_and_repair(raw)

log_path = r"C:/Users/prazl/AppData/Local/Temp/claude/r--apex-alg-apex-beta/d91361a4-152b-4fdf-8401-377effdfa068/scratchpad/bad_triangles.log"
text = open(log_path).read()
flagged = set()
for m in re.finditer(r'\[([\d,]+)\]', text):
    for tok in m.group(1).split(','):
        if tok.strip().isdigit():
            flagged.add(int(tok.strip()))

face_mask = np.isin(cleaned.faces, list(flagged)).any(axis=1)
direct_faces = set(np.where(face_mask)[0].tolist())
print(f"direct faces (touching flagged verts): {len(direct_faces)}")

removed_verts = set(cleaned.faces[face_mask].ravel().tolist())
face_mask2 = np.isin(cleaned.faces, list(removed_verts)).any(axis=1)
expanded_faces = set(np.where(face_mask2)[0].tolist())
print(f"expanded faces (1-ring via vertex sharing): {len(expanded_faces)}")

# Are these removed faces spatially compact, or scattered across the model?
removed_face_centroids = cleaned.vertices[cleaned.faces[list(expanded_faces)]].mean(axis=1)
print(f"removed-face centroid bbox: {removed_face_centroids.min(axis=0)} .. {removed_face_centroids.max(axis=0)}")

# Connected components of the removed-face SET itself (via face adjacency) -
# if the removal is spatially scattered, it'll have many disconnected
# blobs; the two giant boundary loops would come from ONE huge connected
# blob that happens to be shaped like a long strip (thin removed regions
# produce big, thin, non-planar holes even if "locally" connected).
face_adj = cleaned.face_adjacency
adj_map = defaultdict(set)
for a, b in face_adj:
    if a in expanded_faces and b in expanded_faces:
        adj_map[a].add(b)
        adj_map[b].add(a)
g = nx.Graph()
g.add_nodes_from(expanded_faces)
for a, nbrs in adj_map.items():
    for b in nbrs:
        g.add_edge(a, b)
comps = list(nx.connected_components(g))
print(f"connected blobs of removed faces: {len(comps)}, sizes: {sorted([len(c) for c in comps], reverse=True)}")

# for the biggest blob(s), check spatial extent (is it compact or a long strip?)
for comp in sorted(comps, key=len, reverse=True)[:3]:
    pts = cleaned.vertices[cleaned.faces[list(comp)]].reshape(-1, 3)
    bbox = pts.max(axis=0) - pts.min(axis=0)
    print(f"  blob size={len(comp)}: bbox extent={bbox}, diag={np.linalg.norm(bbox):.2f}mm")
