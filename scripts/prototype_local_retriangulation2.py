"""Diagnose why fill_holes leaves the mesh non-watertight after local
removal, and try a stronger local hole-filling strategy."""
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

# inspect boundary edges (open edges) before fill
edges = new_mesh.edges_sorted
from collections import Counter
edge_counts = Counter(map(tuple, edges))
boundary_edges = [e for e, c in edge_counts.items() if c == 1]
print(f"boundary (open) edges before fill: {len(boundary_edges)}")

# group boundary edges into loops via trimesh's own facility
try:
    outline = new_mesh.outline()
    print(f"outline entities: {len(outline.entities)}")
except Exception as e:
    print(f"outline failed: {e}")

trimesh.repair.fill_holes(new_mesh)
edges2 = new_mesh.edges_sorted
edge_counts2 = Counter(map(tuple, edges2))
boundary_edges2 = [e for e, c in edge_counts2.items() if c == 1]
print(f"boundary edges after fill_holes: {len(boundary_edges2)}")
print(f"watertight after fill_holes: {new_mesh.is_watertight}")

# try trimesh.repair.broken_faces
broken = trimesh.repair.broken_faces(new_mesh)
print(f"broken_faces count: {len(broken)}")
