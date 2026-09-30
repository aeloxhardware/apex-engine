"""Prototype: locally repair turkish.stl's self-intersecting patches by
removing the faces incident to flagged vertex clusters and re-filling the
resulting holes, leaving the rest of the mesh completely untouched.
"""
import sys
import re
import numpy as np
import trimesh
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1

MODEL = "res/test_models/turkish.stl"
raw = trimesh.load(MODEL)
cleaned, report = s1._clean_and_repair(raw)
print(f"cleaned: faces={len(cleaned.faces)}, verts={len(cleaned.vertices)}, watertight={cleaned.is_watertight}")

log_path = r"C:/Users/prazl/AppData/Local/Temp/claude/r--apex-alg-apex-beta/d91361a4-152b-4fdf-8401-377effdfa068/scratchpad/bad_triangles.log"
text = open(log_path).read()
indices = set()
for m in re.finditer(r'\[([\d,]+)\]', text):
    for tok in m.group(1).split(','):
        if tok.strip().isdigit():
            indices.add(int(tok.strip()))
print(f"flagged vertex indices: {len(indices)}")

# Expand: any face touching a flagged vertex gets removed. Then also
# expand by one more ring (faces touching the REMOVED faces' vertices) to
# make sure we clear out the whole tangled micro-patch, not just its core.
face_mask = np.isin(cleaned.faces, list(indices)).any(axis=1)
print(f"faces directly touching flagged vertices: {face_mask.sum()}")

removed_verts = set(cleaned.faces[face_mask].ravel().tolist())
face_mask2 = np.isin(cleaned.faces, list(removed_verts)).any(axis=1)
print(f"faces after 1-ring expansion: {face_mask2.sum()}")

keep_mask = ~face_mask2
new_mesh = cleaned.copy()
new_mesh.update_faces(keep_mask)
new_mesh.remove_unreferenced_vertices()
print(f"after removal: faces={len(new_mesh.faces)}, watertight={new_mesh.is_watertight}")
print(f"z range: {new_mesh.vertices[:,2].min():.2f} .. {new_mesh.vertices[:,2].max():.2f}")

# fill holes
trimesh.repair.fill_holes(new_mesh)
trimesh.repair.fix_normals(new_mesh)
trimesh.repair.fix_winding(new_mesh)
print(f"after fill_holes: faces={len(new_mesh.faces)}, watertight={new_mesh.is_watertight}")
print(f"z range: {new_mesh.vertices[:,2].min():.2f} .. {new_mesh.vertices[:,2].max():.2f}")
