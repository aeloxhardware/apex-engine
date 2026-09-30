"""Inspect the geometry around one flagged vertex cluster (4520-4537) to
understand WHY these triangles self-intersect - near-duplicate vertices?
degenerate slivers? genuine crossing geometry?"""
import sys
import numpy as np
import trimesh
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1

MODEL = "res/test_models/turkish.stl"
raw = trimesh.load(MODEL)
cleaned, report = s1._clean_and_repair(raw)
print(f"cleaned: faces={len(cleaned.faces)}, verts={len(cleaned.vertices)}")

cluster = list(range(4515, 4540))
pts = cleaned.vertices[cluster]
print("vertex positions in cluster 4515-4539:")
for i, p in zip(cluster, pts):
    print(f"  {i}: {p}")

# pairwise distances within cluster to find near-duplicates
from scipy.spatial.distance import pdist, squareform
d = squareform(pdist(pts))
np.fill_diagonal(d, np.inf)
print("\nclosest pairs (distance < 0.5mm):")
for i in range(len(cluster)):
    for j in range(i+1, len(cluster)):
        if d[i, j] < 0.5:
            print(f"  {cluster[i]} <-> {cluster[j]}: {d[i,j]:.6f}mm")

# faces touching this cluster
face_mask = np.isin(cleaned.faces, cluster).any(axis=1)
touching_faces = np.where(face_mask)[0]
print(f"\n{len(touching_faces)} faces touch this vertex range")
for fi in touching_faces[:20]:
    f = cleaned.faces[fi]
    area = cleaned.area_faces[fi]
    print(f"  face {fi}: verts={f}, area={area:.6f}")
