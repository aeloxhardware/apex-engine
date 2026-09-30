"""Inspect the two failing loops (size 54, 83) to understand why they
self-intersect when flattened, and find a robust fix."""
import sys, re
import numpy as np
import trimesh
import networkx as nx
from collections import Counter
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1

exec(open("scripts/prototype_local_retriangulation4.py").read().split("total_loops = 0")[0])

for comp in components:
    comp_edges = [e for e in boundary_edges if e[0] in comp and e[1] in comp]
    loops = _split_loop_at_pinch_points(comp_edges)
    for loop in loops:
        if len(loop) in (54, 83):
            pts = new_mesh.vertices[loop]
            print(f"loop size {len(loop)}: bbox = {pts.min(axis=0)} .. {pts.max(axis=0)}")
            centroid, axis_u, axis_v = _local_plane_frame(pts)
            pts2d = np.column_stack([(pts - centroid) @ axis_u, (pts - centroid) @ axis_v])
            from shapely.geometry import Polygon
            poly = Polygon(pts2d)
            print(f"  is_valid={poly.is_valid}, is_simple={poly.is_simple}")
            # check planarity: distance of points from best-fit plane
            normal = np.cross(axis_u, axis_v)
            dist = (pts - centroid) @ normal
            print(f"  out-of-plane distance: max={np.abs(dist).max():.4f}, mean={np.abs(dist).mean():.4f}")
