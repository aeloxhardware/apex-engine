import sys, re
import numpy as np
import trimesh
import networkx as nx
from collections import Counter
from scipy.spatial import cKDTree
from shapely.geometry import Polygon
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1

exec(open("scripts/prototype_local_retriangulation13.py").read().split("total_loops = 0")[0])

for comp in components:
    comp_edges = [e for e in boundary_edges if e[0] in comp and e[1] in comp]
    loops = _split_loop_at_pinch_points_v2(comp_edges)
    for loop in loops:
        loop_pts = new_mesh.vertices[loop]
        centroid, axis_u, axis_v = _local_plane_frame(loop_pts)
        pts2d = np.column_stack([(loop_pts - centroid) @ axis_u, (loop_pts - centroid) @ axis_v])
        poly = Polygon(pts2d)
        if not (poly.is_valid and poly.is_simple):
            print(f"loop {loop} len={len(loop)}: valid={poly.is_valid}, simple={poly.is_simple}, area={poly.area:.4f}")
            print(f"  reason: {shapely_explain(poly) if False else 'n/a'}")

def shapely_explain(poly):
    return None
