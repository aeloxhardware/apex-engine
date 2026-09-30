"""Inspect one perfectly-planar (out_of_plane=0.0) but still-failing loop
to see exactly why shapely reports it as non-simple."""
import sys, re
import numpy as np
import trimesh
import networkx as nx
from collections import Counter
from scipy.spatial import cKDTree
from shapely.geometry import Polygon
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1

exec(open("scripts/prototype_local_retriangulation11.py").read().split("total_loops = 0")[0])

for comp in components:
    comp_edges = [e for e in boundary_edges if e[0] in comp and e[1] in comp]
    loops = _split_loop_at_pinch_points(comp_edges)
    for loop in loops:
        if len(loop) in (15, 10):
            loop_pts = new_mesh.vertices[loop]
            centroid, axis_u, axis_v = _local_plane_frame(loop_pts)
            pts2d = np.column_stack([(loop_pts - centroid) @ axis_u, (loop_pts - centroid) @ axis_v])
            poly = Polygon(pts2d)
            print(f"loop size {len(loop)}: is_valid={poly.is_valid}, is_simple={poly.is_simple}")
            print(f"  loop vertex order: {loop}")
            print(f"  2D points:\n{pts2d}")
