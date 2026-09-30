import sys
import numpy as np
import tetgen
import trimesh
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1

MODEL = "res/test_models/turkish.stl"
raw = trimesh.load(MODEL)
cleaned, report = s1._clean_and_repair(raw)
unioned = cleaned.union(cleaned, engine="manifold")
print(f"unioned: faces={len(unioned.faces)}, watertight={unioned.is_watertight}")

hole_points = s1._find_internal_cavities(unioned)
resolution = 1.5
max_tet_volume = resolution ** 3
tet = tetgen.TetGen(unioned.vertices, unioned.faces)
for hp in hole_points:
    tet.add_hole(hp)
print("Calling tetrahedralize on MANIFOLD-UNIONED mesh...")
tet.tetrahedralize(switches=f"pq1.2a{max_tet_volume:.10f}")
print("done, no crash")
node = np.asarray(tet.node)
elem = np.asarray(tet.elem)
print(f"n_node={len(node)}, n_elem={len(elem)}")
print(f"node z range: {node[:,2].min()} .. {node[:,2].max()}")
