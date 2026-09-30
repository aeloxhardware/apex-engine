import sys
import numpy as np
import tetgen
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1
import trimesh

MODEL = "res/test_models/turkish.stl"
raw = trimesh.load(MODEL)
cleaned, report = s1._clean_and_repair(raw)

repaired = s1._repair_with_pymeshfix(cleaned)
print(f"repaired is None: {repaired is None}")
if repaired is not None:
    print(f"repaired watertight={repaired.is_watertight}, n_faces={len(repaired.faces)} (orig {len(cleaned.faces)})")
    print(f"repaired z range: {repaired.vertices[:,2].min()} .. {repaired.vertices[:,2].max()}")

    hole_points = s1._find_internal_cavities(repaired)
    resolution = 1.5
    max_tet_volume = resolution ** 3
    tet = tetgen.TetGen(repaired.vertices, repaired.faces)
    for hp in hole_points:
        tet.add_hole(hp)
    print("Calling tetrahedralize on PYMESHFIX-REPAIRED mesh (no platform refinement)...")
    tet.tetrahedralize(switches=f"pq1.2a{max_tet_volume:.10f}")
    print("done, no crash")
    node = np.asarray(tet.node)
    elem = np.asarray(tet.elem)
    print(f"n_node={len(node)}, n_elem={len(elem)}")
    print(f"node z range: {node[:,2].min()} .. {node[:,2].max()}")
