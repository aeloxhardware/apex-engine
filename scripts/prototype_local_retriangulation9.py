import pickle
import numpy as np
import tetgen
import trimesh

with open(r"C:/Users/prazl/AppData/Local/Temp/claude/r--apex-alg-apex-beta/d91361a4-152b-4fdf-8401-377effdfa068/scratchpad/final_mesh.pkl", "rb") as fh:
    verts, faces = pickle.load(fh)

mesh = trimesh.Trimesh(vertices=verts, faces=faces, process=False)
print(f"watertight={mesh.is_watertight}, faces={len(mesh.faces)}")
broken = trimesh.repair.broken_faces(mesh)
print(f"broken faces: {len(broken)}")
for fi in broken[:10]:
    f = mesh.faces[fi]
    area = mesh.area_faces[fi]
    print(f"  face {fi}: verts={f}, area={area:.6f}, pts={mesh.vertices[f]}")

import sys
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1
hole_points = s1._find_internal_cavities(mesh)
resolution = 1.5
max_tet_volume = resolution ** 3
tet = tetgen.TetGen(mesh.vertices, mesh.faces)
for hp in hole_points:
    tet.add_hole(hp)
print("Calling tetrahedralize...")
tet.tetrahedralize(switches=f"pq1.2a{max_tet_volume:.10f}")
print("SUCCESS, no crash")
node = np.asarray(tet.node)
print(f"n_node={len(node)}, node z range: {node[:,2].min():.2f} .. {node[:,2].max():.2f}")
