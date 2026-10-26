"""Tranche une pièce avec un maillage tétraédrique fourni (TetraNet) à la place de TetGen.

Usage : python slice_with_tetranet.py <stl> <prefix_tetranet> <sortie.pkl>
  <prefix_tetranet> : fichiers <prefix>.nodes (float64 x3) et <prefix>.tets (int32 x4) du CLI tetranet.
Le moteur est le chemin par défaut de main.py (méthode de la chaleur) ; seule l'étape de
tétraédrisation est remplacée.
"""
import os, pickle, sys, time
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
stl, prefix, out = sys.argv[1:4]

from core.mesh_types import TetMesh
import pipeline.step1_ingest as s1
import main

nodes = np.fromfile(prefix + ".nodes", dtype="<f8").reshape(-1, 3)
tets = np.fromfile(prefix + ".tets", dtype="<i4").reshape(-1, 4).astype(np.int64)
print(f"[tetranet] maillage injecté : {len(nodes)} sommets, {len(tets)} tétraèdres", flush=True)


def injected(mesh, *a, **k):
    # le moteur a pu translater le modèle (pose sur le plateau) : on aligne les sommets
    dz = float(mesh.bounds[0][2]) - float(nodes[:, 2].min())
    tm = TetMesh(tet_vertices=nodes + np.array([0.0, 0.0, dz]), tetrahedra=tets)
    tm.tet_adjacency = s1._build_tet_adjacency(tets, len(nodes))
    print(f"[tetranet] décalage z appliqué : {dz:.4f} mm", flush=True)
    return tm


s1._tetrahedralize_with_translation_retry = injected
s1._tet_mesh_unusable = lambda tm: None

t0 = time.time()
params = main.PipelineParams()
res = main.run_pipeline(stl, params)
print(f"[tetranet] pipeline : {time.time() - t0:.1f} s", flush=True)
layers = res.layers
tp = res.toolpath.points
data = dict(
    nodes=np.asarray(res.tet_mesh.tet_vertices), phi=np.asarray(res.scalar_field.values),
    isovalues=np.array([l.isovalue for l in layers]),
    layer_ids=np.array([l.id for l in layers]),
    diag_layer=np.array([d.layer_id for d in res.step8_diagnostics]),
    contours=[[np.asarray(c, dtype=np.float32) for c in l.contour_3d] for l in layers],
    n_regions=len(res.region_dependency_graph.regions),
    pass_layer=np.array([seg[0] for p in res.schedule.passes for seg in p.segments]),
    pos=np.array([p.position for p in tp], dtype=np.float32),
    angles=np.array([p.axis_angles for p in tp], dtype=np.float32),
    extr=np.array([p.extrusion for p in tp], dtype=np.float32),
    diag_pos=np.array([d.position for d in res.step8_diagnostics], dtype=np.float32) if res.step8_diagnostics else np.zeros((0,3),np.float32),
    diag_err=np.array([d.residual_angle_error_deg for d in res.step8_diagnostics], dtype=np.float32),
    collisions=np.array([tp[e.point_index].position for e in res.collision_report.events], dtype=np.float32) if res.collision_report.events else np.zeros((0,3),np.float32),
)
with open(out, "wb") as f:
    pickle.dump(data, f)
print("[tetranet] sortie écrite :", out, flush=True)
