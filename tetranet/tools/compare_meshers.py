"""Tranche une pièce avec un mailleur donné et écrit une ligne JSON de mesures.
Usage : python compare_meshers.py <stl> <tetgen|tetranet> [resolution]
Un processus par mesure : aucun cache partagé entre les deux mailleurs."""
import json, os, sys, time
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
import main
from pipeline.step9b_support_check import SupportCheckParams

stl, mesher = sys.argv[1], sys.argv[2]
resolution = float(sys.argv[3]) if len(sys.argv) > 3 else 1.5
params = main.PipelineParams(resolution=resolution, mesher=mesher, support_check=SupportCheckParams())
t0 = time.time()
row = dict(model=os.path.basename(stl), mesher=mesher, resolution=resolution)
try:
    r = main.run_pipeline(stl, params)
except Exception as e:
    row.update(error=f"{type(e).__name__}: {str(e)[:200]}", seconds=round(time.time() - t0, 1))
    print(json.dumps(row)); sys.exit(0)
phi, nodes = np.asarray(r.scalar_field.values), np.asarray(r.tet_mesh.tet_vertices)
zc = []
for l in r.layers:
    if l.contour_3d:
        P = np.vstack([np.asarray(c) for c in l.contour_3d]); zc.append((P[:, 2].mean(), np.ptp(P[:, 2])))
zc = np.array(zc) if zc else np.zeros((0, 2))
dz = np.diff(zc[:, 0]) if len(zc) > 1 else np.zeros(0)
n = len(r.toolpath.points)
row.update(
    seconds=round(time.time() - t0, 1), tets=int(r.tet_mesh.tet_count), layers=len(r.layers), points=n,
    phi_max=round(float(phi.max()), 2), phi_z_corr=round(float(np.corrcoef(phi, nodes[:, 2])[0, 1]), 3),
    layer_z_extent_max=round(float(zc[:, 1].max()), 2) if len(zc) else 0.0, z_drops_over_half_mm=int((dz < -0.5).sum()),
    unreachable=len(r.step8_diagnostics), unreachable_pct=round(100 * len(r.step8_diagnostics) / max(n, 1), 1),
    collisions=len(r.collision_report.events), collisions_unresolved=int(r.adjustment_report.unresolved_collisions) if isinstance(r.adjustment_report.unresolved_collisions, int) else len(r.adjustment_report.unresolved_collisions),
    unsupported_pct=round(100 * r.support_check_report.unsupported_fraction, 3),
    ordering_defects=int(r.support_check_report.ordering_defect_count), mesher_report=r.ingest_report.mesher if hasattr(r, "ingest_report") else None,
)
print(json.dumps(row))
