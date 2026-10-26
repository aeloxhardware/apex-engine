"""Phase 0 : mesure TetGen (via step1 du moteur) sur le corpus. Une ligne JSON par pièce.

Usage : python tetranet/tools/baseline_tetgen.py <resolution> <stl> [<stl> ...]
"""
import json, os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from pipeline import step1_ingest

res = float(sys.argv[1])
for path in sys.argv[2:]:
    t0 = time.time()
    try:
        md, tm, rep = step1_ingest.run(path, res, allow_multi_mesh=False, place_on_platform=True,
                                       check_surface_coverage=True)
        out = dict(model=os.path.basename(path), resolution=res, seconds=round(time.time() - t0, 2),
                   faces=int(len(md.faces)), nodes=int(tm.node_count), tets=int(tm.tet_count),
                   report={k: (v if isinstance(v, (int, float, str, bool, type(None))) else str(v))
                           for k, v in vars(rep).items()})
    except Exception as e:
        out = dict(model=os.path.basename(path), resolution=res, seconds=round(time.time() - t0, 2),
                   error=f"{type(e).__name__}: {str(e)[:200]}")
    print(json.dumps(out), flush=True)
