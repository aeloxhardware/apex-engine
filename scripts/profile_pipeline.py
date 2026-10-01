"""
Profil cProfile du pipeline step1->step8 sur un modele, trie par temps propre,
avec un recapitulatif par etape (temps cumule des fonctions `run` de chaque
pipeline/stepN_*.py). Sert a choisir le prochain pas du portage C++ sur une
MESURE, jamais sur une supposition (limite #56).

Usage : APEX_USE_CPP=1 .venv/Scripts/python.exe -u scripts/profile_pipeline.py <modele.stl>
        [--resolution=0.8] [--top=40] [--out=out/profile.prof]
"""
import cProfile
import os
import pstats
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from main import PipelineParams  # noqa: E402
from scripts.export_slicing_preview import run_slicing_only  # noqa: E402


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    model = args[0] if args else os.path.join(ROOT, "res", "test_models", "death.stl")
    params = PipelineParams(resolution=float(flags.get("--resolution", 0.8)))
    out = flags.get("--out", os.path.join(ROOT, "out", "profile.prof"))
    prof = cProfile.Profile()
    t0 = time.time()
    prof.enable()
    run_slicing_only(model, params)
    prof.disable()
    wall = time.time() - t0
    prof.dump_stats(out)
    st = pstats.Stats(prof)
    print(f"\nmur {wall:.1f} s (profile)")
    print("\n== par etape (cumul de run) ==")
    for key, (_cc, _nc, _tt, ct, _callers) in sorted(st.stats.items(), key=lambda kv: kv[0][0]):
        fname, _line, func = key
        base = os.path.basename(fname)
        if func == "run" and base.startswith("step"):
            print(f"  {base:32} {ct:8.2f} s")
    print(f"\n== top {flags.get('--top', 40)} temps propre ==")
    st.sort_stats("tottime").print_stats(int(flags.get("--top", 40)))


if __name__ == "__main__":
    main()
