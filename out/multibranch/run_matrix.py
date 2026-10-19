"""Lance measure_matrix pour une liste (modele, variante), avec delai maximal, resultat dans results.jsonl.
Usage : python run_matrix.py <variante1,variante2> [modeles] [--timeout=1500]
Autonome : un echec ou un delai ne stoppe jamais le lot."""
import json
import os
import subprocess
import sys
import time

ROOT = r"R:\APEX_NEW_VERSION\slicer_engine"
PY = os.path.join(ROOT, ".venv", "Scripts", "python.exe")
OUT = r"D:\apex_multibranch"
MODELS = ["real_overhang", "dino", "benchy", "death", "fig"]


def main():
    flags = {a.split("=")[0]: a.split("=", 1)[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    variants = args[0].split(",")
    models = args[1].split(",") if len(args) > 1 else MODELS
    timeout = int(flags.get("--timeout", 1500))
    env = dict(os.environ, APEX_USE_CPP="1", PYTHONUNBUFFERED="1")
    for model in models:
        for variant in variants:
            t0 = time.time()
            log = os.path.join(OUT, "logs", f"{model}_{variant}.log")
            status = "ok"
            with open(log, "w", encoding="utf-8", errors="replace") as lf:
                try:
                    p = subprocess.Popen(
                        [PY, "-u", "scripts/measure_matrix.py", os.path.join(OUT, "models", model + ".stl"), variant,
                         "--collisions", "--out=" + os.path.join(OUT, "results.jsonl")],
                        cwd=ROOT, env=env, stdout=lf, stderr=subprocess.STDOUT)
                    try:
                        code = p.wait(timeout=timeout)
                        if code != 0:
                            status = f"echec {code}"
                    except subprocess.TimeoutExpired:
                        subprocess.run(["taskkill", "/PID", str(p.pid), "/T", "/F"], capture_output=True)
                        status = f"delai {timeout}s"
                except Exception as exc:
                    status = f"erreur {exc!r}"
            if status != "ok":
                with open(os.path.join(OUT, "results.jsonl"), "a", encoding="utf-8") as f:
                    f.write(json.dumps({"model": model, "variant": variant, "status": status,
                                        "seconds": round(time.time() - t0, 1)}) + "\n")
            print(model, variant, status, round(time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
