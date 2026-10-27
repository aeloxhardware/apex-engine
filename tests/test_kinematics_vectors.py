"""R-008 : les vecteurs publies se rejouent, et le fichier est a jour."""
import json
import os
import subprocess
import sys

import numpy as np

from file_io.gcode_export import firmware_to_machine, machine_to_firmware

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PATH = os.path.join(ROOT, "tests", "vectors", "kinematics_vectors.json")


def _doc():
    with open(PATH, encoding="utf-8") as f:
        return json.load(f)


def test_vectors_replay_forward_and_inverse():
    doc = _doc()
    assert len(doc["vectors"]) >= 20
    for v in doc["vectors"]:
        fw = machine_to_firmware(v["tip"], v["b_apex"], v["pivot_to_tip"], previous_y=v["previous_y"],
                                 singularity_radius=doc["singularity_radius"], tilt_sign=v["tilt_sign"],
                                 z_reference=v["z_reference"], b_sign=v["b_sign"])
        exp = v["firmware"]
        assert np.allclose(fw, [exp["X"], exp["Y"], exp["Z"], exp["B"]], atol=1e-8)
        if v["inverse_exact"]:
            tip = firmware_to_machine(exp["X"], exp["Y"], exp["Z"], exp["B"], v["pivot_to_tip"],
                                      v["tilt_sign"], v["z_reference"], v["b_sign"])
            assert np.linalg.norm(tip - np.array(v["tip"])) < doc["tolerance_mm"]


def test_vectors_file_up_to_date():
    r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "gen_kinematics_vectors.py"), "--check"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
