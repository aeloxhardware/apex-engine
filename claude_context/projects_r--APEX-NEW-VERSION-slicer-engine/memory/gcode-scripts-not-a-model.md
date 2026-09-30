---
name: gcode-scripts-not-a-model
description: The untracked scripts/gen_*_gcode.py must NOT inspire file_io/gcode_export.py — they generate fixed shapes only
metadata:
  type: feedback
---

When writing Apex's G-code export (`file_io/gcode_export.py`), do not base it on `scripts/gen_l_overhang_gcode.py`, `gen_elbow_tube_gcode.py`, `preview_gcode.py` or `export_machine_animation.py`.

**Why:** the user (2026-09-29) said these scripts were built to generate and export FIXED shapes; their structure is not a model for a general exporter.

**How to apply:** design the exporter from `ToolpathMachine` / `MachineProfile` (generic N axes) and the firmware command format; reading those scripts is at most for firmware facts (axis names, units), never for code structure. Related: [[machine-kinematics]], [[maker-faire-rome-deadline]].
