---
name: safe-long-runs
description: Heavy research computations crashed the user's PC (47 GB RAM) — kill Python children explicitly, batch big queries, test one layer first
metadata:
  type: feedback
---

On 2026-09-27 a carving run (uniform subdivision + one giant trimesh closest_point call per layer) saturated the user's 47 GB of RAM and crashed the computer. TaskStop had killed only the bash shell; the Python child kept running unseen.

**Why:** on this Windows machine, stopping a background task does not stop the Python process it spawned; memory-heavy numpy/trimesh calls grow with points × candidates.
**How to apply:** after stopping a task, find and kill leftover python.exe by CommandLine (Get-CimInstance Win32_Process + Stop-Process); batch large distance/proximity queries (e.g. ≤ 50k points per call); try one layer or one model part first and measure time and memory before a full run; check free memory before launching parallel heavy jobs. Related: [[bird-automation-goal]].
