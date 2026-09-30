---
name: maker-faire-rome-deadline
description: Maker Faire Rome end of October 2026 — show printed parts + live printing; G-code export and collision-checked toolpaths come first
metadata:
  node_type: memory
  type: project
  originSessionId: df56fa55-73de-4689-8817-bf5ea9359353
  modified: 2026-09-27T15:16:35.993Z
---

Maker Faire Rome is less than a month after 2026-09-27 (end of October 2026). The user wants to show **pre-printed parts** and **live printing on the stand** — no live slicing.

**Why:** the demo needs a real print, so real G-code and a toolpath verified against the real toolhead (a collision in front of the public is the worst case). Compute speed is a commercial requirement (the user called 10 min per phase unacceptable) but not a prerequisite for the fair.

**How to apply:** prioritize (1) checking/writing the G-code export — status unknown, first read the untracked scripts `scripts/gen_*_gcode.py`, `preview_gcode.py`; (2) step-9 collision check with aetheris-4.glb. Defer P4's 10 surface patterns, multi-model validation, C++ port of the Bird path and ironing until after the fair if time runs short. Roadmap written in CLAUDE.md "ROADMAP — Maker Faire Rome". Related: [[bird-automation-goal]], [[safe-long-runs]].
