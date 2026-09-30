---
name: research-code-stays-local
description: User decision — research/S4_Slicer (R&D code) stays local and unversioned; do not offer again to move it into git
metadata:
  node_type: memory
  type: project
  originSessionId: df56fa55-73de-4689-8817-bf5ea9359353
  modified: 2026-09-26T20:09:39.348Z
---

Decided by the user on 2026-09-26: the research and development code in `research/S4_Slicer/` (bird_loop, fold_judge, bird_slice_deformed, etc.) stays local, outside git (the folder is in .gitignore). Only fixes to Apex itself (pipeline/, tests/, CLAUDE.md, LAST_SESSION.md) get committed.

**Why:** the user's own choice. They were offered to copy the key scripts into a versioned folder and declined.
**How to apply:** don't offer again to version these scripts, and don't `git add` anything under research/. When a research result becomes part of Apex, port it into pipeline/ with tests, as was done for step6 (limite #57). Related: [[bird-automation-goal]].
