---
name: disk-r-small
description: "R: (project drive) is only 8 GB — heavy outputs fill it; archive to H: (489 GB free), never build big repos on R:"
metadata:
  node_type: memory
  type: project
  originSessionId: df56fa55-73de-4689-8817-bf5ea9359353
  modified: 2026-09-30T10:02:52.102Z
---

The project lives on R:, an 8 GB volume (measured 2026-09-30: ~1.7 GB free). A full-snapshot git repo attempt filled it to 0 bytes mid-session. Each research model chain writes ~150-300 MB into `research/S4_Slicer/output_models/`.

**Why:** a full disk silently breaks background computations and git operations.

**How to apply:** check `df -h /r` before long chains; put large temporary repos/builds on H: (large, free); obsolete research outputs were archived (moved, not deleted) to `H:/apex_archive_2026-09-30/` with `LISTE_ARCHIVE.txt` — the user has not yet confirmed deletion. Related: [[safe-long-runs]], [[research-code-stays-local]].
