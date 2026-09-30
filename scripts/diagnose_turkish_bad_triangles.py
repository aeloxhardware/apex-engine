"""Identify the exact self-intersecting triangles TetGen refuses on
turkish.stl, using TetGen's own detection (switches='d') - the only
reliable judge per CLAUDE.md limit #43."""
import sys
import numpy as np
import tetgen
import trimesh
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1

MODEL = "res/test_models/turkish.stl"
raw = trimesh.load(MODEL)
cleaned, report = s1._clean_and_repair(raw)
print(f"cleaned: faces={len(cleaned.faces)}, watertight={cleaned.is_watertight}")

tet = tetgen.TetGen(cleaned.vertices, cleaned.faces)
try:
    tet.tetrahedralize(switches="d")
except Exception as e:
    print(f"detection mode raised: {type(e).__name__}: {e}")
print("done with detection pass")
