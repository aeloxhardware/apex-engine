"""Try manifold3d boolean union as a self-intersection repair for
turkish.stl, avoiding pymeshfix's destructive collapse."""
import sys
import numpy as np
import trimesh
sys.path.insert(0, ".")
from pipeline import step1_ingest as s1

MODEL = "res/test_models/turkish.stl"
raw = trimesh.load(MODEL)
cleaned, report = s1._clean_and_repair(raw)
print(f"cleaned: watertight={cleaned.is_watertight}, faces={len(cleaned.faces)}, z=[{cleaned.vertices[:,2].min()},{cleaned.vertices[:,2].max()}]")

try:
    unioned = cleaned.union(cleaned, engine="manifold")
    print(f"union result: faces={len(unioned.faces)}, watertight={unioned.is_watertight}, z=[{unioned.vertices[:,2].min()},{unioned.vertices[:,2].max()}]")
except Exception as e:
    print(f"union FAILED: {type(e).__name__}: {e}")

# try tetgen directly on the union result if it worked
