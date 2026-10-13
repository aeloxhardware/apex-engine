#!/bin/bash
cd /r/APEX_NEW_VERSION/slicer_engine
export APEX_USE_CPP=1
PY=.venv/Scripts/python.exe
for m in l_shape_overhang l-triple-axis-overhanged-shape death turkish; do
  b=$(.venv/Scripts/python.exe -c "import trimesh;m=trimesh.load('res/test_models/$m.stl',force='mesh');print(m.bounds.tolist())")
  echo "== $m bounds $b"
  for off in 0 60 100; do
    f=out/offset/${m}_x${off}.stl
    $PY - <<PY2
import trimesh, numpy as np
m = trimesh.load('res/test_models/$m.stl', force='mesh')
c = (m.bounds[0] + m.bounds[1]) / 2
m.apply_translation([-c[0] + $off, -c[1], 0])   # centre sur (off, 0) ; z inchange
m.export('$f')
PY2
    echo "-- $m offset $off $(date +%T)"
    $PY scripts/measure_matrix.py $f ref --collisions --out=out/offset/measure_offset.jsonl 2>&1 | tail -3
  done
done
echo DONE
