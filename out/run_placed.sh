#!/bin/bash
cd /r/APEX_NEW_VERSION/slicer_engine
export APEX_USE_CPP=1
OUT=out/measure_placed_2026-10-13.jsonl
for m in l_shape_overhang l-triple-axis-overhanged-shape death turkish; do
  for v in ref placed; do
    echo "== $m $v $(date +%T)"
    .venv/Scripts/python.exe -u scripts/measure_matrix.py res/test_models/$m.stl $v --collisions --out=$OUT > out/placed_${m}_${v}.log 2>&1
    echo "exit $?"
  done
done
echo DONE
