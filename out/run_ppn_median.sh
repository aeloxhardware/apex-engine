#!/bin/bash
cd /r/APEX_NEW_VERSION/slicer_engine
export APEX_USE_CPP=1
for m in l_shape_overhang l-triple-axis-overhanged-shape death turkish; do
  for v in ppn ppn_median9; do
    echo "== $m $v $(date +%T)"
    .venv/Scripts/python.exe -u scripts/measure_matrix.py res/test_models/$m.stl $v --out=out/measure_ppn_b360_2026-10-06.jsonl
  done
done
echo DONE
