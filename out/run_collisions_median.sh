#!/bin/bash
cd /r/APEX_NEW_VERSION/slicer_engine
export APEX_USE_CPP=1
for m in l_shape_overhang l-triple-axis-overhanged-shape death turkish; do
  for w in 0 9; do
    echo "== $m median=$w $(date +%T)"
    .venv/Scripts/python.exe -u scripts/measure_step9.py res/test_models/$m.stl --median=$w 2>&1 | grep -v "^\[preview\]"
  done
done
echo DONE
