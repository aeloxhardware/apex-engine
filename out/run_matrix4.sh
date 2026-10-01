#!/bin/bash
cd "R:/APEX_NEW_VERSION/slicer_alg"
export APEX_USE_CPP=1
export APEX_CPP_LIB_DIR=cpp/build_dev/lib
for model in l_shape_overhang.stl l-triple-axis-overhanged-shape.stl death.stl turkish.stl; do
  for v in ref radial ppn_radial ppn_radial_median9 ppn_radial_nodal ppn_radial_nodal_median9; do
    echo "=== $(date +%T) $model $v"
    timeout 3600 .venv/Scripts/python.exe -u scripts/measure_matrix.py res/test_models/$model $v --out=out/measure_matrix4.jsonl > "out/matrix4_${model%.stl}_$v.log" 2>&1
    echo "    code $? $(date +%T)"
  done
done
echo FIN
