#!/bin/bash
cd "R:/APEX_NEW_VERSION/slicer_alg"
export APEX_USE_CPP=1
for model in l_shape_overhang.stl l-triple-axis-overhanged-shape.stl death.stl turkish.stl; do
  for v in ref median9 feedback3 feedback_dyn branch005 geo_smooth50; do
    echo "=== $(date +%T) $model $v"
    timeout 3600 .venv/Scripts/python.exe -u scripts/measure_matrix.py res/test_models/$model $v --out=out/measure_matrix.jsonl > "out/matrix_${model%.stl}_$v.log" 2>&1
    echo "    code $? $(date +%T)"
  done
done
echo FIN
