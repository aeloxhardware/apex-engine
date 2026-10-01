#!/bin/bash
cd "R:/APEX_NEW_VERSION/slicer_alg"
export APEX_USE_CPP=1
# feedback3 a epaisseur REELLE egalisee sur la reference (facteur = ep_ref / ep_feedback3, matrice 1)
run() { echo "=== $(date +%T) $1 $2 lh=$3"; timeout 3600 .venv/Scripts/python.exe -u scripts/measure_matrix.py res/test_models/$1 $2 --layer-height=$3 --out=out/measure_matrix3.jsonl > "out/matrix3_${1%.stl}_$2_$3.log" 2>&1; echo "    code $? $(date +%T)"; }
run death.stl feedback3 0.1816
run l_shape_overhang.stl feedback3 0.1984
run l-triple-axis-overhanged-shape.stl feedback3 0.2008
run turkish.stl feedback3 0.2117
echo FIN
