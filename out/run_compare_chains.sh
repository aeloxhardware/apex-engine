#!/bin/bash
cd /r/APEX_NEW_VERSION/slicer_engine
export APEX_USE_CPP=1 APEX_CPP_LIB_DIR=cpp/build_dev/lib
for m in "tilted benchy undertriangulated_repaired" "dino"; do
  echo "== $m $(date +%T)"
  .venv/Scripts/python.exe -u scripts/measure_bird_vs_fem.py "research/S4_Slicer/input_models/$m.stl" --ways=bird \
     --out=out/bird/compare_chains_2026-10-13.jsonl --keep-gcode=/h/apex_bird_compare_2026-10-13 > "out/bird/compare_${m// /_}.log" 2>&1
  echo "exit $? $(date +%T)"
done
echo DONE
