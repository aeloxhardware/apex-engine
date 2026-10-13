#!/bin/bash
cd /r/APEX_NEW_VERSION/slicer_engine
export APEX_USE_CPP=1 APEX_CPP_LIB_DIR=cpp/build_dev/lib
R=research/S4_Slicer
run () {
  echo "== $1 passe $2 $(date +%T)"
  .venv/Scripts/python.exe -u scripts/measure_bird_vs_fem.py "$R/input_models/$1.stl" --ways=bird \
     --from-pass="$R/output_models/loop_${3}_pin_,$2" \
     --out=out/bird/compare_chains_fixed_2026-10-13.jsonl --keep-gcode=/h/apex_bird_compare_2026-10-13/fixed \
     > "out/bird/compare_fixed_${3}.log" 2>&1
  echo "exit $? $(date +%T)"
}
mkdir -p /h/apex_bird_compare_2026-10-13/fixed
run "tilted benchy undertriangulated_repaired" 4 tilted_benchy_undertriangulated_repaired
run "dino" 5 dino
echo DONE
