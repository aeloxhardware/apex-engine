#!/bin/bash
cd /r/APEX_NEW_VERSION/slicer_engine
export APEX_USE_CPP=1 APEX_CPP_LIB_DIR=cpp/build_dev/lib
R=research/S4_Slicer
OUT=out/bird/profile_bird_2026-10-13.jsonl
p () { name=$1; shift; echo "== $name $(date +%T)"; .venv/Scripts/python.exe -u scripts/profile_bird.py "$@" --out=$OUT > "out/bird/profile_$name.log" 2>&1; echo "exit $? $(date +%T)"; }
p ltriple res/test_models/l-triple-axis-overhanged-shape.stl
p benchy "$R/input_models/tilted benchy undertriangulated_repaired.stl"
p dino "$R/input_models/dino.stl" --from-pass="$R/output_models/loop_dino_pin_,5"
p death res/test_models/death.stl
echo DONE
