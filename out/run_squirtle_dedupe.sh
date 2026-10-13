#!/bin/bash
cd "R:/APEX_NEW_VERSION/slicer_engine/research/S4_Slicer"
export APEX_USE_CPP=1
PY=../../.venv/Scripts/python.exe
K=output_models/loop_Squirtle_pin_skin_excluded_final_pass8.json
echo "P5 $(date +%T)"; $PY -u skin_print.py Squirtle 8 _pin --head-tolerance=45 --exclude=$K 2>&1 | tail -6
echo "step89 $(date +%T)"; $PY -u full_step89.py Squirtle 8 _pin --fast --height-clamp --apply --apex-detect 2>&1 | grep -v "^\[" | tail -12
echo "export $(date +%T)"; $PY -u export_gcode.py Squirtle 8 _pin 2>&1 | grep -v "^\[" | tail -12
echo "END $(date +%T)"
