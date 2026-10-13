#!/bin/bash
cd /r/APEX_NEW_VERSION/slicer_engine
export APEX_USE_CPP=1
export PYTHONPATH="C:/Users/test_/AppData/Local/Temp/claude/r--APEX-NEW-VERSION-slicer-engine/52b9fafc-86c0-4106-8a19-7f2e8d68bedf/scratchpad"
PY=.venv/Scripts/python.exe
echo "== l_shape"; $PY -u scripts/export_model_gcode.py res/test_models/l_shape_overhang.stl --out=out/test_l_shape.gcode 2>&1 | grep -v "^\[preview\]"
cd research/S4_Slicer
echo "== benchy"; ../../$PY -u export_gcode.py "tilted benchy undertriangulated_repaired" 4 _pin --strict 2>&1 | grep -v "^\["
echo "== dino"; ../../$PY -u export_gcode.py dino 5 _pin --strict 2>&1 | grep -v "^\["
echo FIN
