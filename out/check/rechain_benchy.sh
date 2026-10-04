cd "R:/APEX_NEW_VERSION/slicer_alg/research/S4_Slicer"
export APEX_USE_CPP=1
PY=../../.venv/Scripts/python.exe
M="tilted benchy undertriangulated_repaired"
t0=$(date +%s)
echo "== P3 (cut-factor 2.0 par defaut)"; $PY -u bird_slice_deformed.py "$M" 4 _pin --skin --first-layer=concentric_conform --skin-exclude=output_models/loop_tilted_benchy_undertriangulated_repaired_pin_skin_excluded_final_pass4.json 2>&1 | grep "retour :\|sortie"
echo "== P5 final (verification collisions du lissage)"; $PY -u skin_print.py "$M" 4 _pin --head-tolerance=45 --exclude=output_models/loop_tilted_benchy_undertriangulated_repaired_pin_skin_excluded_final_pass4.json 2>&1 | tail -6
echo "== step8+9"; $PY -u full_step89.py "$M" 4 _pin --fast --height-clamp --apply 2>&1 | grep -v "^\[" | tail -14
echo "== export"; $PY -u export_gcode.py "$M" 4 _pin 2>&1 | grep -v "^\[" | tail -8
echo "TOTAL $(( $(date +%s) - t0 )) s"
