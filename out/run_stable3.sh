export APEX_USE_CPP=1
for m in death turkish l-triple-axis-overhanged-shape l_shape_overhang; do
  echo "=== $m stable3"; .venv/Scripts/python.exe -u scripts/measure_matrix.py res/test_models/$m.stl stable3 --collisions --out=out/measure_stable3.jsonl 2>&1 | grep -E 'ordre contraint|"(layers|unsupported_fraction|collision_fraction|time_real_h)"|Traceback|Error'
done
