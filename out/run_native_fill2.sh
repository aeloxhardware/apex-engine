export APEX_USE_CPP=1
for m in death turkish l-triple-axis-overhanged-shape l_shape_overhang; do
  for v in ref native_fill; do
    echo "=== $m $v"; .venv/Scripts/python.exe -u scripts/measure_matrix.py res/test_models/$m.stl $v --collisions --out=out/measure_native_fill2.jsonl 2>&1 | grep -E '"(unsupported_fraction|collision_fraction|extruded_mm|time_real_h)"|Traceback|Error'
  done
done
