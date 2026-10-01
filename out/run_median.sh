export APEX_USE_CPP=1
for m in l_shape_overhang l-triple-axis-overhanged-shape death turkish; do
  for v in median9; do
    echo "=== $m $v"; .venv/Scripts/python.exe -u scripts/measure_matrix.py res/test_models/$m.stl $v --collisions --out=out/measure_native_fill2.jsonl 2>&1 | grep -E '"(unsupported_fraction|collision_fraction|time_real_h|b_jumps_gt10|dev_median_deg|dev_over45_frac)"|Traceback|Error'
  done
done
