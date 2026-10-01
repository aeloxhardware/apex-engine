export APEX_USE_CPP=1
for m in death turkish l-triple-axis-overhanged-shape l_shape_overhang; do
  for v in ref corder; do
    echo "=== $m $v"; .venv/Scripts/python.exe -u scripts/measure_matrix.py res/test_models/$m.stl $v --collisions --out=out/measure_corder.jsonl 2>&1 | grep -E 'ordre contraint|"(unsupported_fraction|ordering_defects|collision_fraction|collision_platform|seconds)"|Traceback|Error'
  done
done
