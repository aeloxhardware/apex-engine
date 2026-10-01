export APEX_USE_CPP=1 APEX_CPP_LIB_DIR=cpp/build_dev3/lib
for m in l_shape_overhang l-triple-axis-overhanged-shape death turkish; do
  for v in continuous_native; do
    echo "=== $m $v"; .venv/Scripts/python.exe -u scripts/measure_matrix.py res/test_models/$m.stl $v --out=out/measure_continuous.jsonl | grep -E '"(unsupported_fraction|extruded_mm|ordering_defects|layers|seconds)"'
  done
done
