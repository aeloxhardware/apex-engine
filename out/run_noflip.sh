export APEX_USE_CPP=1
for m in l_shape_overhang l-triple-axis-overhanged-shape death turkish; do
  echo "=== $m noflip"; .venv/Scripts/python.exe -u scripts/measure_matrix.py res/test_models/$m.stl noflip --out=out/measure_noflip.jsonl | grep -E '"(unsupported_fraction|ordering_defects|layers|layer_thickness_measured|extruded_mm)"'
done
