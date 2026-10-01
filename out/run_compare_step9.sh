export APEX_USE_CPP=1 APEX_CPP_LIB_DIR=cpp/build_dev2/lib
for m in l_shape_overhang turkish death; do
  echo "=== $m" ; .venv/Scripts/python.exe -u scripts/compare_step9_modes.py res/test_models/$m.stl
done
