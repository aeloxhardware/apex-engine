#!/bin/bash
cd /r/APEX_NEW_VERSION/slicer_engine
export APEX_CPP_LIB_DIR=cpp/build_vk2/lib
for cfg in "0 0" "1 0" "1 1"; do
  set -- $cfg
  echo "== APEX_USE_CPP=$1 APEX_USE_GPU=$2 $(date +%T)"
  APEX_USE_CPP=$1 APEX_USE_GPU=$2 .venv/Scripts/python.exe -m pytest tests -q -p no:cacheprovider > out/suite_cpp$1_gpu$2.log 2>&1
  echo "exit $?"; tail -3 out/suite_cpp$1_gpu$2.log
done
echo DONE
