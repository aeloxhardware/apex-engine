#!/bin/bash
cd /r/APEX_NEW_VERSION/slicer_engine
for v in 0 1; do
  echo "== APEX_USE_CPP=$v $(date +%T)"
  APEX_USE_CPP=$v .venv/Scripts/python.exe -m pytest tests -q -x -p no:cacheprovider > out/suite_cpp$v.log 2>&1
  echo "exit $?"; tail -5 out/suite_cpp$v.log
done
echo DONE
