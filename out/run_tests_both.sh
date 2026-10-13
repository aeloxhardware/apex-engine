#!/bin/bash
cd /r/APEX_NEW_VERSION/slicer_engine
for m in 0 1; do APEX_USE_CPP=$m .venv/Scripts/python.exe -m pytest -q -p no:cacheprovider 2>&1 | tail -5; done
echo END
