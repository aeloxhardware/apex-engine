cd "R:/APEX_NEW_VERSION/slicer_engine/research/S4_Slicer"
export APEX_USE_CPP=1
# --apex-detect (limite #70) : detection de collision par le step9 rapide natif, x10, trajectoire identique au bit
PY=../../.venv/Scripts/python.exe
run_model () {
  M="$1"; K="$2"
  echo "##### $M passe $K  debut $(date +%H:%M:%S)"
  t0=$(date +%s)
  if ! $PY -u skin_loop.py "$M" $K _pin --first-layer=concentric_conform; then
    echo "  ECHEC de la boucle de lissage : modele abandonne"; return 1
  fi
  echo "  boucle de lissage terminee en $(( $(date +%s) - t0 )) s"
  t1=$(date +%s)
  $PY -u full_step89.py "$M" $K _pin --fast --height-clamp --apply --apex-detect 2>&1 | grep -v "^\[" | tail -25
  echo "  step8+9 en $(( $(date +%s) - t1 )) s"
  t2=$(date +%s)
  $PY -u export_gcode.py "$M" $K _pin 2>&1 | grep -v "^\[" | tail -25
  echo "  export en $(( $(date +%s) - t2 )) s ; TOTAL $(( $(date +%s) - t0 )) s"
}
run_model "tilted benchy undertriangulated_repaired" 4
run_model "dino" 5
