# Scripts maintenus (R-030 point 7, limite dure n° 9)

Un script reste dans `scripts/` s'il est un outil **maintenu** : appelé par la CI ou par un test, ou outil de
mesure / d'export documenté. `tests/test_scripts_registry.py` échoue si un script n'est pas listé ici, ou si un
script listé ne compile plus. Les diagnostics et mesures ponctuelles (54 fichiers, tri du 2026-11-07) sont dans la
branche `archives` : `git show archives:scripts/<nom>.py` (voir `docs/ARCHIVES.md`).

| Script | Rôle | Contrôlé par |
|---|---|---|
| `check_limits.py` | cliquet des limites dures (800 / 80 / 8) | CI |
| `check_hidden_state.py` | pas de `global`, d'environnement `APEX_*` ni de machine en dur | CI |
| `check_palliatives.py` | registre `docs/PALLIATIFS.md` à jour | CI |
| `gen_params.py` | schéma JSON et `params.hpp` générés | CI, `tests/test_config_schema.py` |
| `gen_licences.py` | licences des dépendances | CI |
| `gen_kinematics_vectors.py` | vecteurs de la description de cinématique | `tests/test_nonplanar.py` |
| `generate_test_model.py` | modèle de test `l_shape_overhang.stl` | tests |
| `build_cpp.py` | construit le noyau natif (`--dev` : sans toucher `cpp/lib`) | CI (équivalent cmake) |
| `export_slicing_preview.py` | tranchage seul (step1 à step8) | `tests/test_main_pipeline.py` |
| `export_model_gcode.py` | STL → G-code par le pipeline (`--native`) | outil de production |
| `measure_matrix.py` | mesure sur les 4 modèles de référence (`--native`) | définition de « terminé » |
| `profile_pipeline.py` | profil par étape (`--native`) | outil de mesure |
| `estimate_gcode_time.py` | durée estimée d'un G-code | outil de mesure |
| `preview_gcode_file.py` | contrôle aller-retour d'un G-code | outil de contrôle |
| `sync_from_ui.py` | synchronise les fichiers partagés avec l'interface | contrat |
| `limits_baseline.json` | dette de départ du cliquet | `check_limits.py` |
