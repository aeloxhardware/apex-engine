# Commutateur natif / Python : règle, jalon, état (R-030 points 2 et 3)

## État au 2026-11-07
- Le choix natif / Python n'est plus lu dans l'environnement du processus : il est porté par `core/runtime.py`
  (`RuntimeOptions.native`, `PipelineParams.runtime`, JSON de la requête `params.runtime`). Tests : `pytest` (Python) et
  `pytest --native` ; CI : matrice des deux positions (`.github/workflows/ci.yml`).
- Variables retirées de l'environnement : `APEX_USE_CPP`, `APEX_USE_GPU`, `APEX_CPP_GRID_MODE`, `APEX_TETGEN_WORKER`,
  `APEX_TETGEN_LOG`. `APEX_NATIVE_FILL`, `APEX_NATIVE_LAYER`, `APEX_NATIVE_GRID` étaient déjà des paramètres.
- Reste dans l'environnement : `APEX_CPP_LIB_DIR` (emplacement du module natif lu à l'import ; il désigne un autre binaire,
  il ne change pas le résultat) et `APEX_UI_DIR` (outil `scripts/sync_from_ui.py`, hors moteur).
- Aucun `global` mutable dans le moteur : compteurs de repli, cache du contrôleur de collision et progression vivent dans
  `RuntimeState` (un état par run, `runtime.use(options)`).

## Règle et jalon (décision de l'utilisateur, `apex/integration/RULES.md` §11.3, 2026-11-04)
**Aucune date.** Objectif B : le Python est supprimé. Le jalon est atteint quand **le moteur est 100 % C++** (orchestrateur,
ingestion, toutes les étapes, Element Ardor et Element Aeolis, TetraNet, export, générateur du contrat) **et que le code est
propre, sans palliatif**. D'ici là la double implémentation reste ; chaque pièce portée a son test d'équivalence avant que le
Python correspondant soit supprimé. **R-030 point 2 reste ouvert jusqu'à ce jalon.**
- Replis point par point (reprise Python de l'IK, `step8_machine_toolpath._resolve_all_axis_angles_native`) : la grille est
  déjà portée (`GridSolver`, `Step8Params.native_grid_fallback = True`) ; la reprise Python n'est atteinte qu'avec ce
  paramètre à `False` (oracle). Elle disparaît avec le Python.
- Tout changement de défaut qui modifie le G-code incrémente `APEX_CONTRACT_VERSION` dans le même commit.
- Une « date de bascule » proposée le 2026-11-07 puis confirmée a été retirée le même jour : elle contredisait §11.3.
