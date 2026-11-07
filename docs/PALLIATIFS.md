# Registre des palliatifs (limite dure n° 8)

Un palliatif porte un identifiant et un test qui le caractérise ; sa revue a lieu **à chaque release** (RULES §11.3, plus de
date). Un palliatif sans test fait échouer `scripts/check_palliatives.py`. Un palliatif encore présent au jalon « moteur
100 % C++ » est un défaut du jalon.

| Id | Palliatif | Limite | Test qui le caractérise | Revue | Remplacement visé |
|---|---|---|---|---|---|
| PAL-001 | Lift-and-traverse des trajets à vide (`step8b_travel_routing`) | #49 | `tests/test_step8b_travel_routing.py` | à chaque release | Routage 3D (A* / RRT) dans l'espace libre autour de la pièce imprimée |
| PAL-002 | Reprise Python de l'IK au point que la résolution analytique native ne traite pas | #56 (pas 2) | `tests/test_cpp_equivalence_stages.py::test_native_grid_fallback_matches_python_grid_search` | à chaque release | Supprimée avec le Python (`docs/NATIF_ET_PYTHON.md`) |
| PAL-003 | Repli du tracé de bord natif sur `trimesh.outline()` (bord non-variété) | #56 (pas 14) | `tests/test_cpp_equivalence_mesh.py` (compteur `outline_fallbacks`) | à chaque release | Tracé natif des bords non-variété |
| PAL-004 | Nouvel essai de tétraédrisation sur maillage translaté (fragilité d'ingestion) | #72 | `tests/test_step1_ingest_repair.py::test_translation_retry_recovers_and_moves_nodes_back` | à chaque release | Réparation robuste des auto-intersections (cause : arrondi float32) |
| PAL-005 | Échecs connus PRE-001 exclus de la CI bloquante | — | `tests/known_failures.txt` (job `echecs-connus`) | à chaque release | Cause racine du lissage de surface sur Linux |
