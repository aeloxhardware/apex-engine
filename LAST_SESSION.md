# Dernière session — contexte court

L'historique des sessions (ancien `LAST_SESSION.md`, 1 390 lignes) est dans `docs/journal/` (voir `INDEX.md`).

## Session du 2026-11-07 — R-030 : résorption de la dette du moteur
- Point 3 : plus de variable d'environnement qui change le résultat, plus de `global` mutable. Options explicites
  `core/runtime.py` (`RuntimeOptions`, porté par `PipelineParams.runtime`), état de run `RuntimeState`.
  Tests : `pytest` / `pytest --native`. Contrôle : `scripts/check_hidden_state.py`.
- Point 4 : machine et tête en données (`res/machines/*.json`, `config/machine_profiles.py`) ; plus de
  `PLATFORM_RADIUS` ni de `toolhead_id` par défaut en dur.
- Point 7 : `CLAUDE.md` ramené à ~215 lignes ; limites dans `docs/limites/`, spécification dans `docs/`.
- Points 2 et 9 : règle du commutateur natif (`docs/NATIF_ET_PYTHON.md`), registre des palliatifs (`docs/PALLIATIFS.md`).
- Point 8 (TetraNet) : mesure de la coque, projection des contours (D, opt-in, +18 à +23 % de durée) et protection des
  arêtes (C, non fiable, désactivée). Seuil de 0,1 mm non tenu sur turkish (1,9 % des points). Voir `tetranet/docs/RAPPORT.md`.
- Outil de mesure corrigé : `scripts/export_slicing_preview.py` ignorait `mesher` et copiait `_layers_stage`.
- À vérifier / décider : jalon « moteur 100 % C++ » (RULES §11.3, sans date) pour le point 2, TetraNet (point 8), tri de `scripts/`.
