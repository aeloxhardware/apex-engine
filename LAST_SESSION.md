# Dernière session — contexte court

L'historique des sessions (ancien `LAST_SESSION.md`, 1 390 lignes) est dans `docs/journal/` (voir `INDEX.md`).

## Session du 2026-11-13 — mesure du pipeline complet sur une hélice, durée de la tétraédralisation
Aucun code du dépôt modifié : mesures seulement (scripts et STL laissés hors dépôt, voir « Non conservé »).
- **Modèle** : hélice verticale générée à la demande (moyeu R 8 × 30 mm, 3 pales vrillées 45° → 15°, rayon 52 mm), STL étanche,
  ~30 000 faces. Résolution 1,0 mm, couches 0,2 mm. Un seul modèle synthétique, un seul run par configuration.
- **Durée du pipeline complet** : Python + TetGen 423,8 s ; natif + TetGen 57,5 s ; natif + TetraNet 17,3 s.
  Étapes (Python → natif+TetraNet) : tétraédralisation 69,7 → 2,2 s ; collision 199,4 → 6,3 s ; trajectoire machine 77,8 → 2,1 s.
- **Le noyau natif n'existe pas sur Linux** : `scripts/build_cpp.py` ne gère que Visual Studio. Build direct réussi avec
  `cmake -S cpp -G Ninja -DAPEX_LIB_DIR=<dossier>` (pybind11 du venv, Vulkan désactivé), module chargé par `APEX_CPP_LIB_DIR`,
  `tetranet_mesh` présent. Rien de Vulkan n'est utilisé par une étape : seuls `saxpy_f32/f64.comp` existent.
- **Cause de la lenteur de TetGen** (`element_ardor/step1_tetgen_run.py`, `pq1.2a<résolution³>`, sans `Y`) : TetGen redécoupe la
  surface (35 765 points ajoutés sur les facettes, surface de 28 128 à 125 436 triangles) ; le plafond `résolution³` ne mord pas
  (tétraèdres moyens ~0,016 mm³ pour un plafond de 1 mm³), c'est la finesse de la surface qui fixe le nombre de tétraèdres.
  Sur l'hélice d'origine : 4 passes TetGen (3 échecs `split_subface`, puis 23 s). Avec `Y` (rustine à l'exécution, non versionnée) :
  1 passe, 1 272 736 → 164 915 tétraèdres, ingestion 30,0 → 6,8 s, pipeline 57,5 → 35,9 s. Aucune décision documentée contre `Y`.
- **La sortie change** selon le mailleur (couches : 235 TetGen, 303 TetGen+`Y`, 304 TetraNet ; orientations inatteignables :
  27,6 %, 34,5 %, 40 %). La qualité des couches n'a PAS été mesurée : on ne sait pas si c'est mieux ou moins bien.
- **Défaut relevé** : `generalized_winding_number` (`element_ardor/_self_intersection_geometry.py:166`), chemin Python, alloue un
  tableau points × triangles sans découpage par blocs (17,8 Gio demandés sur un STL de 28 128 faces, ingestion plantée).
  Le chemin natif ne l'a pas. Non corrigé, non enregistré dans `docs/limites/`.
- **Aperçu des couches** : page 3D (cordons cylindriques, curseur de couches, plan de coupe) publiée comme artefact claude.ai
  (privé). Trajets d'extrusion de l'étape 6 relevés en 3D (`_lift_to_3d`), sans déplacements à vide ni ajustement de collision.
  Vérifié avec 110 couches seulement (pas de GPU ici) ; le rendu des 304 couches n'est pas testé.
- **Non fait** : flèche de normales (ligne de courant de ∇φ, affichable ; la modifier demande les points d'attache sur la voie
  chaleur, à refaire) ; tests `pytest` et `pytest --native` ; mesure sur les 4 modèles de référence ; comparaison de qualité
  TetGen / TetraNet ; vérification du mandrin (`chuck_checked=False`).
- **Non conservé** : scripts de génération et de mesure (générateur d'hélice, lanceur chronométré par étape, extraction des
  trajectoires, page de l'aperçu) restés dans le scratchpad de la session, hors dépôt (`scripts/` est trié, R-030 point 9).
- **À décider** : option explicite `Y` pour TetGen (désactivée par défaut, mesurée sur les 4 modèles, `APEX_CONTRACT_VERSION`
  incrémenté si le G-code change) ; script de build Linux du noyau natif ; découpage par blocs de `generalized_winding_number`.

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
