# Session du 2026-09-27 — base figée, lissage non-planaire (recherche), étape 8 corrigée

Suite de la session du 2026-09-26 (méthode Bird automatisée, dont le compte
rendu est dans le commit `2b7edc9`). **Plan approuvé en cours d'exécution :**
`C:\Users\test_\.claude\plans\snug-yawning-lovelace.md` (phases P0 à P6).

CODE SECRET : 7489

## >>> SESSION DU 2026-10-13 (chantier « Bird en C++ puis Vulkan ») — LIRE EN TOUT PREMIER <<<

**Plan approuvé** : `C:\Users\test_\.claude\plans\port-e-tout-bird-hidden-nygaard.md` (étape −1,
phases 0 à 3). Décisions de l'utilisateur : portée = TOUT Bird ; matériel = RTX 3090 seulement ;
précision GPU (fp32/fp64) décidée APRÈS mesure ; démarrage immédiat. Détail technique et
chiffres : **CLAUDE.md #76 et #77** ; suivi : **DT.md**. **EN PAUSE à la demande de l'utilisateur
(avant un /compact)** : reprendre au point « Suite » ci-dessous quand il le dit.

### Git
- Branche `session-2026-10-13` (créée depuis `main`), un commit `fc473c2d` (travail du 06 au
  13/10 + `place_on_platform`). **Tout ce qui suit n'est PAS commité.** Rien n'est poussé.
- **Attention** : après ce commit, l'index contient des artefacts INDEXÉS qui ne viennent pas de
  cette session (`build/.cmake/...`, `build/CMakeFiles/...`, `cpp/build/...`, `.pytest_cache`,
  `_skipped.face/.node`), probablement une autre session (`slicer_ux`). Ne PAS les commiter sans
  vérifier : ajouter les fichiers source un par un, comme pour `fc473c2d`.

### Fait aujourd'hui (tout mesuré)
1. **Liste « moins d'une heure »** : commit ; `place_on_platform` (défaut actif, #76 : death était à
   z = 51,4 mm, turkish à 30,3 mm) ; G-code régénérés (`out/share/*_v6.zip`) ; suite complète.
2. **Étape −1 FAITE** : production Bird basculée sur `main.py --method bird` (#77).
   - `material_travel_check` (contrôle des dégagements) dans le moteur et dans l'export NATIF
     (`TravelCheck`, `PastMaterialIndex::any_past_within`, `apex_core.EXPORT_TRAVEL_CHECK`),
     texte identique ; `main.export_gcode_checked(check_travel=True)`.
   - `pipeline.bird.loop.from_saved_passes` + `run_bird_pipeline(loop_result=)` + option
     `main.py --from-pass=<préfixe>,<k>`.
   - À entrée identique (passes 4 / 5 de la recherche) : mêmes collisions, zones, Y, filament ;
     Benchy 580 → 250 s, dino 1 067 → 383 s. **Production : `out/share/benchy_pass4_v7.zip`,
     `dino_pass5_v7.zip`** (sources sur `H:/apex_bird_compare_2026-10-13/fixed/`).
3. **Phase 0 FAITE** : `scripts/profile_bird.py` (options `--from-pass`, `--gcode`, `--no-skin`),
   compteurs natifs `apex_core.bird_lsq_last_profile()`. Profils dans
   `out/bird/profile_bird_2026-10-13.jsonl` et `out/bird/profile_*.log` : Benchy 300,6 s,
   dino 373 s (depuis la passe 5), death 714 s, l-triple ~40 s, l_shape 29 s.
4. **Gains identiques au bit** (G-code du Benchy comparé OCTET PAR OCTET, avec et sans boucle) :
   - P5 parallèle `SkinCheckIndex.check_zone_parallel` (rang d'insertion = préfixe visible) :
     47,7 → 8,4 s ;
   - solveur LSQ parallèle (pool persistant `cpp/include/apex/internal/parallel_pool.hpp`,
     transposée CSR, `bird_lsq_set_threads(threads, min_parallel)`) : déformation 44,6 → 25,8 s ;
   - densification de `material_travel_check` vectorisée : 9,7 → 0,5 s ;
   - candidats barycentriques (`workers=-1`, lots de 200 000 points toutes couches, fils) :
     20,3 → 16,4 s seulement (numpy partage mal la résolution 3×3) ;
   - **Benchy complet 300,6 → 226,8 s ; depuis la passe 4 : 166,9 s.**
5. **Phase 1 (socle Vulkan) FAITE** : `cpp/CMakeLists.txt` option `APEX_WITH_VULKAN`,
   `scripts/build_cpp.py --vulkan`, volk + VMA dans `cpp/third_party/` (licences MIT, fichiers
   `VERSION.txt` — **jamais un fichier nommé `VERSION`** dans un dossier d'en-têtes : MSVC le
   prend pour `<version>`), `cpp/shaders/*.comp`, `apex/gpu/vk_context.hpp|vk_internal.hpp`,
   `core._cpp.gpu_status() / gpu_enabled() / set_gpu_enabled()` (`APEX_USE_GPU`, défaut éteint).
   RTX 3090, fp64 actif, `saxpy` identique au bit (`tests/test_gpu_runtime.py`). Aucun noyau métier.
6. **Règle du maillage grossier de la boucle Bird DÉSACTIVÉE par défaut** (décision utilisateur,
   `BirdLoopParams.max_mean_tet_volume = None`) : elle faisait passer le dino de 4 263 à 151 695
   tétraèdres.
7. **Module installé** : `cpp/lib` reconstruit avec `--vulkan` (tout ce qui précède dedans).
   Modules de travail : `cpp/build_dev`, `cpp/build_vk2` (jetables).
8. **Tests : 504/504** en `APEX_USE_CPP=0`, `=1`, `=1` + `APEX_USE_GPU=1` (avant le changement
   du point 6, qui ne touche aucun test : 26/26 sur les fichiers Bird et orchestrateur après).

### Constat non corrigé
- Le contrôle des dégagements hérité de la recherche ne teste que 16 points sur 46 mm d'axe pour un
  rayon de 0,2 mm : un cordon entre deux échantillons n'est pas vu. À reprendre (capsule) avec mesure.

### Fait après la reprise (même jour, détail CLAUDE.md #77)
- **Barycentriques natifs** (`apex_core.tet_barycentric_candidates`, tolérance déclarée, ×20) dans
  les trois appelants ; `APEX_NATIVE_BARY=0` pour comparer.
- **Règle explicite d'ex aequo de `locate`** (`LOCATE_TIE = 1e-9`, change les deux chemins).
- **Raffinement par lots** (identique à l'octet). Benchy et dino vérifiés de bout en bout ;
  fichiers sur `H:/apex_bird_compare_2026-10-13/bary2/` (nouvelle référence d'identité :
  `benchy_refine.gcode`, `dino_refine.gcode`).
- Module `cpp/lib` réinstallé (`--vulkan`). **Toujours rien de commité.**

### Suite (ordre du profil, plan approuvé)
1. Instrumenter le tranchage déformé du dino (~70 s non attribuées sur 99,6 s ; P3 probable).
2. P2 « dernier moment accessible » (19 s Benchy, 55 s death ; pool de processus aujourd'hui).
3. Étape 8 Bird en tableaux (`run_from_passes`, ~12 s de boucle d'objets Python).
4. Dijkstra (networkx) + champ initial + topologie en C++ (17 s Benchy, 63 s death) ; reproduire le
   départage exact de networkx.
5. P4 (chaînage `lines_multi`, doublons, fusion des micro-zones).
6. Phase 3 GPU : uniquement là où le C++ multithread laisse plus de 1 s par appel.
- **Méthode de contrôle d'identité** (à refaire après chaque pas exact) :
  `APEX_USE_CPP=1 .venv/Scripts/python.exe -u scripts/profile_bird.py "research/S4_Slicer/input_models/tilted benchy undertriangulated_repaired.stl" --from-pass="research/S4_Slicer/output_models/loop_tilted_benchy_undertriangulated_repaired_pin_,4" --gcode=H:/.../x.gcode`
  puis `cmp` avec `H:/apex_bird_compare_2026-10-13/fixed/tilted benchy undertriangulated_repaired_bird.gcode`
  (et, avec la boucle : sans `--from-pass`, contre `H:/apex_bird_compare_2026-10-13/tilted benchy undertriangulated_repaired_bird.gcode`).
- Rappels : un seul calcul lourd à la fois (machine 12 cœurs ; une compilation complète fausse un
  profil) ; tuer les processus enfants après un arrêt ; sorties lourdes sur H:.

## >>> SUITE DU 2026-10-06 (soir, autonome hors UI) — LIRE D'ABORD <<<

Détail dans CLAUDE.md, limites #72 et #73 ; état courant dans `DT.md`. Rien n'est commité.

- **`config.g` est le fichier RRF actif.** Limites de B : 1667 °/s et 1300 °/s². **Rappeler à
  l'utilisateur, jusqu'à confirmation, de passer Z à `M203 Z8000` / `M201 Z1500`** (400 mm/min actuels :
  Benchy ~1,7 h -> ~9,7 h). Mémoire `nag-config-g-z-limits`.
- Estimation du temps tous axes : `GcodeReport.print_time_fw_s`, `GcodeParams.fw_axis_speeds`.
- Points 6 à 8 du DT mesurés (couverture, décalage sur le plateau, doublons du lissage).
- Doublons du lissage de Squirtle : `skin_paths.dedupe_shared_passes` (recherche), Squirtle rejoué,
  5,29 -> 2,28 %. Ancien état sur H: (`squirtle_pass8_before_dedupe`).
- step1 : repli `_tetrahedralize_with_translation_retry` (turkish décalé échouait selon l'arrondi).
- Performances, pipeline complet death **124 -> 37 s**, résultats identiques : boucle quadratique de
  step10, résiduel par lots dans step8, contrôleur de collision rapide mis en cache.
- **♦️ ouvert** : step10 range chaque diagnostic d'orientation de step8 en « support de secours »
  (115 312 sur death) : décision utilisateur.
- Tests : 473/473 dans les deux positions.

## >>> SESSION DU 2026-10-06 (ancien poste `slicer_engine`, utilisateur présent) — LIRE D'ABORD <<<

**Le suivi à jour est dans `DT.md` (Dev Tracker). Détail technique : CLAUDE.md, limites #69 à #71.**

**Contexte** : retour sur l'ancien poste (`R:\APEX_NEW_VERSION\slicer_engine`, utilisateur `test_`, disque H: présent).
`main` local a été REMPLACÉ par `origin/session-2026-10-01` (historique sans ancêtre commun avec
`origin/main`) ; l'ancien `main` est dans la branche locale `backup-main-2026-10-06`. Les fichiers
non versionnés d'origine (`build/`, `out/`, `.vscode/`, `settings.fr-code`) sont dans
`H:/apex_untracked_backup_2026-10-06/`, avec les anciennes trajectoires Bird du Benchy et du dino
(`bird_outputs/`). **RIEN n'est commité ni poussé.** Les deux `.venv` ont été retirés de l'INDEX git
(`git rm --cached`), pas du disque. R: avait 0 octet libre au départ (partition étendue par l'utilisateur).

**Environnement à recréer sur un nouveau poste** : `python -m venv .venv`, puis
`pip install cmake==4.4.3 ninja==1.13.2 numpy==2.5.2 scipy==1.18.0 trimesh==5.0.0 tetgen==0.8.4
rtree==1.4.1 manifold3d==3.5.2 networkx==3.6.1 shapely==2.1.2 pytest==9.1.1 pymeshfix==0.18.1
pybind11==3.1.0 matplotlib==3.11.2 pillow==12.3.0 libigl==2.6.3 potpourri3d==1.4.0` ;
noyau natif : `.venv/Scripts/python.exe scripts/build_cpp.py`. Toute modification de dataclass de
paramètres exige `scripts/gen_params.py` PUIS un rebuild (le test `test_generated_files_up_to_date`
échoue sinon : arrivé 3 fois dans la session).

**Résultats de la session (tous mesurés, détail en limites #69-#71)**
1. **Signe de B : ERREUR PUIS CORRECTION.** J'ai d'abord posé `b_sign = -1` en lisant « B négatif =
   vers l'avant » comme une inversion à appliquer. Faux : c'est le sens de l'axe. Croquis de
   l'utilisateur + aperçu : B négatif envoie la pointe vers l'EXTÉRIEUR. `b_sign = +1` (B firmware =
   B d'Apex, limites réelles 0-185°). Tout G-code produit avec B négatif est inutilisable.
2. **Normale par point activée par défaut** (#69) : écart tête/normale médian ÷2 à ÷3, temps d'impression
   +65 à +130 %. **Filtre médian de B (9) activé par défaut** : −19 à −31 % de temps, collisions en baisse.
3. **Limites de B** réconciliées avec les fichiers RRF : 360 °/s, 300 °/s² (aucun effet sur le temps).
4. **Chaîne Bird** : `full_step89.py --apex-detect` (détection par le step9 rapide natif) : dino 855 → 74 s,
   Benchy 496 → 50 s, trajectoires identiques AU BIT ; activé dans `out/run_bird_chain.sh`.
5. **Plateau de step9 : surface à z = origine** (elle était à +0,5 mm : disque centré sur l'origine).
6. **Dégagements de trajet à vide** : essai de l'autre axe puis pas double si contact (Benchy 3 → 0, dino 3 → 0).
7. **Export** : `GcodeReport.y_turns`, `GcodeParams.max_y_turns` (Benchy ≈ 32 tours de plateau, l_shape ≈ 240,
   dino ≈ 621), `GcodeParams.travel_clearance_max_factor`.
8. **Contours « invalides » Bird : pas de correction** (doublons de cordons 0,3-1,5 %, sans corrélation).
9. **SDK Vulkan 1.4.363.0** installé (winget) ; aucun code GPU écrit.
10. **`research/S4_Slicer/safe_io.py` recréé** (absent de la branche ; 7 scripts l'importent).

**Fichiers d'essai et livrables** (`out/`) : `test_b_direction.gcode` (À SEC, pointe immobile, sans
chauffe : le seul à lancer en premier sur la machine), `test_l_shape.gcode`, `test_b_preview*.gcode` et
`test_y_ring_preview.gcode` (aperçu UNIQUEMENT : ils chauffent et extrudent, ne PAS les lancer sur la
machine). G-code : `research/S4_Slicer/output_models/loop_tilted_benchy_undertriangulated_repaired_pin_pass4.gcode`,
`loop_dino_pin_pass5.gcode`. Archives Deflate `out/share/*_v5.zip` (les fichiers font ~65 Mo : un canal
qui tronque à 1 Mo les corrompt). Scripts ajoutés : `scripts/gen_b_direction_test.py`
(`--preview`, `--ring`, `--b-max`), `scripts/export_model_gcode.py` (STL → G-code par le pipeline FEM).

**Pièges rencontrés (à ne pas refaire)**
- Un `git checkout -B` échoue si des fichiers non versionnés seraient écrasés : les déplacer d'abord.
- Les calculs lancés en parallèle d'une modification de `main.py` lisent la NOUVELLE valeur (une mesure
  de temps a dû être refaite) : ne rien modifier pendant qu'une mesure tourne.
- `python - <<EOF` avec des chaînes contenant des antislashs : `"\1"` devient un caractère de contrôle
  (CLAUDE.md a dû être réparé) ; utiliser `chr(92)` ou ne pas écrire de chemins Windows en clair.
- `scripts/gen_params.py --help` ÉCRIT les fichiers générés (pas d'aide) ; `--check` seul est inoffensif.
- Les fichiers du projet peuvent être en CRLF : lire avec `newline=''` pour comparer.
- Un échec de test isolé (`test_export_extrusion_and_retraction`, `APEX_USE_CPP=0`, pendant un calcul
  en parallèle) est NON REPRODUCTIBLE (46 exports sous charge, 0 échec) ; cause inconnue.

**En attente de l'utilisateur** : essai à sec de `test_b_direction.gcode` sur la machine (sens de B,
pointe immobile), puis l_shape, puis le Benchy ; validation du Y cumulé côté câblage ; tout ce qui touche
l'UI est réservé à la toute fin (config UI, sélecteur d'angle, conseil « ne pas centrer »).

**Suite de tests : 469/469 dans les deux positions de `APEX_USE_CPP`** (noyau reconstruit).

## >>> RÉCAPITULATIF DU 2026-10-01 (nuit + journée) — LIRE D'ABORD <<<

**Résultats principaux (détail plus bas et dans CLAUDE.md, limites #56 et #59 à #67)** :
1. **Ordre des couches contraint par la tête** (`pipeline/step5b_collision_order.py`, limite #67,
   activé) : collisions de la vraie tête death 41,6 → 0,13 %, turkish 53,6 → 0,21 % ; matière sans
   appui death 0,227 → 0,222 %, turkish 0,70 → 0,17 % ; l_shape et l-triple inchangés. Le pipeline
   FEM devient imprimable sur death et turkish (au pivot près).
2. **Parois trouées depuis l'origine** (limite #65) : 7,7 à 15,7 % de la longueur jamais extrudée ;
   corrigé (`Step6Params.continuous_rings`, activé). Tous les G-code antérieurs sont concernés.
3. **step10 exact ne pouvait jamais réorienter** (filtre de la zone de buse oublié, limite #66) :
   corrigé.
4. **step9 rapide activé par défaut** (contrôlé contre l'exact sur les 4 modèles, 0 erreur) et
   accéléré 20 → 2,9 s ; **`radial` tranché** par la mesure (gagne partout).
5. **step3 : tri stable** (la segmentation dépendait de l'implémentation de numpy).
6. **Portage C++** : pas 11 à 15, 17, 19, 21a-c (step2, step3, step7), 24. Pipeline death complet
   STL → G-code : **57 → 20,6 s** (dont ~7 s de TetGen). Le reste du portage est architectural ;
   l'utilisateur a choisi « salon d'abord ».
7. **Salon** : Benchy et dino re-tranchés par la voie Bird avec tous les correctifs (11 min chacun),
   G-code produits, rapport ci-dessous (« Après-midi du 2026-10-01 »).

**EN ATTENTE DE L'UTILISATEUR** : convention pivot / pointe contre son prévisualiseur (aucun G-code
imprimable avant ; régénération : `export_gcode.py ... --tilt-sign= --z-reference=`) ; fichiers
RRF ; choix de la normale par point (qualité contre temps) ; SDK Vulkan.

**Suite prévue** : traiter les 4 dégagements de trajet à vide en contact (Benchy 3, dino 1) ;
brancher le step9 rapide d'Apex dans la chaîne Bird (4-6 min → secondes) après le salon ; finir le
portage (step10, step5, step4, step1, orchestrateur natif) après le salon.

**Suite de tests : 461/461** dans les deux positions de `APEX_USE_CPP`, `cpp/lib` reconstruit.

## >>> SESSION AUTONOME 2026-10-01 (nuit, utilisateur absent) — JOURNAL <<<

Plan : `~/.claude/plans/ok-je-vis-aller-parallel-tiger.md` (phases 0 à 5). Aucun commit ni push
(l'utilisateur gère git). Décisions prises seul : notées ici, avec leur raison.

**En attente de l'utilisateur (tâches reportées, rien deviné)** :
- **Modèle pivot / pointe** : où sa prévisualisation place la pointe pour les 4 lignes de test
  (Benchy lignes 11, 1838, 67452, 227838), ou sa formule. L'export est désormais paramétré
  (`GcodeParams.tilt_sign`, `z_reference`) : la correction tiendra en une ligne.
- **Fichiers RRF** (`config.g`, homing) : limites de Y (jusqu'à 35 tours cumulés sur le Benchy),
  Z < 0 autorisé ou non, vitesses et accélérations par axe.
- **SDK Vulkan** (droits administrateur) : backend GPU de `MeshDistance`.
- **Essais d'impression** : aucun G-code n'est à lancer tant que la convention pivot n'est pas
  tranchée.
- **Ordre contraint par la tête (step5b, activé)** : rien à arbitrer (aucune dégradation), mais à
  regarder : GIF `out/death_order_step5b.gif`.
- **`radial` par défaut** : TRANCHÉ par la mesure (gagne sur les 4 modèles avec le nouvel ordre).
- **Normale par point** : décision toujours en attente (qualité contre temps d'impression).

### Phase 0 — G-code (FAIT)
- **Pointe au sommet du cordon** (`GcodeParams.tip_at_bead_top`, activé par défaut). Fait
  vérifié, pas un réglage : step5 place la première isovaleur à v_min + δ/2
  (`step5_layers.py:747`), donc les points de trajectoire sont au MILIEU de la couche. La
  pointe est relevée de h/2 le long de l'axe de la tête (h = 0,2 pour le corps comme pour le
  lissage). Première couche : pointe à 0,2 au lieu de 0,04-0,11.
- **Modèle de tête paramétré** : `tilt_sign` (sens du déport radial du pivot) et `z_reference`
  (`"tip"` = référence d'origine, `"pivot"` = Z firmware = hauteur du pivot). Défauts =
  comportement d'origine.
- **Garde-fous dans `GcodeReport`** : `tip_z_min`, `head_z_min`, `tip_below_plate`,
  `head_below_margin`, étendues X/Y/Z firmware, `written`. `GcodeParams.strict` refuse
  d'écrire un fichier dont la pointe passe sous le plateau.
- **Outil** `scripts/gcode_tilt_profile.py <gcode> [--png=]` : tables par tranche de B sous
  4 conventions (apex, pivot-z, sign-, raw). Sur l_shape : la convention « raw » (X/Z lus comme la
  pointe) met 34 % des mouvements sous le plateau, dès B > 45° ; « pivot-z » les y met tous.
- Tests : 12/12 dans `tests/test_gcode_export.py` (3 nouveaux).

### Phase 1 — mesures (en cours, `out/run_matrix.sh`, résultats `out/measure_matrix.jsonl`)
- Script `scripts/measure_matrix.py <stl> <variante>` : une ligne JSON par mesure (step9b,
  épaisseur de couche réelle mesurée, diagnostics step8, sauts de B, trajets).
- **Découverte** : `step8_diagnostics` (orientations inatteignables) est ÉNORME avec la
  cinématique couplée : 19 325 / 20 405 points sur l_shape, 289 020 sur death. Le champ FEM
  demande des inclinaisons vers le centre que B ne sait pas faire (limite #59). À traiter en
  phase 3.1.

### Phase 2 — portage C++ (en cours)
- Pas 11 (schéma des paramètres + `params.hpp`), 13 (marching tetrahedra), 14 (tracé de bord) :
  FAITS, au bit près ; détail dans CLAUDE.md limite #56. step5 sur death : 2,56 → 1,14 s.
- Clipper2 2.0.1 vendoré (`cpp/third_party/clipper2`, BSL-1.0 vérifiée) pour le pas 15 (step6).
- Module de développement : `scripts/build_cpp.py --dev` + `APEX_CPP_LIB_DIR=cpp/build_dev/lib`.
  `cpp/lib` n'est pas reconstruit tant que la matrice de mesure tourne : **à reconstruire en fin
  de session** (`scripts/build_cpp.py`).

- Pas 17 (IK couplée native), pas 12 (`point_in_polygon_2d`), TetGen spéculatif (step1), step6
  parallèle + indices de contours : tous identiques au bit ou à 1e-14 près. **death @ 0,8 : 105 s →
  18 s profilé** ; turkish step6 8,0 → 3,5 s.
- Pas 19 (step9 rapide natif, death 233 → 14 s) et pas 24 (export G-code natif, Benchy 20,9 →
  0,8 s, texte identique ligne à ligne) : FAITS. step1 : TetGen en parallèle de sa sonde, et
  réutilisation de la sonde (death 13,8 → 7,7 s, turkish 20,2 → 12,8 s), identique au bit.
- Piège rencontré : comparer à une version rechargée depuis git exige de partager ses CLASSES
  (enums) avec le module courant, sinon faux écart (`scripts/check_against_git.py` corrigé).

### Phase 3 — fiabilisation
- 3.2 : `MachineProfile.platform_radius` (150 mm pour la machine actuelle), utilisé par step9.
  Refus des modèles impossibles à placer (`ModelTooLargeError`, plus petit cercle englobant de
  l'empreinte, indépendant du placement) dans `main` et `export_slicing_preview`.
- 3.1 : **inclinaisons atteignables ACTIVÉES par défaut** (limite #64 de CLAUDE.md : mesurées sur
  les 4 modèles, matière sans appui en baisse partout, écart tête / normale 16-49° → 3-6°).
- **À DÉCIDER PAR L'UTILISATEUR** : normale par point (`Step8Params.per_point_normal(_nodal)`),
  meilleure orientation (écart divisé par 2) contre +38 à +130 % de temps d'impression.
- `feedback3` : améliore 3 modèles sur 4, dégrade turkish de 6 % → reste opt-in.
- **Collisions (step9 rapide, désormais natif : death 233 → 14 s)** : avec `radial`, l-triple
  10,5 % → 0 %, turkish 61 % → 54 %, **death 31 % → 41 % (dégrade)**. `radial` maintenu par défaut
  (métrique step9b meilleure partout) mais **À VALIDER PAR L'UTILISATEUR**. Surtout : le pipeline
  FEM par défaut n'est PAS imprimable sur death/turkish (défaut d'ordre : la tête touche de la
  matière déposée bien plus tôt, ~9 mm plus haut). La voie Bird reste la voie d'impression.

- **step10 corrigé** (limite #59) : en machine couplée, la réorientation ne propose plus que des
  candidats à C figé ; en mode step9 rapide, step10 revérifie avec le même contrôleur rapide.
  Pipeline complet death @ 0,8 (step1→step10, step9 rapide) : **57 s**, 24 789 collisions,
  `converged=True` mais **24 772 en support de secours** et 17 seulement résolues par
  réorientation. « Convergé » ne veut donc PAS dire imprimable : la réorientation (B seul, C
  imposé) ne peut presque rien contre le défaut d'ordre.
- **step9 rapide contre exact** (`scripts/compare_step9_modes.py`) : l-triple 680 événements
  communs sur 681 (le seul manquant à 0,000 mm de pénétration), écart de pénétration ≤ 0,072 mm
  (borne théorique 0,125), **1 s contre 552 s**. Les 3 autres modèles : `out/compare_step9_3models.log`.

### DÉFAUT MAJEUR CORRIGÉ : parois trouées (limite #65 de CLAUDE.md)
Le découpage en tronçons monotones (actif par défaut) ne traçait jamais le segment entre deux
tronçons : **7,7 % de la longueur perdue sur death, 15,7 % sur l_shape**. Le contour n=0 n'était
pas refermé non plus. `Step6Params.continuous_rings` corrige les deux et est **activé par défaut**
(mesuré sur les 4 modèles : death −15 %, turkish −2,6 %, l-triple +0,8 %, l_shape +3 % de matière
sans appui ; longueur extrudée +8 à +20 %). **Tous les G-code produits avant ce correctif ont des
parois trouées** sur la voie FEM. Voie Bird (recherche) : ses couches sont exactement planes dans
l'espace déformé (écart 2e-16 mm), donc un seul tronçon par anneau ; seul le segment de fermeture du
contour n=0 y manquait (≈ 0,1-0,9 % de la longueur), corrigé aussi dès le prochain tranchage.

### Bug step10 corrigé (limite #66)
En mode exact, step10 oubliait le filtre de la zone de contact de la buse : chaque point testé
était en collision avec son propre cordon, donc aucune réorientation possible. Corrigé + test.

### step8 en tableaux (au bit) : death 1,92 → 0,82 s.

### Fin d'après-midi du 2026-10-01 — corrections après les rendus de l'utilisateur
- **Pivot à 46 mm** (utilisateur ; 45 avant).
- **La pointe quittait sa trajectoire PENDANT les mouvements** (saut de B, interpolation linéaire
  du firmware) : jusqu'à 26 mm sur le Benchy. `GcodeParams.max_tip_deviation=0.05` découpe les
  mouvements ; après : 0 mouvement à plus de 0,1 mm. Test ajouté.
- **Trous du dino** : 83 % dans la zone de lissage du sommet de la tête, retirée en bloc après 3
  tours de rognage. `skin_loop.py` passe à 8 tours par défaut (converge en 4-5) : surface
  découverte dino 1,32 → 0,19 %, Benchy 0,11 % (1 zone retirée au lieu de 5).
- **Bosses du Benchy** : pas de matière hors pièce (0 % au-delà de 0,3 mm) ; ce sont les zones NON
  lissées (escaliers bruts) : bas de coque (collision plateau, 26 %), flancs de cabine
  (inatteignables avec C imposé, 38,7 %). Pistes : excentrer la pièce / la surélever ; l'axe
  linéaire de l'upgrade.
- Les fichiers reçus par l'utilisateur via un autre canal étaient TRONQUÉS à 1 000 000 d'octets ;
  envoi en zip Deflate (`out/share/*_v3.zip`).
- Outils : `scripts/render_gcode_check.py`, `scripts/gcode_surface_coverage.py`.
- Push GitHub BLOQUÉ (connexion demandée par git-credential-manager sur ce poste) ; commit local
  sur la branche `session-2026-10-01`, les correctifs de fin d'après-midi NON commités.

### Après-midi du 2026-10-01 — « salon d'abord » (décision de l'utilisateur)
Chaîne Bird refaite avec les correctifs du jour (`out/run_bird_chain.sh`, journal `out/bird_chain.log`) :

| | Benchy (passe 4) | dino (passe 5) |
|---|---:|---:|
| durée totale (déformation déjà faite → G-code) | **11 min 9 s** | **11 min 5 s** |
| dont boucle de lissage (P2/P4/P5 ×3, P3 final) | 384 s | 285 s |
| dont étapes 8+9 (step9 de recherche, Python) | 259 s | 363 s |
| zones de lissage retirées (collision) | 5 | 3 |
| collisions tête restantes | 12 points du corps, pénétration ≤ 0,38 mm DANS la marge de 1 mm | **0** |
| points du corps à > 45° de la normale | 1,6 % | 0,1 % |
| dégagements de trajet à vide en contact franc | 3 | 1 |
| filament | 4,53 m | 6,26 m |
| temps d'impression (limites de B) | **3,28 h** | **4,08 h** |

Sur ce PC, P2 prend 3-16 s et P3 43-57 s (340-510 s sur l'ancien) : la refonte 3D de P3 n'est plus
prioritaire. Incident : P4 a échoué une fois sur un verrou transitoire de Windows en réécrivant son
.pkl → `research/S4_Slicer/safe_io.py` (réessai des écritures jusqu'à 5 s), importé par les 7
scripts de la chaîne. `export_gcode.py` accepte `--tilt-sign`, `--z-reference`,
`--no-tip-at-bead-top`, `--strict` : régénération en une commande une fois la convention pivot
tranchée. G-code produits : `output_models/loop_*_pin_pass{4,5}.gcode` (NON imprimables avant).

### Matinée du 2026-10-01 (utilisateur présent)
- `radial` TRANCHÉ par la mesure (gagne sur les 4 modèles avec le nouvel ordre).
- Portage : pas 21a (step2), 21b (step3, chemin par défaut), 21c (step7 entier), tous IDENTIQUES au
  Python sur les modèles réels ; MeshDistance rend l'indice du triangle le plus proche.
- step3 : tri des arêtes rendu STABLE (le résultat dépendait de numpy : 11 ou 24 régions sur
  turkish) ; mesuré neutre ou meilleur (turkish 0,198 → 0,172 %).
- step9 rapide 20 → 2,9 s ; step5b 6 → 3 s. **Pipeline complet death : 57 → 20 s.** Suite 461/461.
- Non fait, et pourquoi : reste de step5 (≈ 0,6 s, repères locaux dont les sommes diffèrent entre
  numpy et C++ : le porter casserait l'équivalence aval pour un gain minime).

### RÉSULTAT MAJEUR : step5b, ordre des couches contraint par la tête (limite #67), ACTIVÉ
Si la tête posée sur la couche L touche la couche M, alors L avant M (tri topologique stable).
Avec les contraintes d'appui (M sous L ⇒ M avant L), défaut retenu : collisions de la vraie tête
**death 41,6 → 0,13 %, turkish 53,6 → 0,21 %** ; matière sans appui **death 0,227 → 0,222 %,
turkish 0,70 → 0,198 %** ; l_shape et l-triple inchangés. Aucun modèle ne se dégrade. Suite 450/450.
Pipeline complet step1→10 : death 35 s / 88 points en collision (27 855 avant) ; turkish 75 s /
322 points (~76 000 avant).

### step9 rapide ACTIVÉ par défaut (main.PipelineParams)
Contre l'exact : l-triple et l_shape entiers, turkish (2 000 points) et death (600 points) par
échantillon : 0 faux positif, 0 faux négatif, écart ≤ 0,087 mm. Suite : **445/445** dans les deux
positions, `cpp/lib` reconstruit.

### Mesures diverses (fin de nuit)
- `enforce_upward_growth=False` avec radial : death 0,23 → **11,4 %**, turkish 0,70 → 3,55 % de
  matière sans appui → le redressement reste indispensable (piste « flip conditionnel » close).
- Remplissage : couverture réelle conforme (death 22 % pour 20 % visés) ; pas de défaut.
- Solveur step4 : Jacobi CG (1,25 s) reste le meilleur disponible (SuperLU 110 s, ILU scipy ne
  converge pas).
- Contrôles physiques ajoutés : `scripts/check_wall_coverage.py`, `scripts/compare_step9_sampled.py`,
  `scripts/check_step8_arrays.py`, `scripts/check_native_fill.py`.

### Pas 15 — step6 natif (Clipper2), opt-in
`apex::fill_region` (`cpp/src/core/toolpath2d.cpp`) : parois par offset Clipper2 + remplissage
LINES / GRID / CROSSED_ZIGZAG. Option `APEX_NATIVE_FILL=1`. Écart par couche ≤ 1e-4 sur l_shape et
l-triple ; death : 12 couches > 1 %, dues aux onglets (GEOS prolonge une fente fine en pointe de
4 mm qui ronge la paroi, Clipper2 coupe en carré : le natif est plus juste). step6 sur death
1,81 → 0,51 s. Mesure step9b : `out/measure_native_fill.log`.
`_lift_to_3d` : `np.cross` remplacé par un produit vectoriel scalaire, au bit (0,8 s sur death).

### Phase 1.3 — contours « invalides » des tranches Bird (mesure)
La conclusion du 30/09 (« recouvrements bénins de la projection XY ») était **fausse**. Les couches
Bird sont tranchées dans l'espace DÉFORMÉ, où elles sont planes (`tet_vertices = pk`, champ z
déformé) : les contours invalides y sont de vraies auto-intersections, reçues telles quelles par
step6.
| modèle (passe) | couches touchées | auto-intersections | part de l'aire |
|---|---:|---:|---:|
| death_repaired_pin (3) | 104/217 | 163 | 53,5 % |
| Squirtle (8) | 111/218 | 132 | 27,0 % |
| Benchy (4) | 26/313 | 27 | 4,6 % |
- **Propre à Bird** : sur le pipeline FEM par défaut, AUCUN contour invalide (death 0/581 groupes,
  l_shape 0/166, l-triple 0/170 ; `scripts/check_contour_validity.py`).
- Effet sur step6 (remplissage mesuré / attendu, `bird_slice_deformed.py`, attendu calculé sur
  les polygones réparés) : death valides 1,22 / invalides **1,55** ; Squirtle 1,20 / **1,34**.
  Les couches invalides reçoivent 13-27 % de remplissage en plus. **Pas conclu** : dans un pli, la
  zone recouverte de l'espace déformé correspond à DEUX endroits réels distincts, donc un
  remplissage « en double » y est peut-être juste. Décider exige de mesurer après retour à la
  pièce réelle (doublons de cordons à < w/2 l'un de l'autre). Reporté à la phase 4 (Bird dans
  Apex) ; aucune correction appliquée.

## >>> NOUVEAU POSTE (fin du 2026-09-30) — LIRE D'ABORD <<<

**Migration** : le projet est maintenant dans `R:\APEX_NEW_VERSION\slicer_alg` (utilisateur
`prazl`), et non plus dans `slicer_engine`. Les chemins `C:\Users\test_\…` cités plus bas sont
périmés : les plans sont dans `~/.claude/plans/` (restaurés depuis `claude_context/`), la mémoire
dans `~/.claude/projects/r--APEX-NEW-VERSION-slicer-alg/memory/`. **Il n'y a PAS de disque H:**
sur ce poste : les mentions de sauvegarde et d'archive sur H: ci-dessous ne s'appliquent plus.
- Dépôt : `origin` = `aeloxhardware/apex-dev` (instantané complet : recherche, `.venv`, sorties,
  `claude_context/`). Commits faits par l'utilisateur ; ne pas commiter sans demande.
- `.venv/pyvenv.cfg` repointé vers `C:\Users\prazl\…\Python314` (même 3.14.7).
- Module C++ reconstruit sur ce poste (VS 2022 via `scripts/build_cpp.py`) ; l'ancien cache
  CMake (chemins `slicer_engine`) est mis de côté dans `cpp/build_old_slicer_engine/`, l'ancien
  module dans `cpp/lib/*.pyd.bak1001`.
- GPU de ce poste : **GTX 1050 Ti** (la RTX 3090 citée plus bas était sur l'ancien). Runtime
  Vulkan présent, **SDK Vulkan absent** ; CUDA 12.9 installé.

**Réponses de l'utilisateur (2026-09-30)** :
- zéro Z du firmware = buse qui touche le plateau → pas de décalage de 45 mm, l'export est juste ;
- vitesses validées : 25 mm/s, trajets à vide 80 mm/s, 195 °C ;
- GPU : **Vulkan** (portage prévu sous macOS et cartes AMD/Intel) ;
- archive H: : sans objet (pas de disque H:).

**Reste à faire** (liste tenue avec l'utilisateur) :
- *Avant Maker Faire* : 1) premiers essais d'impression Benchy et dino (plus d'obstacle
  logiciel) ; 2) vérifier les ~100 points par modèle à plus de 45° de la normale ;
  3) real_overhang (127 contacts) : imprimable ou retiré de la démo ; 4) confirmer que les
  « contours invalides » des couches courbes sont bénins.
- *Fiabilisation Apex* : mesurer sur les 4 modèles les mécanismes opt-in du 30/09 ; remesurer
  step9b ; porter la chaîne Bird (en option) ; cinématique couplée dans step3/7/10 ; alternance
  parfaite de B ; rayon du plateau + refus des modèles trop grands.
- *Après le salon* : axe D / axe linéaire ; interface (ne pas centrer, sélecteurs d'angle) ;
  performances (P3 en 3D, **backend Vulkan de `MeshDistance` à écrire** après installation du
  SDK, portage C++ step5-7, parallélisation) ; config/schema, hook step10, path-planning, GIF P6,
  ironing ; chantier de fond « matière sans appui ».

## >>> (précédent) session autonome du 2026-09-30, sur l'ancien poste <<<

**Sauvegardes (fin de session, le poste doit être remplacé)** — PÉRIMÉ, voir ci-dessus :
- GitHub `aeloxhardware/apex`, branche main : tout le code versionné et la documentation.
- NON poussé sur GitHub (envoi vers `apex-dev` bloqué par le contrôle de sécurité du mode Auto) :
  `research/S4_Slicer` (dont `DECISIONS_2026-09-30.md` et les scripts du jour), les sorties,
  `.venv` et `build`. Copie complète locale seulement : `H:\apex_backup_2026-09-30\slicer_engine`
  (27 061 fichiers, 2,35 Go), sur l'ancien ordinateur.

Plan : `~/.claude/plans/snug-yawning-lovelace.md` (T1 à T9). Décisions prises seul et chiffres :
`research/S4_Slicer/DECISIONS_2026-09-30.md`.

**Fait :**
- **T1 — P5** : tous les points sont testés, y compris contre les lignes de la zone en cours.
  - Dino contraint : le contact manqué (zone 23, 0,03 mm) est vu, plus un second jamais vu
    (zone 4, 0,05 mm).
  - Arbres par blocs : 890 s → 239 s, pour un résultat identique.
  - Dans la boucle de rognage, P5 échantillonne tous les 2 mm ; seul le P5 final est exhaustif.
- **T2 — Export G-code** : `file_io/gcode_export.py` et `tests/test_gcode_export.py`, 9 tests.
  - Pivot de 45 mm, Y = atan2, F à vitesse de pointe constante, `feed_plan`.
  - Rétraction de 0,8 mm à 35 mm/s.
  - Dégagement sur un seul axe (Z ou R) selon l'orientation de la tête.
  - Script de recherche : `research/S4_Slicer/export_gcode.py`.
- **T3 — Aperçu** : `scripts/preview_gcode_file.py`, rendu PNG et contrôle aller-retour.
- **T4 — P2** : l'écart tête / normale est jugé au pire des 3 sommets et du centre (non encore
  remesuré sur le Benchy).
- **T5 — Refonte de P3** :
  - la boucle de rognage ne lance P3 qu'une fois, à la fin ; pendant la boucle, P5 teste contre le
    corps non creusé ;
  - creusement parallèle : identique au bit près, environ ×7.
- **T7 — Portage dans Apex (opt-in, désactivé par défaut)** :
  - `CollisionSimParams.fast_kdtree` et la classe `FastCollisionChecker` (step9) ;
  - `clamp_head_height` et `reorient_coupled_tiers` (step10) ;
  - branchement dans `main.py` : `height_clamp`, `coupled_reorientation_tiers` ;
  - tests ajoutés ;
  - CLAUDE.md : limite #59 corrigée, nouvelles limites #61 et #62.
- **T9 — C++** : le profil du creusement désigne la distance point-maillage (67 % du temps).
  - Écrit : `cpp/include/apex/mesh_distance.hpp`, `cpp/src/core/mesh_distance.cpp` et le binding
    `MeshDistance`.
  - Backends : série (oracle), threads, et Vulkan, qui se rabat sur le CPU.
  - Tests d'équivalence face à trimesh.

**Suite de la journée (après la mise en ligne interrompue) :**
- **Exports G-code** du Benchy, du dino et de Squirtle : aller-retour à 1,1e-4 mm sur les trois ;
  aperçus PNG (`output_models/*_preview.png`).
- **Ralentissement dû à B** : diagnostic, puis filtre médian de 9 points avant contrôle
  (`full_step89 --smooth-b=9` ; dans Apex, `Step8Params.free_axis_median_window`, limite #63).
  Squirtle : 4,90 h → 2,97 h, même niveau de contacts.
- **`skin_loop`** : exclusion cumulée jusqu'à zéro contact (défaut trouvé sur Squirtle : 106 mm²
  creusés sans lissage). Squirtle réparé, avec 7 zones exclues.
- **P2** : pré-calcul d'accessibilité mis en cache (il est identique à chaque tour, 291 s) et
  parallélisé.
- **Distance point-maillage native** dans le creusement de P3 : ×3,1, régions identiques.
- **Sorties obsolètes archivées** dans `H:/apex_archive_2026-09-30/` (718 Mo, liste jointe).
  R: ne fait que 8 Go.
- **Validation sur 7 modèles** (chaîne complète, puis G-code ; tableau dans DECISIONS, section
  T6) : zéro contact dans le lissage partout.
  - Corps sans aucun contact : dino, death, l-triple, l_shape.
  - Corps avec des contacts dans la marge : Benchy (9, ≤ 0,22 mm) et Squirtle (9, ~1,1 mm).
  - **real_overhang : 127 contacts non résolus** (tête presque horizontale à mi-hauteur) ; son
    impression avec la tête actuelle est à étudier.
  - Temps réels avec le filtre : Benchy 2,46 h, dino 3,48 h, Squirtle 2,97 h, death 2,39 h.
- **T4 mesuré** (Benchy) : juger aux sommets retire 172 faces du lissage (−3,8 %). P2 complet
  prend 38,6 s, grâce à la parallélisation.
- **« Contours invalides »** des couches déformées (death 108/217, Squirtle 93/218) : c'est la
  projection XY de couches COURBES dans la pièce réelle, qui peut se recouvrir légitimement.
  Déjà présent avant (Benchy 65/284). Probablement sans gravité, à confirmer.

**Questions pour l'utilisateur** (toutes répondues depuis, voir la section « NOUVEAU POSTE ») :
- **apex-dev** : la poussée de tout le répertoire a été bloquée par le contrôle de sécurité ; le
  dépôt privé est vide. Trois options : ajouter une règle d'autorisation, faire l'envoi vous-même,
  ou réduire le périmètre.
- **Archive** `H:/apex_archive_2026-09-30/` : peut-on la supprimer ?
- **GPU** : installer le SDK Vulkan, ou passer par CUDA ?
- **Référence Z du pivot** : le zéro Z du firmware correspond-il à la pointe à B = 0 (hypothèse
  retenue) ou au pivot ? Dans le second cas, il faut ajouter 45 mm.
- **Vitesses provisoires** : impression 25 mm/s, trajets à vide 80 mm/s, 195 °C. À valider.

## >>> (précédent) soir du 2026-09-29, avant compaction <<<

Journée du 2026-09-29 résumée dans les sections qui suivent ; détail complet
dans `research/S4_Slicer/DECISIONS_2026-09-29.md`.

**Acquis du jour (tout mesuré) :**
- step5 corrigé (limite #60, commit `1f00506`, poussé sur
  github.com/aeloxhardware/apex ; `origin` = apex, ancien projet dans la
  branche distante `dev.old`, apex-beta dissocié). 386/386 dans les deux
  positions.
- Benchy et dino refaits avec ce correctif (P3 + P5), 0 collision sur le
  lissage.
- Contrôle de collision rapide de toute la trajectoire (corps + lissage) :
  `research/S4_Slicer/full_step89.py --fast` (arbre k-d, critère de P5,
  validé contre step9, 3 à 4 min par modèle au lieu de plusieurs heures).
- Correction des collisions du corps, les deux mécanismes combinés
  (`--height-clamp --apply`) : borne de B selon la hauteur (vraie tête,
  plateau 300 mm) + réorientation par paliers 45 / 90 / sans limite, lissée
  sur ±1 à 3 mm, vérifiée. Benchy 83 → 9 (tous dans la marge, ~0,8 mm de
  jeu), dino 217 → **0**. Prix : ~100 points par modèle à plus de 45° de la
  normale. Trajectoires corrigées enregistrées :
  `output_models/loop_*_pin_machine_toolpath_pass{4,5}.pkl` (positions
  repère plateau, angles, extrusion, corps / lissage).
- Axe D sans limite (fin de course à roulette) : couvre toutes les
  orientations du corps (mesuré). Y en option : couvre tout aussi.
- Contraintes d'inclinaison Bird (`--reachable-tilt`, `--reachable-normal=w`)
  et décalage sur le plateau (`--offset=dx,dy`) : sélectionnables,
  désactivés par défaut ; compromis mitigé (atteignable ≠ sans collision).

**Ce qui reste, dans l'ordre proposé :**

*Avant Maker Faire (fin octobre)*
1. **Export G-code** (`file_io/gcode_export.py`), écrit de zéro — NE PAS
   s'inspirer de `scripts/gen_*_gcode.py` (formes fixes, décision de
   l'utilisateur). Entrée : les `_machine_toolpath_pass*.pkl`. Commande
   `G00 Y<C> X<R> Z<Z> B<B> F<avance>`, compensation du pivot de 45 mm,
   ralentissement de l'avance aux limites de B (`skin_print.feed_plan`).
   **Question à poser d'abord : rétractions** (distance, vitesse, lever la
   buse pendant les trajets à vide ?).
2. Aperçu du G-code produit (relecture + rendu) avant la machine.
3. Premiers essais d'impression (Benchy, dino) pour fixer les réglages
   provisoires : recouvrement lissage / corps, tolérance tête / normale 45°,
   épaisseur de lissage.

*Pour fiabiliser*
- Vérifier à l'impression la qualité des points à plus de 45° de la normale.
- Porter dans Apex le contrôle rapide + borne + réorientation (à la place
  de la boucle lente de step9).
- Valider sur death, l_shape, l-triple, turkish, Squirtle, real_overhang.
- Défauts connus : points au-delà de 45° jugés au centre des faces ; P3
  couche par couche ; contact dans le lissage vu par `full_step89` mais pas
  par P5 (P5 ne teste pas la zone en cours).
- Accélérer la boucle de rognage (P3 une seule fois, à la fin) — accord de
  l'utilisateur non donné.

*Après le salon*
- Axe D : profil, stratégie de choix entre solutions redondantes, modèle 3D
  de la nouvelle tête, remesure des collisions. Y en option.
- Interface : conseil « ne pas centrer » + contrainte de l'auto-organiseur
  (déjà dans la roadmap de CLAUDE.md), sélecteurs d'angle, refus des modèles
  trop grands.
- Cinématique couplée dans step3, step7, step10.
- GIF de l'ordre d'impression (P6), portage de la recherche dans Apex, C++.

*Tenue du projet*
- CLAUDE.md : corriger la limite #59 (« Bird convient à la machine
  actuelle » est faux en partie) et y reporter les résultats du jour.
- Commit : CLAUDE.md modifié (roadmap « ne pas centrer ») et LAST_SESSION.md
  non commités — à faire sur demande. Code de recherche : local.

## >>> MATINÉE AUTONOME DU 2026-09-29 (utilisateur au lycée) <<<

Détail : `research/S4_Slicer/DECISIONS_2026-09-29.md`.

- **Défaut de step5 corrigé** (limite #60, commit `1f00506`, poussé sur
  github.com/aeloxhardware/apex en écrasant l'ancien projet, à la demande de
  l'utilisateur ; ancien contenu gardé en local : branche
  `backup/apex-old-main`). 386/386 dans les deux positions.
- **Benchy et dino refaits** avec ce correctif (P3 puis P5) : 29 et 30 zones,
  0 collision sur le lissage.
- **Toute la trajectoire dans step8 + step9** (`full_step89.py`, nouveau) :
  - normale du corps = gradient de φ ; écart tête / normale médian 9,1°,
    1,7 % au-delà de 45° ;
  - 16 % des normales penchent vers le centre, 38 % ont une composante
    tangentielle notable. « Bird convient déjà à la machine actuelle » est
    donc faux en partie ;
  - step9 (1 point sur 10, 91 % de la trajectoire en 30 min) : 8 contacts,
    tous dans le corps (6 plateau, tête très inclinée ; 2 pièce), 0 dans le
    lissage. Le débit tombe à 8 pts/s.
  - **Suite (après-midi)** : `full_step89.py --fast` (critère k-d de P5,
    validé contre step9 sur le même échantillon). Contrôle COMPLET : Benchy
    83 points en contact sur 64 793 (69 plateau, 21 pièce), dino 217 sur
    88 453 (134 plateau, 197 pièce), 0 dans le lissage, 3-4 min par modèle.
    Contacts plateau jusqu'à 78 mm sous son plan : plan infini, à trier avec
    le rayon réel du plateau (non renseigné). Prochaine étape logique : la
    contrainte d'inclinaison Bird (paramètre), puis remesurer.
- **Contrainte d'inclinaison Bird FAITE** : `bird_loop.py --reachable-tilt`
  (sélectionnable, désactivée par défaut, sorties `_reach`). Résultats :

  | | Benchy libre | Benchy contraint | dino libre | dino contraint |
  |---|---:|---:|---:|---:|
  | surface en surplomb (plan : 10,9 % / 6,4 %) | 5,6 % | 9,8 % | 0,7 % | 0,3 % |
  | normales vers le centre | 32,4 % | 4,1 % | 4,3 % | 9,2 % |
  | corps > 45° pour la tête actuelle | 1,7 % | 0,3 % | — | — |
  | contacts step9 (points) | 83 | 107 | 217 | 148 (dont 1 dans le lissage) |

  Compromis réel et non universel. Borner la rotation par cellule ne fixe
  pas la pente finale : l'ajustement de la déformation crée ses propres
  pentes (dino). Contact dans le lissage vu par step9 mais pas par P5 :
  écart de critère à examiner. Idée à soumettre : décentrer la pièce sur le plateau pour
  que ses surplombs regardent vers l'extérieur.
- **Plateau 300 mm** (utilisateur) : pris en compte, contacts INCHANGÉS sur
  les 4 cas. Les contacts plateau sont de vraies traversées (tête retournée,
  B > 90°, près du plateau).
- **Décalage sur le plateau** (`bird_loop.py --offset=dx,dy`, Benchy) :

  | | surplomb | vers le centre | > 45° tête actuelle |
  |---|---:|---:|---:|
  | centré, libre | 5,6 % | 32,4 % | 6,2 % |
  | centré, contraint | 9,8 % | 4,1 % | 2,0 % |
  | +70 mm, libre | 3,8 % | 52,1 % | 25,1 % |
  | −70 mm, libre | 4,1 % | 38,2 % | 23,2 % |
  | +70 mm, contraint | 11,3 % | 21,5 % | 5,7 % |
  | −70 mm, contraint | 3,6 % | 18,6 % | 7,0 % |

  La position compte beaucoup, le sens aussi. Borner le signe de la rotation
  par cellule ne garantit pas des normales atteignables : il faudrait
  contraindre la normale FINALE dans la déformation (à soumettre). Roadmap
  CLAUDE.md : conseil UI « ne pas centrer » + contrainte de l'auto-organiseur.
- G-code : ne PAS s'inspirer de `scripts/gen_*_gcode.py` (formes fixes).
- **Contrainte sur la NORMALE FINALE** (`bird_loop.py --reachable-normal=<w>`,
  sélectionnable ; résidus ajoutés à la déformation, bird_core) : domine la
  contrainte sur le signe de rotation. Benchy centré :

  | | surplomb | > 45° tête actuelle | contacts step9 |
  |---|---:|---:|---:|
  | Bird libre | 5,6 % | 6,2 % | 83 (69 plateau) |
  | w = 0,5 | 6,8 % | 2,1 % | 903 (893 plateau) |
  | w = 2 | 8,8 % | 0,9 % | — |
  | w = 10 | 10,4 % | 0 % | — |
  | w = 0,5, décalé −70 mm | 3,6 % | 1,0 % | 1 894 (1 561 pièce) |

  **Leçon** : une orientation atteignable n'est pas une orientation sans
  collision. Pencher les couches près du plateau fait toucher le corps de la
  tête. Piste suivante : borner l'inclinaison selon la hauteur (dégagement de
  la tête), ou mettre le contrôle de collision dans la boucle Bird.
- **Contacts du corps (point 1)** : `full_step89.py --height-clamp --reorient`
  (les deux, comme voulu). Plus aucun contact plateau sur les deux modèles.

  | | départ | après borne + réorientation | non résolus |
  |---|---:|---:|---:|
  | Benchy | 83 | 23 contacts (dans la marge de 1 mm sauf 1) | 22 |
  | dino | 217 | 192 contacts, 97 résolus | 95 (tête ~horizontale, B ≈ 102°) |

  **Suite, réorientation ÉCRITE dans la trajectoire** (`--height-clamp
  --apply` : paliers de tolérance 45 / 90 / sans limite, lissage ±1 à 3 mm,
  vérification complète) :

  | | départ | après correction |
  |---|---:|---:|
  | Benchy | 83 | 9, tous dans la marge (≈ 0,8 mm de jeu réel) |
  | dino | 217 | **0** |

  Prix : une centaine de points déposés à plus de 45° de leur normale (qualité
  locale), aucun en collision. `--apply` enregistre désormais la trajectoire
  corrigée (`_machine_toolpath_pass<k>.pkl`, générée pour le Benchy et le dino), entrée de l'export
  G-code.
- **Variantes de machine** (`variants_reach.py`) : actuelle 2,3 % des
  normales du corps à plus de 45° ; + D (±90°) 1,9 % ; + Y (C libre) 0.
  Pistes : contraindre les rotations Bird aux directions atteignables
  (logiciel, pour Maker Faire) ; Y en option ; D utile surtout avec Y
  (vitesse, choix anti-collision), à condition d'écrire une stratégie de
  choix entre solutions redondantes.
- **Réponses de l'utilisateur** : axe D SANS limite (fin de course à
  roulette) : + D couvre alors toutes les orientations du corps (0 % au-delà
  de 5°). Contrainte « Bird limité aux directions atteignables » : à faire en
  PARAMÈTRE sélectionnable (non fait). Ancien dépôt : branche distante
  `dev.old`. `origin` = aeloxhardware/apex ; apex-beta dissocié.

## >>> APRÈS-MIDI DU 2026-09-28 (utilisateur présent) : décisions P4, CINÉMATIQUE, P5 en cours <<<

**Décisions de l'utilisateur, mises en œuvre**
- **Hilbert.** Avertissement structuré par zone, avec un choix A / B / C
  expliqué aux novices (`sp.HILBERT_OPTIONS`). On choisit avec
  `--hilbert-remedy=` ou `--zone-hilbert-remedy=`. Sans choix, la zone est
  tracée en A et l'avertissement est émis.
  - Mesure (lignes à w ± 5 %, sur un plan 94 %) : B vaut 90 à 93 % partout ;
    C vaut 74 à 88 % ; le texte de C a été corrigé en conséquence.
- **Cordon de bord.** `sp.with_edge_ring` : il suit la distance au bord
  mesurée sur la surface, et il est déposé avant le motif. Le recouvrement
  motif / cordon vaut 25 %, provisoire, à fixer après un essai.
- **Micro-zones.** Elles sont rattachées à une voisine quand la fusion reste
  accessible. Sinon, elles ne sont PAS lissées (statut 6 dans P2-bis) : ni
  creusées, ni tracées.
- **Lissage et ordonnancement** jugés justes. **P5 à programmer** : Apex ne
  sera vendu que fini.

**DÉCOUVERTE CRITIQUE : cinématique (limite #59 de CLAUDE.md).**
- Les faits, confirmés par l'utilisateur :
  - la machine actuelle a R, C, Z et B seulement, donc C est imposé par la
    position du point ;
  - B ne penche que d'un côté, de 0 à 185° ;
  - l'upgrade, après Maker Faire, ajoutera un axe linéaire ;
  - types de machines futurs, à prendre en charge : axe D (seconde rotation
    de tête), plateau coulissant le long de X, ou les deux. C'est noté dans
    la roadmap de CLAUDE.md.
- **Apex corrigé.**
  - `MachineProfile.position_coupled_axes` et `radial_direction` ;
  - étape 8 `_resolve_all_axis_angles_coupled` ;
  - profil par défaut = machine actuelle ;
  - 7 tests ; **384/384 dans les deux positions**.
- Pas encore adaptés : les directions candidates de step3 et step7, et step10.
- **P2 corrigé.** L'IK est couplée à la position, avec une tolérance tête /
  normale de **45°** (provisoire, `--head-tolerance`).

| modèle | surface retenue | inatteignable | zones |
|---|---:|---:|---:|
| Benchy | 32,2 % | 37,9 % | 48 |
| dino | 43,5 % | 11,9 % | 64 |

- P4 relancé : couverture 99,2 % (Benchy) et 99,1 % (dino).

**P5 en cours** : `research/S4_Slicer/skin_print.py`. Chaque zone est insérée
après la dernière couche du corps qui la porte. Les angles viennent de la
cinématique réelle. Chaque point est vérifié contre la vraie tête (plateau, et
matière déjà imprimée à ce moment). Une zone en collision est retirée, et
écrite dans `_skin_excluded_pass<k>.json`. On relance alors P3 avec
`--skin-exclude=<json>` (nouvelle option de `bird_slice_deformed.py`), puis
P5.

- P3 relancé sur le Benchy (commande :
  `bird_slice_deformed.py "tilted benchy undertriangulated_repaired" 4 _pin --skin --first-layer=concentric_conform`).

**Premier passage de P5 sur le Benchy** : 50 zones.
- **Durée** : 31 s. La première version, qui bouclait capsule par capsule
  comme l'étape 9, dépassait 10 min ; elle est remplacée par un arbre k-d sur
  la matière imprimée densifiée à 0,25 mm.
- **Collisions** : 6,2 % des points de lissage touchent la vraie tête, avec
  une marge de 1 mm.
- **Où** : surtout au bord des zones, au pied de parois déjà imprimées. Le
  pré-calcul de P2 teste la tête au centre de chaque face ; P5 la teste le
  long des trajectoires.
- **Règle du plan** (zone entière retirée) : 30 zones sur 50 partent, soit
  86 % de la surface lissée, dont la coque pour 1 point sur 1 717.

**Décision de l'utilisateur : ROGNAGE LOCAL.**
- P5 écrit ses points de collision dans `_skin_collisions_pass<k>.json`.
- `surface_skin.py --skin-trim=<json cumulé> --trim-radius=2` classe en
  « collision objet » les faces situées à moins de 2 mm de ces points.
- `skin_loop.py` enchaîne P2 → P4 → P3 → P5, jusqu'à zéro collision ou
  3 tours au plus ; les zones encore en collision sont alors retirées par P5.
- P2 applique désormais les mêmes règles que l'étape 9 : zone de contact de
  la buse de 3 mm, marge de 1 mm (tête gonflée, plateau relevé).
- Commande :
  `skin_loop.py "tilted benchy undertriangulated_repaired" 4 _pin --first-layer=concentric_conform`
  (environ 9 min par tour, dont P3 : 7 min).

**Fin de P5 (suite) :**
- Limites de B données par l'utilisateur : **750 °/s et 300 °/s²**. Le
  firmware coordonne les axes (réponse de l'utilisateur).
- P5 garde donc les angles EXACTS et calcule le ralentissement
  (`feed_plan`, option `--feed=25` mm/s). Les points au-delà de 45° sont
  rognés comme les collisions.
- Étape 8 d'Apex : troncature désactivée par défaut
  (`clamp_angle_rates`), 1 test ajouté, 28/28. La suite complète est à
  relancer après les boucles.
- Une première boucle, faite avec la troncature, était fausse : la coque
  était retirée pour des retards de B.
- Les boucles sont relancées de zéro sur le Benchy et le dino (journaux
  `output_models/skin_loop_benchy.log` et `skin_loop_dino.log`).

**>>> REPRENDRE ICI (fin d'après-midi du 2026-09-28) <<<**

**Vocabulaire convenu avec l'utilisateur** : **step1 à step10** désignent les
étapes du pipeline d'Apex (code versionné) ; **P0 à P6** désignent les phases
du plan de lissage (recherche).

**Prochaine tâche choisie** : le défaut de **step5** d'Apex. Sur le
tranchage de référence du Benchy, la couche 6 perd 100 % de sa matière et la
couche 174 en perd 23 %. La mesure compare le contour de step5 à la section
réelle des tétraèdres, dans l'espace déformé. Hypothèse non vérifiée : un
plan passant par des sommets. `research/S4_Slicer/check_carve.py` contient
la mesure. L'utilisateur m'a arrêté au moment où je commençais à la lire :
rien n'a été modifié.

**Proposé, en attente de l'accord de l'utilisateur** : accélérer la boucle
de rognage. P3 fait 80 % du temps ; il ne tournerait qu'une fois, à la fin,
et les tours de rognage vérifieraient contre le corps non creusé (plus
prudent de 0,2 mm). Une boucle passerait de 40 à environ 13 min par modèle.
Ensuite viendrait P3 en découpe 3D.

**Gantt réordonnable** : https://claude.ai/artifact/625oV9B7u1cXGXp7yNKXLQ.
Ajout recommandé : vérifier le CORPS (couches Bird) avec la vraie tête, via
step8 et step9 avec la cinématique couplée. Jamais fait : P5 ne contrôle que
le lissage.

**Correspondance entre le code de recherche et les steps d'Apex** (pour le
portage) :
- Bird REMPLACE step3 et step4 : les rotations et la déformation donnent le
  champ φ, c'est-à-dire la hauteur déformée ; P1 va ici.
- step5 et step6 sont réutilisés dans l'espace déformé, suivis du retour
  `DeformedToOriginal`. P3 modifie les contours de step5, P4 s'ajoute à
  step6.
- P2 devient un nouveau step, placé entre step4 et step5.
- P5 va dans step7 (ordre), step8 et step9 (angles et collisions).
- La boucle de rognage correspond à un ajustement de type step10.
- Choix d'architecture laissé à l'utilisateur : Bird remplace step3 et
  step4, ou devient une deuxième source de champ, en option. Ma
  recommandation : l'option, tant qu'il n'est pas mesuré sur les quatre
  modèles.

**Avant de porter le code de recherche dans Apex** (hors C++ et
parallélisme) :
1. défaut de step5 ;
2. corps vérifié avec la vraie tête ;
3. validation sur les autres modèles (death, l_shape, l-triple, turkish,
   Squirtle, real_overhang) ;
4. défauts connus :
   - points au-delà de 45° jugés au centre des faces ;
   - plateau traité comme un plan infini ;
   - P3 couche par couche ;
5. essai d'impression pour fixer les réglages provisoires ;
6. mise en forme pour Apex : étapes `run(entrée, paramètres)`, dataclasses,
   tests pytest, pyvista exclu.

L'export G-code s'écrira directement dans Apex (`file_io/gcode_export.py`),
avec la compensation du pivot de 45 mm et le ralentissement de l'avance.

Machine (voir la mémoire `machine-kinematics`) : R, C, Z et B, C imposé par
la position ; B entre 0 et 185°, 750 °/s et 300 °/s² ; firmware à axes
coordonnés. Plateau en carton, un seul tube coudé réussi à ce jour.

**>>> P5 TERMINÉ (boucles relancées de zéro, angles exacts, limites de B) <<<**

| modèle | tours | zones gardées | surface lissée | collisions | temps des passes de lissage (25 mm/s) |
|---|---|---:|---:|---:|---|
| Benchy | 3 + étape finale (4 zones retirées, 73 mm²) | 29 | 1 955 mm² | 0 / 2 910 points | 3,6 → 7,3 min (×2,01) |
| dino | 3 (convergé sans étape finale) | 30 | 2 932 mm² | 0 | 5,3 → 11,3 min (×2,15) |

- Écart tête / normale : médiane 28 à 31°, maximum 45° (la tolérance).
- Suite Apex : **385/385 dans les deux positions**.
- Images : `output_models/loop_*_skin_paths_final_pass*.png`.
- Ordres d'impression : `output_models/loop_*_print_order_pass*.pkl`.
  Chaque passe de lissage porte ses points, ses normales, ses angles exacts
  et son facteur de vitesse par segment.

**Ancien résultat sur le Benchy (avant les limites de B) :**

| tour | points de collision | zones | surface lissée |
|---|---:|---:|---:|
| 1 | 274 | 50 | 2 875 mm² |
| 2 | 20 | 35 | — |
| 3 | 11 | 33 | 2 200 mm² |

- Tour final : les 4 zones restantes (73 mm²) sont retirées, P3 est refait
  sans leur creusement (`--skin-exclude`), et P5 vérifie. Résultat :
  **0 point sur 3 153, 29 zones, 2 127 mm²**. Les contours du corps sont à
  0,200 mm sous la surface.
- Image :
  `output_models/loop_tilted_benchy_undertriangulated_repaired_pin_skin_paths_final_pass4.png`.
- Ordre d'impression final :
  `output_models/loop_tilted_benchy_undertriangulated_repaired_pin_print_order_pass4.pkl`.
- **À surveiller.** Dans 9 zones, quelques points dépassent la tolérance de
  45° (48° au plus). P2 juge chaque face à son centre, alors que la
  trajectoire suit la normale moyennée aux sommets. Remède possible : un
  contrôle par point dans P2, ou un rognage par P5 comme pour les
  collisions.
- **Non testé de bout en bout** : le tour final a été ajouté dans
  `skin_loop.py` après coup ; je l'ai exécuté à la main sur le Benchy.
- **Reste pour P5 et P6** :
  - le dino ;
  - le lissage temporel de B entre points (`_clamp_rate`) ;
  - le GIF de l'ordre d'impression dans la pièce ;
  - les étapes 7, 8 et 9 sur la trajectoire complète (corps + lissage) ;
  - le G-code.

## >>> MATINÉE AUTONOME DU 2026-09-28 : P4 FAIT (tâches A à E, fini vers 8 h 30) <<<

Plan suivi : section « MATINÉE AUTONOME » de
`C:\Users\test_\.claude\plans\snug-yawning-lovelace.md`. Tout est en recherche
(`research/S4_Slicer/`, local). Apex n'a pas été touché, et il n'y a eu ni
commit ni push. Toutes les décisions prises seul sont dans
`research/S4_Slicer/DECISIONS_2026-09-28.md` (résumé plus bas).

### Ce qui est fait

**Les 10 motifs existent** (`surface_paths.PATTERNS`, mêmes noms que
`surface_patterns`) :
- concentrique conforme ;
- concentrique (cercles) ;
- rectiligne ;
- rectiligne aligné ;
- ligne monotone ;
- concentrique carré ;
- concentrique triangulaire ;
- spirale d'Archimède ;
- spirale octagramme ;
- Hilbert.

**Les 9 motifs tracés directement sur la surface sont conformes dans les 54 cas
du banc** (`test_patterns.py --real`). Le banc croise 6 formes : disque,
demi-cylindre, demi-sphère, coque du Benchy, zone plane du Benchy et zone du
dino. Résultats :
- espacement médian de 0,398 à 0,400 mm, entre les 5e et 95e centiles dans ±5 %
  (coins vifs exclus du critère mais affichés : jusqu'à 0,52 mm dans un coin
  droit du conforme, effet géométrique de tout décalage de contour) ;
- couverture de 95,7 à 99,9 % ;
- aucun point hors zone ;
- 0,0 à 15 s par motif et par zone, les spirales étant les plus lentes (log map).

Le banc 4a est toujours à 4/4.

**Carré, triangle et octagramme** : ils sont construits comme le maximum ou le
minimum de distances signées, mesurées sur la surface, à des lignes qui passent
par le centre. Chaque côté est donc une famille de parallèles exactement
espacées, et les coins restent vifs.

**Spirales** : ligne de niveau de « rayon − w·θ/2π ». Le maillage est ouvert
le long de la coupure de l'angle θ, puis les tours sont recollés.

**Ordre de dépôt revu après la première planche** :
- rectiligne par bandes (coque : 107 → 17 passes) ;
- anneaux au plus proche voisin ;
- spirales dans l'ordre de la spirale.

Deux remarques de l'utilisateur, faites pendant la matinée, ont été prises en
compte :
- la ligne monotone est revenue à l'ordre ligne après ligne ;
- les arcs de la spirale octagramme sont triés sur la valeur exacte du rayon
  (trajets à vide sur le disque : 386 → 147 mm).

**Chaîne complète** : `skin_paths.py <modèle> <passe>`. Options :
`--pattern=`, `--zone-pattern=z:motif,…`, `--line-angle=`,
`--zone-angle=z:deg,…`, `--overlap=`, `--thickness=`, `--concentric-order=`,
`--zones=`. Pour chaque zone, elle trace le motif, descend la trajectoire de
t/2 sous la surface, rééchantillonne à 0,5 mm au plus, donne la normale en
chaque point et un facteur d'extrusion de 1. Sortie :
`loop_<m>_skin_paths_pass<k>.pkl`.

| modèle | motif | zones tracées | couverture pondérée | profondeur t/2 ± 0,05 | durée |
|---|---|---:|---:|---:|---:|
| Benchy passe 4 | rectiligne | 128/128 | 99,4 % | 100 % | 11 s |
| dino passe 5 | rectiligne | 93/93 | 99,6 % | 100 % | 13 s |
| Benchy passe 4 | chacun des 9 autres motifs | 128/128 | 97,0 à 99,4 % | 100 % | 10 à 30 s |

Onze micro-zones du Benchy (moins de 4 mm² chacune) restent sous 90 % de
couverture. Elles sont plus étroites que quelques buses.

### Images à regarder

Toutes sont dans `research/S4_Slicer/output_models/` :
- `patterns_bench.png` : planche de 6 formes × 10 motifs. La couleur donne
  l'ordre de dépôt (violet d'abord, jaune ensuite), les traits rouges sont les
  trajets à vide ;
- `patterns_spirals_hilbert.png` : zoom sur les deux spirales et Hilbert ;
- `loop_tilted_benchy_undertriangulated_repaired_pin_skin_paths_pass4.png` et
  `loop_dino_pin_skin_paths_pass5.png` : toutes les couches de lissage dans la
  pièce, colorées par moment d'impression.

### Hilbert : mesure de l'écart (`test_hilbert.py --real`), DÉCISION DE MIDI

Méthode : dépliage ARAP (libigl), courbe au pas w dans le plan, puis retour sur
la surface. L'espacement est mesuré sur la surface entre passes voisines. Sur
un plan, la mesure elle-même ne dépasse pas 94 % dans ±5 %, à cause des
demi-tours de la courbe : c'est la référence.

| zone | surface étirée > 5 % | espacement réel (5e-95e) | dans ±5 % | couverture |
|---|---:|---:|---:|---:|
| disque plat (référence) | 0 % | 0,400-0,422 | 93,9 % | 95,5 % |
| demi-cylindre (développable) | 0 % | 0,400-0,422 | 94,0 % | 95,2 % |
| zone plane du Benchy | 0 % | 0,400-0,422 | 94,3 % | 94,7 % |
| coque du Benchy | 19,9 % | 0,385-0,425 | 87,6 % | 94,6 % |
| zone du dino | 77,8 % | 0,356-0,446 | 52,6 % | 93,4 % |
| demi-sphère | 100 % | 0,324-0,461 | 26,5 % | 91,1 % |

Ces chiffres ne comportent AUCUN remède. Les options, à choisir par
l'utilisateur :
1. accepter tel quel (bon sur les pièces mécaniques et le Benchy) ;
2. découper les zones trop bombées pour Hilbert (seuil d'étirement) ;
3. corriger le pas localement : grille resserrée là où la surface s'étire, ce
   qui déforme le motif ;
4. basculer ces zones sur un motif tracé directement.

### Décisions prises seul (détail et mesures dans `DECISIONS_2026-09-28.md`)

Toutes sont réversibles.
1. **Graine du rectiligne décalée de 1,2e-4 mm.** Un plan passant exactement
   par une rangée de sommets fragmentait la graine.
2. **Repli fast marching** si la chaleur signée échoue.
3. **Extraction des lignes de niveau vectorisée** : ×9, résultat identique à
   1e-15 mm.
4. **Cercles tracés à partir d'un petit cercle de 0,8 mm plus une distance**,
   au lieu du rayon du log map (celui-ci se dégradait sur la coque).
5. **Graines nettoyées** : un point à moins de 1 µm d'un sommet est ramené sur
   ce sommet, et le repli est étendu aux distances non finies.
6. **Méthode de distance choisie selon la qualité des triangles** : fast
   marching si la qualité médiane est d'au moins 0,5, chaleur signée sinon.
7. **Critères de mesure** :
   - niveaux voisins seulement, à l'intérieur de la zone ;
   - comparaison à la ligne voisine entière ;
   - coins vifs exclus du critère mais affichés.
8. **Ordre de dépôt** :
   - rectiligne par bandes ;
   - anneaux au plus proche voisin ;
   - ligne monotone ligne après ligne (à la demande de l'utilisateur) ;
   - spirales dans l'ordre de la spirale, sur la valeur exacte du rayon.
9. **Rectiligne aligné** : plan à un multiple entier de w de l'origine machine,
   au plus près du centre de la zone.
10. **Carré, triangle et octagramme** construits comme ci-dessus, orientés selon
    l'angle de la zone (défaut 45°) ; octagramme = deux carrés tournés de 45°.
11. **Spirales** : maillage ouvert le long de la coupure, centre (rayon < w)
    laissé vide.
12. **libigl installé** (MPL-2.0) pour l'ARAP, avec un sommet épinglé.
13. **Sommets en nœud papillon dupliqués** dans `ZoneMesh` (zone 19 du Benchy).
14. **`skin_paths.py`** :
    - normale de descente moyennée par l'aire des faces de la zone ;
    - épaisseur par défaut 0,2 mm ;
    - une zone en échec est signalée, jamais ignorée ;
    - l'option « lissage des facettes » n'est pas encore branchée ici.

### Questions de midi

1. **Hilbert** : quel remède (options 1 à 4 ci-dessus) ?
2. **Validation visuelle** des planches de motifs et des deux rendus de pièce.
3. **Anneau de contour au bord des zones** : faut-il tracer un cordon qui suit
   le bord de chaque zone avant le motif, comme la paroi d'un slicer plan ?
   Sans lui, les motifs à grille ou à lignes laissent une frange en dents de
   scie le long du bord (visible sur Hilbert et sur le rectiligne en biais).
4. **Micro-zones** (moins de 4 mm², plus étroites que quelques buses) : faut-il
   les imprimer quand même, les fusionner avec une voisine, ou ne pas les
   lisser ?
5. **Suite** : P5 (insertion dans l'ordre d'impression, puis vérification de
   collision avec la vraie tête à l'étape 9), ou d'abord la priorité Maker
   Faire (export G-code réel, en commençant par relire
   `scripts/gen_*_gcode.py`) ?

### Préparation G-code (lecture seule, rien modifié) — pour la question 5

Relu dans `scripts/gen_l_overhang_gcode.py`, `scripts/gen_elbow_tube_gcode.py`
et `out/elbow_tube_x1.gcode` (27/09, 20 h 31). Ce sont des générateurs
AUTONOMES pour une forme fixe : un L et un tube coudé. Ils n'utilisent pas le
pipeline. En revanche, ils fixent déjà les conventions de la machine :

- **format** : `G00 Y<C> X<R> Z<Z> B<B> [E<e absolu>] F<…>`, avec en tête
  `G21 G90 M82 G92 E0 M109 S195` ;
- **convention de tête**, vérifiée avec l'utilisateur d'après le script :
  - `B > 0` penche la buse vers le centre du plateau (vers −R) ;
  - la tête ne se déplace que dans un demi-plan vertical (R, Z) : C est imposé
    par le point, `C = atan2(y, x)` ;
  - une pièce centrée sur l'axe fait faire à C des tours complets (d'où le
    décalage de pièce `PART_OFFSET_R` dans le tube) ;
- **pivot de tête** : le dernier G-code a été généré avec
  « pivot B → pointe : 45 mm », donc SANS compensation par le firmware.
  L'export du pipeline devra compenser ce décalage. Valeur à confirmer par
  l'utilisateur ;
- **vitesse** : le firmware traite Y (en degrés) comme un axe linéaire. F est
  donc recalculé par segment pour garder une vitesse matière constante
  (`_feed_for`), avec un repli `CONSTANT_FEED` ;
- **à vérifier face à Apex** : que la cinématique de l'étape 8
  (`rotation_axes`, B perpendiculaire à R) donne le même signe de B que cette
  convention. Depuis la limite #58, une couche plate donne `B ≈ 0`, mais le
  SIGNE de B pour une couche penchée vers +R n'a pas été confronté au
  script ;
- **question pour l'utilisateur** : `elbow_tube_x1.gcode` a-t-il été imprimé,
  et avec quel résultat ?

## >>> REPRENDRE ICI (fin de soirée du 2026-09-27) : P4, étape 4b <<<

Modèle : `tilted benchy undertriangulated_repaired`, `_pin`, passe 4 (et dino,
`_pin`, passe 5). Tout est en recherche (`research/S4_Slicer/`, local). Le plan
à jour (P2-bis, P4 détaillé, décisions) est dans
`C:\Users\test_\.claude\plans\snug-yawning-lovelace.md`. Toutes les décisions
P4 de l'utilisateur y sont : ne pas les redemander.

**Fait ce soir :**
- **P3** validé : sous-peau à 0,200 mm, 0,2 % « autre ». Défaut trouvé et
  corrigé : la fermeture `buffer(+e).buffer(-e)` en onglet bouchait de vrais
  trous (couche 157) → `shapely.set_precision(region, 1e-6)`. Les « trous »
  perdus des couches 74 et 140 étaient des îlots séparés par des micro-fentes,
  correctement réunis. `check_carve.py` contrôle maintenant « creuser n'ajoute
  jamais de matière ». Durée 14 min : refonte en découpe 3D prévue (roadmap
  CLAUDE.md). P3 est à relancer après P2-bis, mais P4 n'en dépend pas.
- **Défaut de l'étape 5 d'Apex découvert, NON corrigé** : sur le tranchage de
  référence, la couche 6 perd 100 % de sa matière et la couche 174 23 %
  (contour de l'étape 5 contre section réelle des tétraèdres, dans l'espace
  déformé). Cause inconnue (hypothèse : plan passant par des sommets). À
  traiter avant toute impression.
- **P2-bis fait** (`surface_skin.py` v2 ; v1 conservée dans
  `surface_skin_v1.py` = oracle) :
  - dernier moment accessible par face, placement de tête groupé, IK par
    normale distincte ;
  - nombre MINIMAL de bandes (couverture d'intervalles) ;
  - frontières lisses (coupe exacte le long de phi = isovaleur) ;
  - coupe aux arêtes à 150° entre faces ; `--facet-smoothing`.

  Mesures :
  | modèle | durée | écart de classement face par face / v1 | zones |
  |---|---|---|---|
  | Benchy | 11 min 30 → **57 s** | aucun (hors « zone trop petite ») | 264 → 128 |
  | dino | 8 min → **27 s** | **0 face** | 167 → 93 |

  La coque du Benchy est UNE seule zone (la v1 la coupait en ~15 bandes). La
  collision fait 95 % du temps restant.
- **Mesures préalables P4** :
  - `measure_unfold.py` : déformation du dépliage. Sans coupe aux arêtes,
    jusqu'à 72 % de la surface déformée de plus de 10 % ; avec la coupe, ~0
    sauf sur le dino (8,6 %, formes bombées) ;
  - `measure_seams.py` / `render_seams.py` : les raccords aux arêtes tombent
    à 72 % sur de vraies arêtes vives.
- **potpourri3d 1.4.0 installé dans `.venv` d'Apex** (MIT ; numpy, scipy).
- **4a fait** : `surface_paths.py` (ZoneMesh, distances, log map, lignes de
  niveau) et `test_surface_paths.py` : **4/4 conformes**. Choix mesurés :
  - distance à une courbe (bord, graine fermée ou traversant la zone) =
    chaleur SIGNÉE puis valeur absolue. Position à 0,001-0,026 mm (95e). La
    chaleur simple décale les lignes (0,05 à 0,26 mm) ; le fast marching se
    dégrade sur les triangles obtus ;
  - distance au centre = rayon du log map (0,000-0,011 mm) ;
  - maillage raffiné à 0,2 mm (espacement ±2 %).
- **4b en cours** : `concentric_conform`, `concentric_circles` et `rectilinear`
  (zigzag avec raccords le long du premier anneau) dans `surface_paths.py` ;
  banc `test_patterns_4b.py [--real]` (vraies zones : coque du Benchy, plus
  grande zone du dino, zone la plus plane du Benchy). Premiers résultats :
  - disque : les 3 motifs OK (espacement 0,392-0,407, couverture 96-99 %,
    0 % hors zone) ;
  - demi-cylindre : conforme et cercles OK.

**Problèmes ouverts, à traiter en premier demain :**
1. `rectilinear` sur le demi-cylindre : `ValueError: Solver factorization
   failed` (chaleur signée, graine = coupe par un plan). À diagnostiquer :
   graine touchant le bord ? plusieurs morceaux ? Repli possible : fast
   marching pour cette graine (exact sur maillage correct).
2. **Durée** : cercles sur le demi-cylindre en 24 s (187 000 sommets après
   raffinage à 0,2 mm ; le log map domine). Mesurer sur les vraies zones ;
   pistes : raffinage à 0,3 mm (espacement ±5 % mesuré à 0,38 mm), ou
   raffinage limité aux triangles trop grands.
3. Ensuite : finir le banc 4b (`--real`), faire un rendu (script à écrire,
   pyvista ou matplotlib, venv de la recherche), **point d'arrêt utilisateur**
   (montrer les trois motifs sur les formes et les vraies zones).
   **Question à poser à ce point d'arrêt** : pour les motifs autres que le
   conforme, faut-il un anneau de contour au bord de chaque zone, comme les
   parois d'un slicer ?

**Après P4 (ordre convenu, question de l'utilisateur « et ensuite step 8, step 9 ? ») :**
1. **P5** : étapes 8 (IK, lissage des angles `_clamp_rate`) et 9 (tête réelle)
   sur les couches de lissage ; retrait des zones qui échouent.
2. **Toute la trajectoire Bird** (corps et lissage) dans les étapes 7, 8, 8b,
   9 et 10. D'après les notes, elle n'y est pas encore passée. Remesurer les
   collisions avec la tête corrigée (limite #58 : toutes les collisions
   antérieures sont suspectes). Points connus :
   - étape 8b : 60 % de résidu sur les pièces difficiles ;
   - étape 9 : lente ;
   - étape 10 : support de secours, à éviter.
3. **Export G-code** (`file_io/gcode_export.py` inexistant ; relire d'abord
   `scripts/gen_*_gcode.py`, `preview_gcode.py`), puis un essai d'impression.
4. **Avant toute impression** : corriger le défaut de l'étape 5 (couches 6 et
   174).

---

## Reprise initiale de la session (conservée pour mémoire)

Tests sur **`tilted benchy undertriangulated_repaired`**, base figée (`_pin`),
**passe 4**. Voir §6 : l'ordre exact et les précautions mémoire,
**obligatoires** (l'ordinateur a planté, §5).

---

## 1. Décisions de l'utilisateur (à ne pas redemander)

- Développement **d'abord en recherche** (`research/S4_Slicer/`, local, JAMAIS
  versionné), portage dans Apex une fois validé.
- **Tranchage dynamique** : aucune constante propre à un modèle. Tout seuil est
  relatif (buse 0,4, couche 0,2, épaisseur de lissage, taille de la pièce), et
  chaque phase est mesurée sur plusieurs modèles.
- **Base plane** = tout le volume sous le surplomb le plus bas (au moins une
  couche). Variante « contact + marge » mesurée : aucune ne gagne partout. J'ai
  recommandé la **base figée entière**, l'utilisateur a répondu « Continuer »,
  donc c'est le défaut actuel.
- **Lissage non-planaire** = une VRAIE couche de matière qui épouse chaque
  surface extérieure accessible : dessus, flancs, dessous, y compris une
  surface accessible seulement PENDANT l'impression (intérieur du toit du
  Benchy). Tête **inclinée selon la normale** (buse vers −n). Moment : **au plus
  tôt**. Corps creusé d'autant sous la couche.
- **Épaisseur** : paramètre « épaisseur des couches de lissage non-planaire »
  (`--skin-thickness`, 0,2 mm par défaut).
- **Ironing** : pas maintenant. Il faut en préparer les bases (facteur
  d'extrusion par passe de surface) et l'inscrire dans la feuille de route.
- **Pré-calcul obligatoire** qui écarte les collisions évidentes tête-plateau
  et tête-objet avant l'étape 9, pour ne pas la saturer.
- **Plateau = plan infini** dans le pré-calcul. Correction future inscrite
  dans la feuille de route de `CLAUDE.md` : refuser dès l'entrée un modèle trop
  grand pour la machine.
- **10 motifs de surface dès la v1**, réglables :
  - Concentrique (cercles imbriqués) ;
  - Rectiligne ;
  - Ligne monotone ;
  - Rectiligne aligné ;
  - Hilbert ;
  - Spirale d'Archimède ;
  - Concentrique carré (carrés imbriqués) ;
  - Conforme aux parois (contours de parois imbriqués jusqu'au centre) ;
  - Concentrique triangulaire (triangles imbriqués) ;
  - Spirale octagramme.

  Le **motif de première couche** est un paramètre **sans valeur par défaut**.
  La frange non couverte de Hilbert (~95 % de couverture) est **acceptée telle
  quelle**.
- **Autres zones sans surplomb à plat** (hors base) : à décider après mesure,
  PAS dans ce plan.
- **Aplatir le dessous** des modèles à l'entrée, pour tous les modèles (fait,
  §3). L'aplatissement de TOUTES les zones presque planes a été demandé, mesuré
  NÉGATIF et retiré (§3). Consigne : « une fenêtre, pas une maison à côté »,
  donc pas de grosse dépendance.
- En cas de doute : **POSER LA QUESTION**, ne pas deviner (consigne explicite).
- Le code de recherche reste local (mémoire `research-code-stays-local`).

## 2. Apex (versionné)

- **Commit `651c95b`, NON POUSSÉ** : correction de l'étape 8. La tête était
  retournée sur toutes les couches (référence de l'IK −Z au lieu de +Z, alors
  que step3/step7 et le modèle de tête sont en +Z), avec B = 179,9° pour une
  couche plate. `step8` et `step10` corrigés, 2 tests d'orientation PHYSIQUE
  ajoutés, qui échouaient bien avant la correction. **377/377** dans les deux
  positions de `APEX_USE_CPP`.
- **NON COMMITÉS** (à commiter ensemble quand l'utilisateur le demandera) :
  - `pipeline/step6_toolpath2d.py` + `tests/test_step6_toolpath2d.py` :
    - motif `CROSSED_ZIGZAG` par défaut, un zigzag continu par zone qui
      alterne entre u et v selon le RANG de l'isovaleur ;
    - raccords qui suivent le bord (`_boundary_connector`) ;
    - zone décalée d'une demi-buse, pour ne pas redéposer la paroi ;
    - `GRID` réellement branché.
  - `CLAUDE.md` : suite de la limite #57 (motifs), nouvelle limite #58 (étape
    8), note de feuille de route « refuser un modèle trop grand », compteur
    377/377.
- Référence step9b (0,466578 %) : **obsolète**, à remesurer (étapes 6 et 8
  modifiées).

## 3. Recherche (`research/S4_Slicer/`, local) — ce qui existe

| fichier | rôle | état |
|---|---|---|
| `bird_core.py` | `planar_base_mask(tet, overhang_deg, layer_height, max_height=None)` ; `optimize_rotations(..., fixed_zero_cells)` et `calculate_deformation(..., pinned)` (contraintes DURES, par retrait des inconnues) ; `run_one_iteration(..., planar_base=)` ; `load_and_tetrahedralize(..., flatten_bottom=)` | fait |
| `bird_loop.py` | base figée par défaut (suffixe `_pin`), `--no-pin` (référence), `--base-contact` (`_pinc`, 3 couches) ; dessous aplati (`flatten_bottom=0.2`) ; modèles pris dans `input_models/` puis `res/test_models/` ; sauve `_base_mask.npy` | fait |
| `mesh_flatten.py` | `flatten_bottom` : sommets de la 1re couche ramenés sur le plateau, SAUF si une face entièrement écrasée ne regardait pas vers le bas (> 60° de −Z) ou si une face se retournerait (itératif). Maillage rendu tel quel si rien ne bouge | fait, validé sur 9 modèles |
| `measure_pass.py <prefixe> [passe]` | mesures d'une passe : écart de la face du dessous, couches qui la coupent, plus grand repli, surplomb réel (tête inclinée) | fait, reproduit les chiffres connus |
| `surface_patterns.py` | les 10 motifs en 2D, `generate(pattern, zone, spacing, angle)` | fait, 39/40 cas conformes (`test_surface_patterns.py`) |
| `surface_skin.py <modele> <passe> [variante]` | **P2** : surface subdivisée (arête ≤ 2 mm) ; moment au plus tôt par rang de couche ; IK de l'étape 8 avec buse vers −n ; pré-calcul (plateau z<0, objet = matière déjà imprimée = pièce ∩ φ ≤ φ du moment, `PrintedSolid`) ; zones découpées en bandes de moment si bloquées. Sortie `_skin_pass<k>.npz` + `_skin_map_pass<k>.ply` | fait |
| `render_skin_map.py <ply>` | carte colorée des statuts (vert retenue, orange plateau, rouge objet, gris contact) | fait |
| `skin_carve.py` | **P3** : `skin_psi`, `_sections` (sens canonique des arêtes, obligatoire), `_clip`, `ExactDistance`, `carved_region`, `region_loops`, `first_layer_toolpath` | **À CORRIGER AVANT TOUT USAGE, voir §5** |
| `bird_slice_deformed.py <m> <k> [variante] [--skin] [--skin-thickness=] [--first-layer=<motif>]` | tranchage déformé ; avec `--skin`, contours des couches touchées recalculés par `skin_carve`, sortie `_deformed_slice_skin_pass<k>.pkl` | fait |
| `check_carve.py <m> <k> <variante> [ep.]` | vérification de P3 : points de contour classés surface / sous-peau / autre | fait |
| `render_part_gif.py <pkl> <stl>` | GIF pyvista dans la pièce transparente | fait |
| `test_mesh_flatten.py` | banc de l'aplatissement du dessous (venv d'Apex) | fait |

Modèles réparés par pymeshfix : `input_models/death_repaired.stl` (−2,4 % de
volume) et **`input_models/tilted benchy undertriangulated_repaired.stl`**
(−0,05 %). l_shape et l-triple exigent le maillage fin (`0.5`), sinon ils n'ont
que 12 tétraèdres dans la base.

## 4. Mesures acquises

**P1, surplomb réel** (plus grand repli entre parenthèses, en mm²) :

| modèle | libre | base figée | contact + marge |
|---|---:|---:|---:|
| Benchy | 6,78 % (14) | 6,91 % (14) | = figée |
| dino | 0,31 % (15) | 0,40 % (22) | = figée |
| Squirtle | 6,43 % (30) | 6,31 % (21) | = figée |
| death | 13,33 % (22) | 13,85 % (29) | 13,40 % (22) |
| real_overhang | 4,31 % (10) | 5,17 % (18) | 4,68 % (25) |
| l_shape (0.5) | 4,41 % | 4,35 % | 4,28 % |
| l-triple (0.5) | 7,59 % | 7,47 % | 7,75 % |
| **tilted benchy undertri. (dessous aplati)** | 5,71 % (20), base 3,45 mm / 18 couches | **5,80 % (20), base 0 / 1 couche** | — |

La face du dessous est plate (0,00 mm, 1 couche) sur tous les modèles figés.

**Découvertes :**
- La boucle d'origine laissait la pièce DÉRIVER, faute d'ancrage. C'est
  probablement l'« anomalie du recentrage ». Le critère « ne bouge plus »
  (25 %) avait été calé avec cette dérive, il est à revérifier.

**P2 :**

| modèle | faces retenues | collision plateau | collision objet | zones | durée |
|---|---:|---:|---:|---:|---:|
| Benchy (`_pin`) | 50,2 % | 31,5 % | 13 % | 161 | 13 min |
| dino | 50,7 % | — | — | — | 8 min |
| real_overhang | 62,7 % | — | — | — | — |
| l_shape | 48,1 % | — | — | — | — |

Défauts connus de P2 :
- **lenteur** : `collides` est rappelé récursivement dans `assign`. À optimiser
  en mettant en cache le statut par (face, moment) ;
- **rayures** orange et vert sur les faces verticales basses : c'est du bruit
  de seuil sur une surface facettée. Ne PAS « accepter » par vote majoritaire
  (ce serait dangereux) ; plutôt tester plusieurs points par face.

**P3** (version par interpolation aux nœuds, avant l'exactitude) :

| modèle | sous-peau | autre | remarque |
|---|---|---:|---|
| Benchy | 43,5 % à **0,200 mm** exactement | 14,8 % | max 1,9 mm |
| dino | — | 22 % | maillage grossier |

Le reste « autre » a deux causes :
- l'interpolation de ψ dans de grands tétraèdres ;
- les **parois minces lissées des deux côtés**, dont l'âme disparaît parce que
  tous leurs nœuds sont « dans la peau ».

D'où la version exacte (§5).

## 5. INCIDENT MÉMOIRE (l'ordinateur a planté, 47 Go saturés) — précautions obligatoires

**Cause** : `skin_carve.carved_region` avec `ExactDistance` est actuellement en
version à **subdivision UNIFORME**. Chaque section proche de la peau est
découpée en triangles de 0,2 mm, puis les distances sont calculées en UN
appel `trimesh.proximity.closest_point` par couche, ce qui fait des millions
de points. Lancée sur dino, elle a été « arrêtée » par TaskStop. Sous Windows,
cela tue le shell mais **pas le processus Python enfant**, qui a continué
jusqu'à saturer la mémoire.

**Avant toute reprise :**
1. Vérifier et tuer les Python orphelins : `Get-CimInstance Win32_Process
   -Filter "Name='python.exe'"`, filtrer sur `CommandLine`, puis `Stop-Process`.
2. **Remplacer** la boucle uniforme de `carved_region` (bloc `if fine_tri:`)
   par la version **ADAPTATIVE**, qui avait été écrite puis refusée à cause du
   plantage :
   - `ExactDistance.psi` rend `(psi, points_reels)` ;
   - la fonction `_split4(tri)` remplace `_subdivide` ;
   - boucle sur au plus 12 profondeurs : calculer ψ exact aux sommets et la
     plus grande arête RÉELLE de chaque triangle. Si `min ψ > arête`, garder
     le triangle entier. Si `max ψ < −arête`, l'écarter. Si l'arête déformée
     est ≤ `step`, découper linéairement (`_clip`). Sinon, `_split4`.
   - Justification : la distance est 1-lipschitzienne, donc le travail suit
     la longueur de la frontière et non l'aire.
3. Calculer les distances par **paquets bornés** (par exemple 50 000 points
   par appel à `closest_point`), jamais un seul appel géant.
4. Tester d'abord sur **UNE couche** (option à ajouter, ou appel direct), en
   mesurant temps et mémoire, avant de lancer toutes les couches.
5. Arrêter une tâche = tuer aussi ses Python enfants. Ne pas lancer de calculs
   lourds en parallèle sans vérifier la mémoire libre.

## 6. Prochaines étapes (ordre convenu)

Sur `tilted benchy undertriangulated_repaired`, variante `_pin`, passe **4**
(les fichiers `loop_tilted_benchy_undertriangulated_repaired_pin_*` existent) :

1. Tranchage de référence :
   `APEX_USE_CPP=1 ../../.venv/Scripts/python.exe bird_slice_deformed.py "tilted benchy undertriangulated_repaired" 4 _pin`
2. P2 : `surface_skin.py "tilted benchy undertriangulated_repaired" 4 _pin`,
   puis `render_skin_map.py` sur le `.ply`. Montrer la carte à l'utilisateur.
3. P3 (APRÈS le §5) : `bird_slice_deformed.py ... 4 _pin --skin
   --first-layer=<motif demandé à l'utilisateur>`, puis `check_carve.py`.
   Cible : sous-peau à 0,2 mm ± 5 %, part « autre » réduite, aucun point hors
   de la pièce.
4. Point d'arrêt utilisateur (cartes + chiffres).
5. **P4** : dépliage des zones (harmonique ou LSCM, `scipy.sparse`, découpage
   si distorsion trop forte), les 10 motifs appliqués sur la zone dépliée puis
   ramenés sur la surface à −n·t/2, polylignes redécoupées à 0,5 mm, facteur
   d'extrusion par passe (préparation de l'ironing).
6. **P5** : insertion au plus tôt dans l'ordre d'impression, angles
   (normale), lissage avec `step8._clamp_rate`, vérification étape 9 avec
   `aetheris-4.glb`. Une zone qui échoue est retirée, son creusement annulé,
   et on retranche une fois.
7. **P6** : GIF avec les passes de surface en couleur distincte, feuille de
   route (ironing, autres zones planes, portage Apex), mémoire.

**Autres points ouverts :**
- commit du zigzag et de `CLAUDE.md` ;
- push de `651c95b` (seulement quand l'utilisateur le demande) ;
- remesure step9b ;
- rayon du plateau et refus des modèles trop grands (feuille de route) ;
- Bird n'améliorait pas l-triple en maillage grossier (réglé en fin).
