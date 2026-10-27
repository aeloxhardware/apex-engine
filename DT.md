# Dev Tracker (DT) — Apex Slicer

État au **2026-10-13, soir** (chantier « Bird en C++ puis Vulkan » démarré : étape −1, phase 0, phase 1 ; voir #77). Convention de rédaction : section « Dev Tracker » de `CLAUDE.md`
(🟢 résolu · ⏳ en cours · ⬛ non traité, sans successeur · ♦️ à résoudre avant de continuer).
Détail technique : `CLAUDE.md` (limites #56 à #75) ; récapitulatif de session : `LAST_SESSION.md`.
**Tests : 505/505 dans les positions `APEX_USE_CPP=0`, `=1` et `=1` + `APEX_USE_GPU=1`** (2026-10-13 soir, suite complète).


## 2026-10-27 — Retrait de FEM et de Bird termine ; collisions de tete mesurees
- 🟢 Voie unique = chaleur. Retires : champ FEM directionnel, multibranche, parties planes, `Constraint[]`, options de step3/step5 liees, copie de recherche S4 Slicer (citation et licence gardees). Contrat d'interface passe a 2 (`CONTRACT_CHANGES.md`). Tests : 348 passes en Python pur, suite native verte (Linux, g++).
- 🟢 Decisions de l'utilisateur : essai a sec de B — B+ va vers le centre du plateau (coherent avec `b_sign=+1`). La Gen 0 de l'Aetheris-4 a ete retiree prematurement apres un accident qui l'a detruite ; elle sera remplacee par la Gen 1, sur laquelle l'essai a sec est a refaire. Placement FEM : sans objet.
- ⬛ A refaire sur la voie chaleur : points d'attache (contraintes utilisateur). Code C++ mort restant : melange de directions et terme de gradient dans `field_solver.cpp`.
- 🟢 Collisions de tete : `CollisionSimParams.printed_part_margin` (defaut -0,15 mm dans `main`, matiere deja deposee seulement ; le plateau garde 1 mm). Decision de l'utilisateur : le modele de tete est exact, effleurer ou mordre un peu la couche n'est pas un choc. Mesure : sur deux G-code (cartesien plat, Benchy polaire) tous les contacts a marge 0 recouvrent la matiere de 0,18 mm au plus ; Benchy 187 -> 48 points a -0,15 (1 a -0,18), plat 28 -> 3. death (1,2 mm) : 7 840 collisions -> 741, supports de secours 7 805 -> 679. Valeur a regler par l'utilisateur.
- (historique) Collisions de tete avant ce reglage (death, resolution 1,2, step9 rapide, part des points) : defaut 15,1 % ; avec step5b 24,1 % (pire) ; couches scindees 50,5 % ; flux non projete 27,8 %. Le defaut est donc la meilleure variante testee. Apres step10 : 0 collision non resolue mais 7 805 des 7 840 sont traitees en « support de secours » (35 reorientations seulement). Piste non mesuree : step5b par morceau de couche avec la normale locale par pose, ou reorientation de step10 sur la grille couplee.
- ⬛ Purge du docx confidentiel de l'historique : document deja dans `aeloxhardware/CONFIDENTIAL_INFO` ; purge a faire en dernier (copie d'essai, puis force-push, ids de commit changes).

## 2026-10-25 — Decision : voie principale = CHALEUR ; Bird retire ; FEM directionnel plus par defaut
- 🟢 Bird supprime (pipeline/bird, tests, scripts, noyau C++ bird_lsq/bird_locate/skin_check, artefacts out/). **Sauvegarde integrale : branche `bird_old`** (poussee, commit 79549f08). L'autodocteur de surface (utilise par step1) est extrait dans `pipeline/mesh_doctor.py`.
- 🟢 Remerciement et citation de S4 Slicer (Joshua Bird) conserves : README (Credits) et `CITATION.cff`.
- 🟢 `PipelineParams` par defaut = variante `heat_iso_reach` (chaleur, flux projete, `isovalue_order`, sans step5b).
- ⬛ A FAIRE : retirer le code FEM directionnel restant (step4 `build_static_system`, step4b multibranch, step4c planar parts, options geodesiques/feedback, `Constraint[]` utilisateur qui en dependent) — non fait faute de pouvoir mesurer la non-regression ; les tests step4 couvrent encore ce code.
- Tests : echecs deja presents hors de ce chantier dans cet environnement Linux/Python 3.13 : test_shape_gradients (nb de tetra), 5 tests gcode natif identique, test_branch_indicator.

## 2026-10-26 — Portage C++ de ce qui restait en Python (détail : `CLAUDE.md` #82)
- 🟢 `config.g` modifié par l'utilisateur ; salon passé avec une récompense spéciale pour l'originalité du projet. ⬛ **Le Benchy n'a toujours pas été imprimé.**
- 🟢 Natif (avec test d'équivalence chacun) : méthode de la chaleur (death 6,5 → 2,6 s), CG de Jacobi parallèle reproductible, réorientation par paliers / locale / borne de hauteur (step10), contrôle rapide de collision par lots (step9 + step10 : ~12 → ~3 s), tri topologique de step5b, couches par défaut de step5, remplissage de couches entières avec ponts (step6, en parallèle), assemblage des points / normale par point / filtre médian (step8), en-tête et pied du G-code.
- 🟢 **Remplissage natif (Clipper2) ACTIF par défaut** avec le noyau natif, mesuré sur les 4 modèles : l_shape 3,80 → 4,13 %, l-triple 5,10 → 5,28 %, death 3,31 → 3,43 %, turkish 0,162 → 0,152 % de matière sans appui (bruit, réversible : `APEX_NATIVE_FILL=0`).
- 🟢 Pipeline complet death @ 0,8 : 42 → 33 s (TetGen ~14 s incompressible) ; hors TetGen 28 → 18 s.
- ♦️ **Reste en Python** : l'enchaînement des étapes (`main.py`, dataclasses `Layer`/`ToolpathPoint` : 360 000 objets sur death), step1 (trimesh/pymeshfix/TetGen en sous-processus), repli par grille de l'IK, échantillonnage de la tête (step9), boucle `run` de step10, `OSCILLATING_SPIRAL`. Un orchestrateur natif demande de remplacer le contrat inter-étapes par des tableaux (étapes 8, 8b, 9, 9b, 10, export) : chantier de la bascule finale, à décider.
- 🟢 Worker TetGen 100 % C++ (2026-10-26) : source vendorée, exécutable séparé `apex_tetgen_worker`, maillage identique au binding (death, turkish, l_shape), step1 un peu plus rapide.
- ⬛ Accélérer TetGen : 95 % du temps est dans raffinement (4,6 s) + lissage (5,5 s) sur death ; lissage parallélisable par couleurs, à mesurer sur step9b des 4 modèles avant toute activation.

## ⚠️ RAPPELS PERMANENTS — à répéter à CHAQUE réponse tant que la case n'est pas cochée (demande de l'utilisateur, 2026-10-18)
- 🟢 **[x] Paramètres du firmware (`config.g`, Z) modifiés** (confirmé par l'utilisateur, 2026-10-26). Rappel arrêté. **Salon (MFR 2026) passé : récompense spéciale pour l'originalité du projet.** ⬛ **Le Benchy n'a toujours PAS été imprimé.**
- 🟢 **[x] Champ scalaire multibranche : ESSAYÉ (2026-10-19, par Claude à la demande de l'utilisateur), verdict dans `CLAUDE.md` #78.** Résultat MITIGÉ, comme les mécanismes précédents : `multibranch=True, multibranch_direction_mode="up"` (variante `mbr_up`) divise par ~2 la matière sans appui sur death (−56 %) et benchy (−53 %), neutre sur real_overhang/dino/l_shape/l-triple, mais dégrade fig (+17 %) et deer (+23 %). Reste opt-in. Rappel arrêté.
- Contexte : MFR 2026 commence ~2026-10-23. Le lot nocturne tourne sur la voie **Bird**.

## Ce qui a été développé
- 🟢 Correctifs : parois trouées (#65), étape 8 (tête retournée, #58), étape 5 (couche perdue, #60), étape 10 (réorientation impossible en mode exact, #66), étape 3 (tri stable), trous du Benchy (#68), pointe qui quittait sa trajectoire (`max_tip_deviation`, #61), plateau de step9 à z = origine (+0,5 mm auparavant, #71), ingestion de turkish décalé (`_tetrahedralize_with_translation_retry`, #72).
- 🟢 `step5b` : ordre des couches contraint par la tête (#67). Collisions : death 41,6 → 0,13 %, turkish 53,6 → 0,21 %, sans dégradation de la matière sans appui.
- 🟢 Étape 9 rapide (arbre k-d) activée par défaut, contrôlée contre l'exact (0 faux positif, 0 faux négatif).
- 🟢 Inclinaisons atteignables (`radial`) activées par défaut (#64) : gagne sur les 4 modèles.
- 🟢 Export G-code (#61) : pivot à 46 mm, `b_sign = +1` (vérifié dans le code), pointe au sommet du cordon, dégagements contrôlés (Benchy et dino 3 → 0 en contact), cumul de Y (`y_turns`, `max_y_turns`), estimation du temps tous axes (`print_time_fw_s`, `fw_axis_speeds`).
- 🟢 Normale par point (`per_point_normal = True`) et filtre médian de B (`free_axis_median_window = 9`) activés par défaut (#69-#71). Valeurs vérifiées dans `step8_machine_toolpath.py`.
- 🟢 Limites de B réglées sur `config.g` (fichier RRF actif) : 1667 °/s et 1300 °/s² (`main._default_machine_profile`, vérifié).
- 🟢 Portage C++ : pas 1 à 15, 17, 19, 21a-c et 24.
- 🟢 Pipeline complet death (défauts, résolution 1,5) : **124 → 37 s**, résultats identiques (#73).
- 🟢 Chaîne Bird : `--apex-detect` (étape 9 : dino 855 → 74 s, Benchy 496 → 50 s, trajectoires identiques au bit) ; doublons du lissage de Squirtle 5,29 → 2,28 % (`dedupe_shared_passes`).
- 🟢 G-code régénérés avec le bon signe de B : Benchy, dino, l_shape (`out/share/*_v5.zip`).
- 🟢 **2026-10-07** — L (pivot -> pointe) devient un réglage (`file_io/pivot_length.py`, `config/pivot_to_tip.json`, sinon 46 mm), lu par `export_model_gcode.py`, `gen_b_direction_test.py` et `estimate_gcode_time.py` (#74). Le firmware ne fait aucune cinématique : G-code à régénérer après tout changement de L.
- 🟢 **2026-10-07** — Étape 10 : un écart d'orientation de l'étape 8 n'est plus compté comme support (`Step10Params.orientation_support_min_residual_deg`, défaut `None` = jamais ; `0.0` = ancien comportement). Supports death 115 312 → 209, turkish 185 658 → 1 244 (vraies collisions seules) ; l_shape et l-triple 8 107 / 5 853 → 0 (#74).
- 🟢 **2026-10-07** — Pièce du direct préparée (`out/salon/`) : **l-triple** recommandée (aucune collision, 68,5 tours de Y, 0,22 h avec le Z recommandé, 0,80 h avec le Z actuel) ; l_shape 239 tours de Y.
- 🟢 **2026-10-07** — Échec isolé de `test_export_extrusion_and_retraction` expliqué : une session `slicer_ux` éditait `gcode_export.py` pendant la suite de tests (#74). Aucun défaut de code.
- 🟢 **2026-10-13** — Pose du modèle sur le plateau à l'ingestion (`place_on_platform`, activé par défaut, #76) : death était à z = 51,4 mm et turkish à 30,3 mm du plateau. Mesuré sur les 4 modèles, neutre à l'ordre du bruit près, 0 collision plateau.
- 🟢 **2026-10-13** — G-code régénérés : salon FEM identique au bit ; Benchy et dino par la chaîne Bird → `out/share/*_v6.zip`. Temps réestimés (Z recommandé / `config.g` actuel) : Benchy 1,89 / 9,76 h, dino 2,62 / 14,90 h, l-triple 0,23 / 0,80 h, l_shape 0,20 / 1,15 h.
- 🟢 **2026-10-13** — **Production Bird basculée sur `main.py --method bird`** (étape −1, #77) : à entrée identique (passes 4 et 5 de la recherche), mêmes collisions, mêmes zones exclues, mêmes tours de Y et filament ; Benchy 580 → 250 s, dino 1 067 → 383 s ; contrôle des dégagements porté dans l'export natif (texte identique) ; `--from-pass=<préfixe>,<k>`. G-code de production **`out/share/*_v7.zip`**.
- 🟢 **2026-10-13** — Phase 0 (profil par phase, `scripts/profile_bird.py`) : Benchy 300,6 s, dino 373 s, death 714 s (#77).
- 🟢 **2026-10-13** — Gains identiques au bit (G-code du Benchy comparé octet par octet) : P5 parallèle 47,7 → 8,4 s, solveur LSQ parallèle (déformation 44,6 → 25,8 s), densification de l'export 9,7 → 0,5 s, candidats 20,3 → 16,4 s ; **Benchy complet 300,6 → 226,8 s**.
- 🟢 **2026-10-13** — Socle Vulkan (phase 1) : `--vulkan`, volk + VMA, `gpu_status()`, `APEX_USE_GPU` ; RTX 3090, float64 actif, `saxpy` identique au bit au CPU.
- 🟢 **2026-10-13** — Règle du maillage grossier de la boucle Bird **désactivée par défaut** (décision de l'utilisateur ; `BirdLoopParams.max_mean_tet_volume = None`) : elle retétraédrisait le dino de 4 263 à 151 695 tétraèdres (#77).

## Réponses de l'utilisateur
- 🟢 **2026-10-13** : règle du maillage grossier désactivée par défaut ; candidats barycentriques en C++ acceptés avec tolérance déclarée.
- 🟢 **2026-10-13** : barycentriques natifs (×20, tolérance déclarée, #77) + règle explicite d'ex aequo de `locate` + raffinement par lots (identique à l'octet) ; Benchy et dino vérifiés de bout en bout (mêmes collisions par tour et finales ; dino 158 au lieu de 157 à l'étape 9, toutes résolues). Candidats : Benchy 16,4 → 1,2 s, dino 15,2 → 3,4 s.
- 🟢 `config.g` est le fichier RRF actif ; pivot à 46 mm ; zéro Z = buse au contact du plateau ; le firmware ne fait aucune cinématique ; Y est l'axe de rotation actif ; sens physique de B (croquis : B négatif = pointe vers l'extérieur) ; RRF = imprimante bac à sable ; GPU : Vulkan ; normale par point activée.
- 🟢 **2026-10-07** : supports de l'étape 10 = vraies collisions seulement (recommandation retenue). La méthode de mesure de L est CONFIDENTIELLE (brevetable, matériel pas encore construit) : jamais publique, jamais versionnée.

## À résoudre avant de continuer
- ♦️ **Essai à sec sur la machine** (`out/test_b_direction.gcode` : pointe immobile, sans chauffe ni extrusion). Il doit confirmer que la pointe reste fixe et que B positif penche le corps vers l'extérieur. Il bloque toute impression réelle, puisque le sens physique de B n'est vérifié que sur l'aperçu. Ensuite : `out/test_l_shape.gcode` (0,28 m, 0,26 h), puis le Benchy. Les fichiers `test_b_preview*.gcode` et `test_y_ring_preview.gcode` servent à l'aperçu uniquement (ils chauffent et extrudent).
- ♦️ **Z dans `config.g` : 400 mm/min (6,7 mm/s), 70 mm/s².** Le pivot couple Z à B : avec cette valeur, le Benchy passe de 1,66 h à ~9,7 h. **Recommandé : `M203 Z8000` et `M201 Z1500`. La modification est à faire par l'utilisateur ; rappel à répéter jusqu'à confirmation.** C'est bloquant, car les temps d'impression du salon en dépendent. Pour l'état actuel, régler `fw_axis_speeds = (1666.7, 600.0, 6.67, 1666.7)`.
- ♦️ **Placement : la voie FEM garde les coordonnées du STL** (#75). death est à ~376 mm de l'axe et turkish à ~150 mm, sur un plateau de 150 mm de rayon : leurs chiffres FEM des limites #64 à #74 décrivent une pièce hors plateau. Recentrées : death 0,54 % sans appui, 42 collisions (au lieu de 0,27 % / 226) ; turkish 1,10 %, **30 012 collisions**. Ça bloque le choix de la voie par défaut et toute référence chiffrée. Pour débloquer : poser et placer le modèle à l'ingestion (déjà prévu pour `real_overhang`, posé à z = −15 mm), puis remesurer.

## Sans blocage, du plus facile au plus avancé
1. ⬛ Y cumulé à valider côté câblage : Benchy ≈ 32 tours, l_shape ≈ 240, dino ≈ 621 (`max_y_turns` permet de fixer une limite).
2. ⬛ Mesure de L : se fait à la main tant que le matériel prévu n'existe pas ; noter la valeur avec `file_io.pivot_length.record_pivot_length`. Code de la méthode future : `research/private/pivot_calibration/` (local, ignoré par git).
3. 🟢 Pose automatique du modèle sur le plateau (z min = 0) à l'ingestion : faite le 2026-10-13 (#76).
4. ⬛ Aperçu : rétrécissements sur la ligne à B > 90° (rendu du prévisualiseur, non vérifié).
5. ⬛ `resegmentation_hook` réel de step10 (laissé à `None` dans `main.py`, #23).
6. ⬛ Machines à prendre en charge (axe D, axe linéaire) : étape 8 testée pour les deux. Il manque la traduction G-code (nom et sens des axes : information à fournir par l'utilisateur), les directions de step3/step7 et step10.
7. ⬛ Reste du portage C++ : step10, step5, step4, step1 et l'orchestrateur natif ; bascule `APEX_NATIVE_FILL`. Profil de death (37 s) : TetGen ~12 s, puis des postes éparpillés sous 5 s.
8. ⬛ Backend GPU Vulkan de `MeshDistance` : aucun gain mesurable aujourd'hui (CPU : 2 M de points en 1,0 s). À rouvrir pour des centaines de millions de requêtes (pré-calcul du lissage).
9. ⬛ **Champ scalaire multibranche (FEM), mis en feuille de route par l'utilisateur le 2026-10-07.** Un champ unique ne peut pas croître depuis le plateau ET décroître sous une partie suspendue : `enforce_upward_growth` redresse toute direction descendante, et le retirer coûte death 0,23 → 11,4 % (limites #51 à #54). Idée : une branche de champ par partie portée (pilier, coque suspendue…), chacune croissant depuis SON attache, raccordées aux jonctions. Le « flip conditionnel » est clos (mesuré) ; le poids de direction par région reste une variante plus simple, non mesurée. Utile même avec Bird : sur cette machine Bird laisse la moitié des points de death hors d'atteinte (#59, #75).
10. ⬛ `support_feedback_iterations=3` reste en option (dégrade turkish de 6 %).
11. ⬛ UI, à la toute fin : config UI (pont JSON), sélecteur graphique d'angle, conseil « ne pas centrer », contrainte de l'auto-organiseur, refus d'un modèle trop grand.
12. ⏳ **Lier Bird à l'algorithme C++ et inclure la citation** (demande de l'utilisateur, perdue des documents et rétablie le 2026-10-07). Porter la voie Bird (recherche, `research/S4_Slicer`) dans le cœur C++ d'Apex, accélérée (P2 11 min, P3 10 min). S4 Slicer (Joshua Bird, https://github.com/jyjblrd/S4_Slicer) est sous **GPL-3.0** et fournit un `CITATION.cff` : garder sa notice de licence et son attribution dans le code dérivé, et le citer (README / à propos de l'interface). **Phase 1 faite (#75)** : toute la chaîne Bird portée dans `pipeline/bird/` et branchée dans `main.py` (`--method bird`, `--output-gcode`) ; tranchage, P2, P3, P4, P5 identiques au bit à la recherche, déformation à 0,3 µm près. Bout en bout sur le Benchy : G-code de 1 036 228 mouvements, 1,78 h, 30,4 tours de Y (salon `_v5` : 1 070 769, ~32 tours). **Phase 2 en cours** : moindres carrés de la déformation en C++ (boucle 12 min → 73-111 s, mêmes décisions), Viterbi du retour dans la pièce (P3 123 → 61 s, au bit), contrôle de collision du lissage P5 (161 → 39 s, mêmes décisions). Voie Bird complète sur le Benchy : **633 → 382 s** (export 132 → 7 s en natif, au bit). **Phase 3 mesurée (#75), à placement égal** : Bird imprime 27-53 % plus vite partout, sans appui un peu meilleur (death 0,54 → 0,42 %), 0 collision sur death mais 77/102 sur les deux pièces fines ; turkish : boucle sans aucune passe acceptée, Bird non applicable (et FEM centré : 30 012 points en collision). Ingestion Bird robuste aux auto-intersections (death, turkish). Plan : `~/.claude/plans/bird-into-apex-cpp.md`. ♦️ **À trancher par l'utilisateur** : Bird voie principale, ou cohabitation FEM (pièces fines : 0 collision) / Bird (pièces suspendues : death 0 collision).
13. ⬛ Repoussé après le salon : motifs sur surfaces, ironing, validation sur tous les modèles.
14. ⬛ Sélection manuelle des parties : réservée à la toute fin.

## Feuille de route (récapitulatif complet)

### A. Maker Faire Rome — échéance fin octobre 2026 (~3 semaines)
Choix de l'utilisateur : **pièces imprimées à l'avance et exposées** + **une impression en direct** sur le stand. Pas de tranchage en direct : la vitesse de calcul n'est PAS un prérequis du salon. Deux prérequis seulement : G-code réel pour la machine, et trajectoire vérifiée en collision avec la vraie tête.
1. ⏳ Export G-code réel (`file_io/gcode_export.py`) : 🟢 écrit et testé (#61) ; reste la validation physique : essai à sec de B → `out/test_l_shape.gcode` → Benchy (♦️ ci-dessus).
2. 🟢 Trajectoire vérifiée en collision avec la vraie tête (étape 9 rapide, `aetheris-4.glb`) : l-triple et l_shape sans collision. La limite #58 (tête retournée jusqu'au 2026-09-27) rendait suspectes les collisions antérieures : mesures refaites depuis.
3. ♦️ Z de `config.g` (`M203 Z8000`, `M201 Z1500`), puis réestimer les temps d'impression des pièces du salon.
4. ⬛ Y cumulé validé côté câblage ; L remesuré si la buse change.
5. ⬛ Impression des pièces d'exposition (Benchy, dino Bird ; G-code `_v5` à régénérer si un réglage change).
6. ⬛ Pièce du direct : l-triple proposée (`out/salon/`), à imprimer d'abord une fois à l'essai.
- **Repoussé après le salon si le temps manque** : motifs sur surfaces (P4 complet), validation sur tous les modèles, portage C++ de la voie Bird, ironing.

### B. Temps de calcul — « 10 minutes pour une étape est commercialement inacceptable »
Constat de départ (voie Bird, Benchy, 24 345 tétraèdres, 313 couches) : P2 (zones de lissage) 11 min 30 ; P3 (corps creusé) ~10-11 min ; tranchage de référence 37 s. Règle de travail : mesurer chaque phase, signaler dès qu'elle dépasse la seconde, lancer en fond avec `python -u`.
- **P2 — trois leviers, dans l'ordre** : (1) dernier moment accessible par face, calculé une fois (résultat attendu identique face par face) ; (2) traitement groupé numpy ; (3) parallélisation multi-cœurs. Estimations non mesurées : Python actuel 690 s → Python + leviers 1-2 ~30-60 s → C++ mono-cœur ~10-20 s → C++ 8 cœurs ~2-5 s → GPU ~1-3 s au total. 🟢 État mesuré depuis : P2 62 s (5 tours), P3 61-73 s au bit, P4 54 s (cache par zone), P5 exhaustif 50 s ; chaîne de lissage 274 s ; voie Bird complète du Benchy **633 → 382 s**.
- **P3 — refonte** : découper UNE fois en 3D puis trancher normalement (distance à la peau aux nœuds ; affinage du maillage près de la peau, précision 0,2 mm ; découpe des tétraèdres le long de ψ = 0 par table de cas ; tranchage habituel). La version lente et exacte sert d'ORACLE. Estimations : actuelle ~10-11 min → découpe 3D Python ~5-15 s → C++ multi-cœurs ~1-3 s → + GPU ~1-2 s. **Grande inconnue** : l'affinage — s'il multiplie les tétraèdres par 5-10, le tranchage ralentit d'autant (3 s → 15-30 s) ; première mesure à faire : le nombre de tétraèdres après affinage (s'il passe par TetGen, mono-cœur, il fixe le plancher). ⬛ Non entamée sous cette forme ; P3 a été accélérée par Viterbi natif et creusement parallèle.
- **GPU (Vulkan, décision de l'utilisateur)** : se prêtent au GPU : distance à la peau par nœud, découpe et coupes de couches par tétraèdre, pré-calcul de collision de P2 par face (~30 M de tests point-dans-tétraèdre). Ne s'y prêtent pas : affinage du maillage, chaînage des segments en contours, décalages de polygones de l'étape 6 (Clipper/shapely) — CPU multi-cœurs, parallèle par couche. Réserves : GPU grand public rapides en simple précision seulement (casse l'équivalence au bit sur laquelle repose le portage : porter d'abord, paralléliser ensuite, #56) ; une fois P2 et P3 à quelques secondes le goulot sera ailleurs (boucle Bird ~100 s par passe, étape 9). Contrainte de forme à respecter dès maintenant : calculs indépendants par nœud / tétraèdre / face / couche, sans état partagé. Sur `MeshDistance` : aucun gain mesuré aujourd'hui (⬛ point 8).
- **Plafond réaliste du pipeline FEM** : death 37 s, TetGen ~12 s incompressible ; descendre sous ~18 s exige parallélisation ou changement algorithmique, pas du portage.

### C. Moteur (après le salon, par ordre de priorité indicatif)
1. ♦️ Trancher la voie par défaut Bird / FEM après correction du placement (point 12 ci-dessus).
2. ⏳ **Bird lié au cœur C++ avec citation de S4 Slicer (GPL-3.0)** : phase 1 faite, phase 2 en cours (déformation, Viterbi, P5, export natifs), phase 3 mesurée ; reste le portage de P2/P4/étape 8 Bird en C++, la notice de licence dans le code dérivé, la citation dans le README et l'interface (`CITATION.cff`).
3. ⬛ **Champ scalaire multibranche** (point 9) : seule piste restante contre le défaut de fond du champ unique ; aucune mesure encore. Variante plus simple non mesurée : poids de direction par région.
3 bis. ⏳ **Voie « parties à couches planes »** (décision du 2026-10-19, essai en cours) : chaque partie de step3 tranchée en plans dans sa propre orientation, ordre d'impression entre parties, sans champ scalaire unique. À mesurer sur les 8 modèles (step9b, collisions) avant toute activation.
3 bis-bilan. RESULTAT voie 2 (2026-10-19, `planar`, 8 modeles, step9b) : meilleur que FEM seulement sur l_shape (1,14 vs 2,19 %) et l_triple (1,84 vs 2,54 %) ; 10 a 50 fois pire sur dino/benchy/death/deer/fig (0,18 a 8,7 % contre 0,02 a 0,6 %), collisions 28-71 %. Bogues corriges en route : champs de parties qui partaient tous de 0 (couches fusionnees), step5b qui brouille l'ordre (desactive pour cette voie), detection des parties pendantes (balayage + union-find, persistance 1 mm). Cause residuelle : segmentation step3 trop grossiere (1 seule partie sur real_overhang). SUITE PREVUE (nuit du 2026-10-19) : segmentation par graphe de Reeb de la distance geodesique au plateau (branches = aretes, ordre = arbre depuis le plateau, direction = tangente de branche), puis meme mesure sur 8 modeles. Code non commite : pipeline/step4c_planar_parts.py, options planar_parts, scripts/measure_matrix.py (variantes planar, planar_a<DEG>), export_layer_order_gif.py (--planar-parts, --merge-angle). Resultats : D:/apex_multibranch/results.jsonl, GIF dans D:/apex_multibranch/.
3 bis-reeb. RESULTAT squelette (2026-10-21, `planar_reeb`, `planar_reeb_b6` : graphe de Reeb de la distance geodesique, tranches 4 / 6 mm, coupe aux coudes > 35 deg) : ordre visuellement correct sur real_overhang (pilier, bras, crochet), mais PIRE que la reference FEM sur les 8 modeles, y compris l_shape (3,2 / 4,6 % contre 2,2 %). Meilleur cas des gros modeles : dino 0,30 % (b6) contre 0,019 % ; deer 1,19 % contre 0,16 % ; fig 2,1 % contre 0,30 % ; collisions 4 a 75 % contre 0 a 16 %. Fragmentation : jusqu'a 36 162 couches (fig, 4 mm). Verdict : voie 'parties planes' (regions OU squelette) non viable en l'etat ; opt-in, non active.
3 quater. ⏳ **Champ = distance au plateau (méthode de la chaleur), couches = ses isosurfaces** (2026-10-21/22, `use_heat_geodesic_field`, existant mais jamais mesuré ; limite #79) : ordre visuellement correct sur real_overhang et death (pilier, cabine, coque en dômes depuis son attache, rebords en dernier). Nouvelle mesure : matière FLOTTANTE (morceau de couche qui ne touche ni le plateau ni la matière déjà imprimée, `D:/apex_multibranch/islands.py`). Deux causes de défauts trouvées : (1) le tri historique `(rang de région, isovaleur)` de step5, corrigé par `LayerParams.isovalue_order` (opt-in) : death 0,79 → **0,023 %** flottant, îlots > 1 mm³ 10 → 0 (FEM : 7,1 %) ; (2) step5b (couche entière, une seule normale par couche) : 6,5 % flottant et bandes dans les montants, à ne pas utiliser avec ce champ. Collisions de la tête death 28,9 % (47,9 % avec step5b, 0,5 % pour le FEM). Restent : collisions, orientation (25,7 % des points à plus de 45°), dépôt vers le bas (B > 90°) que step9b compte comme « défaut d'ordre » (85 % de ses 7 319 défauts sur death). `heat_t1` (temps h²) : FAUX succès, trajectoire tronquée. **DÉCISION DE L'UTILISATEUR (2026-10-22) : cette voie est la voie principale** (ordre jugé à l'œil nettement meilleur, temps de calcul 20 à 60 s) ; Bird reste disponible, non abandonné. Bloquant avant impression : collisions de tête (death 28,9 %, ENDGAME à 400 mm de haut : 386 428 événements pour 1,94 M de points extrudés, et 59 % d'orientations inatteignables). Suite prévue, dans l'ordre : (1) projeter le flux de la chaleur sur le demi-plan atteignable (comme #64), (2) réécrire step5b par morceaux de couche ou lancer step10 sur cette voie, (3) mesure sur les 8 modèles puis mode activable.
3 ter. ⬛ **Voie « balayage le long d'un squelette »** (décision du 2026-10-19, ROADMAP, non commencée) : sections perpendiculaires au squelette pour les parties fines (pattes, bois), à combiner avec les parties planes pour les volumes épais. Points durs : coudes, branches, racines multiples (deer : 4 pattes).
4. ⬛ Fin du portage C++ : step10, step5, step4, step1, orchestrateur natif, bascule `APEX_NATIVE_FILL`, retrait de l'échafaudage pybind11 (cible 100 % C++, Python supprimé à la bascule finale).
5. ⬛ Poser automatiquement le modèle sur le plateau à l'ingestion ; refuser dès l'entrée un modèle trop grand (`MachineProfile.build_volume` à renseigner).
6. ⬛ Machines à prendre en charge (`rotation_axes` / `position_coupled_axes`) : actuelle `[C, B]`, `[0]` ; + axe D `[C, B, D]`, `[0]` ; + plateau coulissant `[C, B]`, `[]` ; les deux `[C, B, D]`, `[]`. Pour chaque type : test d'étape 8 (🟢 fait pour D et axe linéaire), directions candidates de step3/step7 (dépendent de la position quand C est couplé), réorientation de step10, traduction G-code (nom et sens des axes). L'upgrade de la machine (axe linéaire) est prévue après Maker Faire.
7. ⬛ `resegmentation_hook` réel de l'étape 10 (re-jouer step3→step8 avec exclusion d'angles, limite #23 ; c'est `main.py` qui doit le fournir).
8. ⬛ Débit de step9 dégressif avec la matière accumulée (mode exact, limites #27/#30) : à chiffrer ; le mode rapide l'a contourné.
9. ⬛ Porter dans l'export la planification d'avance (`feed_plan`) pour la voie Bird si besoin ; réorientation de la recherche contre l'alternance de B.

### D. Interface et produit (à la toute fin)
- ⬛ `config/params.py` + `config/schema.py` : 🟢 contrat JSON et schéma générés (pas 11) ; reste le pont avec l'UI TypeScript.
- ⬛ Sélecteur graphique des angles : arc de cercle avec la valeur en degrés ; convention affichée = angle ENTRE LES FACES (180° = plat), conversion unique et explicite.
- ⬛ « Ne pas centrer les modèles » : conseil dans l'interface (orienter les surplombs vers l'extérieur du plateau, 300 mm de diamètre) + contrainte dans l'auto-organiseur ; effet à MESURER avant d'en faire une règle (le décalage mesuré n'a pas de règle nette, #72).
- ⬛ Refus d'un modèle trop grand dès l'ingestion (plateau traité en plan infini dans le pré-calcul de collision du lissage, tant que le rayon n'est pas dans `MachineProfile`).
- ⬛ Aperçu : rétrécissements à B > 90°.

### E. Tout dernier chantier — sélection manuelle des parties
Réservé : à n'ouvrir QUE lorsque l'algorithme automatique sera au point (jamais comme échappatoire à un défaut). Contenu : sélection de parties par l'utilisateur, paramètres par partie, ordre d'impression choisi. À préserver dès maintenant sans coder : pondérations spatiales PAR TÉTRAÈDRE (jamais un scalaire global) ; partition passée en argument (jamais lue dans `RegionDependencyGraph`) ; ordre d'impression = permutation explicite ; extension via le registre `core/constraints.py`.

## Reste à faire, classé du plus court au plus long (2026-10-13)
Durées = **estimations de l'IA, non mesurées** (une session = quelques heures de travail).

### 1. Minutes — actions de l'utilisateur
- ♦️ `M203 Z8000` et `M201 Z1500` dans `config.g`.
- ♦️ Trancher la voie par défaut : Bird seul, ou cohabitation FEM (pièces fines) / Bird (pièces suspendues).
- ♦️ Décider la purge du docx confidentiel dans l'historique git (`fdf3dbdc`).
- ⬛ Fournir le nom et le sens des axes firmware pour D et l'axe linéaire.
- ⬛ Valider le Y cumulé côté câblage (Benchy ≈ 32, l_shape ≈ 240, dino ≈ 621 tours).
- ⬛ Noter L à la main (`record_pivot_length`) ; remesurer si la buse change.

### 2. Moins d'une heure — 🟢 traité le 2026-10-13
- 🟢 Commit du travail des 06-07/10 (branche `session-2026-10-13`, rien poussé).
- 🟢 Pose du modèle à z = 0 à l'ingestion (#76).
- 🟢 G-code régénérés (`out/share/*_v6.zip`) et temps réestimés ; changer `fw_axis_speeds` si `config.g` change.
- 🟢 Suite complète 495/495 dans les deux positions de `APEX_USE_CPP`.

### 3. Quelques heures
- ♦️ Placement FEM : recentrer à l'ingestion, puis remesurer death et turkish (débloque les chiffres de référence).
- ⬛ Notice GPL-3.0 dans le code dérivé de Bird + citation S4 Slicer (README, `CITATION.cff`, interface).
- ⬛ Refus d'un modèle trop grand (`MachineProfile.build_volume`).
- ⬛ Aperçu : rétrécissements à B > 90° (à diagnostiquer d'abord).
- ⬛ `feed_plan` dans l'export de la voie Bird.
- ⬛ Chiffrer le débit dégressif de step9 en mode exact (#27/#30).
- ⬛ Bascule par défaut de `APEX_NATIVE_FILL`.
- ⬛ Poids de direction par région (variante simple du multibranche, à mesurer sur les 4 modèles).

### 4. Un à quelques jours
- ⬛ Étape 10 en C++.
- ⬛ `resegmentation_hook` réel (re-jouer step3→step8, orchestré par `main.py`).
- ⬛ Machines : directions couplées de step3/step7, réorientation de step10, traduction G-code (axe D, axe linéaire).
- ⬛ Pont JSON vers l'UI TypeScript (le schéma est déjà généré).
- ⬛ Sélecteur graphique d'angle (arc, angle entre les faces).
- ⬛ « Ne pas centrer les modèles » : mesure d'abord (#72 : pas de règle nette), puis conseil UI.
- ⬛ Phase 2 Bird : P2, P4 et étape 8 en C++.
- ⬛ Refonte P3 (découpe 3D unique ; mesurer d'abord le nombre de tétraèdres après affinage).
- ⬛ Réorientation de la recherche contre l'alternance de B.

### 5. Une à plusieurs semaines
- ⬛ Portage C++ de step5, step4 puis step1 (TetGen reste un sous-processus).
- ⬛ Orchestrateur natif, puis retrait de l'échafaudage pybind11 et du Python (cible 100 % C++).
- ⬛ Backend GPU Vulkan : pas de gain mesuré aujourd'hui ; seulement si une étape dépasse la centaine de millions de requêtes.
- ⬛ Auto-organiseur de plateau avec la contrainte « surplombs vers l'extérieur ».

### 6. Incertain, long
- ⬛ **Champ scalaire multibranche** : chantier de recherche, aucune mesure ; l'issue n'est pas garantie (six mécanismes ont échoué avant l'audit du circuit).
- ⬛ Motifs sur surfaces (P4 complet), ironing, validation sur tous les modèles : repoussés après le salon.
- ⬛ Sélection manuelle des parties : tout dernier chantier, après que l'automatique est au point.

## Fait (rappel condensé)
- 🟢 Environnements recréés (venv d'Apex, venv de recherche, noyau C++), `.venv` retirés de l'index git.
- 🟢 SDK Vulkan 1.4.363.0 installé (RTX 3090 vue), aucun code GPU.
- 🟢 Mesures closes sans correction : contours « invalides » Bird (0,3-1,5 % de doublons), couverture dino/Squirtle (0,11 / 0,61 % à > 0,6 mm), décalage sur le plateau (pas de règle nette).
- 🟢 `research/S4_Slicer/safe_io.py` restauré.

## Rappels
- **Le travail du 2026-10-06 et du 2026-10-07 n'est pas commité.** Dernier commit : `fdf3dbdc` (2026-10-04). `main` local suit `origin/session-2026-10-01` ; l'ancien `main` est dans `backup-main-2026-10-06`.
- **Confidentialité** : `Methode_calibration_longueur_pivot_pointe.docx` est retiré de l'index et ignoré par git, mais il figure dans le commit déjà poussé `fdf3dbdc` (dépôt privé `aeloxhardware/apex`, branche `session-2026-10-01`). Le purger de l'historique demande une réécriture et un push forcé : à décider par l'utilisateur.
- Le dépôt suit des `__pycache__` et des fichiers de `cpp/build/`, qui apparaissent comme modifiés.
- Sauvegardes sur H: (`apex_untracked_backup_2026-10-06`).
- Tout G-code antérieur à la régénération `_v5` du 2026-10-06 est obsolète (B de signe opposé).
