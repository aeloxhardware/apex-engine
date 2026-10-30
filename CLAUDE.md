# Apex Slicer — Contexte projet — FICHIER MODIFIABLE PAR L'IA

Moteur de slicing géométrique (Python) pour une imprimante 3D **4 axes** (3 linéaires
+ 1 rotatif). L'interface utilisateur est en **TypeScript**, séparée de ce moteur.
Ce fichier centralise tout le contexte nécessaire pour reprendre le développement
sans redemander les spécifications déjà actées. L'IA peut et doit modifier ce fichier lorsqu'elle opère.

> Note : les trajets à vide sans routage sont désormais corrigés par un
> palliatif lift-and-traverse — voir limite #49 pour le mécanisme, ses
> limites honnêtes (taux de résidu mesuré sur géométrie difficile), et
> pourquoi ce n'est pas un vrai path-planning 3D.

---

## VOIE UNIQUE : LA MÉTHODE DE LA CHALEUR (décision de l'utilisateur, 2026-10-27)

Le champ scalaire EST la distance au plateau à travers la matière, par la méthode de la chaleur
(`pipeline/step4_field_solver.py`) ; les couches sont ses isosurfaces, ordonnées par isovaleur
(limite #79). **Le champ FEM directionnel et la voie Bird sont supprimés**, avec tout ce qui ne
servait qu'eux (champ multibranche, parties planes, contraintes utilisateur `Constraint[]`,
options de step3 et step5 liées). Leur code et leurs mesures restent dans l'historique git et
dans la branche `bird_old` ; ils ne sont plus cités que comme inspiration (S4 Slicer, Joshua
Bird, voir README et `CITATION.cff`). Contrat d'interface passé à la version 2
(`CONTRACT_CHANGES.md`). Les mesures chiffrées des limites ci-dessous qui parlent de « FEM » ou
de mécanismes retirés sont des archives, pas des comportements du code actuel. À refaire sur la
voie chaleur : points d'attache (contraintes utilisateur).

## NUMÉRO DE VERSION ET DÉPÔT COMMUN (décision de l'utilisateur, 2026-10-26)

L'application Apex (moteur + interface) a **un seul numéro de version**, rangé dans le dépôt
commun **`aeloxhardware/apex`** : fichier `VERSION` (source de vérité) et `VERSIONS.md`
(journal). **À chaque changement de version du moteur**, la session met aussi à jour ces deux
fichiers dans `apex` (récupérer d'abord son dernier état : l'interface a pu changer le numéro
entre-temps), et commence ses messages de commit par le numéro complet. **Copie aussi le numéro dans `config/APEX_VERSION`** (lu par `file_io/apex_signature.py`, écrit dans `meta.apex_version` du bloc chiffré des G-code). Format du numéro,
niveaux à monter et règle des conflits : `apex/integration/RULES.md`, §9. Le détail des
changements du moteur reste dans le DT et les limites numérotées de ce fichier.

Le moteur et l'interface se développent séparément, une session par dépôt :
`apex/integration/RULES.md` (règles et contrat entre les deux) et
`apex/integration/REQUESTS.md` (demandes ouvertes, à lire en début de session).

## RAPPELS PERMANENTS À L'UTILISATEUR (demande du 2026-10-18)

**À la fin de CHAQUE réponse, rappeler en une ligne le point ci-dessous, tant que
`DT.md` (section « RAPPELS PERMANENTS ») ne le marque pas fait `[x]`.** Ne cesser un rappel
que lorsque l'utilisateur dit l'avoir fait et que la case est cochée dans `DT.md`.
1. ~~Modifier les paramètres du firmware~~ : **fait** (confirmé par l'utilisateur, 2026-10-26, `DT.md` coché). Rappel arrêté. Salon MFR 2026 passé avec une récompense spéciale pour l'originalité du projet ; **le Benchy n'a toujours pas été imprimé**.
2. ~~Essayer le champ scalaire multibranche~~ : **fait le 2026-10-19** (limite #78, résultat mitigé,
   opt-in). Rappel arrêté.

## COMMENT LIRE CE FICHIER (à lire AVANT le reste)

Ce document est long, dense et écrit avec assurance. **Cette assurance n'est
pas une preuve, et sa structure n'est pas une contrainte.** Il décrit ce qui a
été fait, pas les limites de ce qui peut l'être. Lu naïvement, il induit un
biais mesuré en session réelle : optimiser *à l'intérieur* du découpage
existant au lieu de se demander si ce découpage pose la bonne question.

### Autorisation explicite de l'utilisateur

**Tout est modifiable pour peu que ça marche et que ce soit mesuré** : le
découpage en 10 étapes, les signatures `run(input, params)`, les dataclasses
de `core/mesh_types.py`, l'ordre des étapes, l'existence même d'une étape.
Ajouter une étape, en fusionner deux, faire remonter une donnée d'une étape
tardive vers une étape précoce, remplacer un critère par un autre : rien de
tout cela n'exige de demander la permission. Une fonctionnalité qui impose de
restructurer la moitié du pipeline se restructure. Les seules obligations qui
demeurent : **mesurer avant de conclure**, et **ne jamais activer par défaut
un mécanisme non mesuré sur les quatre modèles** de `res/test_models/`.

### La question à poser en premier quand un indicateur ne bouge pas

Sept mécanismes ont été écrits et mesurés pour réduire la matière déposée sans
appui. Six ont échoué. Tous les six modifiaient une **entrée** du circuit
existant (choix de direction, ancrage, ordre des couches). Aucun n'avait
interrogé le circuit lui-même. Le septième — celui qui a enfin fait bouger la
métrique — est né d'une question de l'utilisateur, pas du raisonnement mené
dans ce cadre : *« l'algorithme ne doit-il pas simplement vérifier, pour
chaque tétraèdre, s'il est en l'air ? »*

La réponse était que ce test existait déjà **deux fois**, et qu'il était
branché au mauvais endroit : `step2` évalue le surplomb contre la verticale
MACHINE, une seule fois, **avant que le champ scalaire n'existe** ; `step9b`
le mesure correctement mais *après tout*, sans rien réinjecter. Le solveur
n'était jamais pénalisé pour les surplombs qu'il créait lui-même.

**Règle de méthode qui en découle, à appliquer avant d'écrire un mécanisme de
plus :**

1. **Après deux ou trois mécanismes mesurés sans effet, arrêter d'en ajouter
   et auditer le circuit.** Un quatrième mécanisme sur une boucle ouverte ne
   la fermera pas.
2. **Pour chaque critère, demander : contre quelle RÉFÉRENCE est-il évalué, et
   est-ce celle qui compte à cet endroit du pipeline ?** Le surplomb mesuré
   contre la verticale machine et le surplomb mesuré contre `∇φ` portent le
   même nom et ne sont pas la même quantité.
3. **Suspecter tout ce qui est calculé UNE FOIS, en amont de ce qu'il
   contraint.** Une valeur figée avant l'existence de l'objet qu'elle régit
   est une boucle ouverte jusqu'à preuve du contraire.
4. **Un indicateur intermédiaire ne prouve rien.** Voir le 4ᵉ enseignement de
   la section « Prochaine étape » : direction, ancrage et ordre ont tous été
   améliorés sans que la métrique de référence (step9b) suive.

### Ce que ce fichier ne dit pas

Une limite documentée n'est une preuve de rien tant que le code n'a pas été
relu — la limite #50 a décrit pendant des sessions entières un mécanisme qui
n'existait pas. Les « pistes écartées » l'ont été dans un contexte donné ; une
piste écartée pour une raison qui ne tient plus se rouvre. Les sections
« ce qui reste ouvert » formulent des questions **déjà cadrées** : ce cadrage
est lui-même à remettre en cause en priorité, puisque c'est exactement là que
le raisonnement s'est enlisé six fois.

## Objectif du projet

Calculer des trajectoires d'impression optimales évitant **tout support**, en
inclinant/courbant les couches selon la géométrie locale, plutôt que d'imprimer
des couches planes horizontales comme un slicer 3-axes classique (Orca, PrusaSlicer).

Fonctionnalités visées, au-delà du slicing de base :
- Couches **courbes** (pas seulement planes par sous-volume), via un champ scalaire
  continu résolu sur un maillage tétraédrique.
- Remplissage de l'intérieur vers l'extérieur, avec gestion dédiée des **ponts internes**
  (plafonds de cavités non détectés par une simple analyse de normale).
- **Points d'attache** manuels : l'utilisateur peut déformer les couches après le
  slicing automatique (esthétique ou solidité), en ajoutant des contraintes ponctuelles
  au champ scalaire plutôt qu'en éditant directement la géométrie des couches.
- **Peinture de couleur** sur le modèle + impression groupée par couleur (repassage
  buse 5 axes) pour limiter les purges : trous temporaires comblés au repassage,
  avec bascule automatique vers une purge classique si le trou n'est pas réalisable
  ou pas rentable (calcul de coût trou vs purge).
- **Simulation de collision** avec un modèle 3D réel de la tête d'impression, avec
  boucle d'ajustement automatique en cas de collision détectée.

Architecture voulue par l'utilisateur : **plusieurs petits fichiers Python**, un par
étape du pipeline + un fichier de paramètres + un orchestrateur — pas un monolithe.
Chaque étape expose une fonction pure `run(input, params) -> output`, sans état caché.

---

## Fonctionnalités produit (ce que l'outil doit permettre à l'utilisateur)

Liste de toutes les fonctionnalités discutées, au-delà du slicing plan classique :

- **Tranchage en couches courbes non-planes**, épousant la géométrie locale plutôt que
  des plans horizontaux — objectif central du projet, rendu possible par le 4ᵉ axe.
- **Remplissage de l'intérieur vers l'extérieur** (contours concentriques, du centre
  vers la paroi externe), pas l'inverse.
- **Détection et traitement dédié des ponts internes** (plafonds de cavités non
  supportés par en dessous), distincts des surplombs classiques : trajectoire
  réorientée pour traverser le pont perpendiculairement à son axe de plus grande
  portée, avec ancrage renforcé sur les bords.
- **Choix automatique de l'inclinaison de couche** par zone du modèle, pour éliminer
  tout besoin de support — le modèle est segmenté en sous-volumes, chacun avec sa
  propre direction d'impression optimale.
- **Édition manuelle de l'orientation/déformation des couches après le slicing
  automatique**, via un système de **points d'attache** que l'utilisateur peut
  déplacer et faire pivoter — pour des raisons esthétiques ou de solidité. Les points
  d'attache déforment le champ scalaire sous-jacent (pas les couches directement),
  ce qui garantit que les couches recalculées restent géométriquement cohérentes.
  Un garde-fou signale (sans forcément bloquer) si une édition manuelle réintroduit
  un surplomb au-delà du seuil critique.
- **Peinture de couleur directement sur le modèle 3D**, zone par zone.
- **Impression groupée par couleur avec repassage** (profite de la liberté du 4/5
  axes pour repasser sur n'importe quel point déjà imprimé) plutôt que de purger à
  chaque changement de couleur à chaque couche.
- **Gestion de trous temporaires** : quand une zone d'une couleur a besoin qu'une
  autre couleur soit posée avant elle mais que ce n'est pas encore le cas, un vide
  temporaire est laissé puis comblé lors d'un repassage ultérieur.
- **Bascule automatique trou vs purge classique** : si combler un trou n'est pas
  géométriquement réalisable (buse ne peut pas atteindre le fond sans collision) ou
  pas rentable (coût de repassage/risque d'adhérence à froid supérieur au coût d'une
  purge), le système purge normalement à la place.
- **Modularité/extensibilité de l'algorithme** : possibilité d'ajouter de nouvelles
  fonctions et contraintes (esthétiques, structurelles, etc.) sans réécrire le cœur
  du pipeline — chaque nouvelle contrainte est un terme de plus injecté dans le
  champ scalaire ou un post-traitement optionnel.
- **Modèle 3D réel de la tête d'impression** (buse + porte-buse + carriage proche)
  utilisé pour la vérification de collision, plutôt qu'une simple pointe ponctuelle.
- **Simulation de collision complète en fin de pipeline**, contre la pièce déjà
  imprimée (état incrémental), le mandrin/mors, et le plateau à son angle courant.
- **Ajustement automatique en cas de collision détectée** : réorientation locale de
  l'angle du plateau en priorité, sinon ré-segmentation ciblée de la zone concernée,
  et en dernier recours ajout d'un support classique avec log explicite pour
  l'utilisateur.

---

## Arborescence du projet

> **Depuis le 2026-10-29, chaque algorithme de tranchage vit dans son dossier `element-xxxx/`** (décision de
> l'utilisateur) : `element-ardor/` (pipeline principal 5 axes ; c'est l'ancien `pipeline/`), `element-aeolis/`
> (non planaire ; ancien `nonplanar/`), `element-aequor/` et `element-arva/` (à venir, README seulement). Un nom à
> tiret n'étant pas importable, **`pipeline/` et `nonplanar/` à la racine ne sont plus que des alias d'import** (un
> `__init__.py` qui redirige `__path__`) : tous les `from pipeline.stepN_... import` et `python -m
> pipeline._tetrahedralize_worker` restent valables. L'arborescence ci-dessous cite encore `pipeline/` : lire
> `element-ardor/`. Le code partagé (`core/`, `file_io/`, `config/`, `cpp/`) reste à la racine.

```
slicer_engine/
├── CLAUDE.md                        # ce fichier
├── main.py                          # FAIT — orchestrateur CLI, enchaîne step1 à step10
├── scripts/
│   ├── generate_test_model.py       # FAIT — génère res/test_models/l_shape_overhang.stl (surplomb réel)
│   ├── export_layer_order_gif.py    # FAIT — GIF animé des couches ajoutées UNE PAR UNE dans
│   │                                 #   leur ordre d'impression réel (même dégradé que le GLB ;
│   │                                 #   couches `support_unresolved` en ROUGE VIF). Rejoue
│   │                                 #   run_slicing_only, ne recalcule rien lui-même. Rend via
│   │                                 #   matplotlib 3D (+ Pillow) — seul usage de ces deux
│   │                                 #   dépendances dans le projet, installées pour ce script.
│   ├── export_slicing_preview.py    # FAIT — GLB du slicing SEUL (step1→step8), court-circuite
│   │                                 #   step9/step10 le temps que la zone de garde buse soit faite
│   └── build_cpp.py                 # FAIT — construit le noyau natif (localise vcvars64.bat via
│                                     #   vswhere, puis cmake+ninja). Une seule commande.
├── cpp/                             # NOYAU NATIF — portage en cours, voir limite #56
│   ├── third_party/clipper2/        # Clipper2 2.0.1 (BSL-1.0), offsets polygonaux (pas 15)
│   ├── third_party/eigen/           # Eigen 3.4.0 (MPL-2.0), header-only, VERSIONNE.
│   │                                 #   Introduit au pas 7 seulement, quand step4 en a eu besoin.
│   ├── CMakeLists.txt               # FAIT — C++20, /O2, /fp:precise (JAMAIS /fp:fast : il
│   │                                 #   casserait l'oracle d'équivalence)
│   ├── include/apex/
│   │   ├── types.hpp                # FAIT — Vec3, AABB, segment_aabb. Réduit au strict
│   │   │                             #   nécessaire : on ne transcrit pas les 40 dataclasses
│   │   │                             #   d'avance, sinon on fige un contrat avant d'avoir porté
│   │   │                             #   le code qui le consomme.
│   │   ├── printed_index.hpp        # FAIT — grille spatiale des capsules imprimées
│   │   ├── kinematics.hpp           # FAIT — IK N axes, SEQUENTIELLE et REPRENABLE
│   │   ├── support.hpp              # FAIT — test de support point/capsules, arrêt anticipé
│   │   ├── tet_topology.hpp         # FAIT — tables topo tétra (faces frontière, arêtes)
│   │   ├── kdtree.hpp               # FAIT — arbre k-d 3D statique (pas 5). Portage
│   │   │                             #   ARCHITECTURAL (cKDTree déjà natif) : retire la
│   │   │                             #   dernière primitive spatiale scipy du cœur.
│   │   ├── bvh.hpp                  # FAIT — BVH de triangles + intersects_location (pas 6).
│   │   │                             #   Remplace rtree comme INDEX, PAS la requête : la
│   │   │                             #   boîte de rayon reste celle de trimesh (ray_bounds),
│   │   │                             #   sans quoi l'auto-intersection réapparaît (limite #56)
│   │   ├── field_solver.hpp         # FAIT — pas 7 : blend/gradients/CSR/solveur (step4)
│   │   ├── mesh_repair.hpp          # FAIT — pas 8 : SAT triangle-triangle (step1_ingest),
│   │   │                             #   BVH interne sur AABB exactes, broad-phase non
│   │   │                             #   contraint (contrairement au pas 6 : narrow-phase
│   │   │                             #   déterministe des deux côtés)
│   │   ├── mesh_intersect.hpp       # FAIT — pas 9 : segments d'intersection + winding
│   │   │                             #   number (_self_intersection_resolve). Tolérance
│   │   │                             #   explicite (solve 3×3 par cofacteurs) ; ordre des
│   │   │                             #   segments par face NON reproductible (rtree),
│   │   │                             #   divergence topologique mineure assumée (×11,5)
│   │   ├── interpolate.hpp          # FAIT — pas 10 : pondération inverse-distance de
│   │   │                             #   l'offset hors-plan (step6), ×2,8 — DERNIER poste
│   │   │                             #   Python isolé visible au profil hors TetGen
│   │   ├── toolpath2d.hpp           # FAIT — pas 15 : parois + remplissage d'une zone (step6),
│   │   │                             #   Clipper2, tolérance déclarée, opt-in APEX_NATIVE_FILL=1
│   │   ├── (isosurface, geometry2d, past_material, gcode_export, params : pas 11-24, limite #56)
│   │   └── internal/
│   │       └── triangle_aabb_bvh.hpp  # FAIT — BVH interne partagé (pas 8+9), sans
│   │                                    #   contrat d'équivalence (structure privée)
│   ├── src/core/printed_index.cpp   # FAIT — portage de _PrintedGeometryIndex (2 modes de
│   │                                 #   requête, pour mesurer séparément le gain du langage
│   │                                 #   et celui de l'algorithme)
│   ├── src/core/kinematics.cpp      # FAIT — chaîne analytique de geometry_utils + boucle
│   │                                 #   step8. Repli grille NON porté : rend la main pour
│   │                                 #   le point fautif puis REPREND (voir limite #56)
│   ├── src/core/support.cpp         # FAIT — pas 3
│   ├── src/core/tet_topology.cpp    # FAIT — pas 4
│   ├── src/core/kdtree.cpp          # FAIT — pas 5 : median-split, query_nn (1-NN batch,
│   │                                 #   distance euclidienne bit-exacte / scipy) + query_radius.
│   │                                 #   Tie-break plus petit indice. 9 sites .query() câblés.
│   ├── src/core/bvh.cpp             # FAIT — pas 6 : ×38 à ×221 selon la densité du maillage,
│   │                                 #   impacts bit-identiques. narrow phase TRANSCRITE de
│   │                                 #   trimesh (jamais un Möller-Trumbore : autres flottants)
│   ├── src/core/field_solver.cpp    # FAIT — pas 7 : blend des directions (×3,2, gain
│   │                                 #   ALGORITHMIQUE par élagage AABB), gradients de forme (×30),
│   │                                 #   assemblage FEM direct en CSR (×2,04, motif partagé)
│   ├── src/core/linear_solver.cpp   # FAIT mais NON CÂBLÉ — CG Eigen mesuré à ×1,01 : aucun
│   │                                 #   gain, et coûterait la bit-exactitude (limite #56, pas 7)
│   ├── src/core/mesh_repair.cpp     # FAIT — pas 8 : triangles_cross (SAT 11 axes) +
│   │                                 #   find_self_intersecting_faces (broad+narrow phase),
│   │                                 #   n'exerce jamais le chemin nominal death.stl —
│   │                                 #   gain mesurable sur turkish.stl (cascade de réparation)
│   ├── src/core/mesh_intersect.cpp  # FAIT — pas 9 : triangle_triangle_segment (cofacteurs,
│   │                                 #   tolérance ~1e-14), generalized_winding_number
│   │                                 #   (quasi bit-exact), find_intersection_segments
│   │                                 #   (×11,5 sur turkish.stl, ordre non reproductible)
│   ├── src/core/interpolate.cpp     # FAIT — pas 10 : une seule passe par requête, sans
│   │                                 #   matérialiser la matrice de distances (K,M)
│   ├── bindings/module.cpp          # FAIT — pybind11. ÉCHAFAUDAGE, retiré à la bascule finale.
│   ├── build/                       # (ignoré par git) artefacts cmake/ninja
│   └── lib/                         # (ignoré par git) apex_core.cp314-win_amd64.pyd
├── config/                          # PAS ENCORE ÉCRIT — sera livré en C++ (params.hpp + JSON),
│                                     #   voir limite #56 : c'est ce qui deviendra le contrat UI TS
│   ├── params.py                    #   dataclasses de tous les paramètres globaux
│   └── schema.py                    #   (de)sérialisation JSON <-> dataclasses, pont avec l'UI TS
├── pipeline/
│   ├── __init__.py
│   ├── step1_ingest.py              # FAIT — ingestion mesh + tétraédrisation (+ détection
│   │                                 #   automatique des cavités internes pour tetgen, voir limite #7 ;
│   │                                 #   + raffinement local près du plateau via add_region, limite #36)
│   ├── step2_overhang.py            # FAIT — détection surplombs + ponts internes
│   ├── step3_segmentation.py        # FAIT — tests pytest formalisés (9 tests)
│   ├── step4_field_solver.py        # FAIT — résolution FEM du champ scalaire (14 tests)
│   ├── step5_layers.py              # FAIT — extraction des couches courbes (14 tests)
│   ├── step5b_collision_order.py    # FAIT — ordre des couches contraint par la tête (limite #67)
│   ├── step6_toolpath2d.py          # FAIT — remplissage 2D par couche (11 tests)
│   ├── step7_scheduling.py          # FAIT — couleur + ordonnancement trous/purges (19 tests)
│   ├── step8_machine_toolpath.py    # FAIT — resolution IK generique N-axes (16 tests)
│   ├── step9_collision_sim.py       # FAIT — collision incrementale (capsules), N-axes (18 tests)
│   └── step10_adjustment.py         # FAIT — boucle d'ajustement 3 niveaux (15 tests)
├── core/
│   ├── __init__.py
│   ├── mesh_types.py                # FAIT — toutes les dataclasses partagées (voir détail plus bas)
│   ├── geometry_utils.py            # FAIT — build_face_adjacency, fibonacci_sphere_directions,
│   │                                 #   rotation_axis_circle_directions, composed_rotation_directions
│   │                                 #   (N axes composes, cinematique generique),
│   │                                 #   composed_rotation_angles_and_directions,
│   │                                 #   refine_axis_angles_local_search (IK generique, step8),
│   │                                 #   point_in_polygon_2d
│   ├── _cpp.py                      # FAIT — chargeur du noyau natif + interrupteur APEX_USE_CPP
│   │                                 #   (Python par DÉFAUT ; absence du .pyd jamais une erreur).
│   │                                 #   Échafaudage de migration, retiré à la bascule.
│   └── spatial.py                   # FAIT — KDTree (pas 5) + RayIntersector (pas 6) :
│                                     #   scipy/trimesh OU apex_core selon APEX_USE_CPP, choix
│                                     #   figé à la construction. Sous-ensembles réellement
│                                     #   utilisés : .query() k=1, .query_ball_point(), .data ;
│                                     #   .intersects_location(). Construits par step1
│                                     #   (mesh_data.kdtree/bvh), step3/4/5/7.
├── file_io/                          # renommé depuis `io/` — collision avec le module
│   │                                 #   standard Python `io` (voir limite #20)
│   ├── mesh_loader.py                # FAIT — load_toolhead_mesh (GLB -> Trimesh mm, scene-space)
│   ├── pivot_length.py               # FAIT — longueur pivot-pointe L courante (config/pivot_to_tip.json),
│   │                                 #   lue par l'export (limite #74)
│   ├── visualization_export.py       # FAIT — export_visualization_glb (scène GLB : mesh + couches ;
│   │                                 #   trajectoire OMISE par défaut, `show_toolpath=True` la rétablit)
│   └── gcode_export.py               # PAS ENCORE ÉCRIT
├── res/
│   └── test_models/
│       └── l_shape_overhang.stl      # FAIT — modèle de test avec surplomb réel (même géométrie que
│                                      #   tests/test_step3_segmentation.py::l_shape_data)
└── tests/
    ├── __init__.py
    ├── test_step1_ingest.py         # FAIT — 11 tests (raffinement local pres du plateau,
    │                                 #   seed garanti a l'interieur via _find_guaranteed_interior_point)
    ├── test_step2_overhang.py       # FAIT — 6 tests, tous passants et stables
    ├── test_step3_segmentation.py   # FAIT — 13 tests (dont 4 sur max_region_merge_angle_deg, modele en T)
    ├── test_step4_field_solver.py   # FAIT — 14 tests (melange spatial de directions)
    ├── test_step5_layers.py         # FAIT — 27 tests (clustering spatial + tri DFS de continuite,
    │                                 #   modeles tripode/portique)
    ├── test_step6_toolpath2d.py     # FAIT — 45 tests (offset hors-plan interpolé, limite #31 ;
    │                                 #   motifs de remplissage, limite #57)
    ├── test_step7_scheduling.py     # FAIT — 19 tests
    ├── test_geometry_utils.py       # FAIT — 7 tests (composed_rotation_directions,
    │                                 #   _rodrigues_rotate_single)
    ├── test_step8_machine_toolpath.py  # FAIT — 20 tests (dont orientation physique de la tête, limite #58)
    ├── test_step9_collision_sim.py  # FAIT — 18 tests
    ├── test_step10_adjustment.py    # FAIT — 15 tests
    └── test_cpp_equivalence.py      # FAIT — 26 tests, ORACLE du portage C++ : le natif doit
                                      #   rendre exactement ce que rend le Python (pas 1-5,
                                      #   dont k-d tree : indices ET distances bit-exactes vs
                                      #   scipy). Se SAUTE proprement si cpp/lib/ est vide
                                      #   (pas de compilateur requis pour lancer la suite).
```

**État global (2026-10-01, session autonome) : 461/461 dans les deux positions de
`APEX_USE_CPP`**, aucun test sauté (module `cpp/lib` reconstruit avec les pas 11 à 24). Le
2026-09-30 : 406/406. Auparavant : 402 dans les deux positions, les 3
tests `MeshDistance` se sautant tant que le module n'était pas reconstruit. Avant : 386/386 tests passants dans les deux positions de `APEX_USE_CPP` (step1 à step10, plus les étapes intermédiaires step1b/step8b/step9b, core/geometry_utils.py, et les 56 tests d'équivalence du portage C++). Ce compteur a longtemps affiché 166, puis 272, puis 318, puis 331, puis 340, puis 344, puis 351, puis 357, puis 362, puis 365, puis 370, puis 375, puis 385, alors que la suite en comptait déjà bien davantage : le remettre à jour à chaque session, ou ne pas s'y fier.

**Un noyau natif C++ est en cours d'intégration** (`cpp/`, voir limite #56) :
`_PrintedGeometryIndex` et la cinématique inverse N axes (step8) y sont déjà
portés — pipeline complet **×3,07** sur `death.stl`, à métrique identique.
Il est **désactivé par défaut** —
il faut `APEX_USE_CPP=1` pour l'emprunter, et l'absence du module compilé
n'est jamais une erreur. La suite de tests doit donc passer dans les DEUX
positions de l'interrupteur, et c'est le critère de recette de chaque
portage.
**Les 10 étapes du cœur algorithmique sont écrites, stables et testées, et
`main.py` les enchaîne bout en bout sur un fichier STL/3MF réel** (vérifié
sur `res/test_models/l_shape_overhang.stl`, un modèle avec un vrai
surplomb). Restent hors cœur : `config/` et `file_io/gcode_export.py`
(export G-code final).

**Environnement Python du projet** : venv local (`.venv/`, ignoré par git)
sur **Python 3.14** — pas 3.12, qui a cessé de fonctionner sur cette
machine (interpréteur plantant avec `STATUS_DLL_INIT_FAILED` dès `python -c
"print(1+1)"`, indépendant du code du projet). Toutes les dépendances
listées ci-dessous s'installent sans compilation sur 3.14 (wheels
précompilées disponibles pour chacune, y compris `tetgen` dont la wheel est
`cp312-abi3` mais compatible ABI3 avec 3.14). Recréer l'environnement :
`python -m venv .venv` puis `.venv\Scripts\python.exe -m pip install
trimesh scipy tetgen rtree manifold3d networkx shapely pytest`.

**Dépendances Python requises** (installer via pip) :
```
trimesh scipy tetgen rtree manifold3d networkx shapely pytest
```
**Dépendance OPTIONNELLE** : `pymeshfix` (wheel `cp312-abi3`, compatible
3.14 comme `tetgen`). Importée paresseusement par `step1_ingest`, uniquement
quand TetGen refuse un maillage pour auto-intersections — un environnement
sans elle fonctionne exactement comme avant sur tout maillage déjà valide.
Sans elle, `res/test_models/death.stl` ne passe pas (voir limite #43).
**`death.stl` est un Benchy imprimé à l'envers** (capot/coque en porte-à-faux
au-dessus d'un pilier central et d'une base circulaire — voir capture
d'écran fournie par l'utilisateur en session) : 31×62×50mm, watertight,
2418 faces, ~22% de la surface en surplomb franc (>45°). Modèle de
référence pour mesurer le taux de matière déposée sans support — la
géométrie même du capot garantit de vrais surplombs à résoudre par
inclinaison 4 axes, pas de simples défauts d'ordre de dépôt.
`pip install <package> --break-system-packages` si environnement géré (Debian/Ubuntu récent).
`shapely` ajouté à l'étape 6 pour l'offset polygonal (`buffer`) et le
clipping d'infill — géométrie 2D robuste (trous, auto-intersections,
polygones multi-parties après un offset qui sépare une forme en deux),
hors de portée d'une implémentation maison contrairement au reste du
pipeline géométrique (point-in-polygon, marching tetrahedra, etc.).

---

## Paramètres globaux (contrat UI TypeScript → moteur Python)

```
MachineProfile {
  // Cinématique décrite comme une liste ORDONNÉE et PARAMÉTRABLE d'axes de
  // rotation — PAS un modèle figé sur une seule machine. Chaque machine
  // future (un seul DOF, deux DOF différemment agencés, plus de deux DOF)
  // se décrit par sa propre liste rotation_axes, sans changer le contrat.
  // Voir limite de conception #18 (composition de N rotations, la machine
  // actuelle n'est qu'un cas d'usage parmi d'autres) et #19 (singularité au
  // centre du premier axe de rotation, si applicable).
  rotation_axes: RotationAxisSpec[] {
    axis: Vec3                        // unitaire, dans le repère du DOF parent
    angle_range_deg: {min: float, max: float}  // (0,360) tour complet, plage
                                                 // restreinte, ou (v,v) fige
  }
  // Exemple d'usage — machine actuelle (commande firmware, gabarit
  // symbolique) : G00 Y<C> X<R> Z<Z> B<B> F<feedrate>
  //   R (firmware "X") = rayon du plateau, mm ; C (firmware "Y") = rotation
  //   du plateau, degrés ; Z (firmware "Z") = hauteur ; B (firmware "B") =
  //   rotation de la TÊTE, perpendiculaire à R (B=0° tête vers le bas,
  //   B=90° horizontale, B=180° retournée vers le haut, jamais négatif).
  //   **Plage réelle mesurée : B ∈ [0°, 185°]**, légèrement AU-DELÀ de la
  //   verticale retournée (corrigé dans `main._default_machine_profile`, qui
  //   figeait 180°). Conséquence : déposer vers le bas EST physiquement
  //   possible sur cette machine — toute règle du moteur qui suppose
  //   l'inverse est fausse par construction (voir limite #51,
  //   `enforce_upward_growth`).
  //   rotation_axes = [ {axis: (0,0,1), angle_range_deg: (0,360)},   // C
  //                      {axis: perpendiculaire à R, angle_range_deg: (0,185)} ]  // B
  platform_origin: Vec3             // centre de rotation du premier axe (si applicable) —
                                     // aussi le point de singularité éventuel (limite #19)
  nozzle_model: Mesh                                 // maillage 3D tête (buse + porte-buse + carriage proche)
  nozzle_diameter: float
  max_angular_velocity: float                        // rad/s
  max_angular_acceleration: float
  build_volume: AABB
  z_clearance_min: float
}

SlicingParams {
  layer_height_base: float
  overhang_angle_threshold: float                    // ex 45°
  bridge_detection_ray_density: float
  nozzle_width: float
  infill_density: float
  infill_pattern: enum
  wall_count: int
}

FieldSolverParams {
  resolution: float                                   // taille de tétra pour la résolution FEM
  max_solver_iterations: int
  convergence_tolerance: float
}

ColorMap {
  per_face_color: Map<face_id, color_id>
  merge_tolerance: float
}

SchedulingParams {
  hole_cost_model: {reheat_penalty: float, adhesion_risk_weight: float}
  purge_cost_model: {volume_per_purge: float, time_per_purge: float}
  max_open_holes: int
}

CollisionSimParams {
  path_sampling_step: float
  safety_margin: float
  max_adjustment_iterations: int
}
```

---

Voici la spécification complète, étape par étape, avec les données précises manipulées à chaque phase et tous les paramètres remontés depuis l'UI TypeScript. Je vais rester au niveau "structures de données + opérations", pas pseudo-code exécutable, mais suffisamment précis pour être directement implémentable.

---

## Paramètres globaux reçus de l'interface (contrat UI → moteur)

```
MachineProfile {
  // Cinématique décrite par rotation_axes (liste paramétrable de N axes) —
  // voir le premier bloc MachineProfile ci-dessus pour le détail complet
  // (gabarit générique, exemple machine actuelle, limites #18/#19). Rappel synthétique :
  rotation_axes: RotationAxisSpec[]   // N axes ordonnés, chacun {axis, angle_range_deg}
  platform_origin: Vec3             // centre du premier axe (si applicable) — singularité éventuelle
  nozzle_model: Mesh                                 // maillage 3D de la tête (buse + porte-buse + carriage proche)
  nozzle_diameter: float
  max_angular_velocity: float                        // rad/s, limite plateau
  max_angular_acceleration: float
  build_volume: AABB
  z_clearance_min: float                             // garde mini buse/pièce
}

SlicingParams {
  layer_height_base: float
  overhang_angle_threshold: float                    // ex 45°
  bridge_detection_ray_density: float                // résolution du ray-cast pour ponts internes
  nozzle_width: float
  infill_density: float
  infill_pattern: enum
  wall_count: int
}

FieldSolverParams {
  resolution: float                                   // taille de voxel/tétra pour la résolution FEM
  max_solver_iterations: int
  convergence_tolerance: float
}

ColorMap {
  per_face_color: Map<face_id, color_id>              // peinture utilisateur
  merge_tolerance: float                               // lissage/regroupement avant ordonnancement
}

SchedulingParams {
  hole_cost_model: {reheat_penalty: float, adhesion_risk_weight: float}
  purge_cost_model: {volume_per_purge: float, time_per_purge: float}
  max_open_holes: int                                  // borne sur les trous simultanés en attente
}

CollisionSimParams {
  path_sampling_step: float                            // pas d'échantillonnage trajectoire pour test collision
  safety_margin: float
  max_adjustment_iterations: int
}
```

---

Les 10 étapes du pipeline

## Étape 1 — Ingestion et préparation du mesh

**Ce que le code fait précisément :**
- Charge le mesh (STL/3MF) via `trimesh`, vérifie la watertightness (`mesh.is_watertight`), répare si besoin (`trimesh.repair.fill_holes`, `fix_normals`).
- Calcule les normales par face (`mesh.face_normals`) et par sommet (moyenne pondérée par aire).
- Construit une structure d'accélération spatiale : BVH (`trimesh.ray.ray_pyembree` ou `rtree`) pour les ray-cast ultérieurs, et un KDTree (`scipy.spatial.cKDTree`) sur les centres de faces.
- Voxelise ou tétraédrise le volume à la résolution `FieldSolverParams.resolution` (réutilisation directe de ton pipeline CuPy tétraédrique existant) — c'est le domaine sur lequel le champ scalaire sera résolu.

**Sortie :** `Mesh` nettoyé + `TetMesh` (ou grille voxel) + structures d'accélération.


## Étape 2 — Détection des surplombs et des ponts internes

**Ce que le code fait précisément :**
- Pour chaque face : angle entre sa normale et la verticale locale de référence (au départ, l'axe Z machine avant toute rotation) ; marque `is_overhang = angle > overhang_angle_threshold`.
- Pour les ponts internes : pour chaque face quasi-horizontale non marquée surplomb, lance un ray-cast vertical vers le bas depuis un échantillonnage de points sur la face (densité = `bridge_detection_ray_density`) ; si le rayon ne touche aucune autre face du mesh avant une distance seuil, marque la zone `is_bridge = true`.
- Construit une carte scalaire par face : `risk_score = f(is_overhang, is_bridge, angle)`, utilisée à l'étape suivante comme poids de contrainte.

**Sortie :** `OverhangMap` (face_id → risk_score, is_bridge bool).


## Étape 3 — Segmentation en sous-volumes candidats

**Ce que le code fait précisément :**
- Construit un graphe d'adjacence de faces (dual du mesh).
- Applique un algorithme de croissance de région (region growing) ou un graph-cut multi-labels (`scikit-learn` / implémentation maison via `networkx` + coupe min) où :
  - le coût d'assigner une face à une direction candidate `d` = pénalité si l'angle normale/d dépasse le seuil, pondérée par `risk_score`
  - le coût de frontière entre deux régions voisines de directions différentes = pénalise la fragmentation excessive (paramètre de régularisation, à exposer côté UI si besoin)
- Résultat : chaque face reçoit un `region_id` et chaque région une direction d'impression candidate (vecteur).
- Construit le graphe de dépendances entre régions sous forme d'**arborescence de croissance depuis le plateau** : parcours en largeur multi-sources depuis les régions au contact du plateau à travers l'adjacence entre régions, puis orientation de chaque arête de la profondeur faible vers la profondeur forte (deux régions de même profondeur poussent en parallèle et ne portent aucune arête). L'ordre suit « par où la matière rejoint le plateau », jamais une comparaison de hauteur — une pièce pendante est plus basse que ce à quoi elle est suspendue mais doit être imprimée après. Acyclique par construction (la profondeur croît strictement le long de chaque arête), donc le tri topologique aval (`networkx.topological_sort`, step5) est toujours défini. `RegionDependencyGraph.region_print_depth` expose l'arborescence.

**Sortie :** `region_id` par face + `RegionDependencyGraph`.


## Étape 4 — Champ scalaire : distance au plateau (méthode de la chaleur)

**Ce que le code fait précisément (`pipeline/step4_field_solver.py`) :**
- Sur le `TetMesh`, assemble la rigidité FEM P1 `L` et la masse concentrée `M`.
- Diffuse la chaleur depuis les nœuds du plateau : `(M + t L) u = source`, `t = h² ·
  heat_method_time_factor` (h = arête médiane).
- Normalise le flux `X = −∇u/‖∇u‖` par tétraèdre ; sur la machine actuelle (C imposé par la
  position, B vers +R seulement), `X` est projeté sur le demi-plan atteignable.
- Intègre `L φ = div X`, recale `φ` à 0 sur le plateau : `φ` est la distance au plateau à travers
  la matière, `‖∇φ‖ ≈ 1`. Les couches sont ses isosurfaces ; l'ordre d'impression par isovaleur
  croissante est la progression de la matière (une pièce suspendue vient après ce qui la porte).

**Sortie :** `ScalarField` (valeur par nœud, gradient constant par tétraèdre).


## Étape 5 — Génération des couches courbes (isosurfaces)

**Ce que le code fait précisément :**
- Choisit une série de valeurs d'isosurface espacées selon `layer_height_base` (converti en espacement de champ via la norme moyenne du gradient local, pour que l'épaisseur physique de couche reste ~constante malgré la courbure).
- Pour chaque valeur, extrait l'isosurface (marching tetrahedra sur le TetMesh — variante 3D du marching cubes, adaptée aux maillages tétraédriques).
- Intersecte chaque isosurface avec le mesh original pour obtenir les contours fermés réels de la couche (`trimesh.intersections` ou opération booléenne locale).
- Ordonne les couches selon le `RegionDependencyGraph` (une couche peut appartenir à plusieurs régions si elle les traverse — dans ce cas, elle hérite de la contrainte de dépendance la plus stricte).

**Sortie :** liste ordonnée de `Layer { contour_3d: Polyline[], region_ids: [...] }`.


## Étape 6 — Remplissage par couche (contours + infill + ponts)

**Ce que le code fait précisément (`pipeline/step6_toolpath2d.py`, FAIT) :**
- Pour chaque couche, projette chaque boucle de `contour_3d` dans son plan
  local (`local_frame` de l'étape 5, déjà calculé — pas recalculé ici).
- Classe les boucles projetées en contours "matière" / trous par
  **imbrication géométrique** (`point_in_polygon_2d`, `core/geometry_utils.py`),
  jamais par signe d'aire (limite #11) — parité généralisée à un nombre
  arbitraire de boucles imbriquées, pas seulement 2.
- Génère les périmètres par offset polygonal successif vers l'intérieur
  (`shapely.buffer(-nozzle_width * n)`), pour `wall_count` parois ; s'arrête
  proprement (sans erreur) si une zone plus fine que `nozzle_width *
  wall_count` fait disparaître l'offset avant d'atteindre le nombre nominal
  de parois — voir limite #13.
- Génère l'infill (balayage de lignes parallèles, `infill_density`/`infill_pattern`)
  à l'intérieur du dernier offset valide.
- Ordonne la trajectoire du centre vers l'extérieur (infill, puis parois
  n=wall_count...1, puis le contour original n=0 en dernier).
- Pour les boucles touchant une zone `is_bridge` (étape 2, détection de
  PRÉSENCE pas de majorité — voir limite #14) : ajoute
  `bridge_anchor_extra_walls` parois supplémentaires resserrées sur le
  périmètre (ancrage = renfort de matière sur les bords, pas densité
  d'infill augmentée — lecture directe de CLAUDE.md), puis réoriente
  l'infill pour balayer perpendiculairement à `bridge_span_direction`
  (étape 2, projeté dans le plan local de la couche).

**Sortie :** `Layer.toolpath_2d` (liste de polylignes ordonnées, repère
local **2D** — pas 3D, l'étape 8 se chargera de relever au repère machine
via `local_frame`).


## Étape 7 — Application de la carte de couleur et planification d'ordonnancement

**Ce que le code fait précisément (`pipeline/step7_scheduling.py`, FAIT) :**
- Projette `ColorMap.per_face_color` sur chaque **arête** de trajectoire :
  relève l'arête du repère local 2D vers la 3D (`step6._lift_to_3d`), puis
  interroge `trimesh.proximity.ProximityQuery.on_surface()` — **pas**
  `mesh_data.kdtree`, dont la distance n'est pas une distance à la surface
  (limite #15). Une arête à plus de `surface_color_distance` de la surface
  est de la matière **interne/invisible** et reçoit `color_id = None` :
  imprimable avec n'importe quelle couleur déjà chargée, ce qui est le
  principal levier de réduction des purges.
- Scinde les polylignes aux transitions de couleur (invariant « une
  polyligne = une couleur », dont dépend tout l'ordonnancement aval).
- Lisse la carte de couleur selon `ColorMap.merge_tolerance` : filtre
  morphologique **1D le long de la trajectoire** (paramétrée en longueur
  d'arc), tout run plus court que la tolérance étant absorbé par son voisin
  le plus long. Lissage le long du chemin et non sur le maillage, car le coût
  réel à borner est le nombre de **changements de couleur** sur le trajet —
  lecture directe : « ne jamais changer de couleur pour moins de X mm ».
- Construit un graphe `(region, color) → dépendances physiques`, en croisant
  `RegionDependencyGraph` (arêtes héritées) et le recouvrement entre couches
  successives (arêtes de support), à la granularité (couche, couleur).
- Calcule la **deadline d'accessibilité** de chaque trou candidat : rayons
  d'échappement dans les directions atteignables en 4 axes
  (`rotation_axis_circle_directions`), scellés par la première couche
  **future** rencontrée (limite #16) ; la deadline est le max sur les
  directions (la meilleure route survit le plus longtemps), moins 1.
- Bascule en purge classique si le trou n'est pas **réalisable**
  (`deadline < couche de comblement`) ou pas **rentable**
  (`hole_cost > purge_cost`), ou si `max_open_holes` serait dépassé.
- Ordonnance en earliest-deadline-first sur les trous en attente, entrelacé
  avec l'ordre topologique géométrique de base.

**Sortie :** `(list[Layer] recolorées/re-segmentées, PrintSchedule)`. Une
`PrintPass` pointe vers des **segments** `(layer_id, index)` — pas seulement
une plage de couches : une passe de couleur ne couvre qu'une partie de chaque
couche traversée. L'ordre d'exécution est exactement l'ordre de
`PrintSchedule.passes` ; l'étape 8 le consomme séquentiellement sans jamais
re-dériver la logique d'ordonnancement.


## Étape 8 — Génération de trajectoire machine 4 axes réelle

**Ce que le code fait précisément :**
- Pour chaque point de chaque `toolpath_2d`, transforme du repère local de la couche vers le repère machine global : calcule l'angle de rotation plateau nécessaire (à partir de l'orientation locale du champ à ce point) + position XYZ résultante après rotation.
- Lisse la séquence d'angles de rotation dans le temps (filtre passe-bas ou spline sur la série temporelle d'angles) pour respecter `max_angular_velocity`/`max_angular_acceleration`.
- Réordonne/interpole les points si nécessaire pour garantir une extrusion continue (pas de saut de position sans rétraction explicite programmée).

**Sortie :** `ToolpathMachine[]` — liste de `{position: Vec3, rotation_angle: float, extrusion: float, feedrate: float, color: id}`.


## Étape 9 — Simulation de collision avec le modèle 3D de la tête

**Ce que le code fait précisément (`pipeline/step9_collision_sim.py`, FAIT) :**
- `MachineProfile.nozzle_model` doit déjà être un `trimesh.Trimesh` chargé
  (ex: `file_io.mesh_loader.load_toolhead_mesh`) — step9 ne charge rien
  lui-même. Absent/invalide : `CollisionReport(nozzle_model_checked=False)`,
  jamais d'exception.
- Échantillonne `ToolpathMachine.points` selon `path_sampling_step`
  (distance cumulée, pas un point sur N).
- Pour chaque point échantillonné, positionne une copie des points de test
  du modèle (sommets **et** points le long des arêtes, voir plus bas)
  selon `position` + `axis_angles`, avec la MÊME convention de composition
  N-axes que step3/step7/step8 (rotations en ordre inverse de la liste,
  translation après) — jamais une nouvelle logique de rotation.
- Fait grandir un index spatial de la matière déjà imprimée **au fur et à
  mesure**, sous forme de capsules (segment + rayon `nozzle_diameter/2`),
  jamais une union booléenne (coût prohibitif sur plusieurs milliers de
  points). Une capsule n'est jamais ajoutée en avance sur le point testé,
  et le cordon en cours de dépôt lui-même est exclu (la buse est par
  définition confondue avec le cordon qu'elle vient de déposer — sans
  cette exclusion, chaque point en extrusion se signalerait comme collision
  contre sa propre matière, faux positif structurel mesuré concrètement).
- Teste contre cet index (distance point-capsule, formule fermée) et contre
  le plateau (`ProximityQuery`, jamais `mesh_data.kdtree` — limite #15),
  avec un broad-phase AABB commun avant tout calcul de distance précis.
- Le mandrin/mors n'a aucune géométrie disponible à ce stade : jamais testé
  silencieusement — `CollisionReport.chuck_checked` reste `False`, signal
  explicite et structuré plutôt qu'un skip invisible.

**Sortie :** `CollisionReport` — `events: list[CollisionEvent]` (point
concerné, type d'obstacle, profondeur de pénétration), plus
`chuck_checked`/`nozzle_model_checked` (booléens explicites, jamais un
silence sur ce qui n'a pas été vérifié).


### Modèle 3D de tête d'impression reçu (`res/3d/toolheads/`)

Un premier modèle réel a été reçu pour `MachineProfile.nozzle_model`,
utilisé par l'étape 9 ci-dessus. Chargement via
`file_io/mesh_loader.load_toolhead_mesh(path)` (FAIT) : `Scene` GLB ->
`Trimesh` unique en espace scène (`scene.to_geometry()`), conversion m -> mm
(×1000) — les deux pièges documentés ci-dessous sont gérés à l'intérieur de
cette fonction, aucun appelant n'a besoin de les reproduire. step9 ne charge
pas lui-même le fichier (ne connaît pas le système de fichiers) : c'est
l'appelant de `pipeline/step9_collision_sim.run()` qui doit peupler
`MachineProfile.nozzle_model` via ce loader avant l'appel.

- **Emplacement** : `res/3d/toolheads/aetheris-4/aetheris-4.glb` — un
  sous-dossier par modèle de tête (convention pour accueillir plusieurs
  têtes interchangeables à l'avenir, cohérent avec la généricité déjà
  actée pour `MachineProfile.rotation_axes`).
- **Modèle physique** : BIQU H2 V2S Lite, buse longue (dimensions
  constructeur : 72.1 × 37.9 × 80.6 mm, avec ventilateur).
- **Format** : GLB (glTF binaire). `trimesh.load(path)` retourne une
  `Scene`, **pas** un `Trimesh` direct — la géométrie est nommée
  `"toolhead"` dans `scene.geometry`. Le mesh brut (`scene.geometry["toolhead"].vertices`)
  est dans l'espace LOCAL du nœud, pas dans l'espace de la scène : ce
  fichier porte une rotation + mise à l'échelle composées sur la
  hiérarchie de nœuds (échelle nette ×0.01 entre le mesh et la racine de
  la scène, vérifié directement sur ce fichier — extentes locales brutes
  3.5×7.7×7.6 contre extentes de scène 0.076×0.077×0.035 m). **Ne jamais
  lire les sommets bruts sans appliquer la transformation de nœud** :
  utiliser `scene.dump(concatenate=True)` (ou `scene.to_geometry()`) pour
  obtenir directement un seul `Trimesh` déjà exprimé dans l'espace de la
  scène, plutôt que de reconstruire la transformation à la main
  (`scene.dump(concatenate=True)` est dépréciée depuis avril 2025 — le
  loader utilise `scene.to_geometry()`).
- **Unité : mètres dans le fichier source** (la spec glTF impose le mètre
  comme unité de base ; confirmé sur ce fichier, aucune extension
  déclarant une unité alternative). **Conversion ×1000 obligatoire** vers
  millimètres avant toute utilisation dans le pipeline, qui travaille
  entièrement en mm (même piège que `platform_z`/`merge_tolerance` déjà
  documenté ailleurs : une unité qui diverge silencieusement entre deux
  parties du système). Extentes résultantes après transform de scène et
  conversion : environ 76 × 77 × 35 mm.
- **Non étanche (`is_watertight=False`), accepté tel quel — pas de
  réparation nécessaire.** Confirmé par l'implémentation de l'étape 9 :
  la détection ne teste jamais l'inclusion/le volume, seulement des
  distances (point-capsule pour la pièce imprimée, `ProximityQuery` pour
  le plateau), aucune des deux n'exige l'étanchéité.
- **Convention de repère du modèle (à respecter par l'étape 9 lors du
  positionnement/orientation de `nozzle_model` le long de la trajectoire)** :
  - **Origine** = point de dépôt de matière (pointe de la buse) — vérifié
    géométriquement sur ce fichier : le sommet le plus proche de l'origine
    de la scène est à 0.1 mm, situé exactement à une extrémité du
    maillage, cohérent avec la convention.
  - **Axe X (rX)** = axe de rotation de la tête.
  - **Axe Y** = axe de déplacement de la tête.
  - **Axe Z** = sens de la tête ; le haut de Z pointe vers le haut
    physique de la tête.
  - Cette convention est celle du FICHIER (repère local du modèle de
    tête), à distinguer du repère `(R, C, Z, B)` de la machine — l'étape 9
    devra composer les deux (orientation de `nozzle_model` selon
    `axis_angles` du point de trajectoire courant, puis translation à sa
    position) au moment de son implémentation, non faite ici.


## Étape 10 — Ajustement en cas de collision détectée

**Ce que le code fait précisément (`pipeline/step10_adjustment.py`, FAIT — boucle jusqu'à `max_adjustment_iterations`) :**
- Consomme `CollisionReport` (étape 9) **et** `list[Step8Diagnostic]` (étape 8).
  Pour chaque problème, tente dans l'ordre :
  1. **Réorientation locale** (`_try_local_reorientation`) : candidats issus de
     `composed_rotation_angles_and_directions` (la MÊME grille que l'IK de
     step8, jamais une nouvelle méthode), essayés par écart angulaire
     **croissant** vis-à-vis de la direction de dépôt d'origine, et
     **arrêtés** dès que l'écart dépasse `max_direction_deviation_deg` —
     c'est exactement la condition « tout en respectant encore le seuil de
     surplomb ». Un candidat qui violerait `max_angular_velocity` vis-à-vis
     de ses voisins est rejeté (voir limite #24).
  2. **Ré-segmentation ciblée** : émet une `ResegmentationRequest` (point,
     plage d'angles à exclure par axe) vers un **hook injecté par
     l'orchestrateur** — l'étape 10 ne re-joue PAS elle-même step3→step8
     (voir limite #23).
  3. **Support de secours** : `AdjustmentLogEntry(method="fallback_support")`,
     log explicite.
- Après une réorientation, ne re-vérifie QUE les points modifiés
  (`_collisions_at_points`) : les capsules de l'étape 9 dépendent de
  `position` seule, jamais de `axis_angles`, donc une réorientation laisse
  la matière imprimée strictement identique et aucun point en aval ne peut
  changer d'état. Après une ré-segmentation, au contraire, passe complète
  (la trajectoire entière est invalidée).
- **Critère d'arrêt** : convergence = zéro collision restante. À
  l'épuisement des itérations avec des collisions ouvertes, celles-ci
  basculent **de force** en support de secours et `converged=False` —
  jamais d'exception (même discipline que `Step8Diagnostic` et
  `CollisionReport.chuck_checked` : une donnée explicite, pas un blocage
  qui priverait l'utilisateur de toute sortie).

**Sortie finale :** `(ToolpathMachine ajustée, AdjustmentReport)` —
`AdjustmentReport` porte `entries` (le journal), `converged`,
`unresolved_collisions` et `iterations_used`, transmis à l'UI TypeScript
pour affichage/validation avant export G-code.

---

## Décisions techniques actées et limites connues

### Dépendances / corrections d'API
- Tétraédrisation via `tetgen` (binding PyPI), accepté explicitement par l'utilisateur.
- **trimesh 5.0** a changé son API par rapport à ce qui est documenté/connu :
  - `remove_degenerate_faces()` / `remove_duplicate_faces()` n'existent plus
    → utiliser `mesh.update_faces(mesh.nondegenerate_faces())` et `mesh.update_faces(mesh.unique_faces())`
  - `trimesh.proximity.spatial.cKDTree` n'existe plus → `scipy.spatial.cKDTree` directement
  - `tetgen.tetrahedralize()` retourne un nombre de valeurs variable selon version
    → lire `tet.node` / `tet.elem` sur l'objet plutôt que la valeur de retour
- Opérations booléennes trimesh (`union`, `difference`) nécessitent `manifold3d` installé
  (`engine='manifold'`).
- `tetgen.TetGen` ne détecte PAS automatiquement les cavités internes (ex: boîte creuse) :
  sans appel explicite à `tet.add_hole(point_interieur)` avant `tetrahedralize()`, tetgen
  peut générer des tétraèdres qui traversent le vide (mesuré : ~6% des tétraèdres sur une
  boîte creuse simple). `step1_ingest._find_internal_cavities()` détecte automatiquement
  ces cavités (coques connexes à volume signé négatif via `mesh.split()`) et fournit les
  points de trou — voir limite #7.
- `scipy.sparse.linalg.cg` (scipy ≥1.9, venv du projet en 1.18.0) attend `rtol=`, pas
  `tol=` (renommé) — utilisé par `step4_field_solver._solve_linear_system`.

### Limites de conception importantes (à garder en tête pour les étapes suivantes)
1. **Le plateau d'impression n'est connu qu'à partir de l'étape 3.** `step2_overhang.py`
   seul ne peut pas distinguer une face posée sur le plateau d'une face réellement en
   surplomb dans le vide — il n'a aucune connaissance de la position du plateau.
   `step3_segmentation.py` résout ça via `platform_z` + `platform_contact_tolerance` :
   les faces au contact du plateau ont leur risque forcé à 0 et deviennent des racines
   du graphe de dépendances.

2. **`is_bridge` et `is_overhang` ne sont PAS mutuellement exclusifs.** Un plafond de
   cavité quasi-horizontal (normale vers le bas, angle≈0°) peut légitimement être les
   deux à la fois. Seule une face dont la normale pointe vers le **haut** (angle≈180°,
   dessus d'un solide plein) est exclue de la détection de pont — elle n'a par
   construction pas besoin d'être supportée par en dessous.

3. **`bridge_max_gap_distance` doit être calibré relativement à l'échelle du modèle**
   (idéalement en fonction de `layer_height_base`), jamais une valeur absolue universelle.
   Un défaut mal calibré peut soit tout marquer comme pont, soit ne rien détecter.

4. **La compatibilité face/direction (step3) se teste comme un surplomb, pas un
   alignement.** Une face est compatible avec une direction candidate `d` si elle
   n'est PAS en surplomb selon `-d` (même critère que step2, généralisé à une direction
   quelconque) — PAS si sa normale est proche de `d` lui-même. Erreur initiale corrigée :
   tester l'alignement direct bloquait à tort toutes les faces de murs verticaux
   (perpendiculaires à `d`), empêchant toute fusion de région même sur un cube trivial.

5. **Le graphe de dépendances entre régions doit se baser sur une référence globale
   (`up_vector` machine), jamais sur une moyenne locale par paire de régions.**
   Une référence locale par paire n'est pas transitive et produit des cycles
   artificiels dès qu'un ensemble de régions forme un anneau (typiquement une sphère
   segmentée en tranches autour d'un axe) — chaque paire semblait localement cohérente,
   mais l'ensemble ne l'était pas.

6. **Directions candidates 4 axes vs 5 axes** (`core/geometry_utils.py`) : une machine
   à 1 seul DOF de rotation n'atteint qu'un **cercle** de directions autour de
   son axe de rotation (`rotation_axis_circle_directions`), pas une sphère complète
   (`fibonacci_sphere_directions`, réservée à une orientation libre type 5 axes).
   `SegmentationParams.rotation_axes` (une liste de `RotationAxisSpec`, voir
   limite #18) généralise ce choix à N axes composés, chacun avec sa propre
   plage — le cas à 1 DOF ci-dessus reste un cas particulier valide (liste à
   1 élément, plage complète), pas remplacé.

8. **`TetMesh.tet_adjacency` (adjacence nœud-nœud) ne doit PAS servir de poids au
   Laplacien discret**, malgré ce que suggérait son commentaire d'origine dans
   `step1_ingest.py`. Un Laplacien de graphe (poids uniformes) ignore la géométrie
   réelle des tétraèdres — l'étape 4 assemble la vraie matrice de rigidité FEM P1
   (élément par élément, `step4_field_solver._assemble_stiffness_matrix`).
   `tet_adjacency` reste utile comme outil topologique (propagation d'influence de
   surplomb par sauts d'adjacence), pas comme poids de la forme bilinéaire.

9. **Un système Laplacien pur est singulier** (le champ constant est toujours solution)
   et **un terme de direction limité aux zones à risque laisse le reste du champ
   indéterminé** (donc nul, une pièce sans surplomb n'aurait aucune couche). L'étape 4
   résout les deux : ancrage automatique du plateau (pénalité molle de poids élevé sur
   les nœuds au contact du plateau, `platform_anchor_weight`) + terme de direction de
   gradient appliqué à **tous** les tétraèdres avec un poids de base
   (`gradient_baseline_weight`), renforcé près des surplombs
   (`gradient_overhang_weight`). Piège écarté : déduire les nœuds du plateau de la
   topologie du graphe de dépendances ("régions sans arête entrante") est faux dès que
   la segmentation ne produit qu'une seule région (tout objet simple) — cette région
   unique serait alors à la fois racine et couvrirait tout le modèle. La bonne source
   de vérité est `RegionDependencyGraph.platform_contact_faces` (booléen par face,
   calculé par step3, qui seul connaît la position du plateau).

10. **Marching tetrahedra (étape 5) doit orienter ses triangles à partir des valeurs de
    champ, jamais de l'ordre local des sommets du tétraèdre** (`tetgen` ne garantit
    aucune convention de handedness cohérente d'un tétraèdre à l'autre — même raison
    que le `abs(volume)` de l'étape 4). Convention retenue : normale du triangle
    pointant du côté "en-dessous de l'isovaleur" vers le côté "au-dessus", dérivée des
    sommets de valeur min/max du tétraèdre — garantie correcte car le champ est
    linéaire par tétraèdre (P1).

11. **Le sens de parcours (winding) d'un contour extérieur et d'un trou intérieur
    N'EST PAS opposé** via `trimesh.Trimesh.outline()`, contrairement à l'intuition
    initiale (théorème de Stokes mal appliqué). Vérifié directement sur un anneau plat
    synthétique à normales cohérentes : les deux boucles donnent le même signe d'aire
    projetée. La distinction extérieur/trou doit passer par un test d'imbrication
    géométrique (point-in-polygon), pas par le signe d'une aire signée.

12. **L'espacement des isovaleurs (étape 5) ne peut être qu'approximativement
    constant.** `‖∇φ‖` varie spatialement (le champ résout un compromis
    lissage/direction/contraintes, pas une équation eikonale `‖∇φ‖=1`) ; une seule
    valeur d'espacement est dérivée de la **médiane pondérée par volume** de
    `‖∇φ‖` sur tout le `TetMesh` (pas la moyenne brute, faussée par le gradient
    localement élevé près de l'ancrage plateau à poids élevé).

13. **Une paroi plus fine que `nozzle_width * wall_count` ne peut pas recevoir
    `wall_count` parois pleines** — géométriquement impossible, forcer le nombre
    nominal produirait des offsets dégénérés/auto-intersectants. L'étape 6
    (`_offset_walls`) s'arrête proprement dès qu'un offset devient vide ou
    résiduel (aire `< degenerate_area_ratio * aire_initiale`), sans lever
    d'erreur ni log — le nombre de parois réellement généré pour cette zone
    est simplement `< wall_count`, cohérent avec la convention déjà en place
    (pas d'objet placeholder, moins d'éléments dans la liste retournée).

16. **Un trou temporaire n'est scellé que par la matière déposée APRÈS lui**
    (étape 7). Le calcul de la deadline d'accessibilité ne doit compter comme
    bloquants que les échantillons d'indice de couche **strictement supérieur**
    à celui du trou. Sans cette restriction, un trou de **surface** — le cas
    d'usage même, la couleur étant une propriété de surface — voit son rayon
    d'échappement quitter le solide immédiatement, tous les échantillons
    retombent sur sa propre couche, et le report est refusé quel que soit le
    paramétrage : le mécanisme de trou devient structurellement inatteignable.

17. **Une `ColorMap` par face ne peut pas représenter une frontière plus fine
    que le maillage.** `per_face_color` associe UNE couleur à une face
    entière : sur un cube brut (12 triangles, faces latérales couvrant toute
    la hauteur), une peinture « z<0 vs z>0 » produit en réalité une découpe
    diagonale. Toute fixture de test (ou tout modèle utilisateur) dont la
    peinture suit une frontière géométrique doit être suffisamment subdivisé
    — sinon on teste/affiche une frontière que la donnée n'encode pas.

18. **Le système doit décrire une cinématique de rotation ARBITRAIRE (N axes
    composés), jamais une seule machine figée en dur.** Corrige une hypothèse
    posée dès l'étape 3 (à l'origine : un seul DOF de rotation supposé,
    `rotation_axis_circle_directions` comme unique alternative à
    `fibonacci_sphere_directions`). `core/mesh_types.RotationAxisSpec`
    (axe + plage angulaire) et `core/geometry_utils.composed_rotation_directions`
    généralisent à N axes ordonnés, chacun avec sa propre plage — SANS
    retirer ni remplacer les deux fonctions existantes, qui restent des cas
    particuliers valides et testés (0 axe → sphère libre, 1 axe à plage
    complète → dégénère exactement vers le cercle à 1 DOF).
    `SegmentationParams`/`Step7Params` exposent `rotation_axes: list[RotationAxisSpec]`
    ; `MachineProfile.rotation_axes` décrit la cinématique réelle de la
    machine utilisée. **Résultat mesuré (pas supposé) pour la machine
    actuelle en particulier**, un cas d'usage parmi d'autres possibles, pas
    la référence du type : `C` balaie un tour complet et `B ∈ [0°, 180°]`
    (jamais négatif), donc l'ensemble des directions de dépôt atteignables
    couvre la **sphère complète** — exactement comme la paramétrisation
    standard d'une sphère par (longitude, colatitude) où la colatitude ne
    parcourt que `[0°, 180°]`. Une cinématique à 2 axes tous deux restreints
    (ex: chacun `[0°,90°]`) ne couvre en général PAS la sphère complète —
    vérifié par test (`tests/test_geometry_utils.py`), pour ne jamais
    retomber sur l'hypothèse cachée « 2 axes ⇒ sphère complète ».

19. **Singularité géométrique au centre du plateau (`R=0`) : `C` est
    indéterminé**, comme au pôle d'un système de coordonnées polaires
    (`atan2(0,0)` n'est pas défini). L'étape 8 doit hériter le `C` du point
    de trajectoire précédent plutôt que de recalculer un `atan2` instable —
    cohérent avec le lissage temporel des angles déjà prévu pour respecter
    `max_angular_velocity`/`max_angular_acceleration`. Ne dégrade rien
    géométriquement : un point exactement centré ne contraint par
    construction aucune orientation de matière réelle.

23. **Une ré-segmentation ne peut PAS être locale, contrairement à ce que
    laissait entendre la formulation initiale de l'étape 10 (« relance la
    segmentation pour cette zone uniquement »).** Changer la direction
    d'impression d'une région change son `target_direction` (step3), donc le
    terme de gradient du système **statique** de step4 — pas la partie
    contrainte, la partie structurelle — donc le champ scalaire, résolu
    **globalement** par FEM, donc les couches, les trajectoires,
    l'ordonnancement et les angles. Aucune ré-exécution partielle n'est
    géométriquement cohérente. Conséquence actée : `step10` n'appelle PAS
    lui-même step3→step8 (il lui faudrait `mesh_data`, `tet_mesh`,
    `overhang_map`, le graphe de régions, le champ, les couches, la
    `ColorMap` et les params de 7 étapes — il deviendrait l'orchestrateur,
    rôle qui revient à `main.py`). Il émet une `ResegmentationRequest` vers
    un **hook injecté** ; l'appelant décide comment re-jouer la chaîne et
    rend une nouvelle `ToolpathMachine`. `step10` reste une fonction pure et
    testable avec un hook synthétique.

24. **Une réorientation ponctuelle peut casser le lissage temporel que
    l'étape 8 venait d'appliquer.** Modifier l'angle d'un seul point
    ré-introduit exactement le saut de vitesse angulaire que
    `step8._clamp_rate` avait éliminé — « résoudre » une collision en
    produisant une trajectoire mécaniquement inexécutable. `step10._violates_rate_limits`
    rejette donc tout candidat violant `max_angular_velocity` vis-à-vis de
    ses voisins immédiats. **Limite assumée non résolue** : seule la vitesse
    est re-vérifiée, pas l'accélération (qui exigerait de raisonner sur trois
    points consécutifs et de re-lisser toute la série, donc de rejouer
    step8).

25. **Corollaire de la limite #21 pour la conception des FIXTURES de test.**
    Comme l'étape 9 n'échantillonne que sommets et arêtes, un obstacle qui
    traverse l'**intérieur d'une face** peut être écarté par une simple
    rotation, sans que la géométrie ait « vraiment » bougé. Mesuré en
    écrivant les tests de step10 : une tête large et symétrique censée
    rendre la collision irrésolvable voyait sa collision disparaître à 30°.
    Une fixture qui doit rendre une collision **irrésolvable par rotation**
    doit être invariante par cette rotation — typiquement un obstacle porté
    par l'axe de rotation lui-même (`tests/test_step10_adjustment.py::_axis_straddling_toolhead`).

27. **La boucle de distance point-capsule de `_check_printed_part` était en
    pur Python (`for segment: for vertex`), un problème de performance
    invisible sur les petits volumes de test unitaires mais dominant sur
    une vraie trajectoire.** Avec ~1500 points d'échantillonnage de tête
    (limite #26) et plusieurs dizaines de capsules candidates par point de
    trajectoire testé, la boucle Python interprétée devenait le coût
    dominant de l'étape entière. Corrigé par `_points_to_segment_distances`
    (nouvelle fonction), qui vectorise le calcul sur TOUS les points de
    test d'un candidat en une seule opération numpy au lieu d'une boucle
    Python par sommet — même formule mathématique que
    `_point_to_segment_distance` (conservée, toujours utilisée ponctuellement),
    appliquée en bloc. Gain mesuré sur un cas synthétique de 1000 points de
    trajectoire : le run ne terminait pas sous 139s de CPU (profilage tué
    avant complétion) avant le fix, **10.86s** après. Sur le modèle de test
    réel (`res/test_models/l_shape_overhang.stl`, 9391 points de trajectoire,
    3482 échantillonnés par `path_sampling_step=2.0`), le débit initial
    mesuré est de ~15 points de trajectoire testés par seconde.
    **CORRECTION (mesure ultérieure plus longue)** : la première rédaction
    de cette limite concluait « stable, sans dégradation, comportement
    linéaire confirmé par mesure » — c'était **faux**, conclusion tirée
    d'un échantillon trop court (les premières centaines de points). Une
    mesure poussée jusqu'au 1750ᵉ point montre une dégradation nette et
    continue : **16.5 → 10.3 → 7.7 → 6.0 → 4.9 → 4.2 → 3.7 pts/s** à
    mesure que l'index passe de quelques centaines à ~2700 capsules. La
    grille spatiale de `_PrintedGeometryIndex` réduit le nombre de
    candidats mais ne le borne pas. Point ouvert, à traiter après la zone
    de garde buse (voir « Prochaine étape »).

35. **Le mécanisme `bgmesh`/`metric` de TetGen (maillage de fond portant un
    champ de taille cible par nœud) est INUTILISABLE dans la version
    installée du binding (`tetgen==0.8.4`) — testé et abandonné, pas
    seulement contourné.** Investigation complète avant d'écarter :
    - Le contournement de la dépendance manquante `pyvista` (le paramètre
      `bgmesh=` du binding exige un objet `pyvista.UnstructuredGrid`,
      absent du projet) fonctionne bien : `TetGen._tetgen.load_bgmesh_from_arrays(v,
      t, target_size)` sur le binding compilé accepte des numpy arrays
      purs.
    - **Mais** la docstring du binding précise explicitement que
      `switches` (une string non vide, utilisée dans ce projet pour `-a`/
      `-q`) *"overrides all keyword arguments except for bgmesh
      associated keyword arguments"* — `metric=True` est donc
      SILENCIEUSEMENT IGNORÉ dès qu'un `switches=` non vide est passé. Une
      première mesure "positive" (773 tétraèdres au lieu de 307) était un
      FAUX POSITIF : l'écart venait entièrement d'une valeur `-a`
      différente entre les deux essais comparés, pas du bgmesh — confirmé
      en reproduisant exactement 773 tétraèdres avec le même `-a` et AUCUN
      bgmesh chargé.
    - Forcer `-m` malgré tout (littéralement dans la string switches, à
      n'importe quelle position ; ou via kwargs purs sans string switches
      du tout) **fait planter TetGen nativement (segfault) de façon
      reproductible** sur `l_shape_overhang.stl`, testé sous au moins 6
      variantes (avec/sans qualité, avec/sans `add_hole`) — échec
      systématique à chaque fois.
    - **Aucun contournement trouvé** : ce mécanisme reste hors d'usage pour
      ce projet avec cette version de tetgen. Le raffinement local du
      maillage près du plateau utilise `TetGen.add_region` (deux régions
      distinctes) à la place, mécanisme alternatif fonctionnel et retenu.

41. **Repère visuel demandé pour le GLB de visualisation** : une flèche
    verte indiquant le "haut" d'impression (du plateau vers le haut),
    ajoutée dans `file_io/visualization_export.py`
    (`_up_arrow(mesh_data)`). Orientée vers +Z (convention déjà actée dans
    tout ce fichier : le plateau est à `z=platform_z`, l'impression
    progresse en Z croissant) — pas de nouveau champ sur `PipelineResult`
    pour porter `up_vector` explicitement, cette convention suffit. Taille
    proportionnelle à l'étendue du mesh source (pas une taille fixe,
    resterait illisible ou disproportionnée selon l'échelle du modèle),
    positionnée à côté de la bounding box pour ne pas se superposer
    visuellement au modèle.

49. **Fix de la limite #46 (trajets à vide sans routage) — palliatif
    lift-and-traverse validé par collision, PAS un vrai path-planning 3D,
    avec un taux de résidu non négligeable sur géométrie difficile mesuré
    et documenté honnêtement plutôt que masqué.**

    ### Mécanisme retenu
    Décidé explicitement avec l'utilisateur (deux options soumises : un vrai
    path-planning 3D type A*/RRT dans l'espace libre, écarté comme chantier
    disproportionné pour cette session, vs un palliatif lift-and-traverse
    plus simple mais documenté comme limite à améliorer plus tard —
    l'utilisateur a choisi le second, avec la condition explicite que ce
    soit marqué comme limite du slicer, pas présenté comme une solution
    complète). Nouveau module
    `pipeline/step8b_travel_routing.py::route_long_travels(toolpath,
    machine_profile, sim_params, params) -> (ToolpathMachine,
    list[TravelRoutingDiagnostic])`, positionné dans l'orchestration ENTRE
    step8 et step9 (jamais dans step8 lui-même) : il a besoin de la
    trajectoire ENTIÈRE et de `MachineProfile.nozzle_model` pour valider un
    lift par vraie détection de collision contre la matière déjà imprimée,
    alors que step8 produit la trajectoire point par point sans connaître
    l'état de la matière en avance.

    Pour chaque trajet à vide (`extrusion=0.0`) de longueur >
    `TravelRoutingParams.length_threshold` (défaut 5.0mm), remplace la
    ligne droite d'origine par un chemin en 3 segments — lever
    verticalement de `lift_height` (défaut 3.0mm) au-dessus du point de
    départ, traverser horizontalement à cette hauteur (`traverse_samples`
    points intermédiaires, défaut 1 — un seul point médian), redescendre —
    validé par les MÊMES primitives de détection que step9
    (`_positioned_toolhead_vertices`, `_PrintedGeometryIndex`,
    `_grow_printed_index_to`), jamais une seconde implémentation. Si le
    lift standard collisionne encore, la hauteur double successivement
    jusqu'à `lift_height * max_lift_multiplier` (défaut 4.0, donc plafond
    12mm) avant d'abandonner. Nouvelle dataclass
    `TravelRoutingDiagnostic` (`core/mesh_types.py` : `point_index`,
    `travel_length`, `max_lift_tested`, `still_colliding`) — signal
    explicite quand même la hauteur plafond reste en collision, jamais un
    succès silencieux faux (même discipline que `Step8Diagnostic`,
    `CollisionReport.chuck_checked`).

    Opt-in : `PipelineParams.travel_routing: Optional[TravelRoutingParams] =
    None` par défaut — `None` = comportement historique inchangé (les
    trajets à vide restent des lignes droites non routées), même discipline
    que `junction_refinement_resolution`. **`main.py::run_pipeline` ET
    `scripts/export_slicing_preview.py::run_slicing_only` orchestrent TOUS
    LES DEUX l'appel à `route_long_travels`** entre step8 et step9/la fin du
    slicing — synchronisés dès l'implémentation initiale cette fois (piège
    déjà rencontré pour `junction_refinement_resolution`, limite #47/#48,
    évité ici par vigilance directe plutôt que découvert après coup).

    ### Diagnostic de performance rencontré ET résolu pendant l'implémentation
    Version initiale (`traverse_samples=3`, `toolhead_sample_spacing`
    hérité de step9 = 5.0mm, `grid_cell_size_multiplier` hérité de step9 =
    4.0) : **ne terminait PAS** en plusieurs minutes sur
    `l-triple-axis-overhanged-shape.stl` (15766 points de trajectoire, 127
    trajets à vide longs détectés) — plusieurs process tués manuellement
    après timeout, y compris sous `cProfile` (le surcoût de profilage
    aggravait encore le problème).

    Diagnostic par script instrumenté jetable (mesure directe minute par
    minute, jamais supposée) : le temps était dominé à >99% par
    `_path_collides`, et **corrélé directement à la taille de l'index de
    matière déjà accumulée, pas au nombre de trajets lui-même** — mesuré
    3.06s cumulés au 20ᵉ trajet testé (index à 5005 capsules), 25s cumulés
    au 80ᵉ trajet (index à 12089 capsules), une dégradation continue
    cohérente avec la même pathologie déjà documentée aux limites #27/#30
    pour step9/step10. Cause racine identifiée par lecture directe de
    `_PrintedGeometryIndex.candidates_near` : sa grille spatiale (taille de
    cellule = `nozzle_diameter * grid_cell_size_multiplier` = 0.4×4.0 =
    1.6mm par défaut) est dimensionnée pour des requêtes PONCTUELLES autour
    d'un seul point de trajectoire (l'usage de step9/step10) — une requête
    de routage porte sur l'AABB d'un CHEMIN ENTIER (plusieurs dizaines de
    mm avec le lift), dilatant la recherche sur des centaines de cellules
    et retournant la quasi-totalité de l'index accumulé comme candidats à
    chaque appel.

    **Fix, mesuré à chaque étape** :
    - `TravelRoutingParams.routing_grid_cell_size_multiplier = 40.0`
      (dédié, dix fois plus grossier que le défaut de step9 —
      toujours un sur-ensemble EXACT des vrais candidats, jamais de faux
      négatif géométrique, seulement moins de filtrage en amont).
    - Vectorisation de `_path_collides` : un seul passage
      broad-phase+boucle de candidats sur l'UNION de tous les points
      positionnés du chemin candidat, au lieu d'un appel
      `_check_printed_part` complet (avec sa propre requête
      `candidates_near`) par point individuel du chemin.
    - `traverse_samples` réduit de 3 à 1 (un seul point médian de
      traverse — les deux coudes suffisent à couvrir la grande majorité
      des cas, chaque point supplémentaire multipliant directement le
      coût).
    - `routing_toolhead_sample_spacing = 15.0` dédié (plus grossier que le
      5.0mm hérité de step9 — compromis délibéré : le routage teste bien
      plus de candidats au total que la simulation ponctuelle de step9,
      donc chaque appel individuel doit rester bon marché).

    Résultat isolé mesuré (script instrumenté, même modèle) : **127
    trajets routés en 38.22s**, contre non-terminaison avant (aucune limite
    de temps atteinte en plusieurs minutes).

    ### Mesure de bout en bout (`scripts/measure_travel_routing_effect.py`, conservé)
    Sur `l-triple-axis-overhanged-shape.stl` (15766 points, 127 trajets à
    vide longs), avec les défauts finaux (`traverse_samples=1`) :

    | | avant routage | après routage |
    |---|---:|---:|
    | points de trajectoire | 15766 | 16020 |
    | plus grand saut à vide (segment consécutif) | 24.33mm | 22.16mm |
    | trajets encore en collision (`TravelRoutingDiagnostic`) | — | **76/127 (60%)** |
    | temps step1→8 (+ routage) | 3.5s | 64.4s |

    ### Limite honnête à documenter, pas à masquer
    Le taux de résidu (76/127, soit 60% des trajets longs encore signalés
    en collision même à la hauteur de lift plafond de 12mm) est
    **significatif et attendu au vu de la conception délibérément simple du
    mécanisme** — ce modèle a des surplombs sur 3 axes non coplanaires
    (géométrie tourmentée par construction), et un simple lift vertical ne
    peut pas contourner une pièce qui referme l'espace libre au-dessus de
    la trajectoire sur une hauteur supérieure au plafond testé. C'est
    exactement le cas que `TravelRoutingDiagnostic` existe pour signaler
    plutôt que masquer (même discipline que `Step8Diagnostic`,
    `CollisionReport.chuck_checked` : une donnée explicite sur ce qui n'a
    pas pu être vérifié/résolu, jamais un silence). Le mécanisme reste un
    progrès mesuré (plus long trajet passe de 24.33mm à 22.16mm sur ce
    modèle, et une fraction des trajets EST effectivement routée sans
    collision), mais **ce n'est ni un vrai path-planning 3D, ni une
    garantie d'absence de collision sur toute géométrie** — limite du
    slicer explicitement actée avec l'utilisateur, à améliorer plus tard
    (piste déjà identifiée et écartée pour cette session : A*/RRT dans
    l'espace libre autour de la pièce imprimée).

    ### Tests
    `tests/test_step8b_travel_routing.py` (7 tests nouveaux :
    `_travel_boundaries`, trajet court laissé inchangé, lift sans modèle de
    tête, lift qui évite une capsule imprimée réelle, diagnostic émis
    quand même le lift plafond collisionne, points en extrusion jamais
    modifiés) + 2 tests d'orchestration dans `tests/test_main_pipeline.py`
    (`travel_routing=None` ne déclenche jamais `route_long_travels`,
    `travel_routing` actif route bien un trajet long avant toute
    simulation de collision) — **206/206 tests passants** (197 avant cette
    session + 9 nouveaux). Non-régression : `travel_routing=None` reste le
    défaut de `PipelineParams`, chemin de code jamais emprunté sans
    activation explicite.


59. **Cinématique : sur la machine ACTUELLE, C est IMPOSÉ par la position du
    point (constaté le 2026-09-28, en préparant P5 du lissage non-planaire).**
    La machine n'a que R, C, Z et B, sans axe Y (confirmé par l'utilisateur).
    La tête ne balaie donc qu'un demi-plan vertical : pour atteindre un point,
    le plateau doit tourner de C = angle polaire du point. B reste le SEUL
    degré d'orientation. Il penche dans le plan radial, d'un seul côté
    (B ∈ [0°, 185°], normale de couche vers +R).

    Jusqu'ici, l'étape 8 calculait C pour l'ORIENTATION, sans tenir compte
    de la position, ce que la machine ne peut pas exécuter. La limite #18
    (« sphère complète couverte ») n'est vraie que si C est libre. Cela
    suppose un axe linéaire de plus, prévu à l'upgrade après Maker Faire.

    Mesuré sur les couches de lissage non-planaire (recherche), avec C
    imposé, part des surfaces lissées où la tête reste à moins de 5° de la
    normale :
    - Benchy : 7 % ;
    - dino : 11 %.

    - **Constat mesuré** : « orientation atteignable » ≠ « sans collision ».
      Sur le Benchy centré, avec w = 0,5, le nombre de contacts passe de 83
      à 903 (plateau).
    - **Correctif retenu** (recherche, `full_step89.py`) : borner B selon la
      hauteur au-dessus du plateau (point le plus bas de la vraie tête),
      puis réorienter les points en contact par paliers de tolérance
      (45 / 90 / 181°).
      - Benchy : de 83 à 9 contacts (pénétration ≤ 0,22 mm, dans la marge).
      - dino : de 217 à 0.
      - En échange, environ 100 points par modèle sont déposés à plus de
        45° de leur normale.

    **Correction** :
    - `MachineProfile.position_coupled_axes` : `[0]` pour la machine
      actuelle, `[]` quand un axe linéaire libère C ;
    - `MachineProfile.radial_direction` : direction +R à angles nuls,
      `(0, −1, 0)` pour la machine actuelle, puisque B > 0 autour de X
      penche la normale vers −Y.
    - L'étape 8 (`_resolve_all_axis_angles_coupled`) :
      - tire l'angle couplé de la position (`coupled_axis_angle`) ;
      - résout les axes restants sur la direction ramenée dans le repère de
        cet axe ;
      - émet un `Step8Diagnostic` quand la cible sort de ce qu'ils
        atteignent ;
      - ne lisse jamais l'angle couplé, car le lisser déplacerait le point
        atteint.
    - Ce chemin était en Python dans les deux positions de `APEX_USE_CPP`
      (84 s sur death) ; il passe par l'IK native depuis le pas 17 (limite #56).
    - `main._default_machine_profile` décrit désormais la machine actuelle.

    **Reste non traité, à reprendre** : les directions candidates de step3
    et step7 parcourent encore toute la sphère sans tenir compte de la
    position. **Step10 corrigé le 2026-10-01** : en machine couplée,
    `_try_local_reorientation` ne propose plus que des candidats à C FIGÉ
    (`_coupled_candidates`, grille 4× plus fine sur B seul) — la grille
    produit complète proposait des C inexécutables. Et quand step9 est en
    mode rapide (`fast_kdtree`), step10 revérifie avec le même contrôleur
    (`_fast_checker`) au lieu de l'index exact par capsules (incohérence qui
    rendait step10 interminable sur death.stl). Côté step4, les directions
    cibles sont désormais projetées sur le demi-plan atteignable (limite #64).

    **Limites de la tête et troncature des angles (même jour).**
    L'utilisateur a donné les limites de B : 750 °/s et 300 °/s². Elles sont
    dans `main._default_machine_profile`. Il a aussi confirmé que son
    firmware COORDONNE les axes : quand un axe atteint sa limite, tout le
    mouvement ralentit et l'angle reste exact.

    Le `_clamp_rate` de l'étape 8 modélisait l'inverse : il tronquait
    l'angle, donc la tête prenait du retard. Mesuré en recherche : jusqu'à
    62° d'écart à la normale sur la coque du Benchy.

    La troncature est donc désactivée par défaut
    (`Step8Params.clamp_angle_rates = False` ; `True` rend l'ancien
    comportement). Le bon modèle est un **ralentissement de l'avance**. Il
    est déjà écrit en recherche (`skin_print.feed_plan` : vitesse angulaire
    plafonnée par passes avant et arrière, remise à zéro à chaque changement
    de sens ; ×1,45 sur le temps des passes de lissage du Benchy à 25 mm/s).
    Il est à porter dans l'export G-code.

    - **Cause** : un nœud de la SURFACE à 2,6e-5 et 9,2e-5 mm du plan, sous la
      tolérance de soudure (`weld_tolerance_factor` × diagonale ≈ 8e-5 mm).
      Les points de coupe de ses arêtes fusionnent ; on obtient des triangles
      dégénérés, dont les arêtes répétées (multiplicité 6 à 22) cassent
      `outline()`. La boucle sort ouverte et elle est jetée.
    - **Un nœud intérieur ne casse rien** : seul le bord se rompt.
      L'hypothèse initiale, « plan passant exactement par un sommet », était
      fausse (aucun nœud à moins de 1e-6 mm).
    - **Correctif** : `_extract_closed_contours` retire les triangles
      dégénérés avant le tracé de bord. Sans aire, ils ne portent aucun bord.
    - **Test** : `test_node_within_weld_tolerance_of_plane_keeps_whole_contour`,
      contrôlé comme ÉCHOUANT sans le correctif.
    - Benchy : 0 couche en écart sur 313 après correction. Les sorties de
      recherche produites avant (P3, P5) sont à refaire.

61. **Export G-code (`file_io/gcode_export.py`, écrit le 2026-09-30).** Écrit
    de zéro, sans reprendre les générateurs de formes fixes `scripts/gen_*`.
    Seuls les faits machine qu'ils portaient sont repris, vérifiés avec
    l'utilisateur.
    - **Axes** :
      - firmware Y = `atan2(y, x)` = `C − 90` avec la convention d'Apex,
        recalculé depuis la position ; l'accord avec `axis_angles[0]` est
        mesuré (`GcodeReport.c_mismatch_deg`) ;
      - B firmware = B d'Apex ;
      - Y est déroulé, et hérité au centre.
    - **Pivot** : le firmware commande le pivot de B, à 45 mm de la pointe.
      On a X = R + L sin B et Z = z + L (cos B − 1), avec comme référence la
      pointe à B = 0. **Confirmé par l'utilisateur (2026-09-30)** : le zéro Z
      du firmware est atteint quand la buse touche le plateau, donc aucun
      décalage de 45 mm à ajouter. Le signe est vérifié par un test contre
      `step9._positioned_toolhead_vertices`.
    - **Vitesses validées par l'utilisateur (2026-09-30)** : impression
      25 mm/s, trajets à vide 80 mm/s, 195 °C.
    - **GPU : Vulkan** (décision de l'utilisateur, 2026-09-30), pour couvrir
      aussi macOS et les cartes AMD/Intel.
    - **F** est recalculé par segment pour garder la vitesse de la pointe (le
      firmware traite Y et B comme linéaires), puis ralenti par `feed_plan`
      aux limites de vitesse et d'accélération de B.
    - **Trajets à vide** :
      - rétraction de 0,8 mm à 35 mm/s (décision de l'utilisateur) ;
      - dégagement sur UN seul axe, choisi selon l'orientation de la tête :
        Z (+ pour une tête vers le bas, − pour une tête retournée) ou R dans
        le sens du corps ;
      - un dégagement qui descendrait la tête sous la marge du plateau
        bascule sur l'autre axe ;
      - un contrôle optionnel (`travel_check`) compte les dégagements en
        contact.
    - **Signe de B (2026-10-06) — CORRIGÉ le jour même.** `GcodeParams.b_sign = +1`
      (B firmware = B d'Apex) : B POSITIF penche le corps de la tête vers +R (la buse
      vers le centre), B NÉGATIF envoie la pointe vers l'EXTÉRIEUR (croquis de
      l'utilisateur et aperçu de `slicer_ux`). Une première version posait `−1` : c'était
      un CONTRESENS sur la phrase « B négatif = vers l'avant », qui décrit le sens de
      l'axe et non une inversion à appliquer ; les limites mesurées de B sont [0°, 185°]
      (jamais négatif), cohérentes avec +1. Tous les G-code produits avec `−1` (B
      négatif) sont inutilisables. Géométrie (X, Z, collisions) inchangée, en B d'Apex.
      Présent dans les chemins Python et natif ; `firmware_to_machine` l'inverse. Axe de
      rotation actif à l'impression : **Y** (pas C). Les fichiers RRF
      (`research/printers/aetheris-4/gen0/v0.0.4`) décrivent une imprimante bac à sable
      aux limites volontairement larges.
    - `parse_gcode` et `firmware_to_machine` servent à l'aperçu
      (`scripts/preview_gcode_file.py`) et au contrôle aller-retour.
    - Tests : `tests/test_gcode_export.py`.

    **Suite du 2026-10-01 (session autonome)** :
    - **Le firmware RRF ne fait AUCUNE cinématique** (utilisateur) : ni
      polaire, ni pivot. Tout est dans le slicer, pour qu'une nouvelle machine
      ne demande jamais de mise à jour d'imprimante.
    - **Divergence NON RÉSOLUE avec la prévisualisation UX de l'utilisateur**
      (qui fait référence, conçue pour la machine) : elle montre des
      trajectoires sous le plateau et qui dépassent sur le Benchy, là où le
      modèle ci-dessus place la pointe à 0,04 mm au plus bas et la vraie tête
      jamais sous 1 mm. Lire X/Z comme la pointe (aucune compensation) reproduit
      exactement ses symptômes (outil `scripts/gcode_tilt_profile.py`). **Aucun
      G-code ne doit être imprimé avant que la convention soit tranchée.**
    - Modèle de tête **paramétré** en attendant : `GcodeParams.tilt_sign`
      (sens du déport radial du pivot) et `z_reference` (`"tip"` : défaut
      d'origine ; `"pivot"` : Z firmware = hauteur du pivot).
    - **Pointe au sommet du cordon** (`tip_at_bead_top`, défaut vrai) : step5
      place la première isovaleur à v_min + δ/2, donc chaque point de
      trajectoire est au MILIEU de sa couche ; la pointe est relevée de h/2 le
      long de l'axe de la tête. Fait géométrique, pas un réglage.
    - **Pivot à 46 mm** (mesure de l'utilisateur, 2026-10-01 ; 45 avant).
    - **La pointe quittait sa trajectoire PENDANT les mouvements** (trouvé le
      2026-10-01 grâce aux rendus de l'utilisateur) : le firmware interpole X, Z
      et B linéairement ; avec le pivot à 46 mm, un saut de B entre deux points
      proches (0 → 58°, ou 3,8 ↔ 74° en alternance sur le dino) fait balayer à
      la pointe un arc de plusieurs mm EN EXTRUDANT. Mesuré au milieu de chaque
      mouvement : Benchy 1 953 mouvements extrudés à plus de 1 mm (max 26 mm),
      dino 3 852 (max 16,6 mm). Les contrôles précédents ne regardaient que les
      EXTRÉMITÉS des mouvements. **Correctif** : `GcodeParams.max_tip_deviation`
      (0,05 mm, défaut) découpe chaque mouvement par bisection jusqu'à ce que le
      milieu interpolé par le firmware reste sur la droite visée, l'extrusion
      étant répartie au prorata. Après : 0 mouvement à plus de 0,1 mm sur les
      deux pièces, +9 % de lignes. Le chemin natif (pas 24) ne découpe pas encore :
      il n'est emprunté que si `max_tip_deviation=None`. Test :
      `test_moves_are_split_so_the_tip_follows_the_line`. Reste en amont :
      l'alternance de B elle-même (réorientation de la recherche), qui ralentit
      l'impression.
    - **Garde-fous** dans `GcodeReport` (pointe et tête les plus basses,
      étendues X/Y/Z firmware — Y cumule jusqu'à 35 tours sur le Benchy) ;
      `strict` refuse d'écrire un fichier dont la pointe passe sous 0.

62. **step9 : contrôle rapide par arbre k-d (`CollisionSimParams.fast_kdtree`,
    opt-in, 2026-09-30).** Il est porté de la recherche (`full_step89.py
    --fast`, où il ramenait 30 min à environ 1 min sur le Benchy).
    - **Principe** : la matière est densifiée le long des cordons
      (`fast_densify_step` = 0,25 mm) dans UN arbre construit une fois.
      Chaque point porte l'indice du point de trajectoire qui le dépose, et un
      point testé ne voit que la matière antérieure (même invariant que
      `_grow_printed_index_to`).
    - **Requête** : les K plus proches voisins, plus une requête exacte par
      boule quand la liste est saturée de matière future (la tête survole
      souvent ce qui sera imprimé plus tard).
    - **Écart assumé** : la distance est surestimée d'au plus un demi-pas
      (0,125 mm), à couvrir par la marge de 1 mm.
    - Désactivé par défaut tant qu'il n'est pas mesuré sur les quatre modèles
      de `res/test_models/`.
    - **Mesures contre l'exact (2026-10-01)** : l-triple, trajectoire entière,
      680 événements communs sur 681 (le manquant à 0,000 mm), 1 s contre 552 s ;
      l_shape, 0 événement des deux côtés ; turkish, par ÉCHANTILLON
      (`scripts/compare_step9_sampled.py` : 1 000 points signalés + 1 000 non
      signalés revérifiés en exact), **0 faux positif, 0 faux négatif**, écart
      de pénétration ≤ 0,087 mm (borne 0,125), 41 s pour la trajectoire entière
      contre 6 684 s pour les seuls 2 000 points exacts. L'exact sur la
      trajectoire entière de turkish ou death prend plusieurs heures.
      death, par échantillon (300 + 300) : 0 faux positif, 0 faux négatif,
      écart ≤ 0,073 mm, 16 s. **Activé par défaut dans `main.PipelineParams`**
      (le défaut de la dataclass `CollisionSimParams` reste `False`).
    - **Pas encore portés** : la borne de B selon la hauteur et la
      réorientation par paliers de la recherche. Leur place est l'étape 10,
      puisqu'elles modifient la trajectoire.

63. **Étape 8 : filtre médian des axes d'orientation libres
    (`Step8Params.free_axis_median_window`, opt-in, 2026-09-30).** La normale
    de couche est constante par tétraèdre. Près de l'axe du plateau, C étant
    imposé, B oscille d'un point à l'autre entre 0° (butée) et environ 14°.
    Mesuré sur la trajectoire Squirtle de recherche :
    - 5,9 % des segments du corps font sauter B de plus de 10° ;
    - le firmware, qui coordonne les axes, ralentit tout le mouvement à
      chaque inversion : le temps du corps passe de 1,27 h à 4,67 h ;
    - une hystérésis n'y change presque rien (−8 % à 0,5°) ;
    - **un filtre médian de 9 points par polyligne gagne 41 %**.

    Deux règles :
    - **Appliquer AVANT tout contrôle de collision.** Appliqué après, il défait
      les réorientations anti-collision (écart mesuré jusqu'à 115°).
    - L'écart résiduel à la direction visée est recalculé, donc les
      `Step8Diagnostic` restent justes.

    **Limite connue** : une alternance PARFAITE (0/14/0/14) résiste à un
    médian. Les données réelles montrent surtout des pics isolés.

    **Mesuré sur le pipeline FEM par défaut (2026-10-01)** : le médian 9 ne
    change presque rien (l_shape 708 → 671 sauts de B > 10°, temps réel
    0,157 → 0,154 h ; death, turkish, l-triple identiques). Raison : 597 des
    708 sauts de l_shape (84 %) sont aux TRAJETS À VIDE (changement de couche),
    pas des oscillations ; les 111 restants sont près de l'axe du plateau
    (rayon médian 1,4 mm), où C imposé tourne vite. C'est le cas « ne pas
    centrer les modèles » de la feuille de route, pas un défaut de filtre.
    Reste opt-in.

64. **Nouvelle référence step9b, mécanismes opt-in tranchés, et les
    inclinaisons ATTEIGNABLES activées par défaut (session autonome du
    2026-10-01).** Outil : `scripts/measure_matrix.py <stl> <variante>`, une
    ligne JSON par mesure (`out/measure_matrix*.jsonl`) : step9b, épaisseur de
    couche RÉELLE mesurée, diagnostics step8, sauts de B, **écart tête /
    normale locale** (cinématique directe des angles résolus contre +∇φ du
    tétraèdre le plus proche : métrique PHYSIQUE, indépendante de la cible que
    step8 s'est donnée) et **temps réel estimé** (`feed_plan`, limites de B).

    **Référence step9b @ 0,8, après les limites #57/#58** (l'ancienne,
    0,466578 %, est obsolète) : l_shape 3,337 %, l-triple 4,821 %, death
    0,346 %, turkish 1,033 %.

    **Découverte : le défaut d'orientation du pipeline par défaut.** Avec la
    cinématique couplée (limite #59), la tête déposait en moyenne loin de la
    normale de ses propres couches : écart médian 28,7° (l_shape), 49,1°
    (l-triple), 16,4° (death), 28,2° (turkish), jusqu'à 55 % des points au-delà
    de 45°. Deux causes, mesurées séparément :
    - step4 visait des directions qui penchent en TANGENTIEL ou vers le
      CENTRE, inatteignables : B ne penche que dans le plan radial, vers +R ;
    - step8 visait UNE normale par couche (plan moyen), qui n'est dans le plan
      radial que de quelques points d'une couche qui fait le tour du plateau.

    **Mécanismes et verdicts (4 modèles)** :

    | variante | matière sans appui | verdict |
    |---|---|---|
    | `free_axis_median_window=9` seul | inchangée (orientation seule) | opt-in |
    | `support_feedback_iterations=3` à épaisseur égalisée | death −53 %, l_shape −3 %, l-triple −2 %, **turkish +6 %** | opt-in |
    | `support_feedback_iterations=50` | dégrade 3 modèles | opt-in |
    | `branch_cut_threshold=0.05` | **death ×1,94** | opt-in (retrait possible) |
    | `use_geodesic_directions` + lissage 50 | **death ×2,3** | opt-in (confirme #52) |
    | **`radial_tilt_origin`** (projection sur le demi-plan radial atteignable) | l_shape −9 %, l-triple −45 %, death −23 %, turkish −30 % | **ACTIVÉ** |

    `radial` améliore la matière sans appui sur les QUATRE modèles, alors même
    que ses couches sortent plus fines (biais défavorable : turkish 0,29 →
    0,24 mm, plus proche de la cible 0,2), et ramène l'écart médian tête /
    normale à 3-6° (au-delà de 45° : 0 à 1,7 %). Coût : temps réel de turkish
    +15 % (plus de couches) ; les trois autres modèles sont neutres ou plus
    rapides. Activation : `PipelineParams.reachable_tilt_from_machine = True`
    (défaut) ; `main._effective_field_solver_params` renseigne
    `radial_tilt_origin` depuis le profil machine quand `position_coupled_axes
    == [0]`. Sans effet sur une machine à C libre.

    **Collisions avec la vraie tête (step9 rapide, part des points
    échantillonnés)**, mesurées APRÈS la décision ci-dessus (critère oublié
    dans la première matrice) :

    | modèle | sans `radial` | avec `radial` |
    |---|---:|---:|
    | l_shape | 0 % | 0 % |
    | l-triple | 10,5 % | **0 %** |
    | death | 31,0 % | **41,4 % (dégrade)** |
    | turkish | 61,3 % (+ 36 contacts plateau) | 53,8 % (0 plateau) |

    **Normale par point (step8, opt-in, ARBITRAGE À TRANCHER PAR
    L'UTILISATEUR)** : `Step8Params.per_point_normal` (gradient local du
    champ, `gradient_direction_field`) et `per_point_normal_nodal` (variante
    continue : gradient moyenné aux nœuds, interpolé en barycentrique).
    Écart médian supplémentaire divisé par 2 (death 6,2 → 3,5°, turkish 3,1 →
    1,4°), mais temps réel +38 à +130 % : le gradient P1 varie d'un
    tétraèdre à l'autre et fait osciller B, que le firmware ralentit. Meilleur
    compromis mesuré : nodal + médian 9 (death 1,76 → 2,43 h, turkish 4,75 →
    5,59 h). Qualité de dépôt contre temps : aucune mesure ne tranche.

65. **Étape 6 : des trous dans TOUTES les parois, depuis l'activation du
    découpage en tronçons (corrigé le 2026-10-01, `Step6Params.continuous_rings`,
    ACTIVÉ par défaut).** Trouvé en portant step6 en C++ (pas 15), parce que la
    longueur totale dépendait du point de départ des anneaux, ce qui ne devrait
    pas arriver.
    - **Cause** : `_split_ring_into_runs` découpe chaque anneau en tronçons
      monotones à indices DISJOINTS. Le segment qui relie deux tronçons
      n'appartient à aucun : step8 en fait un trajet à vide. De plus, le
      contour d'origine n=0 et ses trous sortent OUVERTS de la
      classification, donc leur segment de fermeture n'est jamais tracé.
    - **Mesure** (longueur des polylignes, découpées contre non découpées) :
      **death 7,7 %**, **l_shape 15,7 %** de la longueur perdue au découpage ;
      fermetures : 0,1 % et 0,85 %.
    - **Correctif** : chaque tronçon reprend le dernier sommet du précédent
      (le segment de jonction prend le sens du tronçon suivant, qui reste
      monotone), et les contours n=0 reçoivent leur point de fermeture.
      Vérifié : longueur découpée = longueur non découpée + fermetures, à
      0,000 mm près. Le test historique « aucun point ajouté au découpage »
      reste vrai sans l'option ; avec elle, un sommet par jonction est
      volontairement répété.
    - **step9b, 4 modèles** (`out/measure_continuous.jsonl`) :

      | modèle | longueur extrudée | matière sans appui |
      |---|---:|---:|
      | l_shape | 6 946 → 8 324 mm (+20 %) | 3,05 → 3,14 % |
      | l-triple | 10 819 → 12 433 mm (+15 %) | 2,64 → 2,66 % |
      | death | 122 920 → 133 441 mm (+8,6 %) | **0,266 → 0,227 %** |
      | turkish | 325 746 → 352 788 mm (+8,3 %) | **0,721 → 0,702 %** |

      Activé par défaut malgré deux légères hausses : ce n'est pas un réglage
      mais un défaut physique (parois trouées), et la métrique varie peu.
    - **Contrôle physique** (`scripts/check_wall_coverage.py [--legacy]`) :
      part de la bande de parois couverte par les cordons réellement extrudés.
      l_shape : médiane **84,0 → 98,8 %**, aire non couverte 247 → 19 mm² ;
      death : **89,1 → 98,5 %**, 2 378 → 347 mm². Le reliquat (~1,5 %) tient aux
      extrémités plates des cordons et aux zones plus fines que deux parois.
      Le même script mesure le remplissage (zone intérieure hors recouvrement des
      parois) : death médiane 22,2 % pour 20 % visés (raccords du zigzag), 2
      couches sur 174 sous la moitié de la visée ; l_shape 30,5 % (zones
      étroites). Pas de défaut de remplissage.
    - Tests : `test_split_ring_shared_vertices_extrude_every_segment`
      (contrôlé : il échoue sans partage), `test_continuous_rings_closes_original_contour`.

66. **Étape 10 : en mode exact, aucune réorientation ne pouvait réussir
    (corrigé le 2026-10-01).** `_collisions_at_points` (revérification) et
    `_build_reorientation_context` (évaluation des candidats) échantillonnaient
    la tête SANS `_filter_nozzle_contact_zone`, contrairement à step9. La
    pointe touchait donc toujours le cordon qu'elle dépose, et tout point testé
    était déclaré en collision. Conséquence : en mode exact (le mode par
    défaut de `CollisionSimParams`), chaque collision finissait en support de
    secours. Le mode rapide (`_fast_checker`) appliquait déjà le filtre.
    Trouvé en comparant step9 rapide et exact par échantillon
    (`scripts/compare_step9_sampled.py`) : l'exact via step10 signalait 300
    points sur 300 là où step9 n'en voyait aucun. Test :
    `test_exact_recheck_matches_step9_exact_events` (échoue sans le correctif).
    Les journaux de step10 antérieurs (mode exact) sont à relire : leurs
    « support de secours » n'en étaient pas forcément.

67. **Étape 5b : ORDRE DES COUCHES CONTRAINT PAR LA TÊTE
    (`pipeline/step5b_collision_order.py`, 2026-10-01, ACTIVÉ par défaut).
    Collisions de la vraie tête divisées par 35 (death) et 90 (turkish) : le
    pipeline FEM devient imprimable sur ces deux modèles, ce qu'aucun des
    mécanismes des limites #50 à #54 n'avait obtenu.**
    - **Le circuit ouvert** (même diagnostic que la limite #53) : step5 ordonne
      les couches par isovaleur et continuité spatiale, sans jamais connaître la
      tête ; la collision n'est mesurée qu'à l'étape 9, après tout. Sur death,
      la tête touchait (limite #64) de la matière déposée ~55 000 points plus tôt
      et ~9 mm plus haut.
    - **Mécanisme** : la tête est posée (IK de step8, mêmes conventions) sur des
      points du contour de chaque couche L, tous les 3 mm ; si elle passe à moins
      de 1 mm de la matière d'une couche M (contour densifié + grille intérieure,
      un seul arbre k-d pour toutes les couches), alors **L avant M**. Ces
      contraintes s'ajoutent à l'ordre de step5 par un tri topologique STABLE
      (sans contrainte, ordre inchangé). Cycle : on place la couche au plus petit
      nombre de contraintes entrantes restantes (heuristique d'arcs de retour) ;
      le premier essai (plus petit indice) laissait 13 105 contraintes violées
      sur turkish et ne gagnait presque rien (53,6 → 47,9 %).
    - **Mesures** (`measure_matrix.py <stl> corder --collisions`,
      `out/measure_corder.jsonl`) :

      | modèle | contraintes | couches déplacées | violées | collisions tête | matière sans appui |
      |---|---:|---:|---:|---:|---:|
      | death | 25 523 | 127 | 0 | **41,6 → 1,2 %** | 0,227 → 0,311 % |
      | turkish | 52 416 | 345 | 216 | **53,6 → 0,6 %** | **0,70 → 0,25 %** |
      | l-triple | 5 224 | 0 | 0 | 0 → 0 | identique |
      | l_shape | 5 725 | 0 | 0 | 0 → 0 | identique |

      Pipeline complet death (step1→10) : 57 → **35 s**, points en collision
      27 855 → 827. Coût de l'étape : +3,5 s (death), +12 s (turkish).

      **Avec les contraintes d'APPUI en plus** (`support_depth=0.3`, défaut :
      M sous un point de L ⇒ M avant L, sauf si la tête impose l'inverse) :

      | modèle | contraintes | violées | collisions tête | matière sans appui |
      |---|---:|---:|---:|---:|
      | death | 26 436 | 161 | **41,6 → 0,13 %** | **0,227 → 0,222 %** |
      | turkish | 53 276 | 277 | **53,6 → 0,21 %** | **0,70 → 0,198 %** |
      | l-triple | 5 583 | 4 | 0 → 0 | identique (aucune couche déplacée) |
      | l_shape | 6 068 | 2 | 0 → 0 | identique (aucune couche déplacée) |

      **Plus aucun arbitrage** : aucun des quatre modèles ne se dégrade sur
      aucune des deux métriques. C'est le défaut retenu.
      Pipeline complet (step1→10, défauts) : **death 35 s, 88 points en
      collision** (27 855 avant) ; **turkish 75 s, 322 points** (~76 000 avant).
    - Sans les contraintes d'appui, death payait +37 % de matière sans appui
      (0,311 %) : l'arbitrage a disparu en les ajoutant, sans le trancher.
    - **Leçon, la même qu'en #53** : sept mécanismes jouaient sur les ENTRÉES de
      l'ordre (direction, ancrage, géodésique) ; celui-ci branche le critère qui
      compte (la tête passe-t-elle ?) à l'endroit où l'ordre se décide.
    - **Vérifié image par image** (`out/death_order_step5b.gif`) : pilier, pont,
      cabine, puis la coque en porte-à-faux EN DERNIER. C'est l'objectif énoncé
      par l'utilisateur en limite #54 (« imprimer les couches suspendues en
      dernier »), atteint ici pour +37 % de matière sans appui, là où l'ordre
      géodésique le payait ×4. Restent quelques petits fragments isolés en l'air
      à mi-parcours (pointe de la proue), à examiner.
    - Tests : `tests/test_step5b_collision_order.py` (tri stable, cycle signalé).
    - **Reste** : 216 contraintes violées sur turkish (cycles réels entre
      couches en anneau ; découper les couches en morceaux les lèverait) ; la
      granularité est la couche entière.

69. **Normale par point ACTIVÉE par défaut (2026-10-06, décision de
    l'utilisateur), réglable : `Step8Params.per_point_normal` (défaut `True`).**
    La cible de l'IK est le gradient LOCAL du champ en chaque point, et non la
    normale du plan moyen de la couche. Mesuré sur les 4 modèles de
    `res/test_models/` (`scripts/measure_matrix.py`, variantes `ref` / `ppn`,
    résultats `out/measure_ppn_2026-10-06.jsonl`) :

    | | l_shape | l-triple | death | turkish |
    |---|---:|---:|---:|---:|
    | écart médian tête/normale | 5,8 → 2,8° | 4,1 → 1,9° | 6,1 → 2,6° | 4,7 → 1,7° |
    | points à plus de 15° | 24,1 → 3,9 % | 10,3 → 0,6 % | 12,0 → 2,6 % | 14,7 → 4,5 % |
    | temps d'impression | 0,16 → 0,26 h | 0,18 → 0,30 h | 1,93 → 4,40 h | 5,01 → 8,40 h |
    | `step8_diagnostics` | 15 300 → 7 000 | 8 560 → 5 330 | 0 → 111 400 | 24 400 → 195 100 |
    | sauts de B > 10° | 708 → 2 291 | 7 → 1 430 | 0 → 19 730 | 1 → 44 120 |

    Matière sans appui (step9b) strictement inchangée : la géométrie des
    couches ne change pas, seule l'orientation de la tête. **Coût** : temps
    d'impression +65 à +130 %, car le firmware coordonne les axes et ralentit
    à chaque saut de B ; sur death et turkish les orientations inatteignables
    augmentent beaucoup. **Non mesuré** : `ppn` combiné au filtre médian
    (limite #63, `free_axis_median_window`, variante `ppn_radial_median9`),
    remède identifié pour les sauts de B. Suite : 466/466 dans les deux
    positions de `APEX_USE_CPP`.

    **Mesure du filtre médian (2026-10-06)**, `ppn` contre `ppn_median9`
    (`out/measure_ppn_b360_2026-10-06.jsonl`) : temps d'impression −19 / −20 /
    −31 / −21 % (l_shape / l-triple / death / turkish), sauts de B > 10° en
    baisse partout (death 19 730 → 10 720, turkish 44 120 → 26 870) ; coût :
    écart médian +0,2 à +0,4°, points à > 15° +0,2 à +1,2 point,
    `step8_diagnostics` +8 à +16 %. Matière sans appui inchangée. **Non
    mesuré : les collisions de la tête avec les angles filtrés** (le filtre doit
    s'appliquer AVANT le contrôle de collision, limite #63) : le filtre reste
    donc opt-in (`free_axis_median_window = 0`).

    **Limites de B réconciliées (2026-10-06)** : `main._default_machine_profile`
    passe de 750 à 360 °/s (vitesse) et reste à 300 °/s² (accélération) : minimum
    des limites de l'utilisateur et des fichiers RRF (`4axisinput.g` : 360 °/s,
    2000 °/s² ; `config.g` : 1667 °/s, 1300 °/s², unités supposées en degrés).
    Effet mesuré sur le temps d'impression estimé : nul (4,40 h et 8,40 h
    identiques à 750 °/s) : l'accélération est le facteur limitant.

    **SDK Vulkan 1.4.363.0 installé (2026-10-06)** via winget
    (`KhronosGroup.VulkanSDK`, hachage vérifié), `VULKAN_SDK=C:\VulkanSDK\1.4.363.0`,
    GPU vu : NVIDIA GeForce RTX 3090. Aucun code GPU écrit.
    **Filtre médian de B et collisions** (`scripts/measure_step9.py --median=N`,
    step9 rapide natif, normale par point active, limite #69) : événements de
    collision `printed_part`, sans filtre → avec filtre de 9 points : l_shape 0 → 0,
    l-triple 0 → 0, death 137 → 118 (−14 %), turkish 1 613 → 1 406 (−13 %) ;
    pénétration maximale inchangée (≈ 1,17-1,19 mm). Le filtre ne dégrade donc
    pas les collisions : il les réduit. Il reste opt-in
    (`free_axis_median_window = 0`) : à activer sur décision de l'utilisateur
    (−19 à −31 % de temps d'impression, +0,2 à +0,4° d'écart médian, +8 à +16 %
    d'orientations inatteignables).

    **Écart avec l'ancien poste** : 94 points de corps (0,011 %) du dino ont un B
    différent (jusqu'à 10,8°) de la trajectoire produite sur l'autre poste ; la
    chaîne Python actuelle donne la même différence : elle vient de l'évolution du
    code depuis, pas de `--apex-detect`.

    **ATTENTION : les valeurs de B ci-dessous (−129,5° à 0°, −114° à 0°, −64° à 0°) viennent de
    la version erronée `b_sign = −1`, corrigée le jour même (limite #61) : ces fichiers sont à
    REGÉNÉRER avec `b_sign = +1` avant toute impression.**

    **G-code régénérés (2026-10-06)** avec B inversé, pivot à 46 mm et dégagements
    contrôlés (`export_gcode.py --strict`) : Benchy 1 070 769 mouvements, 4,57 m,
    1,64 h nominal / 3,32 h avec les limites de B, B firmware −129,5° à 0°,
    pointe la plus basse 0,139 mm, 0 dégagement en contact ; dino 973 269
    mouvements, 6,40 m, 2,24 h / 4,28 h, B −114° à 0°, 0,176 mm, 0 en contact.
    **Y cumulé : Benchy ≈ 32 tours, dino ≈ 621 tours** (à valider côté câblage).
    Fichiers d'essai : `scripts/gen_b_direction_test.py` (essai à SEC du sens de B :
    la pointe reste immobile, sans extrusion) et `scripts/export_model_gcode.py`
    (STL → G-code par le pipeline FEM, avec comptage des collisions de la vraie
    tête) ; `out/test_l_shape.gcode` : 0,28 m, 0,26 h, B −64° à 0°, aucune
    collision. Archives dans `out/share/*_v4.zip` (Deflate, les fichiers font
    ≈ 65 Mo : un canal qui tronque à 1 Mo les corrompt). **Non fait : le sens
    physique de B n'est pas encore confirmé sur la machine.**

    **Fichier manquant** : `research/S4_Slicer/safe_io.py` (réessai des écritures
    sur verrou Windows) est absent de la branche alors que 7 scripts l'importent ;
    ces scripts ne démarrent pas sans lui. Les essais de cette session ont utilisé
    un équivalent minimal hors dépôt (dossier temporaire de session).

71. **Suite autonome du 2026-10-06 (hors UI) : sept points traités.**

    - **Cumul de rotation de Y** : `GcodeReport.y_turns` (étendue de Y déroulé, en tours) et
      `GcodeParams.max_y_turns` (`None` = aucune limite) ; avec `strict`, un dépassement
      refuse d'écrire le fichier. Affiché par `export_gcode.py` et `export_model_gcode.py`.
      Valeurs : Benchy ≈ 32 tours, l_shape ≈ 240, dino ≈ 621. Test :
      `test_y_turns_reported_and_limit_refuses_write` (chemins Python et natif).
    - **Échec de test isolé** (`test_export_extrusion_and_retraction`, `APEX_USE_CPP=0`, une
      seule fois, pendant un calcul en parallèle) : **non reproductible** (46 exports sous 8
      fils de charge CPU, 0 échec ; cinq passages complets propres depuis). Cause inconnue.
    - **Filtre médian de B ACTIVÉ par défaut** (`Step8Params.free_axis_median_window = 9`) :
      mesuré sur les 4 modèles (limite #69/#70 : temps −19 à −31 %, collisions en baisse,
      matière sans appui inchangée, écart médian +0,2 à +0,4°). Réversible : valeur 0.
    - **Non traités (UI)** : config UI (`config/params.py`, `config/schema.py`, pont JSON),
      sélecteur graphique d'angle, conseil « ne pas centrer les modèles », rendu de l'aperçu
      (rétrécissements à B > 90°).

72. **Points 6 a 9 du DT (2026-10-06, `config.g` confirme actif par l'utilisateur).**

    - **Limites de B** (`main._default_machine_profile`) : 1667 deg/s et 1300 deg/s2 (M203 B100000
      /min, M201 B1300 de `config.g`). Les 360 / 300 de la limite #69 etaient un minimum
      prudent entre deux fichiers dont un n'est pas actif.
    - **Z et Y ne sont PAS modelises par `feed_plan`** (seul B l'est). `config.g` donne Z =
      400 mm/min (6,7 mm/s) et 70 mm/s2 : le pivot de B couple Z a B (Z = z + 46 (cos B - 1)), donc
      chaque saut de B est un saut de Z. Mesure (vitesses seules, B 1667, Y 600, accelerations
      ignorees : borne basse) : Benchy 1,66 h nominal -> 9,7 h a Z = 6,7 mm/s, 4,1 h a 20, 2,5 h a 50,
      2,0 h a 100, 1,9 h a 133 mm/s ; test_l_shape 0,13 -> 1,35 / 0,50 / 0,27 / 0,21 / 0,20 h. Vitesse de Z
      demandee (Benchy, ponderee par le temps) : p50 4,8 / p90 80 / p99 553 mm/s. **Recommande a
      l'utilisateur : `M203 Z8000` et `M201 Z1500` dans `config.g` (valeurs de `4axisinput.g`) ; a faire
      par lui, rappel a repeter jusqu'a confirmation** (le 400 ressemble aux 400 pas/mm de M92).
      `GcodeParams.fw_axis_speeds` (X, Y, Z, B ; defaut = valeurs recommandees) alimente
      `GcodeReport.print_time_fw_s` (`firmware_time_estimate`) ; mettre (1666.7, 600.0, 6.67, 1666.7)
      pour l'etat actuel. Outil : `scripts/estimate_gcode_time.py`. Test :
      `test_firmware_time_estimate_slows_on_z_limit`.
    - **Couverture de la peau** (surface de l'original echantillonnee contre les cordons densifies a
      0,15 mm) : dino 0,69 % a plus de 0,3 mm et 0,11 % a plus de 0,6 mm (un amas de 3,5 mm pres de
      (-7, -16, 50)) ; Squirtle 3,3 % / 0,61 % (deux plaques de ~6 mm vers z = 31, a 0,7 mm,
      faces tournees vers le haut). **Le « 6,2 % » de la limite #68 ne se reproduit pas** sur les
      trajectoires actuelles. Les 4 plaques sont a 0,7-2,2 mm de points de coupe du lissage
      (`skin_trim`) : hypothese, non demontree.
    - **Decalage du modele sur le plateau** (`out/run_offset.sh`, `out/offset/measure_offset.jsonl`, defauts
      actuels, modele centre sur x = 0 / 60 / 100 mm) :

      | modele | orientations inatteignables | matiere sans appui | collisions |
      |---|---|---|---|
      | l_shape | 3344 -> ~950 | 2,2 -> 3,1 % | 0 |
      | l-triple | 3415 -> ~1950 | 2,5 -> 4,5-4,7 % | 0 |
      | death | 138k -> 128-131k | 0,54 -> 0,16-0,23 % | 1,5 -> 0,4 % |
      | turkish | 222k -> 202k (60 mm) / 226k (100 mm) | 1,03 -> 0,59 / 0,52 % | 28 % -> 28 / 26 % |

      **Pas de regle nette** : le decalage reduit toujours les orientations inatteignables, aide death et
      turkish, degrade l_shape et l-triple. Ne pas en faire une consigne a l'utilisateur sur cette base.
      Mesure geometrique grossiere (`scripts/measure_offset_reach.py`, normale de surface prise pour la
      direction de depot) : centrer donne le MOINS de flancs inatteignables a plus de 45 deg (47 % contre
      72 % decale, Benchy) ; proxy faible (un flanc vertical se depose verticalement).
    - **Fragilite d'ingestion (trouvee en decalant turkish)** : `turkish.stl` echoue a l'ingestion
      (auto-intersections, `pymeshfix` effondre) pour 3 decalages sur 4 en x (20, 40, 70 mm ; 50 et la
      position d'origine passent), a la precision flottante pres : ce n'est PAS la position mais
      l'arrondi float32 des sommets qui fait basculer un maillage deja marginal. Arrondir les sommets au
      micron repare 2 cas sur 3 (pas 40 mm) : pas un correctif fiable. Cause localisee : selon l'arrondi,
      `resolve_self_intersections` laisse 8 (origine) ou 20 (decale) auto-intersections residuelles, et
      `_repair_local_self_intersections` abandonne dans le second cas. **Correctif** :
      `_tetrahedralize_with_translation_retry` (step1) rejoue la cascade, APRES un echec seulement, sur le
      maillage translate EN MEMOIRE (centre de la boite en x, y = 0, puis +0,37/+0,21 mm) et ramene les
      noeuds a leur place. Les trois cas en echec passent (20, 40, 70 mm). Pas de repli quand des zones de
      jonction sont fournies (coordonnees absolues). Test :
      `test_translation_retry_recovers_and_moves_nodes_back`.
    - **Backend Vulkan de `MeshDistance` : pas de gain mesurable aujourd'hui.** CPU parallele : 2 M de points
      en 1,03 s sur turkish (10 180 faces ; serie 7,9 s). La couleur de step7 en demande moins de 1 M par
      modele : au mieux ~1 s a gagner. Non construit ; a rouvrir si une etape demande des centaines de
      millions de requetes.

73. **Pipeline complet de death : 124 s -> 37 s (2026-10-06), resultats identiques.** Mesure `main.run_pipeline`, defauts (resolution 1,5), `APEX_USE_CPP=1`. Le
    « 20 s » cite jusqu'ici dans le DT etait celui de step1->8 @ 0,8 d'une session anterieure, AVANT la
    normale par point et le filtre median : il ne valait plus pour le pipeline complet.
    - **step10 : 81 -> 6 s.** `point_index in diagnostic_indices` sur une LISTE de ~115 000 diagnostics
      d'etape 8, une fois par point a traiter : boucle quadratique (73 s de temps propre). Ensembles a la
      place ; journal identique (115 312 `fallback_support`, 17 `local_reorientation`).
    - **step8 : 17 -> 6,6 s.** `_median_filter_free_axes` recalculait l'ecart residuel point par point
      (Rodrigues scalaire, 718 000 appels). Par lots, memes operations dans le meme ordre ; 115 192
      diagnostics identiques. Test : `test_median_filter_residuals_match_pointwise_rodrigues`.
    - **Controleur de collision rapide construit une fois (step9 + step10 : 5 constructions -> 1).**
      `FastCollisionChecker` construisait AUSSI l'arbre scipy, inutile en mode natif (arbre desormais
      paresseux), et step10 le reconstruisait a chaque iteration alors que la matiere ne depend que des
      positions. `step9.fast_checker_for` le met en cache (une entree, comparee par valeur). Test :
      `test_fast_checker_is_reused_while_material_is_unchanged`. step10 6 -> 1,7 s, step9 7 -> 5 s.
    - Repartition finale (51 s -> 37 s apres le cache) : step1 12,3 s (TetGen ~12 s d'attente de
      sous-processus), step9 5,1, step8 4,7, step5b 3,9, step4 3,7, step6 3,6, step10 1,7, step5 1,4.
      Plus aucun poste Python isole au-dessus de ~2,5 s : le reste est eparpille.
    - **Constat, non tranche** : step10 range CHAQUE diagnostic d'etape 8 en « support de secours »
      (115 312 entrees sur death). Avec la normale par point, ces diagnostics sont des ecarts
      d'orientation, pas des collisions : le rapport annonce donc ~115 000 points « a supporter » qui
      n'en ont pas besoin. A decider avec l'utilisateur (seuil de diagnostic, ou ne pas basculer les
      diagnostics d'orientation en support).

74. **Session autonome du 2026-10-07 (sans essai machine) : supports de l'étape 10,
    calibrage de L, pièce du direct, échec de test expliqué.**

    - **`Step10Params.orientation_support_min_residual_deg`** (défaut : voir la décision
      ci-dessous). Un diagnostic d'étape 8 (écart tête / direction visée, SANS
      collision) dont l'écart est sous ce seuil est journalisé `orientation_deviation` et
      n'est plus compté comme support ; au-dessus, il suit la voie historique. Une collision
      réelle au même point suit toujours la voie normale (test
      `test_orientation_threshold_never_hides_a_real_collision`). Mesure (défauts, rejeu de
      step10 sur les trajectoires du pipeline complet) :

      | modèle | diagnostics | > 45° | > 60° | supports : `None` / 45 / 60 / `inf` |
      |---|---:|---:|---:|---|
      | l_shape | 8 107 | 0 | 0 | 8 107 / 0 / 0 / 0 |
      | l-triple | 5 853 | 0 | 0 | 5 853 / 0 / 0 / 0 |
      | death | 115 192 | 187 | 28 | 115 312 / 396 / 237 / 209 |
      | turkish | 184 932 | 510 | 100 | 185 658 / 1 735 / 1 332 / 1 244 |

      Colonne `inf` = « jamais ». Il ne reste alors que les vraies collisions non
      réorientables (death 209, turkish 1 244) ; aucun diagnostic ne dépasse 90° sauf 1 sur
      turkish. step10 : death 2,3 → 0,9 s. **Décision de l'utilisateur (2026-10-07) : jamais.
      Nouveau défaut `None` = jamais ; `0.0` rend le comportement historique** (l'infini n'est
      pas représentable dans le JSON de l'interface).
    - **L (pivot -> pointe) réglable** (`file_io/pivot_length.py`, `config/pivot_to_tip.json`) :
      dernière valeur enregistrée, sinon 46 mm. Le firmware ne fait aucune cinématique (#61) :
      L est appliqué dans le G-code, qui est à RÉGÉNÉRER après tout changement de L.
      `export_model_gcode.py`, `gen_b_direction_test.py` et `estimate_gcode_time.py` le lisent
      (ou `--pivot-to-tip=`). Sans valeur enregistrée, `out/test_b_direction.gcode` régénéré est
      identique. Toute méthode de mesure de L est CONFIDENTIELLE et reste hors du dépôt
      (décision de l'utilisateur du 2026-10-07). Tests : `tests/test_pivot_length.py`.
    - **Pièce du direct** (`out/salon/`, pipeline FEM, résolution 0,8, step9 rapide) :

      | pièce | collisions tête | Y cumulé | temps : Z recommandé / `config.g` actuel |
      |---|---|---:|---|
      | l-triple | aucune | 68,5 tours | 0,22-0,23 h / 0,80 h |
      | l_shape | aucune | 239 tours | 0,19-0,20 h / 1,15 h |
      | real_overhang | 1 701 plateau | — | refusé (`strict`) |

      Recommandée : **l-triple** (4 fois moins de tours de Y). `real_overhang.stl` a son
      minimum à z = −15 mm : le pipeline ne pose pas le modèle sur le plateau, et l'export
      refuse à raison. À traiter : poser automatiquement le modèle à z = 0 à l'ingestion.
    - **Échec isolé de `test_export_extrusion_and_retraction` (2026-10-06) : cause trouvée.**
      La session `slicer_ux` a édité `file_io/gcode_export.py` à 10:29:20 (marqueur
      `apex_sliced` inséré en tête de l'en-tête, déplacé à 10:29:44) ; la suite
      `APEX_USE_CPP=0` lancée à 10:29:29 a importé la version intermédiaire, d'où l'assertion
      d'en-tête en échec. Ce n'est pas un défaut du code. **Leçon : d'autres sessions
      (interface) modifient les fichiers du moteur ; ne pas lancer de mesure ou de suite de
      tests pendant qu'une autre session édite.**

76. **Pose du modèle sur le plateau à l'ingestion (`place_on_platform`, 2026-10-13,
    ACTIVÉ par défaut) et régénération des G-code.**

    - **Constat** : la voie FEM gardait les coordonnées du STL. `platform_z` était
      bien détecté (min z), mais le plateau de la machine et de step9 est à z = 0 :
      `res/test_models/death.stl` a son minimum à **z = 51,4 mm**, `turkish.stl` à
      **30,3 mm**, `real_overhang.stl` à −15 mm (refusé par l'export). Les deux
      premiers flottaient au-dessus du plateau : collision plateau jamais exercée,
      et un G-code qui aurait démarré en l'air.
    - **Mesure, 4 modèles** (`out/measure_placed_2026-10-13.jsonl`, résolution 0,8,
      variante `placed` de `scripts/measure_matrix.py`) : l_shape et l-triple, déjà
      à z = 0, **identiques au bit** ; death 0,22 → 0,21 % de matière sans appui,
      collisions de tête 118 → 116 ; turkish 0,17 → 0,19 % (défauts d'ordre 132 →
      223, 313 → 312 couches), collisions 1 406 → 1 431 ; **0 collision plateau**
      sur les deux. Mixte et de l'ordre du bruit (même géométrie, autre origine :
      les isovaleurs et TetGen s'arrondissent autrement). Activé quand même : c'est
      une normalisation de repère, pas un mécanisme, et sans elle l'impression
      démarre à 30-51 mm du plateau. Réversible : `place_on_platform=False`.
    - **À savoir côté interface** : un modèle posé plus haut dans le STL est
      redescendu ; l'UI doit lire `placement_offset_z` pour ne pas désaligner son
      aperçu.
    - **Temps d'impression réestimés** (`scripts/estimate_gcode_time.py`, borne
      basse, accélérations ignorées) :

      | pièce | nominal | Z recommandé (8000 mm/min) | `config.g` actuel (Z 400) |
      |---|---:|---:|---:|
      | Benchy | 1,64 h | 1,89 h | 9,76 h |
      | dino | 2,24 h | 2,62 h | 14,90 h |
      | l-triple | 0,18 h | 0,23 h | 0,80 h |
      | l_shape | 0,12 h | 0,20 h | 1,15 h |
    - **Suite complète : 495/495 dans les deux positions de `APEX_USE_CPP`**
      (243 s en Python, 226 s en natif). Test : `test_place_on_platform_translates_in_z_only`.

79. **Champ = DISTANCE AU PLATEAU par la méthode de la chaleur (`use_heat_geodesic_field`),
    couches = ses isosurfaces (2026-10-21/22). Premier champ dont l'ORDRE est correct sur une
    pièce suspendue sans compromis visible : pilier, cabine, puis la coque en dômes qui
    s'étendent depuis leur attache, rebords en dernier. Idée de l'utilisateur (« des couches
    perpendiculaires aux chemins de Dijkstra »). Opt-in, non activé.**

    - **Nouvelle métrique, indépendante de la gravité : la matière FLOTTANTE**
      (`D:/apex_multibranch/islands.py`, hors dépôt). Chaque tétraèdre reçoit l'instant
      d'impression de la couche d'isovaleur la plus proche ; un morceau de couche (composante
      par faces) est flottant s'il ne touche ni le plateau ni un tétraèdre imprimé avant lui.
      C'est ce que montrent les captures (îlots, bandes) ; step9b ne le mesure pas.
    - **Cause n°1 : le tri historique `(rang de région, isovaleur)` de step5.** Avec ce champ,
      φ EST la progression, et le rang de step3 la casse : sur death, les couches de rang 2
      (φ 23,6 à 39,6) passaient avant celles de rang 3 (φ 8,3 à 69,9). Corrigé par
      `LayerParams.isovalue_order` (opt-in) ; `frontier_order` (déjà écrit) donne exactement
      le même ordre ici.
    - **Cause n°2 : step5b.** Il déplace des couches ENTIÈRES (une isosurface de distance
      traverse toute la pièce) et pose la tête avec UNE normale par couche (`_layer_normal`),
      fausse sur un dôme : bandes tardives dans les montants de la cabine.
    - **Mesures death (résolution 0,8)** :

      | variante | flottant | îlots > 1 mm³ | collisions tête | step9b |
      |---|---:|---:|---:|---:|
      | FEM (référence) | 7,10 % | 72 | 0,5 % | 0,58 % |
      | chaleur, tri par rang | 0,79 % | 10 | 58,0 % | 2,15 % |
      | chaleur + step5b | 6,50 % | 84 | 47,9 % | 1,82 % |
      | **chaleur, `isovalue_order`** | **0,023 %** | **0** | 28,9 % | 2,19 % |

      real_overhang : 0 % flottant (FEM 1,5 %), collisions 17,3 % (FEM 0 %).
    - **step9b ne mesure pas la même chose ici** : sur death, 85 % de ses 7 319 « défauts
      d'ordre » sont là où la couche pousse vers le BAS, 47 % tête retournée (B > 90°), et
      pour 94 % le « support manquant » est à moins de 0,3 mm sous le point, imprimé ~3 points
      plus tard (la suite du même cordon). C'est le test VERTICAL appliqué à un dépôt vers le
      bas, pas un défaut d'ordre. Sa tenue physique (dépôt sous une couche, tête vers le haut)
      n'est PAS vérifiée : à essayer sur machine avant de compter dessus.
    - **Restent ouverts** : collisions (28,9 % sur death), orientation (25,7 % des points à
      plus de 45° de la normale : le flux de la chaleur n'est pas projeté sur le demi-plan
      atteignable, contrairement au FEM, #64), 5 minima locaux de φ (61 micro-îlots).
    - **Orientation atteignable (`FieldSolverParams.heat_project_reachable`, opt-in, 2026-10-22)** :
      le flux de la chaleur est projeté sur le demi-plan (+radial, vertical) avant
      `L*phi = div(X)` (même projection que le FEM, #64 ; `radial_tilt_origin` était déjà
      renseigné par `main` mais ignoré par la branche chaleur). Variante `heat_iso_reach`.

      | | real_overhang avant → après | death avant → après |
      |---|---:|---:|
      | points à plus de 45° de la normale | 14,7 → **0,5 %** | 25,7 → **16,8 %** |
      | points à plus de 15° | 26,3 → **2,8 %** | 56,7 → **21,3 %** |
      | orientations inatteignables (step8) | 15 118 → **3 305** | 315 440 → **99 691** |
      | collisions de la tête | 17,3 → 13,7 % | 28,9 → 26,6 % |
      | matière flottante | 0 → 0,66 % | 0,023 → 0,29 % |
      | step9b sans appui | 4,19 → 3,70 % | 2,19 → 2,88 % |

      L'orientation s'améliore beaucoup, mais les collisions ne bougent presque pas : leur cause
      n'est donc PAS l'orientation seule (restent l'ordre et la forme des couches, `step5b`).
      Contrepartie : un peu de matière flottante revient (flux projeté = plus un gradient).
      Couverture de la surface vérifiée inchangée (`D:/apex_multibranch/coverage.py`).
    - **Piège payé** : `heat_method_time_factor = 1` (temps h² de Crane) semblait résoudre
      real_overhang (0,67 %, 0 collision) ; FAUX succès : u sous-dépasse le seuil de gradient
      loin du plateau, φ y devient plat et la trajectoire perd le pont et la jambe (8 628 mm
      extrudés contre 12 331). Toujours vérifier la matière extrudée avant de conclure.
    - Tests : `test_isovalue_order_ignores_region_rank`, `test_isovalue_order_disabled_by_default`.
      Variantes `heat*` de `scripts/measure_matrix.py`, visualisation à curseur
      `scripts/plot_toolpath_slider.py`.

80. **Étape 6 : un anneau était coupé à CHAQUE oscillation numérique du champ
    (corrigé le 2026-10-24, `Step6Params.run_split_min_amplitude`, 0,05 mm dans
    `main`).** Trouvé en slicant un cylindre plein (Ø30 × 10) : 9 076 courses
    extrudées pour 50 couches, le périmètre extérieur en plus de 100 morceaux.
    `_monotonic_runs` coupait au moindre changement de sens de l'offset hors-plan
    (seuil 1e-12 mm) ; sur une couche plate l'anneau extérieur (257 points) a une
    étendue d'offset de 3,8 µm et 122 changements de sens. Remède : hystérésis
    d'amplitude (la coupure n'a lieu que si l'offset recule de plus que le seuil
    depuis l'extremum, placée juste après l'extremum) ; défaut des fonctions
    1e-12 = ancien comportement. Cylindre : 9 076 → 200 courses (4 par couche),
    60 924 mouvements au lieu de 85 892, 0 collision.
    **Mesure (variantes `legacy_runs` / `ref` de `measure_matrix.py`, Python pur,
    résolution 0,8)** : l_shape matière sans appui 3,042 → 3,038 %, trajet à vide
    −11 % ; l-triple 2,639 → 2,681 % (défauts d'ordre 78 → 82), trajet à vide
    −14 % ; matière extrudée identique au mm près ; 0 collision. death : voir
    `out/measure_matrix` (mesure longue en Python, résultat non consigné ici
    tant qu'elle n'a pas fini). turkish NON mesuré. Les 5 tests de step6 verts,
    2 ajoutés.

81. **`InfillPattern.OSCILLATING_SPIRAL` (2026-10-24, opt-in, jamais activé par
    défaut, NON mesuré sur les 4 modèles).** Idée de l'utilisateur : spirale
    d'Archimède oscillante, alternée. Contrainte posée par l'utilisateur : le
    plateau ne doit jamais reculer (mouvement continu, pics = mouvements radiaux,
    éventuellement « qui n'avancent pas »).
    - **Construction** (`_oscillating_spiral_polyline`) : en coordonnées polaires
      (θ, r) autour de l'AXE DU PLATEAU, r = pas·θ/2π + déviation ; points de
      contrôle dont θ ne décroît jamais ; deux points au même θ = pic purement
      radial (jamais deux de suite : la tête repasserait sur la même ligne, la
      polyligne cessait d'être simple) ; entre deux points, droite (pic net) ou
      cosinus (courbe), toujours monotone en θ.
    - **Réglages** : `spiral_alternation_threshold` (« Seuil d'Alternance », écart
      radial des pics, plafonné à (pas − buse)/2 : deux tours voisins restent à
      ≥ 1 buse), `spiral_peak_spacing`, `spiral_randomization_threshold`
      (« Seuil de Randomisation » 0-1 : amplitude, espacement, forme et pics
      radiaux tirés au hasard), `spiral_seed` (graine fixe + rang d'isovaleur :
      reproductible, tirages différents d'une couche à l'autre).
    - **Alternance** (corrigée le 2026-10-24 sur demande de l'utilisateur : la
      première version décalait les couches impaires d'un demi-pas, contresens sur
      son dessin) : TOUTES les couches suivent la MÊME ligne médiane et se
      superposent ; seules les aspérités diffèrent (rang impair : première
      déviation de signe opposé, tirages différents). Même sens de rotation du
      plateau à toutes les couches.
    - **Centre** : `spiral_axis_origin`, dérivé de `MachineProfile.platform_origin`
      par `main._effective_toolpath2d_params` (appelé par `run_pipeline` ET
      `run_slicing_only`, synchronisés d'emblée) ; l'axe est projeté dans le plan
      de la couche (`_plate_axis_in_local_2d`). Sans lui : centre de la zone.
    - **Densité** : le pas est corrigé d'après la longueur réellement déposée
      (amplitude et longueur d'onde figées avant, sinon la boucle ne converge pas) :
      disque Ø26, cible 0,20 → 0,201 à 0,207 pour une randomisation de 0 à 1.
    - **Raccords** : entre morceaux clippés, le long du bord, REFUSÉS s'ils feraient
      reculer θ (chaîne coupée, donc trajet à vide) : `_angle_nondecreasing`.
    - **Mesuré sur le cylindre** (randomisation 0,5) : 0 collision, 0 orientation
      inatteignable ; dans le G-code, 0 recul de l'axe Y du firmware sur 83 786
      mouvements de remplissage extrudés. Les parois, elles, tournent dans l'autre
      sens (121 courses sur 933).
    - **Limites assumées** : le rayon de départ vaut un demi-pas (trou central
      d'environ 0,8 mm de rayon, C tourne très vite près de l'axe) ; la monotonie
      est celle de l'angle LOCAL autour de l'axe, pas de l'angle 3D sur une couche
      très inclinée ; pièce hors axe : la spirale devient des arcs et les raccords
      peuvent imposer des trajets à vide ; chemin natif (`APEX_NATIVE_FILL`) non
      porté (repli Python). Aucun gain mécanique mesuré : motif disponible, non
      recommandé par défaut. Tests : 8 (`tests/test_step6_toolpath2d.py`).

82. **Portage natif, session autonome du 2026-10-26 (poste Linux, g++ 13, Python 3.13) :
    méthode de la chaleur, remplissage, étapes 5/5b/8/9/10 et en-tête G-code. `config.g` modifié
    par l'utilisateur (rappel arrêté) ; salon passé avec une récompense spéciale pour
    l'originalité du projet ; le Benchy n'a toujours pas été imprimé.**

    Demande : « porter en C++ ce qui reste en Python, en commençant par le plus facile ». Chaque
    pièce a son test d'équivalence (`tests/test_cpp_equivalence.py`), avec une garde anti-faux-vert,
    et un repli Python à grain fin (jamais tout ou rien).

    | pièce | où | équivalence | mesure |
    |---|---|---|---|
    | **méthode de la chaleur** (`heat_method_distance`) | `heat_method.cpp` | tolérance 7,7e-5 relatif sur death (tolérance du CG, pas un écart de calcul) | scipy 6,5 s → **2,6 s** sur death (4 cœurs) |
    | **CG de Jacobi parallèle et reproductible** (remplace le CG d'Eigen, mono-fil) | `linear_solver.cpp` | produits scalaires par blocs fixes de 4096 : indépendant du nombre de fils | sert aussi tout step4 itératif |
    | réorientation par paliers (`reorient_coupled_tiers`) | `reorientation.cpp` | mêmes angles et mêmes paliers que Python | — |
    | réorientation locale couplée (`_try_local_reorientation`) et borne de hauteur (`clamp_head_height`) | idem | mêmes entrées de journal et mêmes angles | `run` de step10 identique |
    | **contrôle rapide de collision par lots** (`fast_penetrations_batch`) | idem | mêmes événements, profondeurs à 1e-9 | death : step9 + re-vérification step10 **~12 s → ~3 s** |
    | tri topologique stable de step5b | `layer_order.cpp` | identique (entiers) | — |
    | couches par défaut de step5 en un appel (`extract_layer`) + volumes de tétraèdres | `isosurface.cpp` | contours **identiques au bit**, repère local à 1e-12 | — |
    | **remplissage de couches ENTIÈRES** (`fill_layers` : classification, ponts, parois, hors-plan, troncons monotones), en parallèle | `step6_layer.cpp` | sans pont : polylignes identiques à 1e-9 ; avec ponts : même longueur par couche (2e-3), découpage différent | step6 death 2,6 s → 1,4 s |
    | assemblage des points bruts, normale par point, filtre médian (step8) | `step8_assembly.cpp`, `step8_filters.cpp` | positions et indices au bit ; normales à 1e-15 (ddot de BLAS) ; médiane exacte | — |
    | en-tête et pied du G-code | `gcode_export.cpp` | texte identique ; `None` hors de [1e-4, 1e16) (le repr de Python diverge) → repli Python | — |

    **Remplissage natif ACTIF PAR DÉFAUT avec le noyau natif** (`APEX_NATIVE_FILL=0` rend shapely ;
    `APEX_NATIVE_LAYER=0` rend le chemin couche par couche). Mesuré sur les QUATRE modèles
    (`measure_matrix.py ref` contre `native_fill`, step9b) : l_shape 3,80 → 4,13 % (701 dépôts sans
    appui dans les deux cas, 7 % de points en moins), l-triple 5,10 → 5,28 %, death 3,31 → 3,43 %,
    turkish 0,162 → 0,152 %. Longueur extrudée : death −0,06 %, turkish −0,20 %, jamais plus de 2,7 %
    sur une couche. Mixte, dans le bruit des mesures de #64/#65 ; réversible. Les ponts passent
    désormais par Clipper2 (shapely auparavant).

    **Pipeline complet de death.stl @ 0,8, `APEX_USE_CPP=1`** (profil de `run_pipeline`, mêmes réglages) :
    42,4 s → 33,4 s ; dont **TetGen ~14 s, incompressible** (sous-processus, par choix de licence).
    Hors TetGen : 28 s → 18 s. Python pur : 288 s.

    **Ce qui reste en Python, honnêtement :**
    - **l'enchaînement des étapes** (`main.py`, dataclasses `Layer`/`ToolpathPoint`) : aucun
      orchestrateur natif. Les entrées/sorties inter-étapes sont des objets Python (360 000
      `ToolpathPoint` sur death : ~2 s de construction + conversions en tableaux). Un orchestrateur
      natif suppose de remplacer ce contrat par des tableaux, donc de toucher aux étapes 8, 8b, 9, 9b,
      10 et à l'export : chantier de la bascule finale, pas fait ici ;
    - **step1** : chargement/réparation par trimesh, pymeshfix, manifold ; TetGen en sous-processus
      (`pipeline/_tetrahedralize_worker.py`). Porter le chargement seul (STL, fusion de sommets,
      normales) gagnerait quelques centaines de ms et risquerait de changer l'ordre des faces ; non fait ;
    - **step8** : repli par grille du point que la résolution analytique ne traite pas (0,14 %),
      construction des `ToolpathPoint` ; **step9** : échantillonnage de la surface de la tête (trimesh),
      densification de la matière ; **step10** : boucle `run` et journal ; **step7** : enveloppe Python
      autour de `schedule_layers` ; `OSCILLATING_SPIRAL` (step6) ; `export_gcode` : construction du
      contrôle des dégagements ;
    - FEM directionnel, multibranche, parties planes (step4/4b/4c) : inchangés, voués au retrait (DT.md).

    **Tests** : les échecs hérités d'un poste Linux sont résolus par des tolérances (ULP de libm sur
    `fw_y_range`, 90 tétraèdres minimum au lieu de >100 selon la version de tetgen) ; le test de
    planification step7 utilise des couches shapely (avec Clipper2, le milieu d'une arête tombe à 7e-8
    de la grille des faces et deux triangles de couleurs différentes sont à égale distance : frontière
    déclarée de `apex/scheduling.hpp`). Reste en échec, hors périmètre :
    `test_branch_indicator_is_flat_on_a_straight_bar_and_rises_on_a_branch` (dépend de la version de
    tetgen, code FEM voué au retrait).

83. **Protocole d'appel du tranchage (R-003, contrats C2 et C3) et numéro de contrat (R-002),
    2026-10-26.** `python main.py steps` et `python main.py slice --request <requete.json>`
    (`slice_cli.py`, 19 tests dans `tests/test_slice_cli.py`) ; `main.py` garde son usage historique
    (`mesh_path` en argument) pour tout autre appel.
    - **Requête** : `contract`, `mesh`, `output_dir` obligatoires ; `output_name`, `params`
      (`PipelineParams` partiel, valeurs par défaut conservées), `gcode_params` (`GcodeParams` partiel,
      `strict=True` et L courant par défaut). Un champ inconnu à n'importe quel niveau ou un numéro de
      contrat différent de `config.schema.APEX_CONTRACT_VERSION` (1) : code 2.
    - **Sorties** : stderr = lignes JSON `step` (`fraction` 0 au début et 1 à la fin de chaque étape de
      `main.PIPELINE_STEPS`, plus `export`) et `log`, le journal libre de `run_pipeline` y étant
      redirigé ; stdout = UNE ligne, l'objet résultat. Codes : 0, 2, 3 (`MeshIngestError` dont R-016 et
      `ModelTooLargeError`, et refus de l'export `strict`), 1 plantage. G-code écrit en `.part` puis
      renommé. `run_pipeline(..., progress=)` émet ces événements ; `PipelineResult.ingest_report`
      porte `placement_offset_z`.
    - **Limite assumée** : aucune fraction INTERMÉDIAIRE dans une étape (aucune ne rend la main avant
      d'avoir fini) : la barre avance par étape, sans estimation de durée (décision de l'utilisateur,
      RULES §5.3). La commande reste `python main.py` tant que l'exécutable natif sans Python n'existe
      pas (R-014).
    - **Contrat** : `x-apex-contract` écrit dans `config/params.schema.json` par `scripts/gen_params.py` ;
      `CONTRACT_CHANGES.md` à la racine.
    - **Message de refus d'export** (`export_refusal_message`) construit avec les chiffres du rapport,
      accord en nombre vérifié par test. Constaté en passant : le Benchy cassé tranché avec TetraNet
      refuse l'export (4 mouvements dont la pointe passe jusqu'à 0,17 mm sous le plateau) ; cause non
      examinée.
    - **Test corrigé** : `test_native_grid_fallback_matches_python_grid_search` supposait que +Z ET -Z
      font échouer l'analytique. Mesuré : +Z se replie 43 fois sur 43, -Z 0 sur 43 (B = 180° est dans
      [0°, 185°]). Garde remplacée par l'égalité exacte « replis = points +Z ».
    - **Ports abandonnés après mesure** (profil complet de death @ 1,5, 15 s au total dont ~6 s
      d'attente de TetGen) : aucun poste Python propre au-dessus de ~0,7 s. L'échantillonnage de la
      tête de l'étape 9 (`__init__` 0,63 s), la boucle de l'étape 10 (0,2 s), l'enveloppe de l'étape 7 et
      la construction du contrôle des dégagements (0,5 s sur le Benchy) ne justifient aucun port pour la
      vitesse. Ils ne se justifient plus que par l'orchestrateur natif (bascule finale, R-014).

84. **Second algorithme « Element Aeolis » (nom choisi par l'utilisateur le 2026-10-29) : tranchage non planaire « surface en hauteur » (dossier `element-aeolis/`, alias d'import `nonplanar`, v0).** Le paquet et l'identifiant `nonplanar-heightfield/1` gardent leur nom technique ; « Element Aeolis » est le nom affiché (`meta.algorithm_name`, `generator` = `Apex engine (Element Aeolis)`).
    Indépendant du pipeline 5 axes : couches `z = S_k(x,y)`, buse verticale, sans tétraèdres.
    Modèles de test de l'utilisateur dans `res/test_models/non_planar/` (5 plaques 50 × 50 mm,
    toutes des surfaces en hauteur) + sortie IceSL de la plaque 5 pour comparaison. Détail,
    méthode et limites mesurées (épaisseur minimale non bornée, pente du dessus, écart avec
    IceSL) : `nonplanar/README.md`. Tests : `tests/test_nonplanar.py`. Non branché dans
    `main.py` ni dans le contrat (C1/C2) : à ouvrir en demande dans `apex/integration/REQUESTS.md`
    quand l'utilisateur décide de l'exposer.
    **Progression et vitesse (2026-10-30, `alpha-0.2.3.0.4.0.0.221`).** Étapes : `import`, `heightfield`, `layer_surfaces`, `toolpath2d`, `collision`, `export` (`slicer.AEOLIS_STEPS`, `main.py steps --algorithm element-aeolis`), chacune avec des fractions RÉELLES (`core/progress.py`) ; le tranchage tourne en trois phases pour cela. Plaques de test : 14-37 s → 3-9 s (tableau dans `element-aeolis/README.md`) : colonnes interpolées une fois (`LayerStack.sample`), rastérisation vectorisée (identique au bit), union de couverture pour l'empreinte (écart ~1e-7 mm² avec l'ancienne : dernier chiffre de `E` parfois différent), test d'intersection préparé dans le zigzag. Reste au profil : overlay shapely du remplissage (~2 s sur plate_6), écriture du texte (~1,5 s). **Appelable par `main.py slice --request` depuis `alpha-0.2.3.0.5.0.0.221`** (champ `algorithm: "element-aeolis"`, params = `AeolisSliceParams`, refus `not_a_heightfield` code 3 ; `slice_cli._slice_aeolis`, 7 tests dans `tests/test_slice_cli.py`) ; étapes pour l'interface : R-026.

85. **Famille de noms des algorithmes de tranchage (décision de l'utilisateur, 2026-10-29) : toujours préfixés
    « Element », un dossier `element-xxxx/` chacun.** **Element Ardor** : pipeline principal 5 axes par la méthode de la
    chaleur (le feu ; dossier `element-ardor/`, fait). **Element Aeolis** : non planaire, couches `z = S_k(x,y)` (fait, limite #84). **Element Aequor** :
    hybride planaire / non planaire, pour un maximum d'imprimantes compatibles (futur, pas commencé).
    **Element Arva** : planaire (futur, pas commencé). « Apex » reste le nom de l'application complète (Apex Slicer), qui contient les quatre Elements (air, eau, terre, feu ; Aetheris-4, la machine, est l'éther). Noms choisis, MARQUES NON
    VÉRIFIÉES : chercher les antériorités avant tout usage public. Les identifiants techniques (`nonplanar/`,
    `nonplanar-heightfield/1`) ne changent pas ; le nom affiché voyage dans `meta.algorithm_name` du bloc chiffré.

## Détail de `core/mesh_types.py` (dataclasses partagées)

Aucune logique de calcul dans ce fichier — uniquement des contrats de données.

- `Vec3` (alias `np.ndarray`), `AABB`
- `MeshData` (vertices, faces, normales, watertight, bvh, kdtree) — sortie étape 1
- `TetMesh` (tet_vertices, tetrahedra, tet_adjacency) — sortie étape 1
- `OverhangMap` (risk_score, is_overhang, is_bridge, bridge_span_direction) — sortie étape 2
- `Region`, `RegionDependencyGraph` (regions, dependencies,
  **`platform_contact_faces`** — booléen par face, ajouté pour l'étape 4, seule
  source de vérité sur le contact plateau) — sortie étape 3
- `ScalarField` (values par nœud, **`gradient` shape (T,3)** — un par tétraèdre,
  pas par nœud : champ P1 = gradient constant par élément) — sortie étape 4
- `ConstraintType` (enum), `ConstraintFalloff` (enum), `Constraint` — contraintes utilisateur
- `Layer` (contour_3d, region_ids, local_frame, **`toolpath_2d`** — liste de
  polylignes `(K,2)` dans le repère local 2D de la couche, pas 3D, remplies
  par l'étape 6 ; **`isovalue`** — l'isovaleur dont la couche EST
  l'isosurface, renseignée par l'étape 5, requise par l'étape 7 pour
  convertir `φ` en indice de couche ; **`segment_colors`** — une couleur par
  polyligne de `toolpath_2d`, `None` = matière interne/couleur libre,
  remplie par l'étape 7 ; color_id = couleur dominante, affichage UI)
  — étapes 5/6/7
- `ColorMap`, `PrintPass` (**`segments: list[(layer_id, index)]`** — pointeur
  vers la géométrie, jamais une copie ; `layer_range` seul ne suffit pas car
  une passe ne couvre qu'une partie de chaque couche), `PrintSchedule`
  (l'ordre de `passes` EST l'ordre d'exécution) — étape 7
- `RotationAxisSpec` (axe + plage angulaire, un DOF de rotation dans le
  repère de son parent) ; `MachineProfile` (**`rotation_axes: list[RotationAxisSpec]`**
  — cinématique décrite génériquement, N axes ordonnés et paramétrables,
  jamais figée sur une seule machine, voir limite #18)
- `ToolpathPoint` (**`axis_angles: list[float]`** — un angle par élément de
  `MachineProfile.rotation_axes`, même ordre ; remplace les anciens
  `radius`/`platform_angle`/`head_angle`, qui figeaient la structure sur
  exactement 2 DOF nommés "plateau"/"tête" — hypothèse invalidée par la
  généralisation à N axes composés, limite #18. La traduction vers une
  commande firmware précise (X/Y/Z/B) est délibérément hors de ce type,
  reportée à un futur `file_io/gcode_export.py`), `Step8Diagnostic`
  (**donnée pure** — position, layer_id, segment_index, erreur angulaire
  résiduelle en degrés — signalant une orientation désirée inatteignable
  dans les plages d'axes disponibles ; l'étape 8 ne lève JAMAIS d'exception
  pour ce cas, elle accumule ces diagnostics pour l'étape 10),
  `ToolpathMachine` — étape 8
- `CollisionEvent`, `CollisionReport` — étape 9
- `AdjustmentLogEntry` (point, méthode, détail) ; `ResegmentationRequest`
  (**contrat du hook de l'étape 10** — point, position, plage d'angles
  exclue par axe, raison ; l'étape 10 décide QUOI, l'appelant décide
  COMMENT re-jouer step3→step8, voir limite #23) ; `AdjustmentReport`
  (**`converged` explicite**, jamais un silence — plus `entries`,
  `unresolved_collisions`, `iterations_used`) — étape 10

## Détail de `core/geometry_utils.py`

- `build_face_adjacency(faces, face_count) -> csr_matrix` — graphe dual du mesh,
  utilisé par step2 (composantes connexes de ponts) et step3 (region growing).
- `fibonacci_sphere_directions(count) -> ndarray` — répartition uniforme sur sphère.
- `rotation_axis_circle_directions(rotation_axis_direction, count) -> ndarray` —
  répartition sur le cercle atteignable par un seul DOF de rotation. Reste un
  cas particulier valide et testé, PAS remplacé par la fonction générique
  ci-dessous (voir limite #18).
- `composed_rotation_directions(axes, reference_direction, samples_per_axis) -> ndarray` —
  généralise à N axes de rotation composés (`RotationAxisSpec`), chacun
  échantillonné sur sa propre plage angulaire (grille produit, rotations de
  Rodrigues successives). Dégénère exactement vers `rotation_axis_circle_directions`
  pour un seul axe à plage complète (vérifié par test), et vers la direction
  de référence seule pour une liste vide. Utilisé par step3
  (`_generate_candidate_directions`) et step7 (`_escape_directions`) — même
  dispatch générique dans les deux, aucun des deux ne connaît le nombre ni
  la nature des axes de la machine réellement utilisée.
- `point_in_polygon_2d(point, polygon) -> bool` — ray-casting standard, utilisé
  par step5 (test de cohérence) et step6 (classification contour/trou) pour
  l'imbrication géométrique — voir limite #11.
- `composed_rotation_angles_and_directions(axes, reference_direction, samples_per_axis) -> (ndarray, ndarray)` —
  même grille produit que `composed_rotation_directions`, mais retourne
  aussi les angles par axe (degrés) associés à chaque direction. Requis par
  step8 : une résolution cinématique inverse doit repartir d'angles
  concrets une fois le plus proche voisin trouvé, pas seulement de la
  direction résultante. `composed_rotation_directions` est désormais un
  wrapper mince autour de cette fonction (comportement inchangé, vérifié
  par la suite de tests existante).
- `refine_axis_angles_local_search(axes, reference_direction, initial_angles_deg, target_direction, iterations, initial_step_deg) -> (ndarray, float)` —
  affinement générique par coordinate-descent d'un jeu d'angles vers une
  direction cible, utilisé par step8 après la recherche de plus proche
  voisin. Seule primitive par-axe : `_rodrigues_rotate_single` (voir
  ci-dessous) — généralise à N axes sans jacobien analytique par-N. Le pas
  initial DOIT être dérivé de l'écart de la grille NN
  (`360/samples_per_axis`), pas une constante fixe, sous peine de rester
  piégé près du candidat NN (mesuré : résidu ~15-20° avec un pas fixe de 2°
  contre <0.1° avec un pas dérivé). Axes unitaires précalculés une seule
  fois par recherche (pas à chaque nudge) ; `_angular_error` en scalaire
  Python pur (`math.acos`), pas `np.clip`/`np.dot` — dispatch numpy
  générique mesuré comme dominant le coût de step8 sur des vecteurs 3D
  triviaux.
- `_rodrigues_rotate_single(vector, axis, angle_rad) -> ndarray` — même
  formule que `_rodrigues_rotate` (Rodrigues), mais pour UN SEUL vecteur
  `(3,)`, en scalaire pur (pas de `np.cross`/`broadcast_to`) — réservée aux
  appelants qui tournent un vecteur à la fois
  (`refine_axis_angles_local_search`, `step10_adjustment._direction_for_angles`).
  `_rodrigues_rotate` (chemin batch, `(N,3)`) reste strictement inchangée
  pour ses appelants batch légitimes (`composed_rotation_angles_and_directions`,
  step9). ~15× plus rapide sur ce chemin scalaire, correction du principal
  goulot de performance mesuré sur step8.

## Détail de `core/spatial.py`

Enveloppe `KDTree` autour du sous-ensemble de `scipy.spatial.cKDTree` que le
pipeline utilise vraiment — rien de plus (contrat verrouillé par lecture des
9 sites d'appel, voir limite #56, pas 5) :

- construction depuis un `(N,3)` ;
- `.query(points)` — plus proche voisin unique (`k=1`, la seule variante
  appelée) ; entrée `(3,)` → `(distance, index)` scalaires comme scipy ;
- `.query_ball_point(point, radius)` — indices dans un rayon, un seul
  appelant (`core/constraints.py`), re-filtré derrière donc ordre libre ;
- l'attribut `.data`.

Dispatch scipy ↔ `apex_core.KDTree` selon `core._cpp.enabled()`, **figé à la
construction** (même discipline que `_PrintedGeometryIndex._native`). Portage
architectural, pas une optimisation : `cKDTree` est déjà natif, gain de
vitesse nul attendu — le but est de retirer la dernière primitive spatiale
scipy du cœur. Construite par step1 (`mesh_data.kdtree`, partagée avec
step4/5/6), step3, step4, step5, step7.

## Détail de `pipeline/step4_field_solver.py`

Voie unique (décision de l'utilisateur, 2026-10-27) : le champ EST la distance au plateau à
travers la matière, par la méthode de la chaleur (Crane, Weischedel, Wardetzky). `run(tet_mesh,
mesh_data, overhang_map, region_dependency_graph, params) -> ScalarField`.
- Matrices FEM P1 (rigidité, masse concentrée) : `_assemble_stiffness_matrix`, `_FemPattern`
  (motif CSR partagé), `_assemble_lumped_mass_matrix`.
- `_heat_method_distance` : diffusion `(M + t L) u = source` avec `t = heat_method_time_factor ·
  h²` (h = arête médiane), flux normalisé `X = −∇u/‖∇u‖`, projeté sur le demi-plan atteignable
  (`_project_to_reachable_radial_tilt`) quand `radial_tilt_origin` est renseigné, puis
  `L φ = div X` et recalage à 0 sur le plateau. Natif : `apex_core.heat_method_distance`.
- Sources = nœuds du plateau (`_identify_platform_nodes`, depuis
  `RegionDependencyGraph.platform_contact_faces` de step3, seule étape qui connaît le plateau).
- `_solve_linear_system` : direct (`spsolve`) sous `direct_solve_node_threshold`, sinon CG de
  Jacobi (ILU en option) ; repli unique à 10× plus d'itérations pour la diffusion.
- Retirés le 2026-10-27 : le champ FEM directionnel (terme de surplomb, mélange de directions,
  rétroaction, directions géodésiques, ancrage), le champ multibranche, les parties planes et les
  contraintes utilisateur (`Constraint[]`, points d'attache — à refaire sur la voie chaleur).

## Détail de `pipeline/step1_ingest.py`

- `_tetrahedralize` : `resolution` fixe le plafond global (`-a` de tetgen).
  `near_platform_resolution` (opt-in, défaut `None`) active un raffinement
  LOCAL près du plateau via `TetGen.add_region` (deux régions : une tranche
  fine près du plateau, le reste à `resolution`) — voir limites #35/#36
  pour le mécanisme retenu et celui explicitement écarté (`bgmesh`/`metric`,
  cassé dans cette version de tetgen).
- `_split_for_platform_refinement` : découpe le mesh en deux coques
  ouvertes (`trimesh.slice_plane(..., cap=False)`) de part et d'autre de
  `platform_z + near_platform_margin`, calcule un point-seed par coque,
  garanti à l'intérieur via `_find_guaranteed_interior_point` (limite #42 —
  un seed hors-coque casse silencieusement `TetGen.add_region`). Utilisée
  uniquement quand `near_platform_resolution is not None`. Le contour de
  coupe est triangulé UNE SEULE FOIS (`_shared_cut_cap`) et n'apparaît
  qu'une fois dans le PLC combiné final, comme facette interne partagée —
  PAS deux fois via deux appels indépendants à `slice_plane(cap=True)`
  (limite #44 : deux triangulations indépendantes du même contour peuvent
  diverger sur un polygone non convexe riche, ce que TetGen refuse comme
  auto-intersection).
- `_shared_cut_cap` : triangule le contour de bord créé par une coupe
  `cap=False` (dédup des sommets, extraction des arêtes de bord sur le
  plan, `triangulate_polygon`) — reproduit la logique interne du paramètre
  `cap=True` de `trimesh.intersections.slice_mesh_plane`, appliquée une
  seule fois plutôt que deux (limite #44).
- `_find_guaranteed_interior_point` : rejection sampling déterministe
  (`trimesh.sample.volume_mesh`, seed fixe) en repli quand le centroïde
  d'une coque tombe hors d'elle-même (forme non convexe) — partagée entre
  `_split_for_platform_refinement` et `_find_internal_cavities` (limite #42).

## Détail de `pipeline/step6_toolpath2d.py`

- `_classify_contours` : classification trou/matière par parité du nombre de
  boucles englobantes (imbrication `point_in_polygon_2d`), généralisée
  au-delà de 2 boucles (ex: anneau dans anneau) ; chaque trou est rattaché à
  son parent immédiat (le plus petit contour matière qui le contient), pas
  au plus grand contour de la couche.
- `_offset_walls` : offset `shapely.buffer` progressif, arrêt propre sur zone
  dégénérée (limite #13) ; retourne aussi le polygone du dernier offset
  valide (frontière de départ de l'infill).
- `_bridge_direction_for_loop` : pont surface/volume par kdtree
  (`mesh_data.kdtree`, même technique que step4/step5), seuil de PRÉSENCE
  pas de majorité (limite #14).
- `_interpolate_out_of_plane` : pondération inverse-distance (`1/d**power`)
  des offsets hors-plan réels du contour, pour donner une position 3D
  fidèle aux points générés par `shapely` (limite #31, suite). Sources =
  contour extérieur **et** trous du même groupe. `LinearNDInterpolator`
  écarté (enveloppe convexe seulement → `NaN` sur contour non convexe).
- `_reorient_bridge_toolpath` : parois normales + `bridge_anchor_extra_walls`
  parois supplémentaires resserrées sur le périmètre (ancrage = matière sur
  les bords, jamais densité d'infill augmentée — ne pas réutiliser
  `infill_density` pour ça, ça casserait le contrat UI où ce paramètre
  représente la densité réelle) + infill balayé perpendiculairement à la
  portée du pont.
- Limite connue non résolue (hors périmètre de l'étape 6, cohérent avec les
  autres limites documentées) : un pont qui s'étend sur plusieurs couches
  consécutives est reclassé indépendamment couche par couche, sans lissage
  temporel de la direction de traversée entre couches voisines.
- `_out_of_plane_offset` / `_to_toolpath_3col` : `Layer.toolpath_2d` porte
  `(u, v, offset)` depuis le fix du débordement (limite #31) — `offset`
  EXACT (mesuré, jamais interpolé) pour le contour n=0 (paroi extérieure,
  identique au contour d'origine), et INTERPOLÉ pour les parois n≥1 et
  l'infill via `_interpolate_out_of_plane` quand
  `Step6Params.out_of_plane_interpolation_power` est renseigné (défaut
  `None` = nul, comportement historique ; `main.py` l'active à `2.0`).
  `_project_to_local_2d` reste strictement 2D (u,v),
  utilisée par toute la géométrie interne (`_classify_contours`,
  `_offset_walls`, `_generate_infill`) qui ne doit jamais voir la 3e
  colonne. `_lift_to_3d` accepte (K,2) OU (K,3) — reconstruit exactement le
  point d'origine avec offset, retombe sur le plan moyen sans (limite #31).

## Détail de `pipeline/step7_scheduling.py`

- `run(layers, mesh_data, tet_mesh, color_map, region_dependency_graph,
  scalar_field, params) -> (list[Layer], PrintSchedule)`. `tet_mesh` est
  requis **en plus** de `scalar_field` : `ScalarField.values` est indexé par
  nœud du TetMesh, donc illisible en un point 3D sans la géométrie des nœuds.
- `_edge_colors_for_polyline` : couleur par **arête** (pas par point —
  l'arête est l'unité extrudée, et raisonner par arête évite les
  sous-polylignes dégénérées à 1 point au découpage). Utilise
  `ProximityQuery.on_surface`, jamais `mesh_data.kdtree` (limite #15).
- `_smooth_edge_colors` : filtre 1D en longueur d'arc (`merge_tolerance`).
  Ce paramètre vit sur `ColorMap` (contrat UI) et n'est **pas** dupliqué dans
  `Step7Params` — deux sources de vérité qui divergent en silence est le
  piège déjà rencontré avec `platform_z` à l'étape 4.
- `compute_hole_deadline` / `_blocking_layer_along_ray` : seule la matière
  **future** scelle un trou (limite #16). Deux approximations assumées et
  documentées : `φ` lu au nœud TetMesh le plus proche (pas d'interpolation
  barycentrique), et test d'accessibilité géométrique — le test rigoureux
  avec le vrai modèle 3D de la tête et l'état incrémental de la pièce est
  explicitement le travail de l'étape 9, sur laquelle l'étape 10 boucle.
- Cas `color_map is None` / vide : une passe par couche, `color_id=None`,
  aucun trou. Le pipeline sans peinture est le cas d'usage par défaut et doit
  rester pleinement fonctionnel.
- Limite connue non résolue : l'ordonnancement est un glouton (maintien de la
  couleur active aussi longtemps que possible) + EDF sur les trous, pas un
  optimum global du nombre de purges.

## Détail de `pipeline/step8_machine_toolpath.py`

- `run(layers, schedule, machine_profile, params) -> (ToolpathMachine,
  list[Step8Diagnostic])`. `layers` doit être les couches **retournées par
  step7** (`colored_layers`), pas celles d'entrée de step7 : `PrintPass.segments`
  pointe vers `toolpath_2d` re-découpé par couleur.
- `_resolve_axis_angles` : pas de cinématique inverse analytique (impossible
  à écrire génériquement pour N axes arbitraires). Recherche du plus proche
  voisin sur `composed_rotation_angles_and_directions`, puis affinement par
  `refine_axis_angles_local_search` — le pas initial de l'affinement est
  dérivé de l'écart de la grille NN (`Step8Params.ik_refine_step_multiplier
  * 360/ik_samples_per_axis`), pas une constante fixe (voir détail
  `core/geometry_utils.py` ci-dessus).
- `_detect_singular_axes` / `_inherit_angle_at_singularity` : généralise
  « R≈0 ⇒ C indéterminé » (limite #19) à **n'importe quel axe** de la
  chaîne, pas seulement le premier — pour chaque axe, distance
  perpendiculaire de la position courante à sa droite de rotation propre
  (`_axis_chain_world_lines`, toutes ancrées à `MachineProfile.platform_origin`).
  Un axe signalé singulier hérite l'angle du point précédent, **avant** le
  lissage temporel, pour que le filtre de lissage ne voie jamais un saut
  artificiel à corriger.
- Aucune exception pour une orientation inatteignable : `Step8Diagnostic`
  (position, layer_id, segment_index, erreur résiduelle en degrés) est une
  donnée pure accumulée par `run`, consommée par l'étape 10 sans connaître
  la résolution IK interne. `ToolpathResolutionError` reste réservée aux
  entrées structurellement invalides (`toolpath_2d` manquant, segment hors
  limites) — même distinction que `SchedulingError` à l'étape 7.
- `_clamp_rate` : limiteur de taux causal (vitesse puis accélération),
  préféré à une spline car les limites machine sont des contraintes dures
  qu'un clamp respecte par construction. `_unwrap_angle_sequence`/
  `_rewrap_angle_sequence` gèrent le wraparound uniquement pour un axe à
  tour complet (`angle_range_deg` couvrant ≥360°) — une plage bornée
  strictement ne peut par construction pas boucler, le dépliage y serait
  une discontinuité inventée.
- `_iter_schedule_points` : détecte les frontières de trajet à vide par
  changement de clé `(layer_id, index)` entre segments consécutifs de
  `PrintPass.segments`, y compris entre deux passes non contiguës —
  `toolpath_2d` ne porte aucun marqueur de trajet à vide par construction
  (décision actée à l'étape 6), c'est entièrement la responsabilité de
  l'étape 8.
- `ToolpathPoint.axis_angles` est la seule sortie de step8 — délibérément
  agnostique de tout firmware. La traduction vers une commande machine
  précise (ordre X/Y/Z/B, noms d'axes) est repoussée à un futur
  `file_io/gcode_export.py`, pas codée ici.

## Détail de `pipeline/step9_collision_sim.py`

- `run(toolpath, machine_profile, params) -> CollisionReport`. Ne charge
  jamais `nozzle_model` lui-même (pas de dépendance au système de
  fichiers) — l'appelant doit le peupler au préalable, typiquement via
  `file_io.mesh_loader.load_toolhead_mesh`.
- `_positioned_toolhead_vertices` : même convention de composition N-axes
  que `composed_rotation_angles_and_directions`/`step8._axis_chain_world_lines`
  (rotations en ordre inverse de la liste, translation après) — appliquée
  ici à un tableau de sommets complet plutôt qu'à un seul vecteur
  direction, mais sans inventer de nouvelle convention.
- `_surface_sample_points` : sommets **+** points échantillonnés le long
  des arêtes (limite #21) — jamais les sommets seuls, sous peine de
  manquer une collision avec un obstacle fin traversant une face.
- `_PrintedGeometryIndex`/`_grow_printed_index_to` : capsules (segment +
  rayon), jamais une union booléenne (coût prohibitif sur plusieurs
  milliers de points). Grille spatiale maison (dict cellule -> indices),
  insertion O(1), pas de reconstruction périodique. Croissance strictement
  incrémentale, jamais en avance sur le point testé, et le cordon en cours
  de dépôt est explicitement exclu (limite #22).
- `_check_printed_part` (distance point-segment fermée, sans bibliothèque)
  vs `_check_static_obstacle` (`ProximityQuery`, jamais `mesh_data.kdtree`
  — limite #15) : deux méthodes de distance pour deux natures d'obstacle
  (primitif capsule vs mesh réel), un même seuil `safety_margin`. Broad-phase
  AABB commun aux deux avant tout calcul précis.
- `_platform_mesh_at_angle` : approxime le plateau par un disque (aucune
  géométrie réelle fournie à ce stade, contrairement à `nozzle_model`) ;
  retourne `None` si `platform_origin`/`rotation_axes` sont absents — le
  test est alors simplement sauté, pas de faux positif ni d'exception.
- Le mandrin/mors n'a délibérément **pas** de nouveau champ dans
  `MachineProfile` (aucune géométrie reçue, contrairement au modèle de
  tête) — `CollisionReport.chuck_checked` reste `False` en toutes
  circonstances, signal explicite plutôt qu'un champ spéculatif ou un skip
  silencieux.
- Aucune nouvelle dépendance (`python-fcl`/`CollisionManager` écartés) :
  AABB (numpy) + distance point-segment (formule fermée) + `ProximityQuery`
  (déjà utilisé par step7) couvrent le besoin.

---

## Détail de `pipeline/step10_adjustment.py`

- `run(toolpath, collision_report, machine_profile, sim_params, params,
  step8_diagnostics=None, resegmentation_hook=None) -> (ToolpathMachine, AdjustmentReport)`.
  Fonction pure : l'entrée n'est jamais mutée (vérifié par test).
- `_try_local_reorientation` réutilise `composed_rotation_angles_and_directions`
  (grille de recherche de step8), candidats triés par écart angulaire
  croissant et **coupés** au-delà de `max_direction_deviation_deg` — le
  seuil de surplomb est donc une condition d'admissibilité, pas un
  post-contrôle.
- `_collisions_at_points` réutilise les primitives de step9
  (`_PrintedGeometryIndex`, `_grow_printed_index_to`, `_evaluate_point`,
  `_surface_sample_points`), jamais une seconde implémentation de la
  détection. La re-vérification est réduite aux seuls points réorientés :
  les capsules ne dépendent que de `position`, jamais de `axis_angles`.
- `_violates_rate_limits` : garde-fou contre une réorientation
  mécaniquement inexécutable (limite #24).
- Les `Step8Diagnostic` **sautent le niveau 1** : l'étape 8 a déjà épuisé
  la recherche IK locale pour ces points, la relancer à l'identique ne peut
  rien produire de nouveau.
- Aucune exception pour une collision non résolue — bascule forcée en
  support de secours + `converged=False`.
- Limite connue non résolue : quand le hook rend une nouvelle trajectoire,
  la re-vérification est **complète** (une ré-segmentation invalide tout le
  champ scalaire, limite #23) — seule la voie « réorientation » bénéficie
  de la re-vérification partielle.

## Détail de `main.py` et `file_io/visualization_export.py`

- `main.run_pipeline(mesh_path, params) -> PipelineResult` : relais pur
  entre les 10 étapes, aucune logique métier propre. `PipelineParams` un
  champ par étape (mêmes défauts que chaque `Params` déjà acté) ;
  `PipelineResult` conserve TOUTES les sorties intermédiaires (pas
  seulement la trajectoire finale), pour qu'un export ou un futur
  appelant UI n'ait pas besoin de rejouer le pipeline pour retrouver une
  sortie d'étape.
- Dépendances multi-étapes respectées explicitement (pas seulement
  chaînées séquentiellement) : step4/5/7 reçoivent `mesh_data`/`tet_mesh`
  de step1 directement ; step8 consomme `colored_layers` (sortie de
  step7), jamais les `layers` d'entrée de step7 ; step9/step10 partagent
  la même instance de `MachineProfile` et de `CollisionSimParams`.
- `resegmentation_hook` de l'étape 10 est laissé à `None` dans `main.py` —
  l'implémenter réellement (re-jouer step3→step8 avec exclusion d'angles,
  limite #23) reste un travail à part, pas fait ici. Sans hook, toute
  collision non résolue par réorientation locale bascule directement en
  `fallback_support`, un comportement déjà explicite et testé, pas un trou.
- `_default_machine_profile()` charge `aetheris-4.glb` si présent (repli
  sur `nozzle_model=None`, log, jamais un crash) — cinématique par défaut
  = la machine actuelle documentée plus haut, un exemple d'usage parmi
  d'autres (limite #18), pas la référence figée.
- `file_io.visualization_export.export_visualization_glb(result, output_path)` :
  assemble une scène GLB (mesh source semi-transparent + couches en
  dégradé par ordre d'impression + trajectoire finale colorée par
  extrusion + marqueurs rouges aux collisions non résolues), ne fait
  AUCUN calcul de pipeline — relit uniquement ce que `PipelineResult`
  contient déjà. Choisi en GLB pour cohérence avec le modèle de tête déjà
  au même format, ouvrable dans n'importe quel visualiseur 3D classique.
- Deux bugs de performance réels trouvés en exécutant `main.py` de bout
  en bout pour la première fois sur le modèle de test (jamais visibles
  sur les mesh synthétiques ~10mm des tests unitaires) — voir limites
  #26 et #27 : le pas d'échantillonnage du modèle de tête confondu avec
  `nozzle_diameter`, et une boucle de distance point-capsule non
  vectorisée. Les deux corrigés dans `pipeline/step9_collision_sim.py`.

## Dev Tracker (convention de rédaction du suivi de développement)

Le « Dev Tracker » (« DT ») est la liste de suivi de ce qui a été développé et de
ce qui reste à faire. **À rédiger dans deux cas seulement** : quand c'est
nécessaire (fin d'un chantier, changement d'état important, reprise après une
interruption, décision de l'utilisateur qui déplace une priorité), ou quand
l'utilisateur le demande (« écris le DT »). Pas à chaque message.

**Marqueurs, un par ligne, toujours en tête :**
- 🟢 **résolu** : fait, ou question tranchée par l'utilisateur ;
- ⏳ **en cours** : tâche commencée et non terminée (au plus une ou deux à la fois) ;
- ⬛ **non traité, sans successeur** : rien ne dépend de ce point, il ne bloque rien ;
- ♦️ **à résoudre avant de continuer** : bloquant, c'est ce qui conditionne la suite.

**Structure (dans cet ordre) :**
1. *Ce qui a été développé* : 🟢 (correctifs, puis nouvelles fonctions) ;
2. *Réponses de l'utilisateur* : 🟢 ce qui est tranché ;
3. *À résoudre avant de continuer* : ♦️ ;
4. *Sans bloquer* : ⬛ ;
5. *Fait* : rappel condensé des 🟢 inchangés.

**Où il vit** : `DT.md` à la racine du projet (version courante, à mettre à jour à chaque
rédaction) ; `LAST_SESSION.md` en garde le contexte de session.

**Règles de rédaction :**
- **Vérifier avant de classer.** Un point n'est 🟢 que si le code ou une mesure le
  confirme, jamais sur la foi d'un document : `LAST_SESSION.md` a déjà annoncé
  « non commités » des correctifs qui l'étaient (vérifié le 2026-10-06).
- **Un point ♦️ dit pourquoi il bloque** et ce qui le débloquera.
- **Actualiser la liste existante** (changer les marqueurs, ajouter les nouveautés)
  plutôt que la réécrire d'un autre style.
- Citer les chiffres mesurés et les limites (#NN) plutôt que des généralités.

## Conventions de code à respecter

- Chaque `pipeline/stepN_*.py` expose `run(input, params) -> output`, fonction pure,
  sans état caché — permet à l'orchestrateur de ré-exécuter une portion du pipeline
  sans tout recalculer (utile pour l'étape 10, qui relance sélectivement 3 et 9).
- Un fichier de paramètres dataclass par étape (`OverhangParams`, `SegmentationParams`, etc.),
  actuellement définis dans chaque fichier `stepN_*.py` lui-même — à terme, `config/params.py`
  devra les centraliser sans casser les imports existants.
- Tests pytest par étape dans `tests/test_stepN_*.py`, mesh synthétiques générés à la
  volée (jamais de fichiers externes committés) via `trimesh.creation.*` et opérations
  booléennes (`union`/`difference`, moteur `manifold3d`).
- Toujours passer par un fichier temporaire + `step1_ingest.run()` dans les fixtures de
  test pour obtenir un `MeshData` réaliste (avec bvh/kdtree construits), plutôt que de
  mocker `MeshData` à la main.

## IDÉE RÉSERVÉE POUR LA TOUTE FIN — NE PAS IMPLÉMENTER MAINTENANT

**Sélection manuelle des parties, avec paramètres et ordre d'impression par
partie.** Idée de l'utilisateur, actée explicitement comme **le tout dernier
chantier**, à n'ouvrir **que lorsque l'algorithme automatique sera au point**.
Ce n'est pas une piste à explorer quand un mécanisme automatique résiste : la
tentation sera précisément d'y recourir comme échappatoire à un défaut non
résolu, et ce serait déplacer le problème sur l'utilisateur au lieu de le
résoudre. **Tant que le pipeline automatique n'est pas satisfaisant, ne pas
commencer ce travail, même partiellement.**

Contenu visé, tel qu'énoncé :
- un **algorithme de sélection** permettant à l'utilisateur de désigner
  lui-même les différentes parties du modèle (pilier, bras, zone suspendue…),
  plutôt que de subir la segmentation automatique ;
- des **paramètres propres à chaque partie** (direction, hauteur de couche,
  et tout autre réglage aujourd'hui global) ;
- un **ordre d'impression choisi par l'utilisateur** entre ces parties.

**Ce qu'il faut préserver dès maintenant pour que ce soit implémentable plus
tard, sans rien coder de la fonctionnalité elle-même** — ces points ne coûtent
rien s'ils sont respectés au fil de l'eau, et coûteraient une refonte s'ils
étaient perdus :

1. **Toute pondération spatiale doit rester définie PAR TÉTRAÈDRE, jamais
   comme un scalaire global.** Un champ `(T,)` calculé automatiquement
   aujourd'hui est exactement la même structure qu'un champ peint à la main
   demain — la sélection utilisateur devient une simple **source alternative**
   du même tableau. Un scalaire global, lui, n'a aucun chemin de migration.
   Vaut pour la modulation de lissage, les pondérations de direction, et tout
   terme ajouté au système de step4.
2. **Ne jamais supposer qu'une partition du modèle provient de step3.**
   Partout où du code consomme un découpage (régions, branches, clusters),
   il doit le recevoir en argument, jamais aller le chercher lui-même dans
   `RegionDependencyGraph`. Une partition fournie par l'utilisateur doit
   pouvoir se substituer à celle de step3 sans modifier l'appelé.
3. **L'ordre d'impression doit rester une donnée, pas une conséquence
   implicite d'un tri.** `_frontier_order` est déjà sur ce modèle : il produit
   une **permutation explicite** des couches. Conserver cette forme — un ordre
   imposé de l'extérieur remplace alors la permutation calculée, sans toucher
   à la géométrie des couches.
4. **`core/constraints.py` est le point d'extension prévu** pour tout apport
   utilisateur au champ scalaire (registre `CONSTRAINT_ASSEMBLERS`). Une
   sélection de partie assortie de paramètres est, formellement, une
   contrainte de plus — ne pas court-circuiter ce registre en câblant des
   termes utilisateur directement dans step4.

Aucune de ces quatre disciplines ne demande d'écrire du code aujourd'hui :
ce sont des **contraintes de forme** sur le code qui sera écrit de toute façon.

56. **PORTAGE C++ — pas 0 (socle de build) et pas 1 (index spatial) FAITS,
    MESURÉS et validés à métrique strictement identique. step9b divisé par
    5,3. Premier chantier du projet où le gain est net, reproductible et
    sans arbitrage : le portage ne change aucun résultat géométrique, il ne
    change que le temps.**

    Décision de l'utilisateur : passer la base de code en C++, cible **100 %
    C++** (Python supprimé à la bascule finale, UI TypeScript dialoguant
    directement avec la bibliothèque), migration **pilotée par le profil**,
    **pytest conservé comme oracle d'équivalence** pendant toute la
    migration, et **aucune bibliothèque dont l'usage commercial ne soit pas
    libre et gratuit**. Le plan complet est dans
    `~/.claude/plans/dapper-watching-balloon.md`.

    ### Le profilage qui a tout réorienté (et deux erreurs de lecture corrigées)

    `death.stl` @ résolution 0.8, cProfile — 467,1 s profilés pour 212,9 s
    réels (surcoût ×2,19) :

    | étape | cumtime | part |
    |---|---:|---:|
    | **step9b** | **323,3 s** | **69 %** |
    | **step8** | **92,8 s** | **20 %** |
    | step1 (dont 16,4 s de TetGen en sous-processus) | 29,7 s | 6 % |
    | step4 | 14,0 s | 3 % |
    | step6 / step5 / step2 / step3 / step7 | 7,1 s | 1,5 % |

    Une seule fonction domine tout : `step9_collision_sim.py:275
    candidates_near`, **119,0 s de temps propre**, plus les **1,41 milliard
    d'appels à `set.add`** (93 s) et 46,7 M d'appels à `dict.get` (8,3 s)
    qu'elle déclenche — **~46 % du temps total du pipeline à elle seule**.

    **Erreur de lecture n°1, corrigée** : une première catégorisation
    concluait « 42 % du temps est du natif déjà compilé, donc hors de portée
    d'un portage ». C'est faux. Les 94 s de `set.add`, 44 s de `numpy.array`
    et 8 s de `dict.get` sont de la manipulation de structures de données
    **provoquée par des boucles Python**, que le C++ supprime avec la
    boucle. Le vrai incompressible est TetGen (16,4 s) plus l'algèbre creuse
    (~7 s), soit **≈ 5 %**. Ne pas classer un appel natif comme
    incompressible sans regarder QUI l'appelle.

    **Erreur de lecture n°2, corrigée** : `step9b` est un diagnostic
    **désactivé par défaut** (`main.py`, `support_check: Optional[...] =
    None`). Les « 12 min » et les 172,7 s des sessions précédentes
    l'incluent parce que les scripts de mesure l'activent. Le pipeline de
    production step1→step8 coûte **89,7 s** (mesuré, pas déduit).

    ### Le socle (pas 0)

    - **MSVC 14.44 était déjà présent** (`D:\Prazlook\Programmes\vs` +
      Build Tools 2019/2022), mais ni CMake, ni Ninja, ni `cl.exe` au PATH.
      `pip install cmake ninja` dans le venv suffit (BSD-3 / Apache-2.0,
      sans droits administrateur).
    - `cpp/CMakeLists.txt` — C++20, `/O2`, **`/fp:precise`**. Ce dernier
      n'est pas un détail de confort : `/fp:fast` réordonnerait les
      opérations flottantes et rendrait l'oracle d'équivalence inutilisable.
    - `scripts/build_cpp.py` — localise l'installation via `vswhere`, capture
      l'environnement de `vcvars64.bat`, lance cmake/ninja dedans. Une seule
      commande, sans terminal développeur préalable.
    - `core/_cpp.py` — interrupteur `APEX_USE_CPP` (**Python par défaut**) et
      `APEX_CPP_GRID_MODE`. **L'absence du module natif n'est jamais une
      erreur** : sans `.pyd` compilé, tout retombe silencieusement sur
      Python. Lancer les tests ne doit jamais exiger un compilateur C++.
    - Le module se construit et se charge sur **Python 3.14 / MSVC
      19.44**, versions récentes où rien ne garantissait que pybind11 3.1
      passerait — vérifié, pas supposé.

    ### Le pas 1 : `_PrintedGeometryIndex`

    Deux défauts **distincts**, qu'il fallait surtout ne pas confondre :

    1. **coût d'interprétation** — une liste de tuples matérialisée par
       requête, un `dict.get` par cellule, un `set.add` par (cellule ×
       capsule) ;
    2. **défaut ALGORITHMIQUE** — une capsule est insérée dans *toutes* les
       cellules que son AABB touche, donc **revisitée une fois par cellule**
       à chaque requête.

    D'où deux modes de requête coexistant dans la classe C++
    (`QueryMode::Faithful`, transcription fidèle du `set` ; `Stamped`,
    dédoublonnage par estampilles) : c'est ce qui permet de mesurer les deux
    gains **séparément**, avec un seul build.

    **Découverte faite en relisant les cinq sites d'appel, absente du plan
    initial** : trois d'entre eux (`step9b._has_support_below`,
    `step9b._has_support_in_any_direction`, `step5._layer_support_mask`)
    reconstruisaient **immédiatement** deux tableaux `(C,3)` par
    compréhension Python sur les objets rendus — 432 candidats en médiane,
    817 653 fois. C'est la source des 41 s de `numpy.array`. Une méthode
    dédiée `candidates_endpoints` leur rend directement les tableaux, ce qui
    supprime d'un coup l'ensemble de dédoublonnage, la liste d'objets et les
    deux compréhensions. Les deux autres consommateurs (`step9`, `step8b`)
    bouclent capsule par capsule et gardent `candidates_near` : ils ne
    changent pas d'une ligne.

    **Les cinq consommateurs ont été relus un par un et vérifiés
    insensibles à l'ordre** (ils cherchent un minimum de distance ou
    l'existence d'un support) — indispensable, puisque l'ordre du Python est
    celui d'itération d'un `set`, donc déjà non significatif. Les tests
    d'équivalence comparent des ensembles triés, jamais des listes.

    ### Résultats mesurés (`death.stl` @ 0.8)

    | configuration | total | step9b seul | métrique step9b |
    |---|---:|---:|---|
    | Python (référence) | 212,9 s | **123,2 s** | 0,466578 % / 1085 défauts d'ordre |
    | C++, **même algorithme** | 126,2 s | **36,5 s** (**×3,38**) | **identique** |
    | C++, **algorithme corrigé** | 112,9 s | **23,2 s** (**×5,31**) | **identique** |

    Part fixe step1→8, mesurée séparément : **89,7 s** — sans elle, les
    totaux auraient laissé croire à un gain de ×1,89 seulement, en diluant
    step9b dans du temps qui n'a pas bougé.

    Décomposition honnête : **le langage seul apporte ×3,38**, la correction
    algorithmique **×1,57 de plus**. Deux chiffres, pas un.

    **Validation croisée la plus solide obtenue sur ce projet** : le
    compteur C++ relève **1 406 655 778 visites** (cellule × capsule) — soit
    *exactement* le nombre d'appels à `set.add` mesuré indépendamment par
    cProfile côté Python. Les deux implémentations parcourent rigoureusement
    la même chose ; l'équivalence n'est pas seulement testée, elle est
    comptée.

    ### Pas 2 : la cinématique inverse (step8) — et le piège du repli « tout ou rien »

    Après le pas 1, le re-profilage place **step8 en tête avec 86,1 s sur
    168,0 s profilées (51 %)**, `step9b` retombant à 21 % et `candidates_near`
    passant de 212 s de temps propre à **3,28 s**. La priorité a donc basculé,
    et suivre le plan plutôt que la mesure aurait été une erreur.

    Portage : toute la chaîne analytique de `core/geometry_utils.py`
    (`_rodrigues_rotate_single`, `_angle_between`, `_perpendicular_basis`,
    `_signed_angle_in_plane`, `_clamp_into_range`, `_forward_direction`,
    `_resolve_single_axis_analytic`, `_resolve_axis_pair_analytic`,
    `resolve_axis_angles_analytic`) plus la boucle de step8 avec
    `_axis_chain_world_lines`, `_detect_singular_axes` et l'héritage d'angle
    à la singularité → `cpp/src/core/kinematics.cpp`.

    **La boucle N'EST PAS parallélisable**, contrairement à ce que supposait
    le plan de portage. `previous_angles` porte une dépendance séquentielle
    réelle : repliement des axes internes (`axes[2:]`), repli en cas de
    singularité locale, héritage d'angle. Le gain vient donc du seul
    changement de langage. Paralléliser exigerait de casser cette chaîne,
    c'est-à-dire de changer le COMPORTEMENT, pas le langage.

    #### Le défaut qui annulait 100 % du gain

    Première version : au premier point que la résolution analytique ne sait
    pas traiter, le C++ rendait la main et **Python rejouait toute la boucle**.
    Mesure : step1→8 à 89,5 s contre 89,7 s avant portage — **rien**.

    Cause, trouvée en instrumentant plutôt qu'en supposant : le repli se
    déclenche **au point 0**. `_resolve_axis_pair_analytic` échoue quand la
    direction visée est **colinéaire à l'axe externe** — le résidu du
    bracketing y est une **tangence** à zéro, jamais un changement de signe,
    donc aucun encadrement possible. Mesuré sur la cinématique réelle
    (C tour complet + B ∈ [0,185]) : **0 échec sur 2000 directions
    aléatoires**, mais échec systématique sur `+Z`. Or une pièce à fond plat
    commence par déposer verticalement.

    **Un point sur 297 994 (0,0003 %) annulait donc l'intégralité du gain.**

    #### Le mécanisme retenu : reprise, pas abandon

    Le C++ remplit tout jusqu'au point fautif, rend son index ; Python résout
    **ce seul point** avec la recherche par grille déjà écrite et testée ;
    le C++ reprend au suivant en recevant les angles via
    `initial_previous_angles`. La chaîne est passée explicitement, donc
    l'équivalence est préservée. Mesure sur `death.stl` : **412 replis sur
    297 994 points (0,14 %)**, le natif traitant tout le reste en 1,43 s
    cumulées — et step1→8 tombant de 89,5 s à **42,8 s**.

    La grille de repli n'est **délibérément pas portée** : elle sert 0,14 %
    des points ; la porter serait du code à maintenir et une occasion de
    diverger pour rien.

    > **Leçon de méthode, valable au-delà de ce cas** : sur un chemin
    > SÉQUENTIEL, une stratégie de repli « tout ou rien » est fragile de façon
    > disproportionnée — la granularité du repli compte plus que la qualité du
    > portage. Vérifier le TAUX de repli réel avant de choisir cette
    > granularité, jamais se fier à un « ce cas n'arrive jamais » écrit dans
    > une docstring (celle de `_resolve_axis_angles` affirmait exactement
    > cela).

    #### Un test qui ne testait rien

    Les quatre premiers tests d'équivalence de l'IK passaient au **vert sans
    jamais exécuter le C++** : le repli silencieux les faisait comparer Python
    à lui-même. Corrigé par `_LAST_NATIVE_FALLBACK_COUNT`, que les tests
    vérifient désormais, plus un test dédié qui **force** le repli sur trois
    points (dont le premier) et exige qu'il ait eu lieu avant de comparer.

    C'est la deuxième fois dans la même session qu'un mécanisme donne
    l'apparence de fonctionner : voir aussi la « erreur de lecture n°1 »
    ci-dessus. Un indicateur vert n'est pas une preuve qu'on mesure ce qu'on
    croit.

    #### Tolérances : ce qui est exact et ce qui ne peut pas l'être

    Les **angles** concordent à **1e-9 degré**. Les **résiduelles** non : elles
    valent `acos(dot(achieved, target))` avec un produit scalaire à ~1−1e-15,
    or `acos` est mal conditionnée en 1 (dérivée ~1e7 ici), si bien qu'un écart
    d'**un ULP** produit ~1e-8 rad. L'écart est du bruit de conditionnement,
    pas une divergence de calcul. La tolérance est donc **déclarée**
    (1e-6 rad) et complétée par l'assertion qui compte réellement : les deux
    implémentations émettent **exactement les mêmes `Step8Diagnostic`**
    (seuil à 5°, quatre ordres de grandeur au-dessus du bruit).

    ### Pas 3 : le test de support — quand l'ALGORITHME compte plus que le langage

    Après le pas 2, le re-profilage remet `step9b` en tête (40,3 s sur 98,2 s
    profilées) et désigne `core/geometry_utils.point_has_support_along_direction`
    comme premier poste PORTABLE : **11,97 s de temps propre pour 273 146
    appels**, plus ~10 s de dispatch numpy qu'elle déclenche (`cross`,
    `normalize_axis_tuple`, `moveaxis`). Seul TetGen la dépasse (16,6 s), et il
    est déjà natif.

    **Or cette fonction était DÉJÀ entièrement vectorisée.** La transcrire ligne
    à ligne n'aurait gagné que le dispatch. Le vrai défaut était ailleurs, et
    c'est le premier des trois portages où il a fallu changer autre chose que
    le langage :

    1. **Aucun arrêt anticipé.** Elle ne rend qu'un BOOLÉEN et exécutait
       pourtant ~15 opérations numpy sur les ~690 candidats à chaque appel.
       `np.any(...)` est vrai dès qu'UN candidat convient : s'arrêter au
       premier donne exactement la même réponse pour un coût proportionnel au
       rang du premier appui, et non à 690.
    2. **Aller-retour inutile.** Le C++ rendait deux tableaux `(C,3)` (3,95 s
       d'allocation mesurées dans `candidates_endpoints`) que Python
       reconsommait aussitôt.

    **Fix : fusionner la requête spatiale et le test géométrique**
    (`_PrintedGeometryIndex.has_support`, `cpp/src/core/support.cpp`). Le natif
    lit les extrémités **directement dans la grille**, ne matérialise jamais
    les candidats côté Python, et court-circuite au premier appui. Un seul
    portage sert les trois sites d'appel (`step9b._has_support_below`,
    `step9b._has_support_in_any_direction`, `step5._layer_support_mask`).

    Résultat : **step9b passe de 29,2 s à 9,6 s**, soit **×14,4 depuis le
    Python d'origine** (138,6 s).

    > **Leçon** : les pas 1 et 2 étaient de la transcription, et le langage y
    > suffisait. Ici non. Avant de porter une fonction DÉJÀ vectorisée,
    > demander ce qu'elle calcule *en trop* — ici, 689 candidats sur 690 dans
    > le cas favorable.

    Deux gardes ajoutées aux tests, toutes deux nécessaires : une qui vérifie
    que le jeu d'essai produit assez de cas **positifs** (sinon un bug rendant
    tout `False` des deux côtés serait vert), et un test du cas limite
    « strictement en dessous » — une capsule exactement à la hauteur du point
    ne le supporte pas (`delta_up > 1e-9`), condition facile à perdre en
    transcrivant et invisible sur des positions aléatoires.

    ### Pas 4 : tables topologiques + granularité d'appel

    Le re-profilage après le pas 3 (66,5 s profilées) désigne trois grappes
    portables, dont **deux sont la même opération** — construire une table de
    sous-simplexes avec un dictionnaire à clés tuple :

    | cible | coût mesuré |
    |---|---|
    | `step4._tet_face_table` | 1,79 s + 1,78 s de générateur sur **12,5 M appels**, plus une part des 2,37 s de `sorted` et 1,43 s de `setdefault` → **~6-7 s** |
    | `step1._build_tet_adjacency` | 2,91 s + 1,35 s de `set.add` → **~4,3 s** |
    | construction de l'index | `_segment_aabb` 1,79 s + `append_segment` 1,01 s + `apex_core.append` 1,19 s → **~4 s** |

    `cpp/src/core/tet_topology.cpp` sert les deux premières : les faces sont
    empaquetées en enregistrements triés une fois, les arêtes en `uint64`
    (`min << 32 | max`) pour que le dédoublonnage soit un `std::unique` plutôt
    qu'un hachage de tuples.

    **La troisième n'était pas un problème de langage mais de GRANULARITÉ
    D'APPEL** : Python franchissait la frontière **744 324 fois** pour insérer
    les capsules une par une. Deux corrections :
    - `append_segments` insère tout le lot en un appel ;
    - en mode natif, **aucun objet `_PrintedSegment` n'est plus créé à
      l'insertion**. Ils sont construits à la demande dans `candidates_near`,
      pour les ~690 candidats réellement rendus et non pour les 744 324
      capsules insérées — dont la quasi-totalité n'est jamais lue.

    > Erreur commise puis corrigée en cours d'écriture : la première version de
    > l'insertion en lot reconstruisait quand même les 744 324 objets avec leur
    > AABB, c'est-à-dire exactement le coût qu'elle prétendait supprimer.
    > Batcher un appel ne sert à rien si l'on garde le travail qu'il portait.

    Un test verrouille le point le plus fragile : **l'ordre d'INSERTION des
    faces frontière**. `_map_boundary_faces_to_surface` interroge ensuite le
    kdtree avec ces faces dans cet ordre, et l'ordre du résultat s'y propage —
    une permutation décalerait des tableaux en aval **sans lever aucune
    exception**.

    ### PARALLÉLISATION : il n'y en a AUCUNE, et c'est délibéré

    Vérifié, pas supposé : aucun `#pragma omp`, aucun `std::thread`, aucun
    `threading`/`multiprocessing`, et `/openmp` n'est même pas activé au build.
    Mesure directe du **rapport CPU/réel = 0,66** :

    | | temps |
    |---|---:|
    | notre code (CPU du processus) | **28,4 s** — mono-thread |
    | attente du sous-processus TetGen | **~15 s** — mono-thread, immuable |
    | total réel | **43,3 s** |

    Jamais au-dessus de 1,0 : aucune bibliothèque tierce ne parallélise dans
    notre dos. En dessous de 1,0 : le déficit EST l'attente de TetGen, dont le
    CPU est compté dans un autre processus. Recoupe exactement les 15,7 s
    relevées au profilage.

    **Ce qui est parallélisable, et ce qui ne l'est pas** :
    - `step1`/TetGen — **non**, mono-thread et conservé : **plancher série dur** ;
    - `step8`/IK — **non en l'état**, chaîne séquentielle via `previous_angles` ;
    - `step9b` — **partiellement** : le contrôle contre l'index complet est
      indépendant par point, la marche incrémentale ne l'est pas ;
    - `step4` (assemblage FEM, matvec CG), `step5`, `step6` — **oui**.

    **Plafond calculé** : en parallélisant parfaitement les 28,4 s sur 8 cœurs,
    ~3,5 s + 15 s de TetGen ≈ **18,5 s**. Ni le portage seul, ni la
    parallélisation seule n'atteindront 10 s.

    **Pourquoi ne pas paralléliser maintenant**, et ce n'est pas un oubli :
    1. les réductions flottantes parallèles **cassent l'équivalence bit-à-bit**,
       c'est-à-dire l'oracle sur lequel repose tout le portage ;
    2. paralléliser un fragment C++ appelé depuis une boucle Python ne rapporte
       presque rien, l'extérieur restant sérialisé par le GIL.

    Séquence retenue : **porter d'abord, paralléliser ensuite** — quand tout
    sera en C++ et que l'oracle Python aura joué son rôle. Le GIL est déjà
    libéré dans les points d'entrée natifs, donc rien ne sera à réécrire.

    ### POINT D'INFLEXION : la phase « point chaud » est TERMINÉE

    Profil après le pas 4 (43,2 s profilées) : **plus aucun poste Python isolé
    ne dépasse ~2,5 s**. Les quatre pas visaient chacun une grappe de 20 à
    200 s ; il n'en reste aucune.

    | poste | temps | nature |
    |---|---:|---|
    | TetGen | 12,49 s | **immuable**, mono-thread |
    | `apex_core.has_support` | 2,49 s | déjà C++ (coût d'appel × 272 868) |
    | `csr_matvec` (CG) | 2,28 s | scipy natif |
    | `_blend_tet_directions` | 2,04 s | Python — exige un KD-tree, et scipy y est déjà en C |
    | `_interpolate_out_of_plane` | 1,50 s | Python, step6 |
    | `apex_core.resolve_axis_angles_batch` | 1,43 s | déjà C++ |

    `step1` seul mesuré à **18,8 s dont 12,4 s d'attente de TetGen** : il ne
    reste que ~6,4 s de code propre, et elles sont éparpillées (2,44 s de pur
    dispatch numpy sur des vecteurs 3D, ~0,6 s d'auto-intersections, le reste
    sous 0,25 s par poste).

    **PIÈGE À NE PAS REFAIRE** : `step2` (318 lignes, 0,1 s) semblait le banc
    d'essai idéal pour porter les structures de données — jusqu'à ce que sa
    lecture montre que sa détection de ponts repose sur
    `mesh_data.bvh.intersects_location`. Le porter imposerait d'écrire un BVH
    complet **pour 0,1 s de gain**. **Choisir un module par sa TAILLE est
    trompeur : ce sont ses DÉPENDANCES qui déterminent le coût du portage.**

    > Suite, au pas 6 : ce BVH a finalement été écrit, et le constat ci-dessus
    > reste exact — `step2` ne gagne que 0,079 s sur `death.stl`. Ce qui a
    > justifié le travail n'est ni sa taille ni son temps propre, mais qu'il
    > **débloque** step2/step7 et que le gain croît avec la densité du
    > maillage (×221 sur `death_hi-res.stl`).

    Conséquence actée : le portage ne peut plus être piloté par le profil. Il
    doit devenir **systématique et guidé par les dépendances** — d'abord
    `types.hpp` réel, puis `core/bvh.cpp` (débloque step1/2/7/9), puis
    `core/kdtree.cpp` (13 usages de `cKDTree`), puis les étapes
    une par une. Les gains de 1-2 s viendront en supplément, sans être
    cherchés.

    ### Le pas 5 : `core/kdtree.cpp` — arbre k-d spatial

    **Portage ARCHITECTURAL, pas une optimisation, et assumé comme tel.**
    `scipy.spatial.cKDTree` est déjà natif (C) : le gain de vitesse attendu
    est nul à légèrement négatif. Le but est de retirer la **dernière
    primitive spatiale de scipy du cœur** — après quoi il ne reste plus, côté
    scipy, que l'algèbre creuse (`scipy.sparse`, CG) et `ProximityQuery`
    (trimesh). C'est la première application de la règle « guidé par les
    dépendances » ci-dessus : aucun poste chaud ne le désignait.

    **nanoflann écarté** au profit d'un **k-d tree maison** (~180 lignes,
    `cpp/src/core/kdtree.cpp`) : median-split sur l'axe le plus étendu, 1-NN
    récursif avec élagage, requête par rayon. Cohérent avec `types.hpp`
    (« aucune dépendance tierce à ce stade ») et suffisant pour le contrat
    réel — voir ci-dessous.

    **Périmètre verrouillé par lecture des 9 sites `.query(...)` du
    pipeline** (step1 `mesh_data.kdtree`, step3, step4 ×3, step5 ×2, step6,
    step7) : **tous** passent un tableau `(M,3)` et demandent le plus proche
    voisin unique (`k=1`, défaut de scipy). Aucun `k>1`, aucune requête
    arbre-contre-arbre. Plus un seul site `query_ball_point` (un point, un
    rayon : `core/constraints.py`, re-filtré derrière donc insensible à
    l'ordre). Le module n'expose donc que `query_nn` (1-NN par lot) et
    `query_radius` — porter davantage serait du code à maintenir pour un cas
    que personne n'appelle.

    **Équivalence.** `query_nn` rend la **distance euclidienne** (`p=2`,
    `sqrt(dx*dx+dy*dy+dz*dz)`, accumulation `x→y→z`, `/fp:precise`) :
    l'oracle confirme **indices identiques ET distances identiques au bit**
    face à `cKDTree.query` sur 2500 requêtes / 6000 points aléatoires. Le
    tie-break (deux points exactement équidistants) est tranché par **plus
    petit indice** ; scipy ne le spécifie pas, mais les données réelles
    (centres de tétraèdres/faces, sommets flottants) ne produisent jamais
    d'égalité exacte — l'oracle utilise des données non dégénérées pour cette
    raison, et le documente.

    **Câblage : `core/spatial.py`**, enveloppe mince (même discipline que
    `_PrintedGeometryIndex._native`). Elle expose exactement le sous-ensemble
    utilisé — `.query()` (k=1, gère l'entrée `(3,)` → scalaires comme
    scipy), `.query_ball_point()`, l'attribut `.data`. Choix scipy/natif
    **figé à la construction**. Les 6 sites qui construisaient
    `scipy.spatial.cKDTree(...)` passent à `core.spatial.KDTree(...)` ; les 9
    sites `.query()` sont **inchangés** (compatibles canard). `import
    scipy.spatial` retiré de step5/step7 (devenu inutilisé).

    **Recette : 326/326 tests dans les DEUX positions de l'interrupteur**,
    + 5 nouveaux tests d'équivalence (pas 5) : `query_nn` bit-exact vs
    scipy, `query_radius` comme ensembles, dispatch du wrapper, `.data`
    préservé, entrée 1-D → scalaires, non-régression interrupteur éteint.
    Métrique géométrique inchangée (le portage ne déplace aucun résultat).
    Gain de vitesse **non mesuré et non attendu** — conforme à l'intention.

    ### Le pas 6 : `core/bvh.cpp` — et LA LEÇON LA PLUS CHÈRE DU PORTAGE

    Cible : `trimesh.ray.ray_triangle.RayMeshIntersector`. Contrairement au
    k-d tree, **la référence Python n'est PAS native** : `pyembree` est absent
    de ce venv (vérifié, pas supposé), donc trimesh tombe sur son
    implémentation pure numpy. Le gain est donc réel, et il est gros.

    Périmètre verrouillé par lecture : **un seul point d'entrée est appelé**,
    `intersects_location` (step2 détection de ponts, step7 deadline
    d'accessibilité). `intersects_first`, `intersects_any` et
    `contains_points` ne sont appelés nulle part — non portés. Les deux
    appelants **ignorent `index_tri`** et réduisent par `.min()`, donc
    l'oracle porte sur l'**ensemble** des couples `(index_ray, location)`
    dédoublonnés, jamais sur l'ordre (même discipline qu'au pas 1).

    #### Le bug : le broad phase N'EST PAS une optimisation

    Première version : BVH descendu par un test rayon/AABB classique
    (tranches), en supposant que le broad phase n'était qu'un filtre
    conservateur — deux sur-ensembles différents, la narrow phase tranchant
    ensuite de la même façon. **Faux.**

    Les micro-mesures disaient pourtant « identique au bit » : 3 000 rayons
    aléatoires, directions verticales / quelconques / non normalisées,
    origines lointaines, rayons ratant tout — tout concordait, sur `death.stl`
    comme sur `turkish.stl`. **C'est la mesure sur le PIPELINE qui a révélé le
    défaut** : `step2` passait de **49 ponts détectés à 0**.

    Cause : `ray_bounds` clippe le paramètre `t` à `buffer_dist = 1e-5` avant
    de construire la boîte du rayon, ce qui **exclut des candidats la face
    d'où part le rayon**. La narrow phase, elle, la laisserait passer :
    `step2` décale son origine de 1e-6 sous la face, l'auto-impact tombe à
    `dot = -1e-6` — exactement le seuil du filtre « vers l'avant » — et le
    calcul flottant réel donne `-9.999999994736442e-07`, strictement
    supérieur au seuil. Autrement dit : **c'est le broad phase de trimesh,
    plus strict que sa propre narrow phase, qui évite l'auto-intersection.**

    > **Leçon, à appliquer à tout portage restant** : une structure
    > d'accélération n'est « qu'une optimisation » seulement si son critère
    > est un sur-ensemble STRICT du critère exact. Quand elle est plus
    > stricte quelque part — ici d'un facteur 10 sur une tolérance —, elle
    > fait partie du comportement observable et doit être transcrite, pas
    > améliorée. Et : **des rayons aléatoires ne testent pas les cas que le
    > pipeline produit réellement.** Aucun des 3 000 ne partait d'une face.

    #### Le mécanisme retenu

    Le BVH remplace **rtree comme INDEX**, pas la **REQUÊTE**. La boîte du
    rayon reste celle de `ray_bounds`, transcrite à l'identique ; le BVH sert
    à trouver vite les triangles dont l'AABB la recoupe — exactement ce que
    faisait `tree.intersection(bounds)`. Les AABB de triangles sont donc
    **exactes et sans marge**, comme celles indexées par
    `triangles.bounds_tree`. Le jeu de candidats est identique **par
    construction**, plus par chance.

    Effet secondaire heureux, mesuré : c'est aussi **plus rapide** que la
    version fautive, la boîte d'un rayon vertical étant une colonne très fine
    que le BVH élague bien.

    La narrow phase est une **transcription**, pas un meilleur algorithme :
    `planes_lines` (ligne infinie, seuil absolu 1e-5), barycentrique par
    Cramer dans le même ordre d'opérations, appartenance à `±tol.zero=1e-13`,
    filtre avant à `-1e-6`, dédoublonnage `unique_rows` quantifié à 8
    décimales (`tol.merge`). Écrire un Möller-Trumbore aurait été plus court
    et plus rapide — et aurait donné d'autres flottants, donc un oracle
    inutilisable.

    #### Mesures

    | modèle | faces | trimesh | natif | gain |
    |---|---:|---:|---:|---:|
    | `death.stl` | 2 418 | 0,057 s | 0,0015 s | **×38** |
    | `turkish.stl` | 10 180 | 0,106 s | 0,0018 s | **×60** |
    | `death_hi-res.stl` | 226 156 | 0,927 s | 0,0042 s | **×221** |

    (1 500 rayons verticaux, impacts **identiques au bit** dans les trois cas.)

    Sur le pipeline, `step2` sur `death.stl` : **0,102 s → 0,023 s**, à
    `bridges=49 / overhangs=359 / risk` identiques. Le gain absolu reste
    petit *sur ce modèle* — c'est le piège documenté au point d'inflexion
    ci-dessus, `step2` ne pesait que 0,1 s — mais il croît avec la densité du
    maillage, et il débloque le portage de step2/step7.

    **Câblage** : `core.spatial.RayIntersector`, même discipline que `KDTree`
    (choix figé à la construction). `step1._build_acceleration_structures` la
    construit, `MeshData.bvh` la porte, step2/step7 l'appellent sans
    changement. Une option passée au chemin natif **lève** au lieu d'être
    ignorée en silence.

    **7 tests d'équivalence** (pas 6), dont deux qui n'existeraient pas sans
    le bug : celui qui verrouille l'exclusion de la face d'origine sur une
    face plate, et celui qui vérifie que sur une face non plate le natif fait
    **le même choix que trimesh** — y compris quand trimesh accepte
    l'auto-impact. L'oracle exige la conformité, pas la propreté.

    ### BILAN CONSOLIDÉ après les pas 1 + 2 + 3 + 4

    | | Python | C++ | gain |
    |---|---:|---:|---:|
    | **pipeline complet** | 237,1 s | **43,3 s** | **×5,5** |

    Métrique identique à chaque pas : 0,466578 %, 1085 défauts d'ordre, 0
    surplomb réel. **326/326 tests verts** dans les deux positions.

    ### Détail des pas 1 + 2 + 3

    | | Python | C++ | gain |
    |---|---:|---:|---:|
    | **pipeline complet** | 237,1 s | **51,2 s** | **×4,6** |
    | step1→8 | 98,5 s | 41,6 s | ×2,4 |
    | **step9b** | 138,6 s | **9,6 s** | **×14,4** |

    Métrique **identique** à chaque étape : 0,466578 %, 1085 défauts d'ordre,
    0 surplomb réel. **321/321 tests verts** dans les deux positions de
    l'interrupteur.

    Le rapport de force s'est inversé : `step9b` ne pèse plus que 19 % du
    temps, `step1→8` en fait 81 %.

    ### Détail des pas 1 + 2 (quatre configurations d'affilée)

    Les quatre mesures ci-dessous ont tourné **dans la même session**, seule
    façon d'obtenir des soustractions fiables vu la variance de la machine
    (voir section suivante). Elles se recoupent : 77,2 = 48,0 + 29,2.

    | | Python | C++ | gain |
    |---|---:|---:|---:|
    | **pipeline complet (step1→8 + step9b)** | **237,1 s** | **77,2 s** | **×3,07** |
    | step1→8 | 98,5 s | 48,0 s | ×2,05 |
    | step9b | 138,6 s | 29,2 s | **×4,75** |

    Métrique **strictement identique** dans les deux cas : `unsupported_fraction`
    = 0,466578 %, 1085 défauts d'ordre, 0 surplomb réel.

    **À garder en tête avant de citer ce ×3,07** : il est obtenu avec
    **1040 lignes de C++** face à ~14 900 lignes de Python encore en place,
    soit ~7 % du code. Les étapes 1 à 7 et 10 tournent toujours intégralement
    en Python. C'est le rendement d'un ordre piloté par le profil, pas celui
    d'un moteur porté — et cela signifie aussi que le gain suivant sera
    beaucoup plus cher, puisque le temps restant est réparti sur beaucoup plus
    de code (step1 ~25 s dont 14 s de TetGen déjà natif et incompressible,
    step4 ~14 s, le reste de step9b).

    **Conséquence pour l'objectif commercial de 10 s** : même en portant
    parfaitement tout le Python restant, l'ordre de grandeur atteignable à
    cette résolution est ~25-35 s, dominé par TetGen et l'algèbre creuse.
    Descendre sous 10 s exigera aussi des changements ALGORITHMIQUES — le
    facteur de duplication de 7,46 sur l'index et les 690 candidats par
    requête en sont deux, déjà chiffrés. À dire maintenant, pas à l'arrivée.

    ### Variance de mesure sur cette machine : ~9-11 %

    Deux exécutions de la MÊME configuration Python ont donné 89,7 s puis
    97,5 s, et 212,9 s puis 237,1 s. **Toute comparaison par soustraction
    entre runs séparés est donc suspecte** — les configurations à comparer
    doivent tourner d'affilée, dans les mêmes conditions.

    Piège concret rencontré : une variante « step1→8 seul » codée avec
    `APEX_USE_CPP=0` restait valable tant que step1→8 était identique dans les
    deux modes, et est devenue fausse dès le portage de step8. Elle affichait
    97,5 s pour step1→8 alors que le pipeline C++ COMPLET faisait 76,3 s —
    une impossibilité qui a révélé le mélange. **Une mesure de référence
    devient obsolète dès que ce qu'elle mesure change.**

    ### Vérification

    - **331/331 tests verts avec `APEX_USE_CPP=0` ET `=1`** (305 d'origine
      + 26 d'équivalence, dont 5 pour le k-d tree du pas 5). Mesuré à
      185,7 s (C++) contre 187,6 s (Python) : les deux positions de
      l'interrupteur passent la suite, à durée équivalente — le pas 5
      n'était pas un portage de vitesse (voir sa section).
    - **Pas 6 (BVH) : 340/340 dans les deux positions** (331 + 9, dont 7
      d'équivalence pour le BVH lui-même) — 181,7 s (C++) contre 183,6 s
      (Python), toujours équivalent en durée sur `death.stl` : le gain du
      pas 6 est diagnostiqué mais pas dominant sur ce modèle (voir sa
      section, « débloque » plutôt que « accélère massivement » ici).
    - **Pas 7 (step4) : 351/351 dans les deux positions** (340 + 7
      d'équivalence pour step4, + 4 sur la réutilisation de la sonde) — 198,8 s (C++) contre 195,5 s (Python). La
      bit-exactitude du fix de la sonde est prouvée par un A/B **en
      processus** (`test_probe_is_reused_when_constraints_add_nothing`,
      `array_equal`), pas par une comparaison de sommes de contrôle : un
      checksum imprimé à 9 décimales est arrondi, donc incapable de
      distinguer l'égalité au bit d'un écart à 1e-15.
    - `tests/test_cpp_equivalence.py` : candidats identiques au
      Python requête par requête (150 requêtes × 2 marges, avec une garde
      contre le faux vert où *tout* serait vide des deux côtés), concordance
      des deux modes, cohérence de `candidates_endpoints` avec
      `candidates_near` **dans les deux implémentations**, index vide, et
      non-régression du défaut (interrupteur éteint ⇒ aucun objet natif).
    - Métrique `step9b` **strictement identique** aux trois configurations.
      Un portage qui déplacerait la métrique serait un bug, pas un progrès.

    ### Ce que la mesure impose pour la suite — LA PRIORITÉ A BASCULÉ

    step9b ne pèse plus que **23,2 s sur 112,9 s** ; les **89,7 s de
    step1→8 sont désormais 79 % du temps**, dont `step8` est de loin le plus
    gros poste. Continuer à optimiser step9b serait suivre le plan au lieu
    de suivre la mesure — exactement ce que ce fichier reproche partout
    ailleurs.

    Deux chiffres cadrent le reste du travail sur l'index, et ils ne disent
    pas la même chose :

    - **facteur de duplication 7,46** — chaque capsule rendue est rencontrée
      7,46 fois. C'est la limite #27/#30, enfin chiffrée. Les estampilles
      suppriment le *coût* du dédoublonnage, pas les visites ; seul un BVH
      (une capsule référencée une seule fois) les supprimerait.
    - **690,6 candidats uniques par requête** — et ce chiffre-là **ne dépend
      pas de la structure spatiale** : un BVH parfait rendrait les mêmes 690
      candidats. C'est un plancher, pas un défaut d'implémentation.

    Conséquence : le levier suivant sur step9b n'est pas un BVH mais
    **l'arrêt anticipé**. `core/geometry_utils.point_has_support_along_direction`
    exécute ~15 opérations numpy sur 690 éléments pour ne rendre qu'un
    **booléen**, 272 868 fois, sans jamais court-circuiter dès le premier
    appui trouvé. À garder pour quand step8 aura été traité.

    ### Le pas 7 : `step4_field_solver.py` — et LE PORTAGE QUI NE PAYE PAS

    **Résultat principal, contre-intuitif et mesuré : `step4` passe de
    9,45 s à 4,60 s (×2,05) — et la moitié de ce gain vient de code PYTHON,
    pas du C++. Le solveur Eigen, seule pièce que l'utilisateur avait
    explicitement demandé d'inclure, se mesure à ×1,01 : il ne rapporte
    RIEN.**

    #### Étape 0 — un système résolu DEUX FOIS (fix Python, bit-exact)

    Trouvé en profilant avant d'écrire du C++ : `build_static_system` résout
    le système une seconde fois **uniquement pour renseigner
    `UpwardGrowthDiagnostic`**, dont le résidu doit être mesuré sur le champ
    réellement résolu. Or `main.py` active `enforce_upward_growth=True`, et
    passe `constraints=None` — donc `solve_with_constraints` re-résout
    **exactement le même système** : `A` et `b` identiques, vérifié.

    Fix : la sonde est conservée (`FieldSolverStaticSystem.probe_values`) et
    réutilisée quand les trois conditions sont réunies — une sonde existe,
    les contraintes ont une contribution strictement nulle (`nnz == 0`,
    testé et non supposé), et aucun `initial_guess` n'est fourni (le
    solveur itératif partant d'ailleurs ne rend pas le même vecteur). Le
    tableau est **copié**, jamais aliasé : `solve_with_constraints` est
    appelable plusieurs fois sur le même système statique, et une mutation
    par l'appelant reviendrait sinon dans le résultat suivant.

    **9,45 s → 7,47 s, champ identique AU BIT** (checksum
    `1257825.311560218` inchangé). Zéro C++.

    > **Pourquoi le faire AVANT le portage, et pas après** : sans ce fix, la
    > réécriture C++ aurait supprimé le solve redondant au passage et se
    > serait attribué ses 2,0 s. C'est exactement l'erreur de mesure que ce
    > fichier documente ailleurs (« une mesure de référence devient obsolète
    > dès que ce qu'elle mesure change »).

    #### Ce qui a été porté, et ce que ça rapporte

    | poste | avant | après | gain |
    |---|---:|---:|---:|
    | `_blend_tet_directions` | 2,35 s | **0,72 s** | **×3,2** |
    | `_compute_shape_gradients` | 0,36 s | **0,012 s** | **×30** |
    | solveur CG (Eigen) | 1,83 s | 1,81 s | **×1,01 — RIEN** |

    **`_blend_tet_directions` : le gain est ALGORITHMIQUE, pas linguistique.**
    Le Python interroge un k-d tree par région avec la totalité des centres
    de tétraèdres — 22 × 387 978 = **8,5 M de distances** sur `death.stl` @
    0,8 — et ce k-d tree est **déjà natif depuis le pas 5**. Or la fenêtre de
    Hann annule toute source au-delà du rayon : mesuré, seules **21,4 %** de
    ces paires sont dans le rayon et **69,2 % sont écartables par un simple
    test d'AABB élargie**. Transcrire la fonction telle quelle n'aurait donc
    rien donné — c'est mot pour mot la leçon du pas 3.

    L'élagage est **exact, pas heuristique** : un tétraèdre hors de la boîte
    élargie est à plus de `radius` sur au moins un axe, donc sa fenêtre est
    nulle et son poids aussi. Le seul usage restant d'une distance hors
    rayon est le repli `nearest_source`, emprunté par les tétraèdres
    qu'aucune source n'atteint : ceux-là sont repris dans une **seconde
    passe sans élagage** qui retrouve l'argmin exact. Un test dédié force ce
    chemin avec un rayon de 1e-3, où presque tous les tétraèdres y tombent.

    #### LE RÉSULTAT NÉGATIF : Eigen ne rapporte rien

    L'utilisateur avait choisi « tout porter y compris Eigen, accepter la
    perte de bit-exactitude » sur la foi d'une estimation de ~0,5-0,8 s de
    gain, tirée du surcoût Python visible au profil (0,39 s dans la boucle
    `cg` de scipy, plus le dispatch `LinearOperator`).

    **La mesure dit ×1,01** : 1,833 s (scipy) contre 1,812 s (Eigen), même
    algorithme, même préconditionneur diagonal, même critère d'arrêt, 978
    itérations, résidus 9,4e-9 et 9,7e-9. La raison est simple une fois
    mesurée : le coût est **dominé par le produit matrice-vecteur creux**
    (1,42 s sur 1,90 s), et `csr_matvec` de scipy est déjà du C aussi rapide
    que celui d'Eigen. Il n'y avait pas de surcoût Python à retirer.

    **Conséquence actée : le solveur Eigen est construit, testé et mesuré,
    mais N'EST PAS CÂBLÉ.** `_solve_linear_system` continue d'appeler scipy.
    Le câbler coûterait la bit-exactitude — donc la garantie que la métrique
    step9b ne bouge pas, invariant tenu sur les pas 1 à 6 — pour un gain
    nul. Même discipline que partout ailleurs dans ce fichier : un mécanisme
    qui ne démontre pas de gain reste disponible, jamais imposé.

    > **Leçon** : le profil disait « 0,39 s de Python dans `cg` ». J'en ai
    > déduit un gain, et c'était un raisonnement, pas une mesure. Le temps
    > Python d'une boucle n'est récupérable que s'il n'est pas déjà masqué
    > par le calcul natif qu'elle pilote.

    #### Équivalence : la tolérance remplace le bit, et c'est déclaré

    Premier portage du projet dont l'équivalence n'est pas bit-à-bit. Deux
    causes, mesurées séparément :
    - `compute_shape_gradients` inverse une 3×3 **par la formule des
      cofacteurs**, là où numpy appelle LAPACK (LU à pivot partiel) : même
      inverse, autres arrondis. Écart **relatif** 7,2e-14 — la précision
      machine. (L'écart absolu de 6,3e-9 n'était pas interprétable : les
      gradients montent à 86 703 sur maillage fin.)
    - le solveur Eigen (non câblé) : ~2,6e-9 relatif.

    Le checksum du champ passe de `1257825.311560218` à
    `1257825.311525613` — 8ᵉ chiffre significatif, cohérent avec 1e-14
    propagé à travers un CG résolu à 1e-8. Les tests d'équivalence
    **déclarent** donc une tolérance, chacune justifiée par sa mesure, au
    lieu d'un `array_equal`.

    #### Périmètre

    Le **chemin par défaut** (celui que `main.py` exerce). Les mécanismes
    opt-in — directions géodésiques (#52), rétroaction de surplomb (#53),
    lissage vectoriel itératif (#55), modulation du lissage, diffusion des
    directions — restent **intégralement en Python**, désactivés par défaut,
    donc jamais exécutés : les porter figerait un contrat pour du code que
    personne n'appelle.

    #### Seconde passe : assemblage CSR direct, et un TROISIÈME calcul redondant

    **4,66 s → 4,04 s.** Deux corrections, encore une fois moitié
    algorithmiques :

    **1. Assemblage direct en CSR (×2,04 sur les deux matrices).** Le chemin
    scipy assemble par `coo_matrix((data,(rows,cols))).tocsr()` : 16 entrées
    par tétraèdre, **6,2 M** sur `death.stl` @ 0,8, à trier puis dédupliquer
    — et il le refait **intégralement pour chacune des deux matrices**,
    alors qu'elles ont *exactement* le même motif de sparsité (celui de
    l'adjacence nœud-nœud, qui ne dépend que de la topologie, jamais des
    poids). Le motif est désormais construit **une fois** (`_FemPattern`) et
    chaque matrice n'a plus qu'à accumuler dedans par recherche
    dichotomique. 0,971 s → 0,476 s, coefficients à 3,9e-16 relatif, second
    membre **exact**.

    Le motif est **passé explicitement** de `build_static_system` aux deux
    assemblages, jamais mémorisé sur le `TetMesh` : un cache invisible
    survivrait à une modification de `tetrahedra` et produirait une matrice
    silencieusement fausse.

    **2. `tet_face_topology` calculé DEUX FOIS** (0,2 s chacun) :
    `_map_boundary_faces_to_surface` en gardait les faces frontière et
    jetait les paires internes, `_build_tet_tet_adjacency` faisait
    l'inverse — même passe, chacun jetant la moitié du résultat. La première
    rend maintenant les paires en sous-produit ; la seconde les accepte.
    Troisième calcul redondant trouvé dans step4 en deux passes de profil.

    **Limite assumée** : le motif natif du terme de direction est un
    **sur-ensemble** (il couvre tous les tétraèdres, là où scipy ne crée
    d'entrées que pour les actifs). Les entrées surnuméraires valent
    exactement zéro et ne coûtent rien en pratique, `build_static_system`
    additionnant ce terme à la rigidité dont le motif est déjà complet. Le
    test d'équivalence vérifie donc l'égalité **mathématique**, pas celle
    des motifs — et le documente plutôt que de la masquer.

    #### Bilan step4 et ce qui reste

    | étape | temps | nature du gain |
    |---|---:|---|
    | départ | 9,45 s | — |
    | fix du solve redondant | 7,47 s | **Python**, bit-exact |
    | blend + gradients de forme | 4,60 s | C++ (blend : algorithmique) |
    | assemblage CSR + topologie partagée | **4,04 s** | C++ + algorithmique |

    **×2,34 au total, dont une part importante sans C++.** Trois calculs
    redondants trouvés en profilant — le solve, la topologie, le motif de
    sparsité —, tous invisibles à la lecture du code.

    **Reste** : le CG à **1,91 s, soit 48 % de step4**, et il est hors de
    portée — Eigen y est mesuré à ×1,01 (voir ci-dessus). Descendre plus bas
    exigerait de changer le PRÉCONDITIONNEUR (le Jacobi actuel laisse 978
    itérations) ou la formulation, c'est-à-dire un changement algorithmique à
    mesurer sur la métrique step9b, pas un portage.
    Mesuré le 2026-10-01 sur le système de death @ 0,8 (87 913 inconnues,
    1,1 M non-zéros) : CG Jacobi **1,25 s** (959 itérations) ; SuperLU direct
    **110 s** ; CG + ILU de scipy (drop 1e-4) ne converge pas en 10 000
    itérations (96 s). Le seul candidat restant est le Cholesky incomplet
    d'Eigen (déjà vendoré), pour ~1 s à gagner sur ~12 : non tenté.

    **Eigen 3.4.0** est vendoré dans `cpp/third_party/eigen/` (header-only,
    MPL-2.0 vérifiée à la source) et intégré au build.

    ### Le pas 8 : SAT triangle-triangle (`step1_ingest.py`)

    Après le pas 7, la question posée était « on porte tout step1 ». La
    lecture du fichier (2144 lignes) a montré que ce n'est PAS un bon
    candidat au sens littéral : ~90 % du fichier appelle directement
    trimesh/pymeshfix (`slice_plane`, opérations booléennes, triangulation
    de polygones, `MeshFix.repair`) — de la logique de réparation de
    maillage sur des bibliothèques tierces, ni portable à coût raisonnable
    ni souhaitable à porter. `_build_tet_adjacency` (l'adjacence
    nœud-nœud du TetMesh) est déjà natif depuis le pas 4
    (`tet_unique_edges`). Le seul point chaud réel restant, identifié par
    lecture et non par profil (le chemin nominal `death.stl` ne le
    déclenche jamais), est le test SAT triangle-triangle
    (`_triangles_cross`) et la boucle qui l'appelle sur toutes les paires
    candidates (`_find_self_intersecting_faces`) — utilisés uniquement sur
    le chemin de réparation d'auto-intersections (déclenché sur
    `turkish.stl` : 50 faces auto-intersectantes détectées par ce test).

    **Portage** : `cpp/include/apex/mesh_repair.hpp` +
    `cpp/src/core/mesh_repair.cpp`. Deux fonctions exposées :
    - `triangles_cross(t0, t1)` — transcription bit-exacte du théorème de
      l'axe séparateur (11 axes : 2 normales de face + 9 produits croisés
      arête×arête), même tolérance `_COPLANAR_EPS = 1e-9`, mêmes axes de
      norme `< 1e-12` ignorés (arêtes parallèles ou face dégénérée).
      Bit-exact par construction : mêmes opérations flottantes, même ordre.
    - `find_self_intersecting_faces(triangles, faces)` — broad-phase +
      narrow-phase complet, transcription de `_find_self_intersecting_faces`
      (exclusion des paires de triangles partageant un sommet, chaque paire
      `(i,j)` avec `j>i` visitée une seule fois).

    **Le broad-phase natif N'A PAS besoin d'être identique à `bounds_tree`
    (rtree) côté Python**, contrairement au BVH de rayons du pas 6. Un BVH
    maison sur AABB exactes (même structure que `TriangleBVH`, réutilisée
    en interne) suffit : le narrow-phase (SAT) est déterministe et
    bit-exact des deux côtés, donc tout sur-ensemble conservateur de
    candidats donne le même résultat final — contrairement au cas des
    rayons, où le broad-phase de trimesh s'est révélé PLUS STRICT que sa
    propre narrow-phase et faisait donc partie du comportement observable
    (voir pas 6). Ici, aucune structure d'index n'est un comportement
    observable : seul le test SAT compte.

    `_triangles_cross` (dans `step1_ingest.py`) dispatche vers
    `apex_core.triangles_cross` quand `APEX_USE_CPP=1`, bénéficiant du même
    coup `_fill_loop_min_weight` (triangulation Klincsek du rebouchage, qui
    appelle `_triangles_cross` en boucle interne) — cette dernière fonction
    elle-même N'EST PAS portée (état Python riche : cache mémoire, DP
    O(n³)), seule la primitive SAT qu'elle appelle l'est.

    **6 tests d'équivalence** (pas 8) : paire qui se croise réellement,
    contact simple sans intersection (bord commun coplanaire, doit rester
    `False`), balayage aléatoire de 400 paires avec garde anti-faux-vert
    (au moins une paire croisée doit être tirée), maillage synthétique avec
    une vraie auto-intersection CSG, maillage propre (ensemble vide des
    deux côtés), et non-régression du branchement dans
    `step1_ingest._find_self_intersecting_faces` elle-même.

    **357/357 tests dans les deux positions de l'interrupteur** (351 + 6) —
    192,57 s (Python) contre 184,09 s (C++) sur la suite complète. Gain non
    mesuré sur le pipeline `death.stl` (chemin nominal, ce code n'y est
    jamais exercé) — conforme à l'attente : ce portage vise `turkish.stl`
    et tout modèle qui déclenche la cascade de réparation, pas le cas
    déjà propre.

    ### Le pas 9 : segments d'intersection + winding number (`_self_intersection_resolve.py`)

    Reprofilage du pipeline complet (`main.run_pipeline`, step1→step8,
    `death.stl`, résolution 0.8, C++ actif pas 1-8) après le pas 8 :
    **40,4 s au total**, et `_finish_geometric_resolution` /
    `resolve_self_intersections` en 2ᵉ position (~13,5 s) — le plus gros
    poste Python restant après TetGen (incompressible, isolé par design).
    Le fichier source (953 lignes) est majoritairement du calcul numpy
    maison (102 refs `np.` contre 9 `trimesh.`) — profil **inverse** de
    `step1_ingest.py` (pas 8), donc un bon candidat contrairement à ce
    dernier.

    **Portées** (`cpp/include/apex/mesh_intersect.hpp` +
    `cpp/src/core/mesh_intersect.cpp`) : `triangle_triangle_segment`
    (15856 appels mesurés, 2,59 s cumulées), `generalized_winding_number`
    (formule fermée, 0,92 s), et le broad+narrow-phase de
    `_intersection_segments` (3,22 s). `_subdivide_face`,
    `_split_keyhole_polygon`, `_weld_line_endpoints`,
    `_close_small_boundary_loops` **NE SONT PAS portées** : dépendance
    forte à shapely (`polygonize`, `unary_union`, `Polygon`) et networkx,
    beaucoup d'état/cas particuliers — même diagnostic que la partie non
    portable de `step1_ingest.py`.

    **`Triangle3` déplacé dans `types.hpp`** (partagé entre
    `mesh_repair.hpp` et `mesh_intersect.hpp`), et le BVH interne sur AABB
    de triangles factorisé dans `apex/internal/triangle_aabb_bvh.hpp`
    (partagé entre `mesh_repair.cpp` et `mesh_intersect.cpp`, structure
    purement interne, aucun contrat d'équivalence).

    #### Équivalence : TOLÉRANCE EXPLICITE — deuxième portage du projet après step4

    `triangle_triangle_segment` résout un système 3×3
    (`np.linalg.solve`, LAPACK) pour localiser le point de référence sur la
    droite d'intersection des deux plans porteurs. Le pas 7 avait déjà
    mesuré qu'une inversion 3×3 par cofacteurs (la seule alternative sans
    nouvelle dépendance) diverge de LAPACK à ~1e-14 relatif — jamais au
    bit. **Décision explicite de l'utilisateur** (posée avant d'écrire le
    code, par analogie avec le pas 7) : cofacteurs + tolérance documentée,
    plutôt que de renoncer à ce portage ou de lier une dépendance
    supplémentaire pour un seul solve 3×3. Mesuré : écart max **9,77e-15**
    sur 562 paires de triangles réellement croisées (2000 tirées) — trois
    ordres de grandeur sous les tolérances propres du module (`_EPS=1e-10`,
    snap à `1e-9 * scale`).

    `generalized_winding_number` est une formule fermée sans solve, donc
    **quasi**-bit-exacte : écart mesuré 1,1e-15 sur 300 points × 1280
    triangles — un ULP de bruit dû à la réduction `np.sum` de numpy
    (potentiellement pairwise), pas à un solve. Documenté comme tel plutôt
    que faussement annoncé bit-exact (erreur de rédaction initiale du
    header, corrigée après la première mesure empirique).

    #### LA VRAIE LIMITE, MESURÉE ET ASSUMÉE : l'ordre des segments par face n'est pas reproductible

    Contrairement au SAT du pas 8 (qui ne rend qu'un booléen, insensible à
    l'ordre des candidats), ici **l'ordre compte** : `_subdivide_face`
    (non portée) itère `per_face[i]` dans l'ordre reçu, et cet ordre
    influence quels points de la chaîne de découpe sont vus en premier,
    donc la triangulation exacte produite.

    Le Python reçoit cet ordre de `tree.intersection(query)`
    (rtree/libspatialindex). **Vérifié, pas supposé** : cet ordre n'est
    **jamais** trié par indice de triangle (mesuré sur les 10180 requêtes
    de `turkish.stl` — 100% non triées) — c'est un ordre spatial interne à
    la structure R-tree, non reproductible sans réimplémenter
    libspatialindex. Un tri par indice croissant côté C++ a été tenté (par
    hypothèse erronée que rtree triait par indice) : ça a réduit mais pas
    éliminé la divergence — conservé malgré tout (rapproche sans coût),
    mais documenté pour ce qu'il est : une approximation, pas une
    transcription.

    **Mesure bout-en-bout** sur un motif synthétique reproduisant
    `turkish.stl` (deux nervures qui se traversent) — seul appelant réel de
    `_intersection_segments`, via `resolve_self_intersections` :

    | | Python (rtree) | C++ (BVH maison) |
    |---|---:|---:|
    | vitesse | 29,98 s | **2,60 s (×11,5)** |
    | faces en sortie | 10270 | 10264 |
    | volume | 49270,855 | 49270,298 (**1,13e-5 relatif**) |
    | `euler_number` | −34 | **−35** |

    L'écart d'`euler_number` est une **vraie divergence topologique
    mineure** (pas seulement du bruit flottant) : signature d'une boucle de
    bord fermée différemment par `_close_small_boundary_loops` en aval (le
    fichier documente lui-même ce mécanisme comme sensible à l'ordre —
    "convergence en 2 passes", "3 boucles de 3 arêtes").

    **Décision explicite de l'utilisateur, posée via AskUserQuestion avant
    d'agir** (trois options : accepter et documenter / ne porter que les
    deux primitives bit-exactes / chercher à reproduire l'ordre rtree) :
    garder le portage complet. Le gain (×11,5, sur le plus gros poste
    Python restant après le pas 8) l'emporte sur une divergence mesurée
    comme mineure, sur un mécanisme déjà **best-effort** par construction
    (`resolve_self_intersections` rend `None` dans de nombreux cas selon sa
    propre docstring — jamais censé produire une réponse canonique unique).
    Chercher à reproduire l'ordre exact de libspatialindex a été explicitement
    écarté comme hors de portée raisonnable pour un écart déjà mineur.

    #### Vérification

    5 tests d'équivalence (pas 9) : balayage aléatoire de
    `triangle_triangle_segment` (500 paires, garde anti-faux-vert sur le
    nombre de paires croisées, tolérance 1e-9), cas limites explicites
    (triangle dégénéré, plans parallèles), `generalized_winding_number`
    (tolérance 1e-10, garde sur un point clairement intérieur ET un
    clairement extérieur), `find_intersection_segments` (mêmes faces
    affectées, même nombre de segments par face — pas le même ordre),
    et `resolve_self_intersections` bout-en-bout (tolérance de volume
    relative 1e-3, `is_watertight` identique).

    **362/362 tests dans les deux positions de l'interrupteur** (357 + 5) —
    195,02 s (Python) contre 168,82 s (C++) sur la suite complète.

    ### Le pas 10 : interpolation inverse-distance (`step6_toolpath2d.py`)

    **DERNIER portage guidé par le profil.** Reprofilage après le pas 9
    (step1→step8, `death.stl`, résolution 0.8, C++ actif) : **43,1 s
    profilées dont 22,4 s d'attente TetGen (52 %)**, et hors TetGen **plus
    aucun poste Python isolé ne dépasse 2,1 s**.
    `_interpolate_out_of_plane` était ce dernier poste (2,10 s de temps
    propre, 2853 appels) — pondération inverse-distance `1/d^power` de
    l'offset hors-plan, dont le coût venait de la matrice de distances
    `(K,M)` que numpy matérialise en entier.

    Porté dans `cpp/include/apex/interpolate.hpp` +
    `cpp/src/core/interpolate.cpp` : une seule passe par requête (plus
    proche voisin, détection de coïncidence et les deux accumulations
    comprises), sans jamais matérialiser la matrice. **×2,8 mesuré** sur
    une taille réaliste (400 requêtes × 300 sources, 300 répétitions :
    1,64 s → 0,58 s).

    **Équivalence : tolérance serrée, pas le bit — mesuré, pas supposé.**
    Ni solve ni BLAS ici (contrairement aux pas 7 et 9), la transcription
    est terme à terme ; mais les deux réductions (`weights.sum(axis=1)`,
    `(weights * offsets).sum(axis=1)`) sont accumulées séquentiellement là
    où numpy peut sommer par paires. Mesure sur 200 tirages aléatoires
    (tailles et coïncidences variées) : écart max **8,88e-16**, et
    **36/200 seulement bit-identiques** — même phénomène qu'au pas 9. Un
    ULP de bruit, sans conséquence sur un offset en millimètres dont les
    seuils d'usage vont de 1e-9 (`_COINCIDENT_EPS`) à 1e-3.

    Les trois cas limites du Python sont transcrits tels quels car ils ne
    sont pas cosmétiques : coïncidence (`0**-power` exploserait → valeur
    source exacte), `total <= 0` (sous-débordement sur de très grandes
    distances → repli sur le plus proche), et `argmin` rendant le **premier**
    minimum en cas d'égalité (comparaison stricte).

    3 tests d'équivalence (pas 10), dont une garde anti-faux-vert vérifiant
    qu'au moins un cas coïncident a bien été tiré.

    **365/365 tests dans les deux positions de l'interrupteur** (362 + 3) —
    185,71 s (Python) contre 171,35 s (C++) sur la suite complète.

    ### Pas 11 à 14 (session autonome du 2026-10-01) — portage systématique

    **Pas 11 : contrat de paramètres.** `config/schema.py` dérive des
    annotations de type l'aller-retour JSON et le schéma JSON de TOUS les
    paramètres. Aucune liste n'est tenue à la main : un champ ajouté apparaît
    partout, et `test_every_field_in_schema` le vérifie. `config/params.py`
    réexporte les dataclasses sans les déplacer (aucun import cassé).
    `scripts/gen_params.py` écrit `config/params.schema.json` (pour l'UI
    TypeScript) et `cpp/include/apex/params.hpp` (structures C++, valeurs par
    défaut comprises, y compris celles propres à `PipelineParams` comme
    `Step6Params(out_of_plane_interpolation_power=2.0)`).
    `test_generated_files_up_to_date` échoue si l'on oublie de relancer le
    générateur. La sonde native `params_defaults_probe` prouve que l'en-tête
    compile et rend les mêmes défauts que Python. Champs `typing.Any`
    (`MachineProfile.nozzle_model`) : non sérialisés, `x-runtime-only`. Un
    champ JSON inconnu est une ERREUR (faute de frappe côté UI jamais ignorée
    en silence). Deux annotations précisées au passage
    (`LayerParams.support_rotation_axes: list[RotationAxisSpec]`,
    `coupled_reorientation_tiers: tuple[float, ...]`).

    **Pas 12** : `point_in_polygon_2d` natif (`apex/geometry2d.hpp`), lancer de
    rayon transcrit terme à terme, au bit près (0,52 s de Python sur death).
    Le reste de `geometry_utils` n'a pas d'appelant chaud mesuré.

    **step6 en parallèle par couche** (`_STEP6_THREADS`, défaut min(8,
    cœurs)) : les couches sont indépendantes, `map` conserve l'ordre, résultat
    identique au bit (`scripts/check_step_parallel.py`). Gain limité par le
    GIL : 2,89 → 1,77 s sur death. Step8 : normale de couche mise en cache
    (29 691 produits vectoriels pour 257 valeurs distinctes).

    **Bilan death @ 0,8 (profilé, même poste)** : 105 s en début de session →
    **18 s**. TetGen 7,2 s (plancher, désormais seul gros poste), step8 3,6 s,
    step4 2,7 s, step6 2,6 s, step5 1,2 s.

    **Pas 13 : marching tetrahedra + soudure** (`apex/isosurface.hpp`,
    `IsosurfaceExtractor`, construit une fois par champ). **Au bit près**, et
    c'est atteignable ici (pas de réduction, pas de solve) : même formule
    d'interpolation, orientation au signe d'un produit mixte, même ORDRE de
    sortie que la version numpy par lots (tous les « caps » puis tous les
    « quads », 2 triangles consécutifs), soudure par `nearbyint` (arrondi au
    pair comme `np.round`) et tri lexicographique STABLE (comme
    `np.unique(return_index)`).

    **Pas 14 : tracé de bord** (`boundary_loops`, remplace
    `trimesh.outline()`). Reproduit exactement outline() quand le bord est
    fait de boucles simples (chaque sommet de degré 2) : le DFS de scipy part
    du plus petit indice, suit son plus petit voisin et se referme, et les
    boucles sortent par indice de départ croissant (vérifié empiriquement).
    Sinon, repli sur trimesh POUR CETTE COUCHE seulement
    (`_NATIVE_OUTLINE_FALLBACKS`).

    Mesure sur modèles réels (`scripts/check_step5_native.py`, même champ,
    step5 rejoué dans les deux modes) : **contours identiques au bit**,
    **0 repli**. step5 : l_shape 0,14 → 0,02 s ; death @ 0,8 **2,56 → 1,14 s**.
    Suite complète : 421/421 dans les deux positions (avant le pas 14).

    **Pas 17 : IK couplée à la position (step8, machine actuelle).** Le
    profil de death @ 0,8 après le pas 14 était sans appel : **step8 = 84 s sur
    105 s**, entièrement dans `_resolve_all_axis_angles_coupled`, écrit en
    Python pur au passage à la cinématique couplée (limite #59) — l'IK native
    du pas 2 ne servait plus du tout. Réorganisé en trois temps
    (`_resolve_all_axis_angles_coupled_native`) : C de tous les points d'un
    coup (sa seule dépendance séquentielle est l'héritage sur l'axe, un
    remplissage vers l'avant) ; cibles ramenées dans le repère de C par
    Rodrigues par lots ; axes internes résolus par l'IK NATIVE du pas 2, sans
    détection de singularité (la boucle Python n'en fait pas pour les axes
    internes), avec la même reprise au point fautif. Aucune ligne de C++
    nouvelle. Équivalence : écart max 2e-14 degré sur 20 000 points (93 %
    bit-identiques ; produits scalaires vectorisés à un ULP près), mêmes
    diagnostics. ×60 sur la fonction ; **pipeline death 105 s → 25,5 s**
    (profilé), step8 84 → 4,3 s. TetGen (12,7 s) redevient la moitié du temps.

    **TetGen en parallèle de lui-même (step1, sans C++).** Le « plancher dur »
    de TetGen cachait deux tétraédrisations complètes EN SÉRIE sur tout modèle
    qui passe par la résolution des auto-intersections (death.stl) : la sonde
    de `_finish_geometric_resolution` (qualité, SANS raffinement plateau,
    6,55 s, résultat jeté) puis la tentative réelle (AVEC raffinement, 6,85 s).
    Ce ne sont pas les mêmes calculs, donc pas de réutilisation possible. Mais
    quand la sonde réussit, cas nominal, elle rend le maillage INTACT, et la
    tentative réelle porte sur ce même objet. Elle est donc lancée EN PARALLÈLE
    de la sonde (deux sous-processus), et son issue n'est retenue que si la
    sonde n'a rien modifié : résultat identique dans tous les cas (vérifié au
    bit sur death, 388 775 tétraèdres ; logique verrouillée par
    `test_speculative_tetgen_uses_parallel_attempt_only_when_probe_left_mesh_intact`).
    step1 sur death : **15,6 s → 8,1-8,6 s**. Interrupteur de module
    `step1_ingest._SPECULATIVE_TETGEN`. Désactivé quand des zones de jonction
    sont raffinées (non mesuré dans ce cas).

    Second doublon, trouvé sur turkish.stl (7 appels TetGen) : quand la sonde
    réussit APRÈS une réparation, la cascade essaie le maillage réparé avec
    raffinement (échec), puis « sans raffinement plateau » — c'est-à-dire
    exactement le calcul de la sonde (même maillage, même résolution ; marge
    et `platform_z` ignorés sans raffinement). `_finish_geometric_resolution`
    rend désormais le TetMesh de sa sonde, réutilisé tel quel (7,26 s).
    Comparaison au bit contre la version git (`step1` complet) : identique sur
    death (388 775 tétraèdres) et turkish (452 144) ; step1 **13,8 → 7,7 s**
    (death) et **20,2 → 12,8 s** (turkish).

    **Pas 19 : contrôle de collision rapide** (`apex/past_material.hpp`,
    `PastMaterialIndex`). L'algorithme Python de `FastCollisionChecker` (K plus
    proches voisins scipy, puis boule exacte si la liste est saturée de
    matière future) calcule EXACTEMENT la distance minimale à la matière déjà
    déposée dans le rayon : un point passé hors des K plus proches est plus
    loin qu'eux, et le cas saturé est rattrapé. Le natif calcule ce même
    minimum avec un arbre k-d dont chaque nœud connaît sa boîte et le PLUS
    PETIT indice de dépôt de son sous-arbre : un sous-arbre entièrement futur
    est élagué d'un coup. Identique au bit (300 requêtes aléatoires contre
    l'algorithme Python), mêmes événements sur death (24 789). step9 rapide :
    **death 233 → 13,7 s (×17)**, turkish ~36 s pour 141 602 points.

    **Pas 21a : étape 2 entière en natif** (`apex/overhang.hpp`,
    `compute_overhang_map`) : angles, score de risque, ponts par rayons
    verticaux (BVH du pas 6), composantes connexes (union-find sur le graphe de
    TOUTES les faces, comme le Python) et axe principal (vecteur propre de la
    covariance 3×3, Eigen). Gain de temps nul (0,02 s) : pièce de
    l'orchestrateur natif. Sur les 4 modèles : surplombs et ponts identiques,
    risque **au bit**, portée à 2e-16. Le signe de la portée est désormais
    CANONIQUE des deux côtés (plus grande composante positive), sans effet
    physique (step6 balaie perpendiculairement). Cas dégénéré assumé : un
    plafond carré a deux variances égales, tout axe du plan est principal, et
    LAPACK et Eigen en choisissent des différents ; le test vérifie donc la
    variance le long des deux axes, pas leur égalité.

    **Pas 21b : étape 3, chemin par défaut, en natif** (`apex/segmentation.hpp`,
    `segment_default`) : compatibilité face/direction, croissance de régions,
    direction par région, profondeurs (BFS multi-sources par composante),
    graphe de dépendances. **IDENTIQUE au Python sur les 4 modèles**, avec et
    sans critère d'angle de fusion, ORDRE DES DÉPENDANCES COMPRIS (transcrit de
    networkx : nœuds croissants, voisins par première apparition). Les options
    `branch_cut_threshold` et `anchor_growth_weight` restent en Python.
    **Découverte préalable** : `_region_growing` triait les arêtes par
    `np.argsort` NON stable, alors que chaque arête figure deux fois au même
    risque. Le résultat dépendait donc de l'implémentation de numpy : sur
    turkish, le tri stable fait passer de 11 à 24 régions. **Tri stable
    adopté** (mesuré : death et les 2 modèles fins inchangés, turkish 0,198 →
    0,172 % de matière sans appui et collisions 0,21 → 0,13 %).

    **Pas 21c : étape 7 entière en natif** (`apex/scheduling.hpp`,
    `schedule_layers`) : couleur par arête (distance à la surface par
    `MeshDistance`, qui rend désormais aussi l'indice du triangle le plus
    proche), lissage par longueur d'arc, découpage « une polyligne = une
    couleur », groupes par couleur, deadline des trous (BVH du pas 6 + arbre
    k-d du champ), arbitrage trou / purge et ordonnancement EDF. **Identique
    au Python** sur 5 scénarios (sans peinture, bicolore, latérale, avec lissage,
    trous gratuits), polylignes, couleurs et passes comprises. Tolérance
    déclarée aux seules frontières : point à égale distance de deux faces de
    couleurs différentes, somme de longueurs comparée à `merge_tolerance`.
    `build_color_dependency_graph`, dont le Python ne lit pas le résultat,
    n'est pas porté.

    **step9 rapide : positionnement de la tête par lots**
    (`_positioned_toolhead_vertices_batch`, matrices de rotation composées dans
    le même ordre, produit matriciel par lots de 1 000 poses). Pose par pose,
    Rodrigues prenait 14,5 s sur 20. death : **step9 20 → 6,3 s**, mêmes 88
    événements (écart de quelques ULP, sans effet sur un contrôle à 0,125 mm).
    Puis : requête de matière passée PAR LOTS et en parallèle
    (`min_past_distance_batch`, identique au bit à l'appel unitaire, testé),
    contrôle du plateau par lots avec raccourci exact (une pose loin en z du
    plateau ne peut pas le toucher), échantillonnage et densification
    vectorisés. **step9 rapide : 20 → 2,9 s** sur death, toujours 88 événements.
    Même technique dans step5b (`_rotation_matrices`) : 6 → 3 s (death),
    12 → 7,4 s (turkish), contraintes identiques.

    **Pas 24 : export G-code** (`apex/gcode_export.hpp`). Le corps du fichier
    est produit en natif, **texte identique ligne à ligne** : 955 226 lignes
    du Benchy, et tous les champs de `GcodeReport` égaux. Trois transcriptions
    obligatoires : le modulo flottant de Python (signe du diviseur, pas
    `fmod`), `f"{v:.4f}"` par `std::to_chars` (arrondi exact comme le dtoa de
    Python), et la norme à QUATRE composantes (X, Y, Z, B) de l'avance.
    `feed_plan` : passes avant/arrière en natif, au bit, à condition d'écrire
    `pow(x, 2.0)` comme numpy (`prev ** 2`) et non `x * x` (27 facteurs sur
    943 000 différaient d'un ULP sinon). En-tête, pied et `travel_check`
    (rappel Python) restent en Python ; `head_zlow` passe en TABLE
    (`head_zlow_from_table`). Benchy : **export 20,9 s → 0,8 s**.

    **Pas 15 : remplissage 2D d'une zone (step6), Clipper2** (`apex/toolpath2d.hpp`,
    `fill_region`). Parois par offset en onglet (limite 5 comme GEOS), puis
    remplissage LINES / GRID / CROSSED_ZIGZAG (balayage pair-impair, raccords qui
    suivent le bord, à défaut raccord droit couvert par la zone). Les ponts
    restent en Python. **Pas au bit, et c'est assumé** : Clipper2 calcule en
    entiers (1e-6 mm ici) et part d'un autre sommet d'anneau. Écart par couche
    ≤ 1e-4 sur l_shape et l-triple (`scripts/check_native_fill.py --continuous`).
    Sur death, 12 couches s'écartent de plus de 1 % : sur une fente fine du
    contour, l'onglet de GEOS s'allonge en pointe de 4 mm qui ronge la paroi, là
    où Clipper2 coupe en carré à la distance d'offset. Le natif est plus proche
    de la vraie géométrie, mais différent. step9b à parois continues :
    l_shape −0,8 %, l-triple −1,1 %, death +0,7 %, turkish −2,6 %. step6 sur
    death **1,81 → 0,51 s**. **Opt-in** (`APEX_NATIVE_FILL=1`, en plus de
    `APEX_USE_CPP=1`) : à basculer par défaut avec le passage au tout-natif.
    C'est en le portant qu'a été trouvé le défaut des parois trouées (limite #65).
    Remesuré avec les défauts actuels (ordre contraint par la tête, parois
    continues) : death +0,8 %, turkish −3 %, l-triple −1,1 %, l_shape −0,8 % de
    matière sans appui, collisions identiques. Reste opt-in par principe (le
    chemin natif ne doit pas déplacer la métrique) jusqu'à la bascule.

    **step8 en tableaux (sans C++, au bit).** `_iter_schedule_points` rend un
    `_RawPoints` (une colonne par champ, assemblé par blocs d'une polyligne)
    au lieu de ~400 000 objets `_RawPoint` ; l'IK couplée concatène ses blocs
    au lieu de convertir chaque ligne en liste ; `_lift_to_3d` calcule sa
    normale par un produit vectoriel scalaire au lieu de `np.cross` (même
    formule, même arrondi, 0,8 s d'appels évités). Trajectoire et diagnostics
    **identiques au bit** contre l'ancienne version
    (`scripts/check_step8_arrays.py`, avec et sans médian / normale par point).
    step8 sur death : **1,92 → 0,82 s**. L'accès point par point (`raw[i]`)
    reste disponible pour les replis et les tests.

    > Leçon : un changement de COMPORTEMENT (ici la cinématique couplée) peut
    > faire sortir silencieusement un chemin du natif. Reprofiler après chaque
    > changement de comportement, pas seulement après chaque portage.

    Outillage : `scripts/build_cpp.py --dev` construit dans
    `cpp/build_dev/lib` sans toucher `cpp/lib`, verrouillé sous Windows tant
    qu'un calcul long l'a chargé ; `APEX_CPP_LIB_DIR` fait charger ce
    module. `scripts/profile_pipeline.py` : profil par étape.

    ### CE QUI RESTE À PORTER — le profil ne guide plus

    Après le pas 10, le paysage mesuré sur `death.stl` @ 0.8 :

    | poste | temps propre | nature |
    |---|---:|---|
    | **TetGen** (`WaitForSingleObject`) | **22,40 s** | sous-processus, **plancher dur** |
    | `argsort` numpy | 1,31 s | déjà natif |
    | `csr_matvec` (CG scipy) | 1,21 s | déjà natif, Eigen mesuré ×1,01 (pas 7) |
    | `ufunc.reduce` | 0,95 s | déjà natif |
    | `apex_core.resolve_axis_angles_batch` | 0,86 s | **déjà C++** |
    | `step8::_smooth_all_axes` | 0,80 s | Python |
    | `core::point_in_polygon_2d` | 0,73 s | Python |
    | `np.linalg.norm` (327 k appels) | 0,72 s | dispatch numpy diffus |

    **TetGen fait à lui seul 52 % du temps** et ne bougera pas. Ce qui reste
    en Python est éparpillé en dizaines de postes à 0,1-0,8 s : le portage
    doit désormais être **systématique** (étapes entières : step5, step6,
    step7), plus guidé par le profil. C'est le même « point d'inflexion »
    que celui documenté après le pas 4, mais cette fois il n'y a plus de
    grappe à viser du tout.

    **Plafond, à dire maintenant plutôt qu'à l'arrivée** : même en portant
    parfaitement tout le Python restant, le plancher sur ce modèle à cette
    résolution est ~25 s (TetGen + algèbre creuse). Descendre plus bas
    exige de la **parallélisation** (plafond calculé ~18,5 s sur 8 cœurs,
    voir la section dédiée) ou un changement algorithmique — pas du
    portage.

    ### Le pas 7 (préparatifs, conservés pour mémoire) : `step4_field_solver.py`

    **Décision de l'utilisateur : porter step4 — assemblage FEM ET
    résolution — EN UN SEUL BLOC, avec Eigen dès le départ.** Pas d'étape
    intermédiaire qui assemblerait en C++ puis repasserait par scipy pour
    résoudre. Choix délibérément à contre-courant du reste du plan (les pas
    1-6 visaient chacun le plus petit incrément mesurable) : c'est le
    morceau le plus dur du pipeline restant — 2522 lignes, 38 fonctions,
    plusieurs mécanismes opt-in (rétroaction de surplomb #53, directions
    géodésiques #52, lissage vectoriel itératif #55) — et l'utilisateur a
    choisi de l'attaquer directement plutôt que par ses bords faciles
    (`step2`, 318 lignes, aurait été le choix « petit pas »).

    **RENDEZ-VOUS EXPLICITE AVANT D'ÉCRIRE LE PREMIER `.cpp` DE STEP4** :
    l'utilisateur veut être notifié (via `AskUserQuestion`) au moment où
    l'implémentation réelle commence, pour pouvoir charger les gros modèles
    (`death_hi-res.stl`, `turkish.stl`) et un effort de raisonnement maximal
    à ce moment précis plutôt que pendant la préparation. **Ne pas
    commencer à écrire `step4.cpp`/`kinematics`-like code pour step4 sans
    reposer cette question d'abord**, même si tout ce qui suit ici est prêt.

    #### Préparatifs déjà faits (ceux-ci NE nécessitaient PAS ce rendez-vous)

    - **Eigen 3.4.0 vendoré** dans `cpp/third_party/eigen/` (header-only,
      dossier `Eigen/` + `COPYING.MPL2` + `COPYING.README` seulement — pas
      le dépôt complet : docs/tests/lapack/blas non copiés, 6,7 Mo). Licence
      **vérifiée à la source** (`COPYING.MPL2`), cohérente avec la pile déjà
      actée dans la limite ci-dessous. PAS ENCORE intégré à
      `cpp/CMakeLists.txt` — volontairement, pour ne pas figer de dépendance
      de build avant que le code qui l'utilise existe (même discipline que
      `types.hpp` réduit au minimum aux pas précédents).
    - **Fumée de compilation validée** : un programme minimal
      (`Eigen::SparseMatrix` + triplets + `SimplicialLDLT`) compile sous
      MSVC 19.44 avec `/std:c++20 /O2` et résout correctement un système 3×3
      de test — vérifié par exécution, pas supposé. Confirme que le
      vendoring est utilisable avant d'investir dans l'assemblage réel.
    - **Aucun contrat `types.hpp` figé pour step4** : le principe déjà
      écrit dans `cpp/include/apex/types.hpp` (« ne pas transcrire de
      dataclass avant d'avoir porté le code qui la consomme ») reste
      respecté. Les dataclasses à transcrire (probablement `TetMesh`,
      `MeshData` en lecture, `OverhangMap`, `RegionDependencyGraph`,
      `ScalarField`, `Constraint`, plus les diagnostics optionnels
      `UpwardGrowthDiagnostic`/`SupportFeedbackDiagnostic`/
      `SmoothnessModulationDiagnostic`) restent à arrêter précisément une
      fois l'implémentation réellement commencée.

    #### Ce qui reste à décider AU moment du rendez-vous, pas avant

    - Périmètre v1 exact : `build_static_system` + `solve_with_constraints`
      + `run` au minimum ; les mécanismes opt-in (limites #52/#53/#55)
      entrent-ils dans le premier portage ou restent-ils Python le temps
      d'un second passage ?
    - Solveur Eigen à utiliser : `SimplicialLDLT`/`SparseLU` en direct
      (équivalent de `spsolve`) et `ConjugateGradient` préconditionné
      (équivalent du CG+ILU actuel, `direct_solve_node_threshold` comme
      bascule) — à confirmer contre le comportement exact de
      `_solve_linear_system`.
    - Stratégie d'équivalence : contrairement aux pas 1-6, un solveur creux
      itératif (CG) ne garantit PAS un résultat bit-identique à scipy même à
      code équivalent (ordre de réduction flottante différent selon
      l'implémentation). L'oracle devra donc probablement basculer d'une
      comparaison bit-exacte à une tolérance numérique explicite pour cette
      seule étape — à trancher et documenter avant d'écrire les tests,
      pas après.

    ### LICENCES — DÉCISION PRISE : AGPL ASSUMÉE, TetGen CONSERVÉ

    **Décision de l'utilisateur, actée : le logiciel sera 100 % gratuit et
    open source, pour tous et pour toujours ; seule l'imprimante associée est
    vendue. Dans ce cadre, TetGen (AGPL v3) et pymeshfix (AGPL v3) sont
    CONSERVÉS.**

    **Clarification ajoutée après une confusion en session (à ne pas
    reproduire) : cette décision n'exclut PAS MIT/BSD/GPL-3.0 — au
    contraire, elle les inclut explicitement.** Le paragraphe plus bas
    (« Ce que le portage avait rendu visible ») mentionne une « pile
    permissive » (Eigen, nanoflann, Clipper2, Geogram, manifold3d) : c'est un
    inventaire de PRÉPARATIFS antérieurs à la décision ci-dessus, choisis
    parce qu'ils sont permissifs et donc sans risque quel que soit l'angle
    retenu ensuite — ce n'est PAS un plafond au-delà duquel une licence
    deviendrait inacceptable. GPL-3.0 et MIT sont toutes deux compatibles
    avec « gratuiciel pour tous et pour toujours » et avec AGPL (déjà
    acceptée pour TetGen/pymeshfix) : ne pas réintroduire un critère
    « permissif uniquement » qui n'a jamais été la décision actée ici.

    Le raisonnement, à ne pas réouvrir sans raison nouvelle :

    - **L'AGPL n'interdit pas l'usage commercial.** Confusion courante :
      « libre » y désigne la liberté, pas la gratuité. Vendre du matériel
      accompagné d'un logiciel AGPL est le schéma standard du secteur —
      **PrusaSlicer, Bambu Studio et OrcaSlicer sont tous sous AGPL-3.0**, et
      les deux premiers sont édités par des fabricants qui vendent des
      imprimantes. Le critère n'est pas « est-ce commercial » mais « le code
      source reste-t-il disponible sous la même licence ».
    - **AGPLv3 et GPLv3 sont compatibles** (chacune a un article 13 autorisant
      la combinaison). Le code du projet peut donc être GPL-3.0 ; l'ensemble
      distribué traîne alors la clause réseau de l'AGPL, **inerte pour un
      logiciel de bureau** (elle ne se déclenche que sur une interaction à
      travers un réseau).
    - **Licence : lier TetGen est permis (précisé le 2026-10-26, sur remarque de
      l'utilisateur).** Le logiciel étant sous GPL-3.0 et TetGen sous AGPL-3.0, la
      combinaison est autorisée (article 13 de chaque licence) ; l'ensemble distribué
      porte la clause réseau de l'AGPL, inerte pour un logiciel de bureau. L'argument
      « programmes distincts » (`pipeline/_tetrahedralize_worker.py` reçoit deux chemins
      de fichiers) est une précaution SUPPLÉMENTAIRE, terrain juridiquement contesté, et
      n'est PAS la raison de garder TetGen séparé. Pas un avis juridique.

    **Conséquences directes sur le plan de portage, toutes des simplifications :**

    - le spike fTetWild est **abandonné** ;
    - le découplage de `step5` est **abandonné** — et avec lui le risque
      d'erreur de surface de ~0,06 mm identifié ci-dessous ;
    - `add_hole` et `add_region`, les deux mécanismes mesurés dont dépend
      step1 (limites #7 et #36/#48), sont **préservés tels quels** ;
    - la géométrie reste **identique au bit près** ;
    - il n'existe donc **aucun endroit du portage qui déplace la métrique**.

    **RÈGLE D'ARCHITECTURE (reformulée le 2026-10-26) : TetGen reste dans un
    exécutable séparé (`apex_tetgen_worker`, sous-processus, chemins de fichiers) pour
    des raisons de ROBUSTESSE, pas de licence.**
    - TetGen plante nativement (segfault) de façon reproductible sur certaines entrées
      (limite #35). En sous-processus, le plantage est une erreur récupérable
      (`_TetgenCrashError`, replis de step1) ; lié dans `apex_core`, il tuerait le
      moteur et l'interface TypeScript avec lui.
    - Le sous-processus permet de lancer deux tétraédrisations en parallèle (step1 de
      death : 15,6 s → 8 s) et garde le mailleur remplaçable.
    - Il ne coûte rien en vitesse : les ~14 s de death sont du calcul TetGen.
    - Lier TetGen reste PERMIS par la licence. Ce serait un choix de robustesse à
      assumer (perte de l'isolation et du parallélisme), à décider explicitement.
    - Un worker 100 % C++ (sans Python) est possible : il remplacerait
      `_tetrahedralize_worker.py`, toujours en exécutable séparé. Il exige la SOURCE de
      TetGen (v1.6.0, `tetgen.cxx` 36 537 lignes + `tetgen.h` + `predicates.cxx`, AGPL-3.0). Le wheel
      n'embarque que le binaire, mais le TARBALL source du paquet PyPI `tetgen` 0.8.4 les contient
      (`src/tetgen.cxx`, `src/tetgen.h`, `src/predicates.cxx`, `src/tetgen-license`) : vérifié le
      2026-10-26. **FAIT le 2026-10-26** : vendoré dans `cpp/third_party/tetgen/` (+ `LICENSE.txt`,
      `VERSION.txt`), exécutable séparé `apex_tetgen_worker` (`cpp/tools/tetgen_worker.cpp`, cible CMake
      jamais liée à `apex_core`, sortie dans `cpp/lib/`). Format binaire en entrée/sortie ; cette variante
      de TetGen lève des `std::runtime_error` dont le texte est relayé tel quel (la cascade de step1 y
      cherche « self-intersection », « recoversubfaces », « split_segment »). `step1_ingest` :
      `_NativeTetgen` (même interface que `tetgen.TetGen`) injecté dans `_tetrahedralize_once`
      (`_tetgen_factory`), utilisé par `_tetrahedralize_once_isolated` quand l'exécutable existe ;
      `APEX_TETGEN_WORKER=0` rend le worker Python. Crash natif → `_TetgenCrashError`, erreur →
      `_TetgenNamedConflictError`, comme avant. **Maillage IDENTIQUE au binding** (nœuds et tétraèdres
      `array_equal`) sur l_shape (sans/avec raffinement plateau), death (318 691 tét.) et turkish
      (328 531 tét., cascade de réparation comprise) ; step1 death 13,3 → 11,3-12,0 s, turkish 20,3 →
      17,6-18,3 s. Test : `test_native_tetgen_worker_matches_binding`.
    - **Temps par phase (death, `pq1.2a3.375`, `APEX_TETGEN_LOG=<fichier>` les journalise)** :
      Delaunay 0,005 s, surface 0,003, récupération du bord 0,04, **raffinement 4,5-4,6 s,
      lissage du maillage 5,3-5,7 s**, amélioration 0,12, sortie 0,03 ; total ~10,2-10,5 s. Turkish
      identique à 10,6 s (bord 0,25 s). **Plus de 95 % du temps est dans le raffinement et le lissage** ;
      Delaunay initial et récupération du bord sont négligeables, donc paralléliser ces deux-là ne
      rapporterait rien. Pistes restantes (non mesurées sur step9b) : lissage de sommets par couleurs
      indépendantes (parallélisable), réglage de `-s`/`-O`, maillage intérieur plus grossier.

    ### Ce que le portage avait rendu visible (conservé pour mémoire)

    Vérifié dans le venv, pas supposé : **TetGen est AGPL v3**
    (`.venv/.../tetgen-0.8.4.dist-info/licenses/src/tetgen-license` ;
    licence commerciale payante disponible auprès du WIAS Berlin), et
    **pymeshfix également**. Le wrapper PyPI `tetgen` est MIT mais *embarque*
    TetGen sous AGPL — la licence du wrapper ne change rien à celle du
    moteur. **Cette contrainte existe déjà aujourd'hui**, indépendamment du
    portage ; celui-ci ne fait que forcer la décision, puisqu'il faut de
    toute façon choisir le tétraédriseur à lier.

    Pile permissive retenue pour la suite (licences vérifiées à la source) :
    Eigen (MPL-2.0), nanoflann (BSD-2), Clipper2 (BSL-1.0), Geogram (BSD-3),
    manifold3d (Apache-2.0). Remplacements : **pymeshfix → manifold3d**
    (déjà une dépendance du projet, auto-union pour résoudre les
    auto-intersections) ; **TetGen → à trancher par un spike mesuré** entre
    Geogram construit sans TetGen et fTetWild (MPL-2.0, mais tire GMP en
    LGPL-3 et peut lier TetGen via `igl::tetgen` — les deux à auditer).

    > **Le changement de tétraédriseur est le SEUL endroit du portage où
    > l'équivalence bit-à-bit est impossible** : un autre mailleur donne un
    > autre maillage, donc un autre champ, donc d'autres couches, donc une
    > autre métrique. Il doit être traité comme un changement
    > ALGORITHMIQUE, mesuré sur les quatre modèles, **isolé** du reste du
    > portage — sinon un bug de portage et un effet du mailleur seront
    > indiscernables.

    1. **Espacement du remplissage dérivé de la taille de la couche.**
       `_generate_infill` calculait `diagonal / density * 0,05` au lieu de
       `nozzle_width / density` : une couche de 35 mm à 20 % recevait une
       ligne tous les 8,75 mm au lieu de 2 mm, et plus une couche était
       grande, plus elle était creuse (mesuré : 43 mm de remplissage au lieu
       d'environ 160 sur une couche de 424 mm²). Désormais
       `nozzle_width / min(density, 1)`, plafonné à 5 000 lignes par couche.
    2. **Les trous ne recevaient aucune paroi.** `_offset_walls` ne traçait
       que le bord extérieur de chaque offset, et seul le contour extérieur
       d'origine était ajouté : ni le contour d'un trou, ni les parois qui
       l'entourent n'existaient. Le remplissage butait contre un bord jamais
       déposé. Le test existant vérifiait « rien DANS la cavité », jamais
       « des parois AUTOUR ».
    **Suite, même jour — motifs de remplissage.** Deux défauts de plus :
    `infill_pattern` était reçu et **jamais lu** (GRID donnait des lignes
    simples), et l'angle de balayage valait 0 à TOUTES les couches (lignes
    empilées sans jamais se croiser). Nouveau motif par défaut, demandé par
    l'utilisateur : `InfillPattern.CROSSED_ZIGZAG` — lignes croisées-alternées,
    un zigzag continu par zone, selon l'axe u une couche sur deux et selon v
    l'autre. Trois points de conception, chacun mesuré ou testé :
    - **raccords qui SUIVENT LE BORD** (`_boundary_connector`) : un raccord
      droit sortait légèrement des bords irréguliers issus du maillage et
      était refusé presque à chaque ligne (médiane : 1 ligne par zigzag).
      Après : 89 % de la longueur balayée dans des zigzags d'au moins 3 lignes ;
    - **zone décalée d'une demi-buse** : le polygone reçu est celui de la paroi
      la plus intérieure, dont le bord EST l'axe de la paroi — des raccords le
      suivant auraient redéposé la paroi.
    GRID produit désormais deux balayages croisés (densité / 2 chacun) ;
    LINES reste le comportement historique. Limite connue : dans une bande
    plus étroite que quelques espacements, le zigzag est surtout fait de
    raccords (dents de scie), l'orientation imposée par l'alternance ne
    suivant pas la bande. 375/375 tests dans les deux positions.

58. **Étape 8 : la tête était RETOURNÉE sur toutes les couches (corrigé le
    2026-09-27).** Deux conventions coexistaient pour la direction à angles
    nuls : step3 et step7 prennent **+Z** (sens d'impression = axe du corps
    de la tête, +Z du modèle de tête), step8 et step10 prenaient **−Z**. Or
    la cible que step8 donne à l'IK est la normale de couche (+∇φ, sens
    d'impression) : avec la référence −Z, une couche plate (cible +Z) donnait
    **B = 179,9°** — tête retournée, corps 40 mm SOUS la buse, dans la pièce —
    et toute orientation était inversée (couche inclinée de 11,5° → B =
    168,5°). Mesuré en positionnant la vraie tête (`aetheris-4.glb`) avec
    `step9._positioned_toolhead_vertices`.

    **Pourquoi personne ne l'a vu** : aucun G-code n'est produit ; les tests
    de step8 construisaient des couches à normale **−Z** (`_simple_frame`,
    commentée « buse vers le bas ») alors que step5 produit **+∇φ**, et ne
    vérifiaient que des RELATIONS entre directions (opposées, égales), jamais
    l'orientation physique. Seule l'étape 9 en voyait les effets : ses
    collisions — et donc les supports de secours de step10 — sont à relire à
    la lumière de ce défaut.

    **Correction** : référence **+Z** dans step8 (et step10 par cohérence :
    step10 ne compare que des écarts entre directions de même référence, son
    comportement ne dépendait pas du signe). `deposition_reversed` (viser
    −normale) redevient exact : tête retournée, corps sous la buse. 2 tests
    ajoutés qui vérifient la position RÉELLE du corps de la tête (au-dessus
    pour une couche plate, en dessous en dépôt inversé), contrôlés comme
    ÉCHOUANT sans la correction.

## Prochaine étape

**Les 10 étapes du cœur algorithmique sont écrites, stables et testées
(272/272 tests).** `MachineProfile` décrit une cinématique générique à N
axes de rotation composés (`rotation_axes: list[RotationAxisSpec]`, limites
#18/#19) ; step3, step7, step8, step9 et step10 dispatchent tous sur la même
famille de fonctions et de conventions génériques
(`composed_rotation_directions`/`composed_rotation_angles_and_directions`,
composition en ordre inverse de la liste) sans hypothèse figée sur le
nombre d'axes de la machine réellement utilisée.

`main.py` enchaîne maintenant step1→step10 sur un fichier STL/3MF réel
(vérifié sur `res/test_models/l_shape_overhang.stl`, généré par
`scripts/generate_test_model.py`), et `file_io/visualization_export.py`
produit une scène GLB inspectable (mesh + couches + trajectoire) — voir
« Détail de `main.py` et `file_io/visualization_export.py` » ci-dessus.

L'orientation des couches près du plateau/des jonctions bénéficie
maintenant de deux fixes complémentaires, activés par défaut : le mélange
spatial de directions (limite #32, `_blend_tet_directions`) et le
raffinement local de maillage près du plateau (limite #36, `add_region`) —
confirmés visuellement sur `out_test.glb`. Le cas jonction multi-région
reste sans raffinement local dédié (limite #36, portée explicitement
limitée au cas plateau) — traitable plus tard avec le même mécanisme
`add_region` si un modèle réel à plusieurs régions l'exige.

**Le diagnostic du modèle en T réel (surplomb 3 axes,
`res/test_models/l-triple-axis-overhanged-shape.stl`) a mené à trois fixes
supplémentaires cette session, tous opt-in (défaut `None`, comportement
historique inchangé) — voir limites #38/#39/#40/#41 pour le détail complet :**
1. **step3** — critère angulaire de fusion de régions
   (`max_region_merge_angle_deg`) : empêche deux sous-volumes
   géométriquement incompatibles (ex: deux bras perpendiculaires) de
   fusionner sur la seule base qu'une direction de compromis reste
   techniquement compatible. 1 région → 5 régions sur le modèle en T.
2. **step5** — clustering spatial des boucles disjointes
   (`split_disjoint_clusters`) : scinde une couche en plusieurs `Layer`
   quand ses contours forment des composantes spatialement disjointes et
   non imbriquées (jonction en T/étoile), chacun avec son propre
   `local_frame` cohérent.
3. **step5** — tri par continuité spatiale (DFS sur un graphe de proximité
   entre clusters d'isovaleurs successives, couplé automatiquement à
   `split_disjoint_clusters`, pas de paramètre séparé) : résout le résidu
   du fix précédent — deux clusters d'une même région step3 mais de
   composantes spatiales disjointes (pied/bras d'un T) ne s'entrelacent
   plus à chaque isovaleur. Grands sauts spatiaux entre clusters réduits de
   42 → 7 (-83%) sur le modèle en T.

Effet de bord mesuré et **non résolu** (limite #40) : le nombre
d'orientations IK inatteignables (`step8_diagnostics`) sur le modèle en T
reste à 528 avec ces trois fixes actifs — point ouvert pour une session
future, indépendant du tri des couches lui-même.

Repère visuel ajouté au GLB (limite #41) : une flèche verte indique le
"haut" d'impression (+Z) dans `file_io/visualization_export.py`.

**Résidu après le fix DFS (limite #40), diagnostiqué et en grande partie
résolu plus avant cette session** : sur le modèle en T, quelques trajets à
vide subsistaient entre branches (7 grands sauts mesurés). Trois pistes
tentées pour les réduire directement en step5 (glouton plus-proche-voisin
sur l'ordre des composantes, fenêtre de recherche élargie du graphe de
continuité, seuil de connexion adaptatif dérivé du gradient local) ont
toutes MESURÉ un résultat pire que l'état déjà commis (156.2 → jusqu'à
203.7 de distance totale de trajet) — code exploratoire écarté, jamais
commité. Cause racine finale identifiée par ce diagnostic : **PAS un
problème d'algorithme de continuité step5**, mais un maillage tétraédrique
bien trop grossier au bout du bras concerné (6 tétraèdres seulement dans
toute cette zone), qui a mené à la découverte du bug #42 (`resolution` sans
effet réel hors zone plateau).

**Confirmé après le fix du bug #42** (aucun changement de code step5
nécessaire) : re-mesuré sur le modèle en T avec les mêmes params
(`max_region_merge_angle_deg=30.0`, `split_disjoint_clusters=3.0`) —
grands sauts **7 → 5**, distance totale de trajet **156.2 → 101.7** (-35%),
zone du bras+X passée de 6 à 34 tétraèdres. Effet de bord positif
supplémentaire mesuré : orientations IK inatteignables (`step8_diagnostics`)
**528 → 426** sur ce même modèle. Les sauts résiduels (5-6 selon le seuil de
mesure) sont désormais **confirmés par mesure comme structurellement
inévitables**, pas seulement supposés tels : augmenter la résolution du
maillage ×30 (1913 → 57992 tétraèdres) ne les réduit pas du tout, et leurs
isovaleurs décroissantes montrent qu'il s'agit de transitions entre branches
disjointes du même rang topologique — voir le dernier point de la limite #40
pour la mesure complète et pourquoi la piste "raffiner les zones de surplomb"
a été écartée. Chantier clos, sauf si une optimisation de l'ORDRE de visite
des branches (type voyageur de commerce sur les extrémités) est un jour
souhaitée — sujet distinct de la continuité de couches.

### Débit de step9 dégressif avec la matière accumulée (limite ouverte)

Le débit de step9 se dégrade avec la matière accumulée (mesuré : 16.5 →
3.7 pts/s entre le 250ᵉ et le 1750ᵉ point échantillonné, index passant de
quelques centaines à ~2700 capsules ; la grille spatiale de
`_PrintedGeometryIndex` limite le nombre de candidats mais pas assez). Sur
`l_shape_overhang.stl`, la trajectoire complète (12105 points, modèle de
tête réel) ne termine pas step9 en moins de 550s alors que les 250 premiers
points seuls prennent 0.8s. Problème de performance distinct de la
correction du signal — non traité, à chiffrer/prioriser séparément.

Ce qui reste, hors cœur algorithmique :

- **Champ scalaire multibranche** (mis en feuille de route par l'utilisateur le
  2026-10-07) : une branche de champ par partie portée, chacune croissant depuis son
  attache, raccordées aux jonctions. Seule piste restante contre le défaut de fond du
  champ unique (surplombs qui doivent croître vers le bas, limites #51 à #54) ; le
  « flip conditionnel » est clos. Aucune mesure encore.

- **Un `resegmentation_hook` réel pour l'étape 10** — `main.py` le laisse
  à `None` pour l'instant (voir détail ci-dessus). L'implémenter
  demanderait de re-jouer step3→step8 avec une contrainte d'angle
  supplémentaire (limite #23) ; c'est `main.py`, seul détenteur de
  `mesh_data`, `tet_mesh`, `overhang_map`, du graphe de régions, de la
  `ColorMap` et des params des 7 étapes concernées, qui devra le fournir.
- **`config/params.py` + `config/schema.py`** — centraliser les dataclasses
  de paramètres (aujourd'hui définies dans chaque `stepN_*.py`) et la
  (dé)sérialisation JSON, pont avec l'UI TypeScript. À faire sans casser
  les imports existants.
- **`file_io/gcode_export.py`** — traduction de `ToolpathPoint.axis_angles`
  (générique, N axes) vers la commande firmware réelle
  (`G00 Y<C> X<R> Z<Z> B<B> F<feedrate>` pour la machine actuelle). Cette
  traduction est **spécifique à une machine** et a été délibérément tenue
  hors de step8/step9/step10 ; un plan dédié devra trancher le format
  exact, la gestion de la rétraction, et comment décrire le mapping
  axe générique → nom d'axe firmware sans re-figer le moteur sur une seule
  machine.
- **Types de machines à prendre en charge** (annoncés par l'utilisateur le
  2026-09-28, voir limite #59). Tous se décrivent par
  `MachineProfile.rotation_axes` et `position_coupled_axes`, sans changer le
  contrat :
  | type | profil |
  |---|---|
  | machine actuelle (R C Z B) | axes `[C, B]`, `position_coupled_axes = [0]` |
  | + « axe D », 2ᵉ rotation de tête, perpendiculaire à la rotation autour de Y | axes `[C, B, D]`, `[0]` : B et D orientent, C reste imposé par la position |
  | + plateau coulissant le long de X (axe linéaire) | axes `[C, B]`, `[]` : C devient libre pour l'orientation |
  | les deux | axes `[C, B, D]`, `[]` : un axe d'orientation redondant |

  L'upgrade de la machine, prévue après Maker Faire, ajoute un axe linéaire.
  Reste à faire pour chaque type :
  - un test d'étape 8 ;
  - les directions candidates de step3 et step7, qui dépendent de la
    position quand C est couplé ;
  - la réorientation de step10 ;
  - la traduction G-code (nom et sens de chaque axe dans le firmware).
- **Interface graphique des paramètres d'angle** (demande de l'utilisateur,
  2026-09-27) : tout paramètre d'angle (ex. seuil d'arête vive du lissage
  non-planaire, défaut 150°) se règle par un **sélecteur graphique** qui
  dessine l'angle comme un **arc de cercle**, avec la valeur en degrés
  inscrite dedans. Convention affichée pour les arêtes : **angle ENTRE LES
  FACES** (angle intérieur, 180° = plat), et non l'angle entre normales des
  slicers classiques (150° entre faces = 30° entre normales). Le moteur peut
  stocker l'une ou l'autre ; la conversion doit être unique et explicite.
- **Ne pas centrer les modèles sur le plateau** (demande de l'utilisateur,
  2026-09-29). Sur la machine actuelle, C est imposé par la position et B ne
  penche que vers l'extérieur (+R) : les inclinaisons atteignables dépendent
  donc de la POSITION de la pièce sur le plateau. Une pièce centrée a des
  surplombs tournés vers le centre, que la tête ne peut pas servir. À faire :
  - un **conseil dans l'interface** : ne pas centrer les modèles, et orienter
    leurs surplombs vers l'extérieur du plateau ;
  - une **contrainte dans l'auto-organiseur de plateau**, qui place les
    modèles en conséquence.

- **Refuser dès l'entrée un modèle trop grand pour la machine** (demande de
  l'utilisateur, 2026-09-27). Le pré-calcul de collision tête-plateau du
  lissage non-planaire (recherche, `surface_skin.py`) traite le plateau
  comme un PLAN INFINI — hypothèse prudente, conservée tant que le rayon
  réel du plateau n'est pas dans `MachineProfile`. Le remède voulu n'est pas
  d'affiner ce test mais de refuser à l'ingestion (step1) un modèle qui
  déborde du volume/plateau utilisable (`MachineProfile.build_volume`,
  aujourd'hui non renseigné), pour qu'aucun modèle trop grand ne traverse le
  pipeline.

### ROADMAP — Maker Faire Rome et temps de calcul (actée le 2026-09-27)

**Échéance : Maker Faire Rome, dans moins d'un mois (fin octobre 2026).**
Ce qui doit y être montré, choisi par l'utilisateur :
- **des pièces imprimées exposées**, imprimées à l'avance ;
- **une impression en direct** sur le stand.

Pas de tranchage en direct : la vitesse de calcul n'est donc PAS un
prérequis du salon. Deux choses le sont :

1. **Export G-code réel pour la machine** (`file_io/gcode_export.py`, jamais
   écrit ; commande `G00 Y<C> X<R> Z<Z> B<B> F<feedrate>`, rétractions à
   trancher). **État inconnu, à vérifier en premier** : des scripts non
   versionnés existent (`scripts/gen_l_overhang_gcode.py`,
   `gen_elbow_tube_gcode.py`, `preview_gcode.py`,
   `export_machine_animation.py`) et produisent peut-être déjà un G-code qui
   imprime. Le relire avant d'écrire quoi que ce soit.
2. **Trajectoire vérifiée en collision avec la vraie tête** (étape 9,
   `aetheris-4.glb`) avant toute impression en direct : une collision devant
   le public est le pire scénario. À relire à la lumière de la limite #58
   (tête retournée jusqu'au 2026-09-27 : toutes les collisions antérieures
   sont suspectes).

#### Temps de calcul : 10 minutes pour une étape est COMMERCIALEMENT INACCEPTABLE

| phase | durée mesurée | cause |
|---|---:|---|
| P2 — zones de lissage (`surface_skin.py`) | **11 min 30** (200 s pré-calcul + 490 s zones) | boucle Python face par face ; faces retestées à chaque découpage en bandes |
| P3 — corps creusé (`skin_carve.py`, couche par couche) | **~10-11 min** | défaut de CONCEPTION : la même découpe refaite 313 fois, par tranches, à coups d'unions shapely de milliers de morceaux |
| tranchage de référence complet | 37 s (dont 3 s de couches) | — |

**Règle de travail qui en découle** : mesurer la durée de chaque phase et
la signaler dès qu'elle dépasse l'ordre de la seconde, au lieu de la
découvrir en fin de calcul. Lancer les calculs en tâche de fond avec
`python -u` (sinon le journal reste vide jusqu'à la fin, stdout étant
tamponné vers un fichier).

**P2 — trois leviers, dans cet ordre** (estimations NON mesurées) :
1. **Dernier moment accessible par face**, calculé une seule fois. La
   matière imprimée ne fait que croître : une face bloquée par l'objet à un
   moment l'est à tous les suivants (la collision plateau ne dépend pas du
   moment). Le découpage en bandes devient une comparaison d'entiers.
   Résultat attendu IDENTIQUE face par face (à vérifier sur la carte).
2. **Traitement groupé** : placer la tête et localiser les points pour des
   milliers de faces en un appel numpy, au lieu de ~10 petits appels par face.
3. **Parallélisation** multi-cœurs du pré-calcul (faces indépendantes).

| version de P2 | estimation |
|---|---:|
| Python actuel | 690 s (mesuré) |
| Python + leviers 1 et 2 | ~30-60 s |
| C++ mono-cœur | ~10-20 s |
| C++ 8 cœurs | ~2-5 s |
| GPU | calcul < 1 s ; total ~1-3 s (parties fixes : chargement, subdivision, voisinages) |

**P3 — refonte : découper UNE fois en 3D, puis trancher normalement.** La
découpe est une propriété de la pièce, pas de chaque couche. (1) distance à
la peau aux nœuds ; (2) affinage du maillage près de la peau (précision
0,2 mm ; tétraèdres actuels de 1 à 5 mm) ; (3) découpe des tétraèdres le long
de ψ = 0 (table de cas, comme les marching tetrahedra) ; (4) tranchage
habituel du corps creusé. La version actuelle, lente mais exacte, sert
d'ORACLE à la version rapide (même méthode que le portage C++).

| version de P3 | estimation |
|---|---:|
| actuelle | ~10-11 min (mesuré) |
| découpe 3D, Python | ~5-15 s |
| découpe 3D, C++ multi-cœurs | ~1-3 s |
| + GPU | ~1-2 s |

**Grande inconnue : l'affinage (2).** S'il multiplie le nombre de
tétraèdres par 5 à 10, le tranchage qui suit ralentit d'autant (3 s →
15-30 s en Python). Première mesure à faire : le nombre de tétraèdres après
affinage. Si l'affinage passe par TetGen (natif, mono-cœur), c'est lui qui
fixe le plancher.

**GPU (question de l'utilisateur, DirectX évoqué)** — ce qui s'y prête et ce
qui ne s'y prête pas :
- **oui** : distance à la peau par nœud, découpe par tétraèdre, coupes des
  couches par tétraèdre, pré-calcul de collision de P2 par face (~30 M de
  tests point-dans-tétraèdre) ;
- **non** : affinage du maillage (topologie changeante), chaînage des
  segments en contours dans une couche, décalages de polygones de l'étape 6
  (Clipper/shapely) — processeur multi-cœurs, parallèle par couche.

Réserves à peser avant le choix d'API :
- **DirectX = Windows seulement** ; Vulkan, WebGPU ou CUDA (NVIDIA
  seulement) si l'UI doit être multiplateforme ;
- les GPU grand public sont rapides en **simple précision** seulement, ce
  qui casse l'équivalence au bit près sur laquelle repose le portage (limite
  #56 : porter d'abord, paralléliser ensuite) ;