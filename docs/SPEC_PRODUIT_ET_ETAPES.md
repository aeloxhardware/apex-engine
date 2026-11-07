# Spécification du produit et des 10 étapes

_Extrait de l'ancien CLAUDE.md (R-030 point 7)._

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

> **Depuis le 2026-10-29, chaque algorithme de tranchage vit dans son paquet `element_xxxx/`** (décision de
> l'utilisateur) : `element_ardor/` (pipeline principal 5 axes ; l'ancien `pipeline/`), `element_aeolis/`
> (non planaire ; l'ancien `nonplanar/`), `element_aequor/` et `element_arva/` (à venir, README seulement). **Depuis le
> 2026-11-07 ce sont de vrais paquets** (R-030 point 7) : les alias d'import `pipeline/` et `nonplanar/` (avec `exec`)
> sont supprimés ; on écrit `from element_ardor.stepN_... import` et `python -m
> element_ardor._tetrahedralize_worker`. L'arborescence ci-dessous cite encore `pipeline/` : lire
> `element_ardor/`. Le code partagé (`core/`, `file_io/`, `config/`, `cpp/`) reste à la racine.

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

**Dépendances déclarées et figées (2026-11-04, R-030 point 5)** : `pyproject.toml` liste ce que le code importe ; `requirements.lock`
(`pip-compile --generate-hashes`) fige les versions avec empreintes (`pip install --require-hashes -r requirements.lock`) ;
licences : `docs/LICENCES_DEPENDANCES.md` (généré par `scripts/gen_licences.py`, vérifié en CI). La liste ci-dessous est historique.

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

**Ce que le code fait précisément (`element_ardor/step4_field_solver.py`) :**
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

**Ce que le code fait précisément (`element_ardor/step6_toolpath2d.py`, FAIT) :**
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

**Ce que le code fait précisément (`element_ardor/step7_scheduling.py`, FAIT) :**
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

**Ce que le code fait précisément (`element_ardor/step9_collision_sim.py`, FAIT) :**
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
l'appelant de `element_ardor/step9_collision_sim.run()` qui doit peupler
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

**Ce que le code fait précisément (`element_ardor/step10_adjustment.py`, FAIT — boucle jusqu'à `max_adjustment_iterations`) :**
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

