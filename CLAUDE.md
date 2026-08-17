# Apex Slicer — Contexte projet — FICHIER MODIFIABLE PAR L'IA

Moteur de slicing géométrique (Python) pour une imprimante 3D **4 axes** (3 linéaires
+ 1 rotatif). L'interface utilisateur est en **TypeScript**, séparée de ce moteur.
Ce fichier centralise tout le contexte nécessaire pour reprendre le développement
sans redemander les spécifications déjà actées. L'IA peut et doit modifier ce fichier lorsqu'elle opère.

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

```
slicer_engine/
├── CLAUDE.md                        # ce fichier
├── main.py                          # FAIT — orchestrateur CLI, enchaîne step1 à step10
├── scripts/
│   ├── generate_test_model.py       # FAIT — génère res/test_models/l_shape_overhang.stl (surplomb réel)
│   └── export_slicing_preview.py    # FAIT — GLB du slicing SEUL (step1→step8), court-circuite
│                                     #   step9/step10 le temps que la zone de garde buse soit faite
├── config/                          # PAS ENCORE ÉCRIT
│   ├── params.py                    #   dataclasses de tous les paramètres globaux
│   └── schema.py                    #   (de)sérialisation JSON <-> dataclasses, pont avec l'UI TS
├── pipeline/
│   ├── __init__.py
│   ├── step1_ingest.py              # FAIT — ingestion mesh + tétraédrisation (+ détection
│   │                                 #   automatique des cavités internes pour tetgen, voir limite #7)
│   ├── step2_overhang.py            # FAIT — détection surplombs + ponts internes
│   ├── step3_segmentation.py        # FAIT — tests pytest formalisés (10 tests)
│   ├── step4_field_solver.py        # FAIT — résolution FEM du champ scalaire (7 tests)
│   ├── step5_layers.py              # FAIT — extraction des couches courbes (12 tests)
│   ├── step6_toolpath2d.py          # FAIT — remplissage 2D par couche (8 tests)
│   ├── step7_scheduling.py          # FAIT — couleur + ordonnancement trous/purges (22 tests)
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
│   └── constraints.py               # FAIT — registre extensible d'assembleurs de contraintes
│                                     #   utilisateur (ATTACH_POINT/AESTHETIC/STRUCTURAL), point
│                                     #   d'extension pour tout nouveau type sans toucher step4
├── file_io/                          # renommé depuis `io/` — collision avec le module
│   │                                 #   standard Python `io` (voir limite #20)
│   ├── mesh_loader.py                # FAIT — load_toolhead_mesh (GLB -> Trimesh mm, scene-space)
│   ├── visualization_export.py       # FAIT — export_visualization_glb (scène GLB : mesh + couches + trajectoire)
│   └── gcode_export.py               # PAS ENCORE ÉCRIT
├── res/
│   └── test_models/
│       └── l_shape_overhang.stl      # FAIT — modèle de test avec surplomb réel (même géométrie que
│                                      #   tests/test_step3_segmentation.py::l_shape_data)
└── tests/
    ├── __init__.py
    ├── test_step2_overhang.py       # FAIT — 6 tests, tous passants et stables
    ├── test_step3_segmentation.py   # FAIT — 10 tests
    ├── test_step4_field_solver.py   # FAIT — 7 tests
    ├── test_step5_layers.py         # FAIT — 12 tests
    ├── test_step6_toolpath2d.py     # FAIT — 8 tests
    ├── test_step7_scheduling.py     # FAIT — 22 tests
    ├── test_geometry_utils.py       # FAIT — 6 tests (composed_rotation_directions)
    ├── test_step8_machine_toolpath.py  # FAIT — 16 tests
    ├── test_step9_collision_sim.py  # FAIT — 18 tests
    └── test_step10_adjustment.py    # FAIT — 15 tests
```

**État global : 120/120 tests passants** (step2 à step10 + core/geometry_utils.py).
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
  //   rotation_axes = [ {axis: (0,0,1), angle_range_deg: (0,360)},   // C
  //                      {axis: perpendiculaire à R, angle_range_deg: (0,180)} ]  // B
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
  laplacian_smoothness_weight: float
  gradient_overhang_weight: float
  resolution: float                                   // taille de tétra pour la résolution FEM
  max_solver_iterations: int
  convergence_tolerance: float
}

Constraint[] {                                        // un par point d'attache / zone
  type: "attach_point" | "aesthetic" | "structural"
  position: Vec3
  target_field_value: float | null
  target_normal: Vec3 | null
  influence_radius: float
  falloff: enum (linear, gaussian, cubic)
  weight: float
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
  laplacian_smoothness_weight: float                 // pondération lissage vs contrainte d'angle
  gradient_overhang_weight: float                    // poids de la contrainte "gradient ~ normale locale"
  resolution: float                                   // taille de voxel/tétra pour la résolution FEM
  max_solver_iterations: int
  convergence_tolerance: float
}

Constraint[] (envoyés dynamiquement, un par point d'attache ou zone) {
  type: "attach_point" | "aesthetic" | "structural"
  position: Vec3
  target_field_value: float | null
  target_normal: Vec3 | null
  influence_radius: float
  falloff: enum (linear, gaussian, cubic)
  weight: float
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
- Construit le graphe de dépendances entre régions : pour chaque paire de régions adjacentes, détermine par ray-cast/test de contact qui doit être imprimée avant qui (résolution par tri topologique, `networkx.topological_sort`).

**Sortie :** `region_id` par face + `RegionDependencyGraph`.


## Étape 4 — Construction et résolution du champ scalaire continu

**Ce que le code fait précisément :**
- Sur le `TetMesh`, assemble l'opérateur Laplacien discret (matrice creuse, `scipy.sparse`).
- Ajoute les contraintes dures/molles au système :
  - Contrainte de surplomb (issue de l'étape 2/3) : le gradient du champ à chaque tétraèdre proche d'une face à risque doit s'aligner avec la direction candidate de sa région, pondéré par `gradient_overhang_weight`.
  - Contraintes utilisateur (`Constraint[]` reçus de l'UI) : pour chaque point d'attache, ajoute une ligne au système forçant `field(position) ≈ target_field_value`, atténuée par `influence_radius` et `falloff`.
  - Lissage global pondéré par `laplacian_smoothness_weight`.
- Résout le système linéaire (`scipy.sparse.linalg.spsolve` ou solveur itératif `cg`/`gmres` si gros système, potentiellement accéléré CuPy comme sur ton pipeline mesh existant) jusqu'à `convergence_tolerance` ou `max_solver_iterations`.

**Sortie :** `ScalarField` (valeur par nœud du TetMesh).


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

7. **Les cavités internes doivent être signalées explicitement à tetgen.** Un mesh
   watertight avec un vide interne (boîte creuse) ne suffit pas : sans
   `tet.add_hole(point)` par cavité, tetgen peut générer des tétraèdres qui traversent
   le vide, corrompant silencieusement le volume résolu par l'étape 4. Bug réel détecté
   tardivement (aucun test antérieur à l'étape 5 n'exerçait la correction volumique du
   `TetMesh` sur une géométrie à cavité — step2/3 ne travaillent que sur la surface,
   step4 n'était testé que sur cube plat/forme en L). Détection automatique et générale
   (pas cas-par-cas) : séparer le mesh en coques connexes étanches
   (`mesh.split(only_watertight=True)`), une coque à volume signé négatif
   (`shell.volume`, propriété trimesh) est une cavité, son centroïde (vérifié via
   `shell.contains()`) sert de point de trou.

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

14. **La détection de pont par couche (étape 6) doit être un seuil de
    PRÉSENCE, jamais de majorité, sur la boucle de contour.** Une boucle de
    contour est le périmètre COMPLET de la silhouette de matière à cette
    hauteur de champ (ex: paroi gauche → dessus → paroi droite → dessous du
    pont → retour) ; la portion réellement `is_bridge` n'en est qu'un petit
    arc, jamais une majorité. Bug réel mesuré : avec un seuil de ratio à 0.5,
    aucune couche d'une fixture pont réaliste (dalle entre deux piliers)
    n'était jamais détectée (ratio réel mesuré entre 0.06 et 0.21 sur toutes
    les couches concernées) — un seuil de majorité ne peut structurellement
    jamais se déclencher, quel que soit sa valeur. Fix : seuil bas de présence
    (`Step6Params.bridge_contour_ratio_threshold`, défaut `0.05`).

15. **`MeshData.kdtree` répond « quelle face est la plus proche », JAMAIS
    « à quelle distance est la surface ».** Il est construit sur les
    **centres de face** (`step1_ingest`) : la distance qu'il retourne est une
    distance à un centre de face, pas à la surface, et l'écart croît avec la
    taille des triangles. Mesure sur un cube brut (12 triangles) : des points
    situés **exactement sur** la surface ressortent à 0.370–3.901 du centre de
    face le plus proche (distance réelle : 0.000) — avec un seuil de 0.5,
    3 points sur 35 seulement étaient classés « en surface ». Sur maillage
    grossier, même l'**indice** devient douteux près des arêtes. Pour toute
    vraie distance (ou tout classement surface/intérieur), passer par
    `trimesh.proximity.ProximityQuery(mesh).on_surface(points)`, qui retourne
    distance réelle **et** triangle correct. Step6 n'était pas affecté car il
    n'utilise que l'indice ; step7 (`_edge_colors_for_polyline`) a dû basculer
    sur `ProximityQuery`. Piège à garder en tête pour les étapes 8/9/10, qui
    manipuleront elles aussi des points hors surface.

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

20. **Le dossier `io/` du projet entre en conflit avec le module standard
    Python `io` (builtin/frozen).** `from io.mesh_loader import ...` échoue
    avec `ModuleNotFoundError: 'io' is not a package` — Python résout
    `io` vers le module standard avant de considérer un package local du
    même nom, quel que soit l'ordre de `sys.path`. Détecté concrètement en
    écrivant l'étape 9 (premier consommateur réel du module). Renommé en
    **`file_io/`** — tous les imports doivent utiliser `file_io.mesh_loader`,
    `file_io.gcode_export`, jamais `io.*`.

21. **Un mesh de tête positionné testé uniquement par ses SOMMETS peut
    manquer une collision réelle.** Un obstacle fin (ex: un mur imprimé)
    peut traverser une FACE du mesh sans passer près d'aucun sommet — faux
    négatif mesuré concrètement en écrivant l'étape 9 (une boîte de test
    traversée en son milieu par un mur n'était détectée par aucun de ses 8
    sommets, tous à plusieurs mm du mur). `step9_collision_sim._surface_sample_points`
    ajoute des points échantillonnés le long de chaque ARÊTE du mesh
    (pas ~`nozzle_diameter`), pas seulement ses sommets — suffisant pour un
    mesh de tête de taille modeste (aetheris-4.glb : ~377 sommets, ~1000
    arêtes), sans recourir à un calcul point-triangle plus coûteux.

22. **La buse est par construction confondue avec le cordon qu'elle vient
    de déposer.** Tester un point en extrusion contre la matière imprimée
    "jusqu'à ce point INCLUS" déclenche une collision à chaque point en
    extrusion (distance quasi nulle à son propre cordon) — faux positif
    structurel, pas un cas limite rare. `step9_collision_sim._grow_printed_index_to`
    exclut la capsule qui se termine AU point testé ; elle ne devient une
    collision possible qu'au(x) point(s) suivant(s), une fois que la buse
    s'en est réellement éloignée.

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

26. **`nozzle_diameter` (largeur de buse) et le pas d'échantillonnage
    géométrique du modèle de tête (limite #21) sont deux échelles
    PHYSIQUEMENT DIFFÉRENTES, confondues à l'origine dans step9.**
    `nozzle_diameter` décrit la taille du CORDON extrudé (0.4mm
    typiquement) ; le pas d'échantillonnage des arêtes du modèle de tête
    doit être à l'échelle de la TÊTE elle-même (plusieurs mm à quelques
    dizaines de mm). Bug de performance réel, découvert en exécutant
    `main.py` de bout en bout pour la première fois (jamais visible sur les
    mesh synthétiques ~10mm des tests unitaires) : avec
    `nozzle_diameter=0.4` comme pas d'échantillonnage sur le vrai modèle de
    tête (`aetheris-4.glb`, arêtes de plusieurs dizaines de mm),
    `_surface_sample_points` générait **17734 points de test par point de
    trajectoire échantillonné** (mesuré), rendant la simulation totale
    inutilisable même sur un modèle trivial. Corrigé par
    `CollisionSimParams.toolhead_sample_spacing` (nouveau champ, défaut
    5.0mm, découplé de `nozzle_diameter`) : **1540 points** avec ce défaut
    sur le même modèle — un facteur ~11.5×, mesuré, pas estimé.

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

28. **La pointe de buse touche la pièce et le plateau PAR CONSTRUCTION —
    la tester comme un obstacle produit une collision à chaque point.**
    L'origine du modèle de tête EST le point de dépôt (convention voulue et
    vérifiée), et `_surface_sample_points` inclut ce sommet parmi les points
    de test. Avec `safety_margin=1.0 mm`, ce seul point suffit à déclencher
    une collision contre un cordon voisin et contre le plateau à chaque
    point de trajectoire. Mesuré sur `l_shape_overhang.stl` : **exactement
    2 collisions par point échantillonné, sans exception** (498 à 250
    points, 998 à 500, 3498 à 1750) — 3588 collisions au total, soit ~100 %
    des points, ce qui rend l'étape 10 à la fois interminable et
    ininterprétable. La limite #22 ne couvre que le cordon en cours de
    dépôt, pas les cordons voisins ni le plateau. Correction prévue
    (`nozzle_contact_clearance`) décrite dans « Prochaine étape ».

29. **L'orientation d'un modèle de tête doit être vérifiée par MESURE, pas
    déduite de la convention documentée.** `aetheris-4.glb` a son corps
    étendu selon **+Y** dans le fichier (axe long mesuré : 77.2 mm en Y,
    contre 76.0 en X et 35.0 en Z), alors que la note de convention de ce
    même document affirmait « Axe Z = sens de la tête ». Conséquence
    concrète avant correction : le mesh chargé allait de `z=-18.9` à
    `z=+16.1`, soit **18.9 mm de tête pendant SOUS la pointe de buse** — à
    une buse posée à `z=1 mm`, **47 % des points d'échantillonnage de la
    tête tombaient sous le plateau** (20 % encore à `z=11 mm`). Corrigé par
    `mesh_loader.TOOLHEAD_BODY_TO_MACHINE_UP` (rotation +90° autour de X,
    appliquée dans le loader donc une seule fois à la source). Si un futur
    modèle respecte réellement la convention Z, cette rotation devra
    devenir une propriété **par modèle**, pas une constante globale.
    Corollaire pour les tests : l'ancien test ne vérifiait que les
    `extents` (invariants à l'orientation à un ordre d'axes près) et
    laissait donc passer le défaut ; il assert désormais la propriété
    physique — la tête ne pend jamais sous sa buse.

30. **Une re-vérification conçue « une fois par point » ne doit jamais être
    appelée depuis une boucle de recherche.** `_collisions_at_points`
    (step10) reconstruit à chaque appel les échantillons de tête ET un
    index de matière imprimée regrossi depuis 0. Correct pour son usage
    d'origine (une passe en fin d'itération), catastrophique appelé depuis
    la boucle de candidats de `_try_local_reorientation` (jusqu'à 144
    candidats × 3588 points) : coût `O(collisions × candidats × position)`.
    Deux pièges successifs mesurés en corrigeant : hoister l'index hors de
    la boucle de **candidats** ne suffisait pas — il restait reconstruit à
    chaque **point** en collision, laissant une dégradation de 2.48 → 0.39
    calls/s. Il a fallu un `_ReorientationContext` partagé sur toute une
    itération externe de `run()`, grossi de façon monotone (même invariant
    que step9) : taux devenu stable, 0.30 → 0.44 calls/s.

31. **Un critère de fusion purement binaire (« au moins une direction
    compatible existe ») fusionne des régions qui ne devraient pas l'être,
    même quand un coût d'angle est ajouté ailleurs.** Diagnostiqué en
    inspectant `out_test.glb` : la base verticale d'un L (aucun surplomb)
    apparaissait courbée dès les premières couches. Racine tracée jusqu'à
    step3, pas step4 : `_region_growing` fusionnait la base et le bras en
    **une seule région** dès qu'**une** direction candidate satisfaisait les
    24 faces réunies (mesuré : 54.8° de la verticale, direction ex-æquo au
    minimum de risque résiduel avec 2 autres candidats sur 27 — le
    départage par `np.argmin` n'était PAS en cause, la verticale pure avait
    un résiduel réellement pire, 3.9281 contre 0.0000). **Un coût de
    déviation ajouté seulement dans `_select_region_directions` ne peut pas
    corriger ça** : cette fonction choisit une direction à l'intérieur
    d'une région déjà figée par `_region_growing`, elle ne peut jamais
    scinder une région existante — vérifié par mesure, resserrer
    `max_angle_from_direction_deg` de 45° à 30° ne changeait rien au nombre
    de régions (`_region_growing` reste un test binaire, insensible à ce
    levier). Fix : `SegmentationParams.direction_deviation_weight` (défaut
    `0.0`, jamais actif par défaut pour ne pas changer le comportement des
    appelants existants) rend la décision de fusion sensible au coût
    (`_best_direction_cost`, fonction partagée avec `_select_region_directions`
    — jamais deux formules divergentes) : fusionner n'est retenu que si son
    coût ne dépasse pas la somme des coûts séparés + `region_count_penalty`.

32. **`region_count_penalty` était présent dans `SegmentationParams` et
    transmis à `_region_growing` depuis le début, mais n'était JAMAIS lu
    dans le corps de la fonction** — seulement mentionné dans un
    commentaire affirmant qu'il « module indirectement » la fusion via
    l'ordre de traitement, ce qui était faux (l'ordre ne dépendait que de
    `risk_score`). Le fix de la limite #31 lui donne son premier usage réel
    ; sa valeur par défaut (`0.1`) n'avait donc jamais été validée
    empiriquement, et s'est révélée **10× trop élevée** : mesuré sur le L
    de test, la scission base/bras se maintient jusqu'à
    `region_count_penalty=0.01`, se rompt dès `0.02` (retombe à 1 seule
    région). Nouveau défaut `0.005` (marge de sécurité, moitié du seuil de
    rupture mesuré), robuste sur un balayage `direction_deviation_weight`
    de 0.10 à 0.50. `direction_deviation_weight=0.2` activé explicitement
    dans `main.py`/`PipelineParams` (pas dans le défaut de
    `SegmentationParams` lui-même) une fois validé par
    `tests/test_step3_segmentation.py::test_l_shape_with_deviation_weight_splits_base_from_arm`.
    Résultat vérifié sur `out_test.glb` régénéré : amplitude Z des couches
    basses passée de ~3-4 mm à ≤0.01 mm (109 couches, contre 142 avant —
    la base ne génère plus de couches artificiellement inclinées).
    Effet de bord observé et non traité ici : le nombre d'orientations
    machine inatteignables signalées par l'étape 8
    (`Step8Diagnostic`) est passé de 0 à 851/10089 points (8.4 %, erreur
    résiduelle médiane 8.0°, max 8.3°), concentrées sur le bras — la
    direction demandée y est désormais plus franche (49.2° au lieu d'être
    diluée dans une moyenne à 54.8° sur toute la pièce), ce qui sollicite
    davantage la grille de recherche IK de step8. Comportement structuré et
    documenté (`Step8Diagnostic`, jamais une exception), pas un défaut de
    ce fix — mais à garder à l'esprit si la qualité d'orientation du bras
    doit être affinée plus tard (densité de grille IK, `ik_refine_step_multiplier`).

33. **Un tri de couches par « rang topologique de région puis isovaleur »
    (au lieu d'isovaleur seule) produit un ordre d'impression physiquement
    aberrant dès qu'une région de rang faible possède un fragment
    géométriquement voisin d'une région de rang élevé.** Découvert en
    inspectant visuellement `out_test.glb` régénéré après le fix de la
    limite #31 (scission base/bras réelle) : de longs trajets à vide en
    éventail au niveau du coude. Diagnostic en 2 temps, la première
    hypothèse formulée puis invalidée par mesure (pas supposée) :
    - Hypothèse initiale : incohérence entre step4 (`nearest_face_per_tet`)
      et step5 (`_map_tets_to_regions`), deux mappings tétraèdre→région
      indépendants basés sur le même KDTree. **Invalidée par comparaison
      directe** (`np.array_equal` sur les deux mappings) : les deux
      fichiers produisent EXACTEMENT le même résultat. Les 57 tétraèdres
      incriminés (18,6 % du volume, jonction base/coude) sont bien et
      correctement classés région "base" (verticale) par proximité de
      surface — leur valeur de champ élevée (9,0-11,0, proche de celle du
      bras voisin) est un effet parfaitement normal du lissage laplacien
      global (continuité spatiale), pas une erreur de classification.
    - Cause réelle : `pipeline/step5_layers.py::run` triait `raw_layers`
      par `(rang_topologique_max_des_régions_traversées, isovaleur)`
      (`_sort_key`, retiré). Avec une région "base" (rang 0) qui possède un
      petit fragment déconnecté à isovaleur élevée (le sommet du coude, vu
      depuis la direction verticale de la base) et une région "bras"
      (rang 1) qui couvre le même volume physique à isovaleur intermédiaire,
      **toute couche touchant uniquement la région de rang 0 était triée
      avant TOUTE couche touchant la région de rang 1**, sans égard à la
      position physique réelle. Mesuré : isovaleur sautant de 4,09 (index
      39, région 1) à 10,5 (index 40, région 1 encore) puis retombant à
      4,20 (index 48, régions 0+1) — deux passages sur la même plage
      4,2-10,7, la seconde fois seulement.
    - Fix : tri par **isovaleur seule** (`raw_layers.sort(key=lambda e: e[0])`).
      Justification : le champ scalaire est résolu globalement (FEM), son
      isovaleur encode déjà par construction la progression spatiale
      cohérente — le rang topologique de région reste correct pour son
      usage d'origine (dépendances physiques entre régions dans
      `RegionDependencyGraph`), mais n'a pas vocation à réordonner des
      couches déjà cohérentes entre elles par leur isovaleur.
      `_compute_region_ranks` conservée (fonction pure, testée), simplement
      plus appelée pour trier. `pipeline/step7_scheduling.py` consomme
      `layers` via `enumerate(layers)` comme ordre d'impression réel — ce
      n'était pas qu'un défaut d'affichage GLB, un mauvais tri ici cassait
      l'ordonnancement en aval.
    - Piège pour les fixtures de test : `tests/test_step5_layers.py::l_shape_data`
      utilise `SegmentationParams()` par défaut (`direction_deviation_weight=0.0`),
      donc ne scinde PAS la pièce — elle n'a jamais exercé ce cas, avant
      comme après ce fix. Nouvelle fixture dédiée
      `l_shape_split_data` (même géométrie, `direction_deviation_weight=0.2`
      explicite) pour réellement tester le scénario à 2 régions.

34. **Une couche peut légitimement contenir des segments de trajectoire
    spatialement DISJOINTS (clusters à des dizaines de mm l'un de l'autre)
    même une fois l'ordre inter-couches corrigé (limite #33) — problème
    DIFFÉRENT, ne pas confondre les deux.** Après le fix de la limite #33,
    `out_test.glb` régénéré montrait encore un faisceau de trajets à vide
    au niveau du coude. Vérifié que ce N'ÉTAIT PAS un reliquat non corrigé
    de la limite #33 : l'isovaleur des couches était déjà strictement
    monotone. Cause distincte, mesurée : `pipeline/step7_scheduling.py`
    construisait chaque `PrintPass` en concaténant les segments d'une
    couche dans leur ordre d'extraction de step6 (paroi/infill par bloc
    géométrique local), sans tenir compte de la POSITION DE LA TÊTE à la
    fin de la passe précédente. Résultat mesuré sur le modèle de test :
    83 sauts de trajectoire > 10 mm sur 10088 segments, un aller-retour
    quasi complet base↔bras à CHAQUE couche de la zone de jonction (chaque
    couche listant systématiquement tout son cluster proche de la base
    puis tout son cluster proche du bras, dans cet ordre fixe).
    Fix : `_order_segments_by_proximity` (nouvelle fonction) réordonne les
    segments d'une passe par plus proche voisin glouton, en partant de la
    position réelle de fin de la passe précédente (`head_position`, état
    maintenu dans `_schedule_passes` comme `active_color` l'est déjà).
    Résultat mesuré : 83 → 33 sauts > 10 mm — **pas zéro, et c'est le
    résultat ATTENDU, pas un échec partiel** : 35 passes sur 109 touchent
    réellement les deux clusters à la fois (mesuré), et réordonner des
    SEGMENTS entre eux (jamais les POINTS à l'intérieur d'un segment,
    scope volontairement limité) ne peut pas éliminer la traversée
    lorsqu'une passe unique doit de toute façon visiter les deux clusters
    — seulement la placer une fois par passe concernée plutôt qu'à répétition
    entre passes consécutives. Les grands sauts restants sont désormais
    strictement internes à une passe (vérifié : le point de fin d'une passe
    et le point de début de la suivante coïncident systématiquement,
    `head_position` fonctionne comme prévu), plus jamais un aller-retour
    entre deux passes déjà bien enchaînées.

    **Distinction explicite avec la limite #33, à ne jamais fusionner** :
    #33 est un défaut d'ORDRE DES COUCHES (une couche entière imprimée au
    mauvais moment dans la séquence globale) — un vrai bug de correction
    physique, résolu. #34 est un choix de TRAJECTOIRE DE DÉPLACEMENT DE
    TÊTE à l'intérieur d'une séquence de couches déjà correctement ordonnée
    — une optimisation de coût (minimiser les trajets à vide), pas une
    correction de validité. Ce fix ne prétend RIEN résoudre côté ordre
    physique inter-couches ; si un futur diagnostic révèle un vrai problème
    d'ordre entre couches (pas seulement de distance de trajet), il devra
    être traité comme la limite #33 (dans step5/step3), jamais confondu
    avec l'optimisation de trajectoire de la limite #34 (step7).

---

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
  voisin. Seule primitive par-axe : `_rodrigues_rotate` (déjà existante) —
  généralise à N axes sans jacobien analytique par-N. Le pas initial DOIT
  être dérivé de l'écart de la grille NN (`360/samples_per_axis`), pas une
  constante fixe, sous peine de rester piégé près du candidat NN (mesuré :
  résidu ~15-20° avec un pas fixe de 2° contre <0.1° avec un pas dérivé).

## Détail de `core/constraints.py`

Point d'extension pour tout nouveau type de contrainte utilisateur — ajouter un
type ne nécessite jamais de modifier `step4_field_solver.py` :

- `ConstraintAssemblyContext` — données partagées (TetMesh, MeshData, kdtree,
  gradients de forme précalculés), sans logique.
- `CONSTRAINT_ASSEMBLERS: dict[ConstraintType, ConstraintAssembler]` — registre ;
  ajouter un type = une fonction d'assemblage + une entrée ici.
- `assemble_all_constraints(constraints, ctx) -> (csr_matrix, ndarray)` — point
  d'entrée unique appelé par step4.

## Détail de `pipeline/step4_field_solver.py`

- Système `A·φ = b` symétrique défini positif, assemblé en trois fonctions
  pures pour permettre une re-résolution rapide après édition d'un point
  d'attache (jamais besoin de relancer step1/2/3, ni de réassembler la partie
  structurelle) :
  - `build_static_system(...) -> FieldSolverStaticSystem` — tout ce qui NE
    dépend PAS des `Constraint[]` (rigidité, direction de gradient, ancrage
    plateau). À appeler une seule fois après step3, résultat à mettre en cache.
  - `solve_with_constraints(static_system, constraints, params) -> ScalarField`
    — chemin rapide, appelé à chaque édition de contrainte.
  - `run(...)` — enveloppe des deux, pour un slicing complet en une passe.
- Solveur : direct (`spsolve`) sous `direct_solve_node_threshold` (défaut 30000
  nœuds), sinon itératif (`cg`, préconditionné) — jamais `gmres` (système
  toujours symétrique par construction).

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

## Prochaine étape

**Les 10 étapes du cœur algorithmique sont écrites, stables et testées
(120/120 tests).** `MachineProfile` décrit une cinématique générique à N
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

### PRIORITÉ IMMÉDIATE — zone de garde autour de la pointe de buse (`nozzle_contact_clearance`)

**L'étape 9 produit actuellement un faux positif structurel à CHAQUE point
de trajectoire** : le pipeline complet est donc inexploitable de bout en
bout, alors que step1→step8 (le slicing proprement dit) fonctionne et donne
un résultat visuellement cohérent.

Mesure directe (`l_shape_overhang.stl`, modèle de tête réel) : **exactement
2 collisions par point échantillonné**, sans aucune exception — 498
collisions à 250 points, 998 à 500, 1498 à 750, 3498 à 1750. Un
`printed_part` + un `platform` à chaque point. Ce n'est pas un motif de
collision réel, c'est une propriété structurelle.

**Cause** : l'origine du modèle de tête EST la pointe de buse (convention
vérifiée et voulue), et `_surface_sample_points` inclut ce sommet à
l'origine parmi les points de test. Or la pointe de buse est **par
définition** en contact avec la matière qu'elle vient de déposer et à
quelques dixièmes de millimètre du plateau aux premières couches. La tester
avec `safety_margin=1.0 mm` garantit une collision à chaque point.

La limite #22 (exclusion du cordon en cours de dépôt) ne suffit pas : elle
n'exclut que la capsule qui se termine AU point testé, pas les cordons
voisins déjà posés sur la même couche, ni le plateau.

**Correction prévue** : nouveau champ `CollisionSimParams.nozzle_contact_clearance`
(mm) — les points d'échantillonnage du modèle de tête situés à moins de
cette distance de l'origine locale (la buse) sont exclus des tests de
collision, la zone de contact buse/matière étant attendue par construction.
Points d'attention pour l'implémentation :
- Le filtrage doit se faire **une fois** sur les points en repère LOCAL
  (dans `_surface_sample_points` ou juste après), jamais par point de
  trajectoire — sinon on repaye le coût à chaque appel.
- Le défaut doit être à l'échelle de la buse et de la marge de sécurité
  (ordre de quelques mm), pas de `nozzle_diameter` — même piège que les
  limites #26 et #21.
- Ne PAS exclure la géométrie au-delà de cette zone : c'est précisément le
  corps de la tête (qui monte 77 mm depuis la buse, voir
  `mesh_loader.TOOLHEAD_BODY_TO_MACHINE_UP`) qui doit rester testé.

**En attendant**, `scripts/export_slicing_preview.py` produit le GLB de
visualisation en court-circuitant explicitement step9/step10 (step1→step8
en ~30 s), avec `CollisionReport(nozzle_model_checked=False)` et
`AdjustmentReport(converged=False)` — jamais un rapport vide qui
prétendrait qu'aucune collision n'existe.

**Second point, indépendant et non résolu** : le débit de step9 se dégrade
avec la matière accumulée (mesuré : 16.5 → 3.7 pts/s entre le 250ᵉ et le
1750ᵉ point échantillonné, index passant de quelques centaines à ~2700
capsules). La grille spatiale de `_PrintedGeometryIndex` limite le nombre de
candidats mais pas assez ; à regarder après la zone de garde, une fois que
le nombre de collisions réelles sera revenu à une valeur plausible.

Ce qui reste, hors cœur algorithmique :

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
