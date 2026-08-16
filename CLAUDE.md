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
├── main.py                          # PAS ENCORE ÉCRIT — orchestrateur, enchaîne les étapes
├── config/                          # PAS ENCORE ÉCRIT
│   ├── params.py                    #   dataclasses de tous les paramètres globaux
│   └── schema.py                    #   (de)sérialisation JSON <-> dataclasses, pont avec l'UI TS
├── pipeline/
│   ├── __init__.py
│   ├── step1_ingest.py              # FAIT — ingestion mesh + tétraédrisation (+ détection
│   │                                 #   automatique des cavités internes pour tetgen, voir limite #7)
│   ├── step2_overhang.py            # FAIT — détection surplombs + ponts internes
│   ├── step3_segmentation.py        # FAIT — tests pytest formalisés (9 tests)
│   ├── step4_field_solver.py        # FAIT — résolution FEM du champ scalaire (7 tests)
│   ├── step5_layers.py              # FAIT — extraction des couches courbes (12 tests)
│   ├── step6_toolpath2d.py          # FAIT — remplissage 2D par couche (8 tests)
│   ├── step7_scheduling.py          # FAIT — couleur + ordonnancement trous/purges (19 tests)
│   ├── step8_machine_toolpath.py    # PAS ENCORE ÉCRIT
│   ├── step9_collision_sim.py       # PAS ENCORE ÉCRIT
│   └── step10_adjustment.py         # PAS ENCORE ÉCRIT
├── core/
│   ├── __init__.py
│   ├── mesh_types.py                # FAIT — toutes les dataclasses partagées (voir détail plus bas)
│   ├── geometry_utils.py            # FAIT — build_face_adjacency, fibonacci_sphere_directions,
│   │                                 #   rotation_axis_circle_directions, composed_rotation_directions
│   │                                 #   (N axes composes, cinematique generique), point_in_polygon_2d
│   └── constraints.py               # FAIT — registre extensible d'assembleurs de contraintes
│                                     #   utilisateur (ATTACH_POINT/AESTHETIC/STRUCTURAL), point
│                                     #   d'extension pour tout nouveau type sans toucher step4
├── io/                               # PAS ENCORE ÉCRIT
│   ├── mesh_loader.py
│   └── gcode_export.py
└── tests/
    ├── __init__.py
    ├── test_step2_overhang.py       # FAIT — 6 tests, tous passants et stables
    ├── test_step3_segmentation.py   # FAIT — 9 tests
    ├── test_step4_field_solver.py   # FAIT — 7 tests
    ├── test_step5_layers.py         # FAIT — 12 tests
    ├── test_step6_toolpath2d.py     # FAIT — 8 tests
    ├── test_step7_scheduling.py     # FAIT — 19 tests
    └── test_geometry_utils.py       # FAIT — 6 tests (composed_rotation_directions)
```

**État global : 67/67 tests passants** (step2 à step7 + core/geometry_utils.py).
Étapes 1 à 7 du pipeline stables et testées ; step8 (trajectoire machine)
est la prochaine à écrire.

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
  // HPI/HPE n'est qu'un cas d'usage parmi d'autres) et #19 (singularité au
  // centre du premier axe de rotation, si applicable).
  rotation_axes: RotationAxisSpec[] {
    axis: Vec3                        // unitaire, dans le repère du DOF parent
    angle_range_deg: {min: float, max: float}  // (0,360) tour complet, plage
                                                 // restreinte, ou (v,v) fige
  }
  // Exemple d'usage — machine HPI/HPE actuelle (commande firmware, gabarit
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
  // (gabarit générique, exemple HPI/HPE, limites #18/#19). Rappel synthétique :
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

**Ce que le code fait précisément :**
- Charge `MachineProfile.nozzle_model` (mesh de la tête complète : buse + porte-buse + partie basse du carriage, dans son repère local outil).
- Pour chaque point échantillonné de `ToolpathMachine` (pas = `path_sampling_step`) :
  - Positionne/oriente une copie transformée de `nozzle_model` selon la position + rotation courante de la trajectoire.
  - Teste l'intersection de ce mesh positionné contre : (a) le mesh de la pièce **déjà imprimée jusqu'à ce point** (pas le mesh final complet — état incrémental), (b) le mandrin/mors de la machine, (c) le plateau lui-même selon son angle courant.
  - Teste avec une marge `safety_margin` (dilate légèrement le mesh de la tête ou réduit la distance seuil d'intersection).
- Chaque collision détectée est enregistrée avec : point de trajectoire concerné, type d'obstacle touché, profondeur de pénétration.

**Sortie :** `CollisionReport[]` (liste de collisions, vide si aucune).


## Étape 10 — Ajustement en cas de collision détectée

**Ce que le code fait précisément (boucle itérative, jusqu'à `max_adjustment_iterations`) :**
- Pour chaque collision reportée, tente dans l'ordre :
  1. **Réorientation locale** : ajuste l'angle de rotation plateau à ce point précis (recherche locale dans une plage d'angles voisins qui élimine la collision tout en respectant encore le seuil de surplomb).
  2. **Ré-ordonnancement de région** : si la réorientation locale est impossible sans casser la contrainte de surplomb, retourne à l'étape 3 avec une contrainte supplémentaire "cette région ne peut pas utiliser cette plage d'angles" et relance la segmentation pour cette zone uniquement (pas tout le modèle).
  3. **Marquage support de secours** : si aucune des deux options ne résout la collision (cas limite), marque la zone comme nécessitant un support classique en dernier recours, avec log explicite pour l'utilisateur.
- Après ajustement, relance l'étape 9 uniquement sur la portion de trajectoire modifiée (pas la simulation complète), jusqu'à convergence ou épuisement des itérations.

**Sortie finale :** `ToolpathMachine[]` validé + `AdjustmentLog[]` (zones modifiées, méthode utilisée, éventuels supports de secours ajoutés) — transmis à l'UI TypeScript pour affichage/validation avant export G-code.

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
    machine utilisée. **Résultat mesuré (pas supposé) pour la machine HPI/HPE
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
- `ToolpathPoint` (**`radius`/`platform_angle`/`head_angle`** — (R, C, B)
  explicites, PAS un `rotation_angle` unique, structurellement insuffisant
  pour deux rotations indépendantes), `ToolpathMachine` — étape 8
- `CollisionEvent`, `CollisionReport` — étape 9
- `AdjustmentLogEntry` — étape 10

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

---

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

Étapes 1 à 7 stables et testées (67/67 tests). `MachineProfile` décrit
maintenant une cinématique générique à N axes de rotation composés
(`rotation_axes: list[RotationAxisSpec]`, voir limites #18/#19) ; step3 et
step7 dispatchent tous deux sur `composed_rotation_directions` sans
hypothèse figée sur le nombre d'axes de la machine réellement utilisée.

`pipeline/step8_machine_toolpath.py` lui-même n'est pas encore écrit, ni
planifié en détail : consommer `PrintSchedule.passes` **dans l'ordre**
(aucune logique d'ordonnancement à re-dériver, gelée à l'étape 7), relever
chaque segment `(layer_id, index)` du repère local 2D vers le cartésien
machine via `Layer.local_frame` (`step6._lift_to_3d` fait déjà la moitié du
travail), PUIS résoudre position/orientation désirée en valeurs concrètes
d'angles pour CHAQUE axe de `MachineProfile.rotation_axes` (pas seulement
2, le cas HPI/HPE), en respectant la plage de chaque axe et la singularité
géométrique au centre du premier axe si applicable (limite #19), lisser la
série temporelle d'angles pour `max_angular_velocity`/`max_angular_acceleration`,
garantir la continuité d'extrusion. Un plan dédié à step8 doit encore
trancher l'algorithme précis de résolution multi-axes (notamment le cas où
la plage d'un axe interne est dépassée : réessai via un axe externe, ou
signal vers l'étape 10) — non fait ici, et à concevoir génériquement plutôt
que spécifiquement pour 2 axes.