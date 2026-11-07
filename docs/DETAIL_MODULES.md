# Détail des modules

_Extrait de l'ancien CLAUDE.md (R-030 point 7)._

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

## Détail de `element_ardor/step4_field_solver.py`

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

## Détail de `element_ardor/step1_ingest.py`

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

## Détail de `element_ardor/step6_toolpath2d.py`

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

## Détail de `element_ardor/step7_scheduling.py`

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

## Détail de `element_ardor/step8_machine_toolpath.py`

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

## Détail de `element_ardor/step9_collision_sim.py`

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

## Détail de `element_ardor/step10_adjustment.py`

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
  vectorisée. Les deux corrigés dans `element_ardor/step9_collision_sim.py`.

