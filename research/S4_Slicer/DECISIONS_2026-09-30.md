# Décisions prises seul — session autonome du 2026-09-30

Plan : `~/.claude/plans/snug-yawning-lovelace.md`. Une ligne par décision : quoi, pourquoi, mesure, réversible ?

## T1 — P5 (skin_print.py)
- `--sampling` vaut désormais 0 par défaut, ce qui teste TOUS les points. Pourquoi : un contact de
  0,02 mm sur le dino tombait entre deux échantillons de 2 mm. Réversible avec `--sampling=2`.
- La tête est aussi testée contre les lignes déjà déposées de la zone en cours et contre le début
  de la ligne en cours (arbres `zone_tree` et `partial_tree`), comme `full_step89 --fast`.
  Réversible.
- Sauvegarde avant le nouveau P5 du dino contraint : `bak_0930/`.
- **Résultat (dino contraint, passe 4)** : sur 83 528 points de lissage, 2 zones en contact.
  - zone 23 : 5 points, pénétration max 0,03 mm. C'est le point vu par `full_step89` et manqué par
    l'ancien P5.
  - zone 4 : 6 points, 0,05 mm. Ce contact n'était vu par AUCUN des deux contrôles (tous deux
    échantillonnent tous les 2 mm).

  Les deux zones sont retirées (1 869 mm² sur 2 993). Ce sont des contacts de quelques centièmes
  de millimètre, dans la marge de 1 mm.
- **Coût** : l'arbre du début de ligne, reconstruit à chaque point, est quadratique (une ligne du
  dino compte 19 566 points), soit 890 s. Avec des arbres par blocs de 256 points : 239 s, et des
  points de collision IDENTIQUES à la version quadratique (vérifié). L'ancien P5 prenait 24 s
  avec un point tous les 2 mm.
  **À faire** : ce surcoût (×10) pèse dans la boucle de rognage. Piste : ne tester que les points
  à moins de X mm d'une matière déposée (arbre des obstacles interrogé une fois par zone).

- Dans `skin_loop`, les P5 de la boucle échantillonnent tous les 2 mm ; seul le P5 final teste
  tous les points. Pourquoi : la boucle cherche où rogner, et c'est la vérification finale qui
  doit être exhaustive. Réversible.
- Attention : la chaîne Squirtle du 2026-09-30 a démarré AVANT ce réglage (P5 exhaustif à chaque
  tour).

## T5 — Refonte de P3
- Pendant la boucle, P5 teste contre le corps NON creusé (`--body-uncarved`), qui est un
  sur-ensemble de la matière, donc plus prudent. P3 ne tourne qu'une fois, à la fin, suivi d'un P5
  final sur le vrai corps creusé. Si ce P5 final trouve encore un contact, on applique le retrait
  de zones existant (P3 `--skin-exclude`, puis P5). L'ancien comportement reste disponible avec
  `--p3-each-turn`.
- Le creusement est parallélisé par couche (`ProcessPoolExecutor`, `--workers`, défaut
  min(10, cœurs − 2)). Chaque processus construit sa propre `ExactDistance`.
- **Oracle (Benchy, passe 4, 4 zones exclues)** : 313 couches, 0 couche de topologie différente,
  écart max des contours **0,00 mm** (identique au bit), par rapport à la sortie série de la
  veille (`bak_0930/benchy_p3_serial.pkl`).
  - Creusement : 286 couches en **68 s** sur 10 processus ; P3 complet en 140 s.
  - Pour comparaison en série : Squirtle, 218 couches en 368 s (1,7 s par couche) ; Benchy, 0,24 s
    par couche. Soit environ **×7** sur le creusement.

- **Profil du creusement** (Benchy, 12 couches, 14,0 s) :
  - `trimesh.proximity.closest_point` = **9,3 s (67 %)**, dont 4,2 s de requêtes rtree ;
  - `shapely.union_all` = 2,7 s (19 %) ;
  - le découpage des sections : moins de 1 s.

  → **Cible de T9 (C++, GPU) : la distance exacte point-maillage** (BVH de triangles, distance
  point-triangle), et non le découpage. L'union reste en Python (shapely).

- **Squirtle, passe 8 (première chaîne complète avec la nouvelle boucle)** : 3 tours de P2, P4 et
  P5 (65 → 26 → 10 points de contact), puis P3 une seule fois (94 s, contre 468 s par tour
  auparavant). La boucle de rognage prend 1 206 s au total.
  **DÉFAUT trouvé** : le P5 final, exhaustif, voit encore 5 zones en contact (44 points), alors
  que leur creusement est fait. On obtient 106 mm² de corps creusé SANS couche de lissage.
  **Correctif** : `skin_loop.final_exclusion_rounds` cumule les exclusions et relance P3 et P5
  jusqu'à zéro contact, 3 tours au plus. Squirtle est réparé avec cette fonction.
- **P2 : le pré-calcul d'accessibilité est identique à chaque tour** (il ne dépend pas des
  rognages) et coûtait 291 s sur les 306 s de P2. Il est mis en cache
  (`_skin_access_cache_pass<k>.npz`, clé = empreinte SHA-256 de toutes ses entrées) et parallélisé
  par paquets de faces (`last_accessible_moment_parallel`). Le gain n'est pas encore mesuré.
- **Contacts du corps restants sur Squirtle** après la borne en hauteur et 4 tours de
  réorientation : 10 points, pénétration de 1,07 à 1,16 mm dans la marge (jeu réel de 0,04 à
  0,13 mm), z médian 3,7 mm, B médian 56°. Ils sont près du plateau, dans le corps : à regarder.
- **Export G-code de Squirtle** : 534 421 lignes, 3,9 m de filament, 1,39 h de temps nominal,
  4,92 h avec les limites de B (×3,5), 18 dégagements en contact franc.
  **À creuser** : ce ×3,5 dû à B est énorme ; ce sont probablement les raccords de lissage, où B
  change vite.

- **Squirtle réparé** : exclusion cumulée des zones [0, 2, 4, 8, 12, 15, 17], puis P3 et P5 →
  **0 contact** dans le lissage.
- **Exports G-code** (le temps « réel » tient compte des limites de B) :

  | modèle | lignes | filament | temps nominal | temps réel | dégagements en contact franc | aller-retour (écart max) |
  |---|---:|---:|---:|---:|---:|---:|
  | Benchy (passe 4) | 955 206 | 4,50 m | 1,63 h | 3,26 h | 6 | 1,1e-4 mm (0 point > 1e-3) |
  | dino (passe 5) | 1 034 994 | 6,40 m | 2,24 h | 4,33 h | 4 | 1,1e-4 mm |
  | Squirtle (passe 8) | 520 962 | 3,90 m | 1,39 h | 4,90 h | 20 | 1,1e-4 mm |

  L'aller-retour valide le pivot, les signes et le déroulement de Y (l'écart est l'arrondi à 4
  décimales). La pointe descend au plus bas à 0,04 à 0,09 mm, la première couche. B va jusqu'à
  111-115°.
- **Ralentissement dû à B, diagnostiqué** : sur Squirtle, le corps passe de 1,27 h à 4,67 h.
  5,9 % des segments du corps font sauter B de plus de 10° : B oscille entre 0° et environ 14°
  d'un point à l'autre, près de l'axe (R médian 12 mm). La normale est constante par tétraèdre et,
  C étant imposé, B bute à 0 quand elle penche vers le centre. L'hystérésis n'y fait presque rien
  (seuil 0,5° : 4,52 h). **Filtre médian de 9 points par polyligne : 2,89 h (−41 %)**, mais
  appliqué APRÈS la réorientation il la défait (écart jusqu'à 115°). D'où
  `full_step89 --smooth-b=9`, appliqué juste après la cinématique inverse, AVANT la borne en
  hauteur, la réorientation et le contrôle, avec l'écart tête / normale recalculé.
  **Mesure sur Squirtle** :
  - temps réel : **4,90 h → 2,97 h (−39 %)** ;
  - contacts : 9 non résolus (contre 10), même pénétration (1,07 à 1,16 mm, dans la marge) ;
  - corps : écart tête / normale médian 6,8°, > 45° pour 1,6 % des points ;
  - le premier tour de réorientation traite 1 488 points (le filtre en déplace beaucoup), et en
    résout 1 479.

  Porté dans Apex : `Step8Params.free_axis_median_window` (limite #63). `--smooth-b=9` est
  appliqué à tous les modèles de la chaîne suivante.
- matplotlib 3.11.2 installé dans le venv d'Apex (aperçu du G-code ; CLAUDE.md le disait déjà
  installé, ce qui n'était plus vrai pour le venv 3.14).

## T6 — Validation sur les autres modèles (chaîne `chain4.sh`, 1 h 25 pour tout)
Tous avec `--smooth-b=9`, borne en hauteur et réorientation, puis export G-code.

| modèle | couches (contours invalides) | lissage : zones gardées, 0 contact ? | corps > 45° | contacts du corps non résolus | temps nominal → réel |
|---|---|---|---:|---:|---|
| Benchy (passe 4) | 313 | 0 contact | 1,7 % | 9 (0,04 / 0,22 mm) | 1,63 → 2,46 h (sans filtre : 3,26 h) |
| dino (passe 5) | 297 | 0 contact | 0,1 % | **0** | 2,24 → 3,48 h (sans filtre : 4,33 h) |
| Squirtle (passe 8) | 218 (93 invalides) | 7 zones retirées, 0 contact | 1,6 % | 9 (jusqu'à 1,16 mm) | 1,39 → 2,97 h |
| death (passe 3) | 217 (**108 invalides**) | 1 zone de plus retirée au P5 final, 0 contact | 1,1 % | **0** | 1,17 → 2,39 h |
| real_overhang (`_a0.5_pin`, passe 11) | 207 (28 invalides) | 0 contact | **8,8 %** | **127** (médiane 0,69, max 1,17 mm ; z ≈ 30, B ≈ 102°) | 0,11 → 0,16 h |
| l-triple (passe 16) | 125 (0) | 0 contact | 0 % | **0** | 0,11 → 0,17 h |
| l_shape (passe 16) | 147 (0) | 0 contact | 1,8 % | **0** | 0,08 → 0,10 h |

- **Cache de P2** : sur death, P2 passe de 95 s (tour 1) à 40 s (tours 2 et 3) ; le reste de P2
  n'est pas en cache.
- **L'exclusion cumulée a servi dès le premier nouveau modèle** : sur death, le P5 final exhaustif
  a vu 4 points manqués par la boucle ; zone 5 retirée, puis 0 contact.
- **À traiter** :
  - real_overhang : 127 contacts non résolus. La tête est presque horizontale (B ≈ 102°) à
    mi-hauteur, et 23,5 % des points sont bornés par la hauteur. Le modèle est-il seulement
    imprimable avec la tête actuelle ?
  - Contours invalides (auto-intersections) sur Squirtle, death et real_overhang : replis Bird
    probables, que l'étape 6 absorbe mal ; à mesurer.
  - Squirtle : 9 contacts du corps à 1,1 mm de pénétration, près du plateau (z ≈ 3,7 mm).

## Rangement (demande de l'utilisateur)
- **apex-dev** : la poussée de TOUT le répertoire (fichiers ignorés compris) a été BLOQUÉE par le
  contrôle de sécurité du mode Auto (« exfiltration »). L'utilisateur a créé le dépôt privé, vide.
  Rien n'a été poussé ; la décision lui revient.
- **R: est un disque de 8 Go**, plein au premier essai de préparation du dépôt (dépôt partiel
  supprimé).
- **Sorties obsolètes ARCHIVÉES** (déplacées, pas supprimées) vers `H:/apex_archive_2026-09-30/`
  (718 Mo, 1 002 fichiers et 2 dossiers ; liste dans `LISTE_ARCHIVE.txt`) :
  - les STL de tentatives de la boucle Bird ;
  - les variantes d'expériences conclues : w = 0,5, 2 et 10, décalages ±70 mm, `pinc`, `a0.5`,
    `folds`, l'ancien `benchy_upsidedown`, les variantes sans `_pin` ;
  - `patterns_bench.pkl` ;
  - mes sauvegardes du jour, le cache CMake périmé et les anciens modules C++.

  Gardés : les `_pin` de chaque modèle, les `_pin_reach` du Benchy et du dino (hors décalages),
  les sorties finales.

## T9 — Distance point-maillage en C++ (`apex_core.MeshDistance`)
- `cpp/include/apex/mesh_distance.hpp` et `src/core/mesh_distance.cpp` : BVH de triangles,
  descente la plus proche d'abord avec élagage par la distance à la boîte, point le plus proche
  d'Ericson. Le cœur ne dépend pas de pybind11 (`MeshView` = pointeurs bruts) ; le binding
  `MeshDistance` n'est qu'un adaptateur.
- **Backends** :
  - 0 : série (oracle) ;
  - 1 : threads ;
  - 2 : Vulkan, qui se rabat sur 1 avec un motif (`mesh_distance_backend`). Le backend Vulkan
    n'est pas écrit.
- **Oracle** : trimesh n'est PAS exact près des sommets (jusqu'à 1,9e-5 mm, 15 points sur 5 000,
  toujours par excès). Le natif coïncide avec la force brute à 1,7e-16 près. Le test compare donc
  à trimesh là où les deux s'accordent, et à la force brute ailleurs.
- **Compilation** : deux obstacles, contournés sans rien supprimer.
  - Le lanceur `Scripts/cmake.exe` du venv échoue en silence. `build_cpp.py` appelle désormais le
    binaire du paquet (`cmake.CMAKE_BIN_DIR`).
  - Le cache CMake datait d'avant le déplacement du projet de G: vers R:. Il est mis de côté dans
    `cpp/build_old_G/`.
  - L'ancien module est gardé : `cpp/lib/*.pyd.bak0930`.
- **Gain sur le creusement (Benchy, 12 couches)** : 13,9 s → 4,5 s (×3,1). Les régions creusées
  sont IDENTIQUES (différence symétrique 0 mm²), grâce à l'arrondi à 1e-6 de `carved_region`.
  Activé par défaut dans `skin_carve.ExactDistance` quand le module le fournit
  (`APEX_P3_NATIVE_DISTANCE=0` rend la voie trimesh), avec 1 thread par processus dans le
  creusement parallèle.
- **Vulkan** :
  - présent : le runtime (`vulkan-1.dll`) et une RTX 3090 ;
  - absent : le SDK Vulkan (en-têtes, `glslc` pour compiler les shaders en SPIR-V). Son
    installation demande probablement une action de l'utilisateur (droits d'administrateur) ;
  - le backend n'est donc pas écrit : `MeshDistance(backend=2)` se rabat sur les threads CPU,
    avec son motif.

  **À décider** :
  - installer le SDK (LunarG, licence Apache-2.0 à vérifier à la source) ;
  - ou passer par CUDA, déjà installé (v13.2 dans le PATH), mais réservé à NVIDIA ;
  - noter que la RTX 3090 calcule en double précision à 1/64 du débit simple : la précision
    simple casserait l'égalité au bit avec le CPU.
- Ce qui reste : l'union shapely (2,7 s sur 12 couches). Elle n'est pas portée : c'est la
  géométrie 2D robuste (limite du projet : shapely/Clipper).

## T2 — Export G-code (file_io/gcode_export.py, Apex)
- **Firmware Y = atan2(y, x)**, soit `C − 90` avec la convention d'Apex. Y est recalculé depuis la
  position et comparé à `axis_angles[0]` (`c_mismatch_deg`). Ce fait machine est relevé dans les
  scripts `gen_*` et vérifié avec l'utilisateur ; aucun code n'en est repris.
- **Pivot** :
  - référence prise : à B = 0, les coordonnées firmware sont celles de la pointe (calibration
    buse verticale) ;
  - on a donc X = R + 45 sin B et Z = z + 45 (cos B − 1) ;
  - **À CONFIRMER par l'utilisateur** : si le zéro Z du firmware est le pivot et non la pointe à
    B = 0, il faut ajouter 45 mm à Z partout. Réversible, une ligne.
- Le signe du pivot est vérifié par un test contre `step9._positioned_toolhead_vertices`, donc la
  même convention que le contrôle de collision.
- **F** = vitesse de pointe × distance firmware / distance de pointe (le firmware traite Y et B
  comme linéaires). Pour une rotation pure, on utilise `rotation_feed` = 3000. Plafond 12000.
  Valeurs provisoires.
- **Dégagement** :
  - Z si |cos B| ≥ |sin B| (+Z pour une tête vers le bas, −Z pour une tête retournée), sinon R dans
    le sens du corps de la tête ;
  - si le mouvement descend sous la marge du plateau (point le plus bas de la vraie tête), on
    prend l'autre axe ;
  - le trajet entre les deux points dégagés, lui, reste multi-axes : « un seul axe à la fois »
    appliqué au mouvement de dégagement.
- **Contrôle des dégagements** : 16 points sur l'axe de la tête (de la pointe au pivot) contre la
  matière déjà déposée, contact franc si la distance est inférieure au rayon du cordon (0,2 mm).
  La tête complète n'est pas testée (coût) : c'est une limite.
- **E** absolu, avec section = hauteur × largeur de ligne (0,2 × 0,4) × flux. Lissage :
  `skin_thickness` = 0,2.
- **Vitesses provisoires** : impression 25 mm/s, trajets à vide 80 mm/s, température 195 °C
  (celle de l'en-tête connu).

## T4 — P2 (surface_skin.py)
- L'écart tête / normale est jugé au PIRE des trois sommets et du centre de la face, et non plus au
  seul centre. Coût : IK × 4.
- **Mesure (Benchy, passe 4, même fichier de rognage)** :
  - faces « inatteignables » (> 45°) : 7 787 → **7 925** (+138) ;
  - faces retenues pour le lissage : 4 488 → **4 316** (−172, −3,8 %) ;
  - le reste est quasi inchangé.

  Les sorties du Benchy ont été restaurées ensuite : elles restent celles du jugement au centre,
  pour la cohérence de la chaîne. Le nouveau jugement s'applique aux prochaines passes de P2
  (Squirtle, death et les autres modèles l'ont déjà eu dans `chain4`).
- **P2 parallélisé, effet mesuré sur le Benchy** : le pré-calcul de collision prend **18,2 s**
  (environ 200 s en série dans la v1), et P2 complet 38,6 s (IK × 4 comprise : 17,8 s).
