# Session autonome du 2026-09-30 (5 h) — export G-code, fiabilisation, refonte de P3

## Contexte

Maker Faire Rome (fin octobre) demande un G-code réel et une trajectoire sans collision. Les
trajectoires corrigées du Benchy et du dino existent (`loop_*_machine_toolpath_pass*.pkl` : position,
axis_angles, extrusion, is_skin), mais aucun exporteur n'existe. La fiabilisation (autres modèles,
défauts connus, portage dans Apex) n'a été mesurée que sur 2 modèles.

**Volume** : environ 16 h de travail estimé, soit ~3,2 fois les 5 h, comme demandé. Liste
ordonnée par priorité : on avance dans l'ordre, et ce qui n'est pas fait reste pour la suite.

## Règles de la session

- Mode Auto, pas de retour au mode Plan. AskUserQuestion seulement pour un point qui régit le cœur
  du slicer. Sinon je tranche seul et je note chaque choix dans
  `research/S4_Slicer/DECISIONS_2026-09-30.md` (quoi, pourquoi, mesure, réversible ou non).
- Aucun commit ni push. Le code de recherche reste local. Apex peut être modifié
  (`file_io/`, `pipeline/`, `tests/`) ; la suite doit passer dans les deux positions de
  `APEX_USE_CPP`.
- Un seul calcul lourd à la fois, `python -u`, `timeout`, `check_memory`. Vérifier en fin de tâche
  qu'aucun Python orphelin ne tourne.
- Ne jamais s'inspirer des scripts `gen_*_gcode.py` (générateurs de formes fixes).
- Pendant qu'un calcul tourne en arrière-plan, j'avance sur une tâche de code qui n'en dépend pas.

## Décisions de l'utilisateur (2026-09-30)

- **Rétraction** : 0,8 mm à 35 mm/s. C'est un paramètre, et le réamorçage est identique.
- **Dégagement sur les trajets à vide** : dynamique, selon l'orientation et la position de la tête.
  On agit sur l'axe le plus approprié, **un seul à la fois**.
- **Refonte de P3** approuvée :
  - parallélisée massivement ;
  - en tout dernier, s'il reste du temps : C++ et Vulkan.
- **P5** : tester tous les points, et aussi contre les lignes déjà déposées de la zone en cours.

## Tâches, par ordre de priorité (estimation réaliste)

### T1 — P5 : les deux correctifs (30 min)
Dans `research/S4_Slicer/skin_print.py` :
- tester chaque point du lissage, et non un point tous les 2 mm ;
- ajouter un arbre de la zone en cours : ses lignes déjà déposées et le début de la ligne en cours,
  comme `full_step89.py`.

Vérification : relancer P5 sur `dino_pin_reach` passe 4. Le point de la zone 23 doit apparaître
(0,02 mm de pénétration), et ce contact doit être le même que celui vu par `full_step89 --fast`.
Mesurer la durée.

### T2 — Export G-code `file_io/gcode_export.py`, dans Apex (3 h)
Écrit de zéro. Entrée : le dictionnaire `_machine_toolpath_pass*.pkl`, avec une fonction
`export_gcode(toolpath, machine, params, path)` et un petit script d'appel en recherche.

1. **`GcodeParams`**, dataclass dont chaque réglage est un paramètre :
   - `nozzle_temp=195` ;
   - `filament_diameter=1.75` ;
   - `layer_height` et `line_width` (pour E) ;
   - `print_feed` et `travel_feed` (mm/s) ;
   - `retract_length=0.8` et `retract_speed=35` ;
   - `travel_clearance=1.0` mm ;
   - `pivot_to_tip=45.0` mm ;
   - `travel_retract_min=2.0` mm.
2. **Traduction des axes** : C = axis_angles[0] (Y firmware, en degrés), R = rayon (X firmware),
   Z, B = axis_angles[1].
   - **Pivot** : le firmware commande le pivot de B, pas la pointe. On a donc
     position du pivot = pointe + 45 mm × axe de la tête, dans le plan (R, Z), avec le signe de la
     convention « B > 0 penche la buse vers −R ».
   - Contrôle : un point imprimé à B = 0 donne Z_pivot = Z + 45 ; B = 90 donne un décalage en R.
   - **Le signe de B est à vérifier contre `step9._positioned_toolhead_vertices`** : c'est la
     même source que le contrôle de collision, donc aucune convention nouvelle.
   - C est déroulé (pas de saut de 360°).
3. **Extrusion** : E absolu (`M82`), E = longueur × section du cordon ÷ section du filament. Pour
   le lissage, on utilise l'épaisseur de lissage (0,2 mm).
4. **Avance** :
   - F est recalculé par segment, car le firmware traite Y (degrés) comme un axe linéaire. On
     garde la vitesse de la pointe constante en calculant F sur la « distance » firmware
     (ΔX, ΔY en degrés, ΔZ, ΔB).
   - Le ralentissement aux limites de B vient de `skin_print.feed_plan`, porté comme fonction dans
     `gcode_export`.
5. **Trajets à vide (dégagement dynamique, un axe à la fois)** :
   1. rétraction ;
   2. mouvement de dégagement sur **un seul axe** : Z si la tête est plutôt verticale (|cos B| ≥
      |sin B|), dans le sens qui éloigne la buse de la matière (+Z pour une tête vers le bas, −Z
      pour une tête retournée), sinon R dans le sens de l'axe de la tête. Il est refusé s'il fait
      descendre la tête sous la marge du plateau (table `zlow` de `full_step89`) ; dans ce cas on
      prend l'autre axe ;
   3. le trajet lui-même (C, B, puis le reste) ;
   4. le mouvement inverse sur le même axe ;
   5. le réamorçage.

   Les points du dégagement sont vérifiés par le contrôle rapide (arbre k-d). Une collision est
   signalée et comptée, jamais cachée.
6. **En-tête et pied de fichier** : `G21 G90 M82 G92 E0 M109 S<t>` ; en fin de fichier, rétraction,
   dégagement, puis `M104 S0`. Chaque changement de passe (corps ou lissage) est commenté.
7. **Tests** (`tests/test_gcode_export.py`) : traduction des axes et du pivot sur des cas connus,
   E cumulatif, rétraction et dégagement sur un seul axe, déroulement de C, fichier bien formé.
   Suite complète dans les deux positions de `APEX_USE_CPP`.
8. Export du Benchy et du dino. Rapport : nombre de lignes, durée estimée (avec `feed_plan`),
   filament utilisé, collisions de dégagement.

### T3 — Aperçu du G-code (1 h)
`scripts/preview_gcode_file.py` (nouveau, dans Apex, sans s'inspirer de `preview_gcode.py` non
versionné sans l'avoir relu) :
- relire le fichier et rejouer la cinématique inverse : pivot vers pointe, puis (R, C, Z) vers x, y, z ;
- rendu des trajectoires, extrusions colorées et trajets à vide en gris, dans le repère de la
  pièce ;
- **contrôle aller-retour** : les positions relues doivent retrouver les positions d'origine à
  1e-3 mm près. C'est le test qui valide le pivot et les signes.

### T4 — Défaut « au-delà de 45° jugé au centre des faces » (1 h + calcul)
Dans P2 (`surface_skin.py`), l'accessibilité est jugée au centre des faces. Il faut la juger aux
sommets, ou au pire point de la face, avec le même critère de tolérance. Vérification sur le
Benchy : faces refusées avant et après, et points au-delà de 45° en P5 qui doivent tomber à 0.

### T5 — Refonte de P3 : corps creusé une seule fois, parallélisé (4 h)
1. Dans `skin_loop.py`, P3 ne tourne plus à chaque tour. La boucle alterne P2, P4 et P5 sur une
   **version rapide du creusement**, et le P3 complet ne tourne qu'**une fois, à la fin**.
   - Le contrôle de P5 n'a besoin que du corps imprimé jusqu'au moment de chaque zone.
   - Il faut donc vérifier dans le code ce que P5 lit du corps creusé. S'il en a besoin, le
     remplacer par le corps NON creusé moins une bande de l'épaisseur du lissage : c'est un
     critère plus prudent, qu'il faut mesurer.
2. P3 lui-même est parallélisé par couche, avec `concurrent.futures.ProcessPoolExecutor` : les
   couches sont indépendantes une fois ψ calculé. Les tâches sont groupées par lots pour limiter
   la sérialisation, et on surveille la mémoire (plantages passés).
3. **Oracle** : l'ancien P3 sur le Benchy et le dino.
   - Les contours doivent être identiques à 1e-6 mm près, couche par couche.
   - Même résultat final de la boucle : zones retenues et zéro collision.
   - Durée avant et après, et facteur par rapport au nombre de cœurs.

### T6 — Valider sur les autres modèles (calcul en arrière-plan, environ 45 min par modèle)
Chaîne complète : tranchage, `skin_loop`, puis `full_step89 --fast --height-clamp --apply`, puis
l'export G-code (T2) quand il existe.
- **Ordre** : Squirtle (il faut relancer, P3 a été tué au tour 3), puis real_overhang,
  l-triple, l_shape et death. turkish seulement si le temps le permet.
- Les modèles sans passe `_pin` utilisent la meilleure variante disponible, notée dans
  DECISIONS.
- **Mesures par modèle** :
  - contours invalides (Squirtle : 93 sur 218, à comprendre : replis Bird ?) ;
  - zones de lissage retenues ;
  - contacts avant et après correction, et paliers de tolérance ;
  - part des points au-delà de 45° ;
  - durée.
- Lancer dès T1 fini, puis enchaîner un modèle après l'autre pendant T2 à T5, avec le nouveau P3
  dès qu'il est validé.

### T7 — Porter dans Apex le contrôle rapide et les deux corrections (3 h)
Porter le contrôle rapide de step9 (arbre k-d sur la matière densifiée, arbre du début de passe),
la borne de B selon la hauteur et la réorientation par paliers. Le tout en opt-in
(`CollisionSimParams.fast_kdtree`, etc.), désactivé par défaut tant qu'il n'est pas mesuré sur les
quatre modèles de `res/test_models/`.
- Réutiliser `step9._positioned_toolhead_vertices`, `_surface_sample_points`,
  `_filter_nozzle_contact_zone` et la logique de `full_step89.py`.
- Tests : même verdict que la boucle lente de step9 sur un petit cas, borne en hauteur, et
  réorientation qui résout un contact synthétique.
- Mettre à jour CLAUDE.md : nouvelle limite, et limite #27 mise à jour.

### T8 — Tenue du projet (30 min, en fin de session quoi qu'il arrive)
- CLAUDE.md :
  - corriger la limite #59 : Bird ne convient qu'en partie à la machine actuelle ;
  - reporter les résultats du 2026-09-29 et du 2026-09-30 (axe D, plateau, contraintes,
    collisions, export) ;
  - mettre à jour le compteur de tests.
- `LAST_SESSION.md` : section « REPRENDRE ICI » avec ce qui est fait, les chiffres, les décisions
  prises seul et les questions pour l'utilisateur.
- Aucun commit sans demande.

### T9 — Si tout le reste est fini : P3 en C++, puis Vulkan (le reste du temps)

**Contrat d'interface (demande de l'utilisateur)** :
- **Pour l'instant**, le module natif prend en entrée les sorties du code Python existant, sous
  forme de tableaux numpy passés par pybind11 :
  - nœuds d'origine et déformés, tétraèdres, ψ aux nœuds, isovaleurs des couches ;
  - en sortie, les segments de coupe par couche.

  Python garde l'union et le chaînage (shapely), puis l'étape 6.
- **Il doit aussi poser les bases pour la suite.** Le cœur C++ ne dépend jamais de pybind11 : il
  reçoit des vues simples (`std::span` ou pointeur plus taille, structure `TetSliceInput` dans
  `cpp/include/apex/`). Le jour où tout sera porté, les étapes C++ amont rempliront directement
  cette structure, et le binding ne sera plus qu'un adaptateur à retirer, comme le reste de
  l'échafaudage (limite #56).
- **Le GPU est accepté dès maintenant, mais optionnel.** Un même point d'entrée
  `slice_tets(input, backend)` propose trois backends :
  - `cpu_serial`, qui sert d'oracle ;
  - `cpu_parallel` (threads) ;
  - `vulkan`.

  Au lancement, on détecte si Vulkan est disponible ; sinon, on se rabat en silence sur le CPU
  (même discipline que `APEX_USE_CPP` : l'absence de GPU n'est jamais une erreur). Le shader de
  compute prend les mêmes tampons que le CPU (nœuds, tétraèdres, ψ, isovaleurs), pour qu'un
  futur appelant C++ n'ait rien à changer.
- **Ordre** :
  1. C++ série, avec des tests d'équivalence contre le Python ;
  2. threads C++ ;
  3. Vulkan : vérifier la licence du SDK et des chargeurs (Apache-2.0 attendu, à vérifier à la
     source), puis la compilation GLSL→SPIR-V au build.

  La précision simple du GPU risque de casser l'équivalence au bit. Il faut la mesurer et la
  déclarer comme tolérance, comme aux pas 7 et 9. Si le matériel ne gère pas le double, le noter.
- Si l'installation du SDK Vulkan demande une action de l'utilisateur (droits d'administrateur),
  ne pas bloquer : garder le backend Vulkan compilable en option (`APEX_WITH_VULKAN`), le noter
  dans DECISIONS et passer à autre chose.

## Enchaînement automatique

Personne ne sera disponible. Tout s'enchaîne sans attendre de validation :
- **Aucun point d'arrêt utilisateur** : une question qui surgit est tranchée par moi, puis notée
  (DECISIONS), avec une solution réversible chaque fois que possible.
- **Chaque calcul lourd est lancé en arrière-plan** par un script de chaîne (`scratchpad/*.sh`)
  qui enchaîne les modèles de T6 un par un. Je reprends la main à chaque notification de fin, et
  je code T2 à T5 et T7 entre-temps.
- **En cas d'échec d'une tâche** (erreur non résolue en 30 min) : noter l'état et la cause, puis
  passer à la tâche suivante. Ne jamais rester bloqué.
- **T8 (compte rendu) est fait avant la fin des 5 h, quoi qu'il arrive.**

## Fichiers critiques

- **Recherche** :
  - `research/S4_Slicer/skin_print.py` (P5, `feed_plan`) ;
  - `skin_loop.py` ;
  - `bird_slice_deformed.py` (P3) ;
  - `surface_skin.py` (P2) ;
  - `full_step89.py` (contrôle rapide, borne, réorientation, `--apply`).
- **Apex** :
  - `file_io/gcode_export.py` (nouveau) ;
  - `tests/test_gcode_export.py` (nouveau) ;
  - `scripts/preview_gcode_file.py` (nouveau) ;
  - `pipeline/step9_collision_sim.py` (T7) ;
  - `main._default_machine_profile`.

## Vérification d'ensemble

- Suite pytest complète verte avec `APEX_USE_CPP=0` puis `=1`, après T2 et après T7.
- Aller-retour du G-code (T3) à 1e-3 mm près sur le Benchy et le dino.
- Oracle de P3 (T5) : contours identiques, et gain de temps mesuré.
- Chaque résultat est chiffré dans `DECISIONS_2026-09-30.md` puis dans `LAST_SESSION.md`.
