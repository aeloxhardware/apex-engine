# DT — Portage de la voie Bird en C++ (Apex Slicer)

État au **2026-10-13, soir** (branche `session-2026-10-13`, commit `6da758cd`).
Convention : 🟢 fait · ⏳ en cours · ⬛ non traité · ♦️ décision/bloquant.
Sources : `CLAUDE.md` #75 à #77, `DT.md`, `cpp/`, `pipeline/bird/`.
Les durées « à venir » sont des **estimations non mesurées**.

## 1. Contexte

- La voie Bird est l'algorithme de déformation du maillage tétraédrique de S4 Slicer
  (Joshua Bird, GPL-3.0). Elle sert à trancher en couches non planaires.
- Elle est intégrée à Apex sous `pipeline/bird/` (phase 1, #75), lancée par
  `main.py --method bird`. La voie FEM reste séparée et inchangée.
- Objectif de l'utilisateur : porter Bird et les modules R&D en C++, les lier au cœur d'Apex,
  puis ajouter un backend Vulkan (RTX 3090 détectée).
- Plan en quatre phases :

| Phase | Contenu | État |
|---|---|---|
| 0 | Mesure (profil par phase) | 🟢 |
| 1 | Python versionné = référence | 🟢 |
| 2 | C++ guidé par le profil | ⏳ |
| 3 | Mesure Bird contre FEM sur les 4 modèles, choix du défaut | ⏳ mesurée, défaut non tranché |

- Étape −1 (2026-10-13) : la production des G-code Bird passe de la chaîne de recherche
  (`research/S4_Slicer`, non versionnée) à `main.py --method bird`.

## 2. Ce qui est porté en C++

Tout est dans `cpp/include/apex/` et `cpp/src/core/`, appelé par `APEX_USE_CPP=1`.
Python reste le défaut et l'oracle. L'équivalence est vérifiée par `tests/test_cpp_equivalence.py`.

| Poste | Fichier C++ | Équivalence | Gain mesuré |
|---|---|---|---|
| Moindres carrés de la déformation et des rotations (`trf_no_bounds`, `lsmr`, QR Householder via Eigen, quartique) | `bird_lsq.cpp` | **Tolérance déclarée** : écart à scipy médian 0,15 µm, max 25 µm pour ~1 mm de déplacement ; coût final à 5e-5 relatif | Passe Benchy 62,6 → 9,9 s (×6,3) ; boucle complète ~12 min → 111,5 s |
| Solveur LSQ **parallèle** (pool persistant, réductions séquentielles) | `bird_lsq.cpp`, `internal/parallel_pool.hpp` | **Au bit** | Déformation Benchy 44,6 → 25,8 s |
| Viterbi de `map_polyline` | `apex_core.viterbi_choice` | **Au bit** | P3 : 123 → 61 s |
| Contrôle P5 (zones exclues) : grille incrémentale, y compris la particularité des blocs de 256 points | `skin_check.cpp` | Mêmes zones, mêmes points ; pénétration à 6e-15 mm | 31,3 → 2,4 s (un point sur 2 mm) ; exhaustif 161 → 39 s |
| P5 **parallèle** (`check_zone_parallel`, matière visible = préfixe de la zone) | `skin_check.cpp` | **Au bit** | Benchy 47,7 → 8,4 s |
| Barycentriques des candidats (LU à pivot partiel, ×20) | `bird_locate.cpp` | **Tolérance déclarée** (accord de l'utilisateur) : score ≤ 4,5e-9, point ≤ 1,2e-10 mm, 0 décision changée aux seuils ±1e-6 | 1,80 → 0,09 s ; sert le tranchage déformé, `locate` (étape 8) et `PrintedSolid.phi_at` (P2) |
| Distance point-maillage (creusement) | `mesh_distance.cpp` | Au bit | Déjà en service (pas T9) |
| Export G-code : découpage `max_tip_deviation` | `gcode_export.cpp` | **Texte identique ligne à ligne** (1 076 351 lignes) | 131,8 → 7,1 s |
| Contrôle des dégagements (`TravelCheck`, `PastMaterialIndex::any_past_within`) | `gcode_export.cpp`, `past_material.cpp` | Texte identique, y compris un cas construit où les dégagements sont bloqués | — |
| Socle Vulkan (volk, VMA, shaders, contexte) | `cpp/src/gpu/`, `cpp/shaders/` | `saxpy` identique au bit en float32 et float64 | Aucun noyau métier |

Gains hors C++, identiques au bit :
- densification du contrôle des dégagements, 9,7 → 0,5 s ;
- `locate` de l'étape 8 en un seul appel, 21,5 → 14,6 s ;
- raffinement du retour Bird par lots, 4 383 appels → 34 sur le dino ;
- cache de P4 par zone ;
- arrêt des paliers de réorientation dès qu'un tour ne résout rien.

Règle d'ex aequo de `locate` (`LOCATE_TIE = 1e-9`) : **changement de comportement
assumé**, valable pour les deux chemins (natif et LAPACK). Elle remplace un « premier maximum
strict » qui dépendait de l'arrondi. Effet mesuré contre `_v7` : 0,23 % des points à |ΔB| > 0,5°
(max 21°), collisions et longueur extrudée identiques.

### Temps de calcul (CPU seul, `main.py --method bird`)

| Modèle | Avant le chantier | Phase 0 (profil) | Après gains |
|---|---:|---:|---:|
| Benchy | 633 s (Python) → 381,7 s (solveur natif) | 300,6 s | **226,8 s** ; ~150 s depuis la passe 4 |
| dino | — | 373 s | **~256 s** depuis la passe 5 |
| death | — | 714 s | non remesuré |
| l_shape / l-triple | — | 29 s / ~40 s | — |

Comparaison à entrée identique (passes 4 et 5 de la recherche) : mêmes collisions par tour de
lissage, mêmes zones exclues, mêmes tours de Y et filament. Benchy 580 → 250 s, dino 1 067 → 383 s.
G-code de production : `out/share/*_v7.zip`.

## 3. Ce qui reste à faire

### 3.1 Portage C++ de la voie Bird (phase 2 : « reste le portage de P2 / P4 / étape 8 Bird »)

Les modules Python de `pipeline/bird/` totalisent 4 760 lignes. Seuls les points chauds ci-dessus
sont natifs. Le reste est du Python (numpy, shapely, scipy).

| # | Poste | Fichier Python | Poids mesuré | Remarque |
|---|---|---|---|---|
| 1 | ⬛ **Tranchage déformé : creusement P3** | `skin/skin_carve.py`, `deformed_slicing.py` | Dino : ~70 s non attribuées sur 99,6 s ; sérialisation vers les processus 8 s (Benchy) | **À instrumenter d'abord.** Creusement parallèle par processus. Turkish y manque de mémoire. |
| 2 | ⬛ P2 : dernier moment accessible par face | `skin/zones.py` | 15 s après le traitement par lots (Benchy 20,7 s ; death 55 s) | Faces indépendantes, bon candidat pour le multi-cœur. Aussi candidat GPU. |
| 3 | ⬛ P4 : chemins de lissage sur surfaces | `skin/surface_paths.py` (1 209 lignes), `surface_patterns.py` | Benchy 54 s (avec cache) | Dépend de shapely : portage coûteux. Clipper2 est déjà vendorisé. |
| 4 | ⬛ Étape 8 Bird : boucle d'objets Python | `machine_toolpath.py` | Benchy 21 s ; death 15,3 s | Après `locate` unique, reste la construction des objets. À passer en tableaux, comme step8 FEM. |
| 5 | ⬛ Dijkstra ×2 + topologie | `deformation.py` | death 63 s ; Benchy 17 s | |
| 6 | ⬛ Champ initial des rotations | `deformation.py` | ~9 s (Benchy), ~27 s (death) | Seule partie des rotations encore en Python. |
| 7 | ⬛ Variante « normale atteignable » du solveur | `deformation.py` | Non utilisée par défaut | À ne porter que si on l'active. |
| 8 | ⬛ Contrôle P5 sur un seul modèle de mémoire | `skin/print_order.py` | death 93 s en séquentiel | Le parallèle est fait. Il reste à le mesurer sur death et dino. |
| 9 | ⬛ `loop.py`, `fold_judge.py`, `chain.py`, `mesh_flatten.py` | — | Sans poids mesuré | Orchestration : à porter avec l'orchestrateur natif, pas avant. |

Ordre suggéré (partout : Python = oracle, port d'abord, parallélisation ensuite, tolérance déclarée
quand le bit est impossible) :
1. instrumenter le creusement P3 ;
2. P2 ;
3. étape 8 Bird en tableaux ;
4. Dijkstra ;
5. P4 en dernier.

Durée estimée : plusieurs jours à quelques semaines, P4 étant le plus long.

### 3.2 Backend Vulkan (phase 1 faite, noyaux non écrits)

- ⬛ Aucun noyau métier. Règle posée : **n'écrire un noyau que si le C++ multithread laisse
  plus de 1 s par appel.**
- Candidats : pré-calcul de collision de P2 (~30 M de tests point-dans-tétraèdre), distance à la
  peau par nœud, coupes de couches par tétraèdre.
- Hors GPU : affinage de maillage, chaînage des contours, décalages de polygones (CPU multi-cœur).
- Piège connu : les GPU grand public sont rapides en simple précision. La RTX 3090 expose le
  float64 (testé), mais plus lentement, et cela casse l'équivalence au bit.
- `MeshDistance` : aucun gain mesuré (2 M de points en 1,0 s sur CPU).

### 3.3 Licence et citation (⏳, demande de l'utilisateur)

- ⬛ Notice GPL-3.0 et attribution à Joshua Bird dans le code dérivé : `pipeline/bird/` l'indique,
  mais **pas le code C++ dérivé** (`bird_lsq`, `bird_locate`, `skin_check`).
- ⬛ `CITATION.cff` de S4 Slicer à conserver, citation dans le README et dans l'interface.
  Il y a un `LICENSE` et un `README.md` à la racine, mais **pas de `CITATION.cff`**.

### 3.4 Défauts connus à reprendre

- ⬛ Le contrôle des dégagements est grossier (hérité de la recherche) : 16 points sur 46 mm
  d'axe pour un rayon de 0,2 mm. Un cordon qui passe entre deux échantillons n'est pas vu. À
  remplacer par un test capsule axe / matière, avec mesure.
- ⬛ **turkish** : « non applicable » en phase 3. Un processus de creusement manque de
  mémoire (OpenBLAS) et le pool se bloque ; 547 304 triangles à la première couche.
- ⬛ Variance de mesure de ~10 % sur cette machine. Comparer en enchaînant les runs, jamais par
  soustraction entre sessions.
- ⬛ La règle du maillage grossier (`max_mean_tet_volume`) est désactivée. À 30 tétraèdres, Bird
  battait pourtant FEM sur l_shape. À revalider si le maillage grossier réapparaît.
- ⬛ Planification d'avance (`feed_plan`) de l'export non portée pour la voie Bird, si nécessaire.
- ⬛ Réorientation contre l'alternance de B (hors Bird au sens strict).

### 3.5 Décisions de l'utilisateur en attente

- ♦️ **Choix de la voie par défaut, Bird ou FEM** (phase 3, à placement égal) :

| Modèle | FEM | Bird | Impression Bird |
|---|---|---|---|
| l_shape | 3,1 % sans appui | 2,58 % | 0,09 h contre 0,2 h |
| l-triple | 5,5 % | 2,41 % | 0,13 h |
| death | 0,54 % | **0,42 %** | **1,50 h**, Y 121 tours |
| turkish | 1,10 %, **30 012 collisions** | — | non applicable |

  Bird imprime 27 à 53 % plus vite partout. Pour la métrique de référence (matière sans appui),
  il est légèrement meilleur.
  Les chiffres FEM de death et turkish des limites #64 à #74 décrivent une pièce hors plateau.
  Seule la mesure recentrée est comparable ; la pose à z = 0 à l'ingestion est faite depuis
  (#76).
- ♦️ Purge du `.docx` confidentiel dans l'historique git (`fdf3dbdc`), qui demande un push forcé.
- ♦️ `M203 Z8000` et `M201 Z1500` dans `config.g` (actions de l'utilisateur) : sans eux, le
  Benchy passe de 1,66 h à ~9,7 h d'impression.

### 3.6 Validation physique (bloque toute impression réelle)

- ♦️ Essai à sec `out/test_b_direction.gcode`, puis `out/test_l_shape.gcode`, puis le Benchy.
  Le sens physique de B n'est vérifié que sur l'aperçu.
- ⬛ Y cumulé à valider côté câblage : Benchy ≈ 31 tours, dino ≈ 620.

## 4. Ce que ce document ne tranche pas

- Les « gains Vulkan » n'ont aucun chiffre mesuré.
- Les temps de death sont ceux de la phase 0 (714 s), antérieurs aux gains parallèles.
- Le plan `~/.claude/plans/port-e-tout-bird-hidden-nygaard.md` cité par `CLAUDE.md` #77 n'est pas
  présent dans ce dossier de plans. Le séquencement détaillé de la phase 2 est donc reconstitué
  d'après les profils, pas d'après le plan.
