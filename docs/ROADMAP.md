# Prochaine étape et feuille de route

_Extrait de l'ancien CLAUDE.md (R-030 point 7)._

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
