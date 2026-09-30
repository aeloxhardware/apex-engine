# Apex Slicer — Bugs résolus

Journal des bugs diagnostiqués puis corrigés, mesures à l'appui. Sorti de
`CLAUDE.md` pour ne pas polluer le contexte actif avec de l'historique clos
— `CLAUDE.md` ne garde que les limites de conception encore ouvertes et les
décisions techniques toujours en vigueur. Un fix déplacé ici reste
recherchable mais n'est plus chargé par défaut dans le contexte de travail
courant.

Convention : chaque entrée garde le format déjà en usage dans CLAUDE.md
(cause mesurée, fix retenu, résultat mesuré) — rien n'est réécrit en la
déplaçant, seulement recopié puis retiré de CLAUDE.md.

---

## `nozzle_contact_clearance` — faux positif systématique en étape 9

**L'étape 9 produisait un faux positif structurel à CHAQUE point de
trajectoire**, rendant le pipeline complet inexploitable de bout en bout
malgré un slicing (step1→step8) fonctionnel et visuellement cohérent.

Mesure AVANT fix (`l_shape_overhang.stl`, modèle de tête réel) :
**exactement 2 collisions par point échantillonné**, sans aucune exception
— 498 collisions à 250 points, 998 à 500, 1498 à 750, 3498 à 1750. Un
`printed_part` + un `platform` à chaque point. Ce n'était pas un motif de
collision réel, c'était une propriété structurelle.

**Cause** : l'origine du modèle de tête EST la pointe de buse (convention
vérifiée et voulue), et `_surface_sample_points` incluait ce sommet à
l'origine parmi les points de test. Or la pointe de buse est **par
définition** en contact avec la matière qu'elle vient de déposer et à
quelques dixièmes de millimètre du plateau aux premières couches. La tester
avec `safety_margin=1.0 mm` garantissait une collision à chaque point. La
limite #22 (exclusion du cordon en cours de dépôt) ne suffisait pas : elle
n'exclut que la capsule qui se termine AU point testé, pas les cordons
voisins déjà posés sur la même couche, ni le plateau.

**Fix retenu** : nouveau champ `CollisionSimParams.nozzle_contact_clearance`
(mm, défaut `3.0`) — les points d'échantillonnage du modèle de tête situés
à moins de cette distance de l'origine locale (la buse) sont exclus des
tests de collision, la zone de contact buse/matière étant attendue par
construction. Filtrage effectué **une seule fois** (nouvelle fonction
`_filter_nozzle_contact_zone`, `pipeline/step9_collision_sim.py`), sur
`toolhead_base_samples` en repère LOCAL juste après
`_surface_sample_points` dans `run()` — jamais recalculé par point de
trajectoire (le mesh de tête et son échantillonnage sont invariants le
long de la trajectoire). À l'échelle de la buse/marge de sécurité (quelques
mm), délibérément PAS dérivé de `nozzle_diameter` (échelle du cordon, même
piège que les limites #21/#26). `clearance <= 0` désactive le filtre
(comportement historique explicitement accessible, testé).

**Garde-fou trouvé et corrigé pendant l'implémentation** : si TOUS les
points d'échantillonnage tombent dans la zone de garde (tête synthétique
plus petite que `clearance`, jamais rencontré sur le modèle de tête réel
dont le corps monte à ~77mm), le filtre produisait un tableau vide, qui
faisait planter `_toolhead_aabb` en aval (`.min(axis=0)` sur un tableau
vide) — et aurait sinon désactivé silencieusement toute détection de
collision pour cette tête. `_filter_nozzle_contact_zone` retombe
explicitement sur le tableau d'origine intact dans ce cas : mieux vaut
re-tester la zone de contact que de ne plus rien tester du tout.

**Résultat mesuré après fix** (`l_shape_overhang.stl`, modèle de tête réel,
250 premiers points de trajectoire échantillonnés) : **498 → 0 collisions**,
`step10` converge proprement (`converged=True`, 0 ajustement nécessaire) —
contre un flot de faux positifs ininterprétable avant. Un test dédié
(`test_obstacle_beyond_clearance_still_detected`,
`tests/test_step9_collision_sim.py`) confirme qu'un vrai obstacle touchant
le CORPS de la tête (au-delà de `clearance`) reste bien détecté — le fix
ne désactive pas la détection au-delà de la zone de garde, seulement la
zone de contact attendue elle-même. Suite complète : **222/222 tests**
passants (217 avant ce fix + 5 nouveaux), zéro régression.

`main.py` appelait déjà step9/step10 réellement (pas de changement
d'orchestration nécessaire) — c'est `scripts/export_slicing_preview.py`
qui continue de les court-circuiter délibérément pour la vitesse
d'itération (`CollisionReport(nozzle_model_checked=False)`,
`AdjustmentReport(converged=False)`, jamais un rapport vide qui
prétendrait qu'aucune collision n'existe), un choix distinct de ce fix,
pas remis en cause ici.

**Point ouvert distinct, non résolu par ce fix** : le débit de step9 se
dégrade avec la matière accumulée (mesuré AVANT ce fix : 16.5 → 3.7 pts/s
entre le 250ᵉ et le 1750ᵉ point échantillonné, index passant de quelques
centaines à ~2700 capsules ; la grille spatiale de `_PrintedGeometryIndex`
limite le nombre de candidats mais pas assez). Confirmé toujours présent
après ce fix : sur `l_shape_overhang.stl`, la trajectoire complète (12105
points, modèle de tête réel) ne termine pas step9 en moins de 550s alors
que les 250 premiers points seuls prennent 0.8s — un problème de
performance distinct, hors périmètre de ce fix (qui portait sur la
CORRECTION du signal, pas sa VITESSE), reste documenté comme limite ouverte
dans `CLAUDE.md`.

## Cavités internes non signalées à tetgen (add_hole manquant)

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

---

## Détection de pont par couche : seuil de majorité jamais atteignable

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

---

## MeshData.kdtree confondu avec une vraie distance à la surface

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

---

## Conflit du dossier io/ avec le module standard Python io

20. **Le dossier `io/` du projet entre en conflit avec le module standard
    Python `io` (builtin/frozen).** `from io.mesh_loader import ...` échoue
    avec `ModuleNotFoundError: 'io' is not a package` — Python résout
    `io` vers le module standard avant de considérer un package local du
    même nom, quel que soit l'ordre de `sys.path`. Détecté concrètement en
    écrivant l'étape 9 (premier consommateur réel du module). Renommé en
    **`file_io/`** — tous les imports doivent utiliser `file_io.mesh_loader`,
    `file_io.gcode_export`, jamais `io.*`.

---

## Collision manquée : mesh de tête testé par sommets seuls

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

---

## Faux positif structurel : la buse confondue avec son propre cordon

22. **La buse est par construction confondue avec le cordon qu'elle vient
    de déposer.** Tester un point en extrusion contre la matière imprimée
    "jusqu'à ce point INCLUS" déclenche une collision à chaque point en
    extrusion (distance quasi nulle à son propre cordon) — faux positif
    structurel, pas un cas limite rare. `step9_collision_sim._grow_printed_index_to`
    exclut la capsule qui se termine AU point testé ; elle ne devient une
    collision possible qu'au(x) point(s) suivant(s), une fois que la buse
    s'en est réellement éloignée.

---

## nozzle_diameter confondu avec l'échelle d'échantillonnage du modèle de tête

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

---

## Orientation du modèle de tête non vérifiée par mesure (tête pendant sous la buse)

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

---

## Re-vérification step10 reconstruite à chaque appel dans la boucle de recherche

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

---

## Débordement du toolpath hors du solide (projection/relevage step6 plan moyen)

31. **La paire projection/relevage de step6 (`_project_to_local_2d` /
    `_lift_to_3d`) aplatit chaque couche sur SON plan local moyen
    (step5 `_compute_local_frame`, Decision 6, déjà documentée
    "approximatif") — perte d'information mesurée concrètement, pas
    théorique.** Sur une couche fortement inclinée et non-plane (mesuré :
    jusqu'à 1.30mm d'écart entre un point de contour réel et le plan
    approximatif, sur `l_shape_overhang.stl` avec segmentation à 1 seule
    région, direction moyenne ~55°), le point 3D reconstruit après
    projection puis relevage ne retombe PAS forcément sur le point
    d'origine — visible concrètement en trajectoire finale : **344 points
    sur 9391** sortaient du volume réel du solide (jusqu'à 0.16mm sous le
    plateau, en pleine extrusion), particulièrement près d'un bord fin
    proche du plateau où l'inclinaison du plan est la plus pénalisante.
    Diagnostic par élimination successive avant fix : `Layer.contour_3d`
    brut (sortie step5) est déjà correct à 100% (0 point hors solide) ;
    `_offset_walls`/`shapely.buffer` est innocenté par test direct (le
    débordement existe même sans aucun offset, sur le simple aller-retour
    projection+relevage) — la cause est exactement la paire
    projection/relevage, rien d'autre.
    **Fix** : `Layer.toolpath_2d` porte désormais 3 nombres par point
    `(u, v, offset_hors_plan)` au lieu de 2 — `step6._out_of_plane_offset`
    calcule cet écart au moment de la projection, `_to_toolpath_3col`
    l'attache à la polyligne, `_lift_to_3d` (accepte (K,2) OU (K,3)) le
    réinjecte au relevage : `origin + u·axis_u + v·axis_v + offset·normal`.
    Portée délibérément limitée au contour n=0 (paroi la plus extérieure,
    identique au contour d'origine — `_offset_walls` l'ajoute tel quel en
    dernier élément) : c'est elle qui doit coller à la vraie surface. Les
    parois intérieures (n=1..wall_count) et l'infill sont des points 2D
    NOUVEAUX générés par `shapely`, sans équivalent 3D d'origine à
    préserver — offset implicitement nul pour eux, comportement inchangé.
    Pas de changement à `_compute_local_frame` (step5) ni à
    `_offset_walls`/`shapely.buffer` (tous deux confirmés hors de cause).
    **Résultat mesuré après fix** : 344/9391 → **1/9391** point encore hors
    solide (résidu sur une paroi intérieure/infill, hors périmètre du fix
    par construction — négligeable, 0.15mm). Confirmé visuellement sur
    `out_test.glb` régénéré.
    **Effet de bord à retenir pour toute future étape consommant
    `toolpath_2d`** : le contrat de données passe de (K,2) à (K,3) —
    toute comparaison géométrique 2D pure (`point_in_polygon_2d`, produit
    matriciel avec une base 2D, direction `p1-p0`) doit explicitement
    tronquer aux 2 premières colonnes (`poly[:, :2]`), et tout
    `reshape(1, 2)` figé doit devenir `reshape(1, -1)` — deux régressions
    de ce type ont été mesurées et corrigées dans
    `step7_scheduling._group_position_3d` et dans les tests existants au
    moment de ce fix, pas hypothétiques.
    **SUITE (session ultérieure) — le résidu ci-dessus est traité : offset
    hors-plan INTERPOLÉ pour les points générés par `shapely`.** Le fix
    initial ne couvrait que le contour n=0 (seul à avoir un équivalent 3D
    d'origine à restaurer) ; les parois n≥1 et l'infill, points 2D nouveaux,
    gardaient un offset nul et restaient donc collés au plan local MOYEN.
    Sur un modèle plus tordu que celui du diagnostic d'origine
    (`l-triple-axis-overhanged-shape.stl`), ce résidu remonte à **74 points
    hors solide** (43 au-delà de 0.05mm, 24 au-delà de 0.1mm, 4 au-delà de
    0.2mm, max **0.390mm**) — avec une signature sans ambiguïté mesurée
    point par point : **0 %** des dépassements viennent du contour n=0
    (vérifié à 0.0000mm de la surface), **100 %** des parois internes et de
    l'infill. Cause : la couche s'écarte de son propre plan moyen de jusqu'à
    **0.845mm** sur ce modèle (médiane 0.079, moyenne 0.217).
    **Fix** : `step6_toolpath2d._interpolate_out_of_plane` — pondération
    inverse-distance (`1/d**power`) des offsets RÉELS connus sur les boucles
    de contour du même groupe (contour extérieur **et** ses trous : les
    ignorer dégraderait l'interpolation à leur voisinage), appliquée à tout
    point généré. Nouveau `Step6Params.out_of_plane_interpolation_power`
    (`Optional[float]`, défaut `None` = désactivé, même discipline opt-in que
    `direction_blend_radius=0.0`), activé à `2.0` dans `main.PipelineParams`
    (même précédent que `field_solver` qui y force `direction_blend_radius=10.0`).
    `scipy.interpolate.LinearNDInterpolator` écarté : la triangulation de
    Delaunay d'un ANNEAU ne couvre que l'enveloppe **convexe**, donc rendrait
    `NaN` à l'intérieur de tout contour non convexe — cas courant ici.
    **Invariant préservé, testé explicitement** : le contour n=0 garde ses
    offsets EXACTS, jamais interpolés.
    **Résultat mesuré** (même modèle, mêmes params) : 74 → **0** point hors
    solide à tous les seuils, max 0.390mm → **0.0021mm**. Balayage de
    l'exposant : `power=1.0` laisse 20 points (max 0.096mm), `power=3.0` en
    laisse 1 (max 0.032mm), `power=2.0` zéro — d'où ce défaut. Géométrie 2D
    strictement inchangée (155 couches, 1016 polylignes, 15178 points de
    trajectoire, identiques avant/après) : le mécanisme ne touche QUE la 3ᵉ
    colonne, jamais le tracé ni l'ordre centre→extérieur.

---

## Transition brusque d'orientation au plus-proche-voisin strict (tet_directions)

32. **Le plus-proche-voisin STRICT utilisé pour `tet_directions` (étape 4)
    produisait une transition brusque d'orientation à la frontière
    plateau/région ET entre régions adjacentes, sans aucune notion de
    distance au plateau.** Diagnostiqué par lecture de code puis confirmé
    empiriquement sur un modèle tripode synthétique (pied + 3 bras non
    coplanaires, régions injectées manuellement) : `region_ids` passe de
    `[0,3]` à `[0,1,2,3]` (les 4 régions simultanément) sur moins de 0.1mm
    d'isovaleur, et `n_loops` par couche saute de 1 à 4 en 3 couches
    consécutives (fragmentation de contour, pas seulement un artefact
    d'orientation). **Fix retenu** : `_blend_tet_directions`
    (`pipeline/step4_field_solver.py`) remplace le plus-proche-voisin par
    une moyenne pondérée par distance inverse sur toutes les sources
    d'influence proches d'un tétraèdre (le plateau, direction `up_vector`,
    distance verticale à `platform_level` ; chaque région, sa
    `target_direction`, distance à SES SEULES faces via un cKDTree dédié).
    Nouveau `FieldSolverParams.direction_blend_radius` (défaut `0.0` =
    DÉSACTIVÉ, comportement historique inchangé — opt-in explicite, même
    discipline que `direction_deviation_weight` de step3) et
    `direction_blend_power` (exposant inverse-distance, défaut `2.0`).
    Fenêtrage LISSE (cosinus/Hann, C1-continu à `radius`) plutôt qu'une
    coupure nette, pour qu'aucune source n'entre/sorte du rayon
    d'influence de façon discontinue. Activé par défaut dans
    `main.py` (`direction_blend_radius=10.0`, `direction_blend_power=2.0`,
    calibré sur `l_shape_overhang.stl`) : angle du gradient réel mesuré
    parfaitement vertical (0.0°) sur les tétraèdres les plus proches du
    plateau, contre jusqu'à 42.7° sans mélange.

---

## Annulation vectorielle de la moyenne pondérée de directions (cas symétrique)

33. **Une moyenne pondérée de directions peut s'annuler vectoriellement si
    deux sources ont des directions quasi-opposées à poids proches**
    (`sum(w_i * direction_i)` de norme quasi nulle) — cas non rare, une
    pièce symétrique (croix, étoile) le produit naturellement en son
    centre. Normaliser un vecteur quasi nul serait numériquement instable
    (direction résultante arbitrairement bruitée). `_blend_tet_directions`
    détecte ce cas (`FieldSolverParams.direction_blend_min_norm`, défaut
    `0.1` — la norme d'une moyenne de vecteurs unitaires vaut au maximum
    1.0 si alignés, tend vers 0 s'ils s'annulent) et retombe sur la
    direction de la SEULE source au poids maximal pour ce tétraèdre précis
    — équivalent local au winner-take-all historique, uniquement là où le
    mélange continu serait mathématiquement dégénéré. Ce repli reste
    ponctuel (mesure nulle en pratique sur un maillage irrégulier), il ne
    réintroduit pas de saut à grande échelle.

---

## Erreur de discrétisation FEM (escalier) près du plateau

34. **Un résidu d'« escalier » dans l'orientation des couches près du
    plateau, visible après le fix #32, N'EST PAS causé par le mélange de
    directions lui-même — c'est une erreur de discrétisation FEM
    classique.** Diagnostiqué par élimination successive, mesure directe à
    chaque étape (pas supposé) :
    - Hypothèse initiale (coupure nette de `direction_blend_radius`)
      infirmée : un fenêtrage lisse (cosinus/Hann, voir limite #32) ne
      change quasi rien à l'écart mesuré sur un tétraèdre témoin
      (32.34°→32.36°).
    - Balayage systématique de tous les poids susceptibles d'expliquer
      l'écart entre la direction CIBLE (`_blend_tet_directions`, vérifiée
      lisse en elle-même) et le gradient RÉELLEMENT résolu par le FEM :
      `platform_anchor_weight` (1e6→1e1, 5 ordres de grandeur : quasi
      aucun effet, 9.77°→7.60° tout en dégradant la planéité du plateau),
      `laplacian_smoothness_weight` (1.0→100.0 : aucun effet), et
      `overhang_influence_depth` (0→3 : aucun effet).
    - Seule la **résolution du maillage** (taille de tétraèdre) a un effet
      net et monotone : mean dev 9.77°→7.07°→5.48°→4.77° à mesure que la
      résolution globale passe de 1.5mm à 0.5mm (307→9088 tétraèdres).
    - **Cause racine** : un tétraèdre P1 proche du plateau peut s'étendre
      de z=0 (ancré à phi=0) jusqu'à plusieurs mm au-dessus, plage sur
      laquelle la direction CIBLE (mélangée par #32) tourne déjà de façon
      continue — mais un tétraèdre P1 ne représente qu'UN SEUL gradient
      constant sur toute son étendue. Erreur de discrétisation classique,
      pas un défaut de pondération — voir limite #36 pour le fix retenu
      (raffinement local de maillage).

---

## Raffinement local du maillage près du plateau (TetGen add_region)

36. **Raffinement LOCAL du maillage tétraédrique près du plateau, via
    `TetGen.add_region` (deux régions, pas un champ continu) — corrige la
    limite #34, en évitant le mécanisme cassé de la limite #35.**
    `pipeline/step1_ingest.py::_split_for_platform_refinement` découpe le
    mesh source en deux coques étanches de part et d'autre de
    `z = platform_z + near_platform_margin` (`trimesh.slice_plane(...,
    cap=True)`), les concatène en un seul PLC (même convention que
    l'exemple officiel `add_region`, `pv.merge([sphere, cube])` — la
    frontière de coupe partagée devient une facette interne), et calcule un
    point-seed par coque (centroïde, `z` explicitement repoussé vers
    l'intérieur de sa région pour ne jamais tomber sur la frontière de
    coupe). `_tetrahedralize` construit alors `TetGen` sur ce PLC combiné
    (pas le mesh source directement), appelle
    `add_region(1, slab_seed, near_platform_resolution**3)` et
    `add_region(2, rest_seed, resolution**3)`, puis tétraédrise avec
    `switches="pq1.2Aa"` (`-A` active les régions, chacune utilisant SON
    PROPRE `max_vol` déclaré via `add_region` plutôt qu'un plafond global
    unique) — **vérifié fonctionnel, aucun crash**, contrairement à la
    limite #35. Nouveaux paramètres opt-in sur `step1_ingest.run()` :
    `near_platform_resolution` (défaut `None` = désactivé, comportement
    historique inchangé — même discipline que `direction_blend_radius=0.0`
    de la limite #32), `near_platform_margin` (défaut `5.0`), `platform_z`
    (défaut `0.0` à l'origine, même convention que
    `SegmentationParams.platform_z` — **devenu `Optional[float] = None`
    avec détection automatique, voir limite #45**).
    Un seed de coque peut légitimement échouer le test `mesh.contains()`
    de trimesh sur une forme non convexe (mesuré : `rest_seed` sur
    `l_shape_overhang.stl`) SANS empêcher TetGen de tétraédriser
    correctement depuis ce seed (propagation topologique interne, pas un
    test d'inclusion géométrique strict) — pas de garde-fou bloquant sur ce
    test, volontairement.
    **Résultat mesuré sur `l_shape_overhang.stl`** (`near_platform_resolution=0.5`,
    `near_platform_margin=5.0`, activé par défaut dans `main.py`) :
    **1677 tétraèdres**, déviation angulaire moyenne du gradient près du
    plateau **4.91°** (médiane 4.21°, max 24.20°) — quasi identique à un
    raffinement GLOBAL complet (`resolution=0.5` partout : 9088
    tétraèdres, 4.77°/4.08°/20.27°), pour **~5.4× moins de tétraèdres**.
    Effet de bord positif mesuré sur le pipeline complet : diagnostics
    d'orientation inatteignable (étape 8) réduits de 531 à **168** sur le
    même modèle (des directions plus précises localement sont en moyenne
    plus faciles à atteindre pour la cinématique 2 axes de la machine).
    Confirmé visuellement par l'utilisateur sur `out_test.glb` régénéré :
    effet d'« escalier » réduit.
    **Portée limitée au cas plateau, PAS aux jonctions multi-région**
    (le second cas mesuré comme problématique par le modèle tripode,
    limite #32) : les frontières de région ne sont connues qu'après
    step3, deux étapes après la tétraédrisation — aucun raccourci
    géométrique équivalent au cas plateau (`platform_z` connu dès
    `mesh.vertices`). Traiter le cas jonction demanderait une
    re-tétraédrisation en deux passes (step1→step3 une première fois pour
    localiser les jonctions, puis re-tétraédriser avec des régions
    informées par ces jonctions, puis rejouer step4+) —
    **fait en session ultérieure, voir limite #48** : mécanisme implémenté
    et validé sur jonction localisée (fixture tripode), mais confirmé
    NON APPLICABLE sur `death.stl` (zone dense trop étendue, deux
    approches de repli testées sans effet mesuré) — limite #48 documente
    le résultat complet, pas une simple fermeture de ce point.

---

## step8 dominé par le coût numpy générique sur des rotations scalaires

37. **step8 (résolution cinématique inverse par point de trajectoire) prenait
    ~98% du temps total du pipeline sur un modèle trivial (24 triangles,
    45,6s au total sur `l_shape_overhang.stl`) — pas step5 comme on aurait
    pu le supposer (0,53s seulement).** Diagnostiqué par chronométrage par
    étape puis `cProfile` (pas supposé) :
    - `_rodrigues_rotate(vectors, axis, angles_rad)` (`core/geometry_utils.py`)
      est écrite et documentée pour un usage BATCH (tourner N vecteurs à la
      fois — utilisée ainsi par `composed_rotation_angles_and_directions`
      pour échantillonner une grille, et par step9 pour tourner tous les
      sommets d'un mesh de tête à la fois, deux usages batch légitimes).
      Mais `refine_axis_angles_local_search._direction_for` l'appelait avec
      **UN SEUL vecteur** à chaque nudge de son coordinate-descent — 342936
      fois pour 10392 points de trajectoire sur ce modèle
      (`ik_refine_iterations=8` × 2 axes × 2 sens × 10392 points).
      `np.cross` sur un tableau à une seule ligne passe par tout le chemin
      générique numpy N-dimensionnel (`moveaxis`, `normalize_axis_tuple`,
      `broadcast_to`) — un surcoût énorme mesuré (57,6s cumulés sur ce seul
      appel) pour un produit vectoriel 3×3 trivial.
    - **Fix** : nouvelle fonction `_rodrigues_rotate_single(vector, axis,
      angle_rad)` — même formule de Rodrigues, en scalaire pur (indexation
      directe `vector[0]`/`[1]`/`[2]`, pas de `np.cross`), vérifiée
      identique par test (`np.allclose`) contre le chemin batch existant.
      `_rodrigues_rotate` elle-même reste **strictement inchangée**
      (signature, comportement, tous ses appelants batch intacts) — le
      risque de régression sur le chemin batch est géré par
      non-modification, pas par un test ajouté après coup (aucun test
      unitaire dédié n'existait pour `_rodrigues_rotate` isolément avant ce
      fix). `refine_axis_angles_local_search._direction_for` ET
      `step10_adjustment._direction_for_angles` (même pattern scalaire,
      même défaut latent bien que moins chaud) basculent sur la nouvelle
      fonction.
    - **Deuxième passe** (même classe de problème, une couche plus bas,
      trouvée par re-profilage après le premier fix) : `np.clip`/`np.dot`/
      `np.linalg.norm` appelés sur des scalaires triviaux dans
      `_angular_error` et dans le clamp d'angle du coordinate-descent
      (~700000 appels chacun) subissaient le même surcoût de dispatch
      générique. Remplacés par du Python scalaire pur
      (`math.acos`, comparaisons directes). Bonus algorithmique trouvé au
      passage : les axes unitaires (`spec.axis / norm`) étaient
      renormalisés à CHAQUE nudge alors qu'ils ne changent jamais pendant
      une recherche — précalculés une seule fois désormais
      (`reversed_axes_unit`).
    - **Résultat mesuré sur `l_shape_overhang.stl`** (24 triangles,
      10392 points de trajectoire) : step8 44,9s → 12,7s (1er fix) → 7,1s
      (2e fix) ; total pipeline step1→step8 45,6s → **7,9s** (~5,8× plus
      rapide). Sur `res/test_models/l-triple-axis-overhanged-shape.stl`
      (30 triangles, modèle 3 axes non coplanaires, 14144 points de
      trajectoire, 0 orientation inatteignable) : step8 9,5s, total
      **10,4s**.

---

## Segmentation : fusion de régions géométriquement incompatibles

38. **La segmentation (étape 3) fusionnait deux régions géométriquement
    incompatibles (ex: deux bras perpendiculaires d'un modèle en T) des
    qu'UNE SEULE direction candidate restait techniquement compatible avec
    les deux — sans jamais regarder À QUEL POINT cette direction de
    compromis s'écarte de l'optimum de chacune.** Diagnostiqué visuellement
    (trajectoire finale en éventail cassé sur un modèle en T réel,
    `res/test_models/l-triple-axis-overhanged-shape.stl`), confirmé par
    mesure : sur ce modèle, la segmentation par défaut ne produisait
    **qu'une seule région** pour tout le solide, y compris ses deux bras à
    90° l'un de l'autre.
    - **Premier critère tenté et ÉCARTÉ** : tolérance sur le RISQUE RÉSIDUEL
      scalaire (`region_count_penalty`, un ratio de dégradation). Mesure
      directe sur le modèle en T : `own_best=shared_best=0.0` pour TOUTES
      les fusions testées, quelle que soit la tolérance (0 à 100%) — avec
      `direction_candidate_count≈26`, il existe presque toujours une
      direction candidate de résidu quasi nul même pour une fusion
      géométriquement incohérente (la sphère de candidats est assez dense
      pour qu'un compromis numérique existe sans être géométriquement
      cohérent). Défaut structurel qui **s'aggrave** avec la densité de
      candidats — dangereux pour une cinématique future dense/continue
      (limite #18).
    - **Fix retenu** : critère ANGULAIRE
      (`SegmentationParams.max_region_merge_angle_deg`, `Optional[float]`,
      défaut `None` = désactivé, même discipline que
      `direction_blend_radius=0.0`/`near_platform_resolution=None`). Une
      fusion n'est acceptée que si la direction optimale de la région
      combinée (argmin du risque résiduel, restreint aux directions déjà
      compatibles avec les deux) reste à moins de
      `max_region_merge_angle_deg` de la direction optimale de CHACUNE des
      deux régions séparément (`_angle_between_deg`, `_region_growing` dans
      `pipeline/step3_segmentation.py`). Un angle entre deux directions ne
      dépend PAS de la densité de candidats, contrairement à un minimum sur
      un grand ensemble — robuste par construction à cette classe de
      problème.
    - **Résultat mesuré sur le modèle en T** (`max_region_merge_angle_deg=30.0`) :
      1 région → **5 régions**, écart hors-plan max des couches (mesure
      indépendante en aval) réduit de 3,86mm à 2,67mm (-31%). Amélioration
      visuelle confirmée par l'utilisateur sur les deux bras (courbes
      fluides), mais résidu signalé à la jonction pied/bras — voir limites
      #39/#40 ci-dessous pour la suite du diagnostic.

---

## step5 : couche non scindée malgré des contours spatialement disjoints

39. **Le clustering par isovaleur (étape 5) devait scinder une couche en
    plusieurs `Layer` distincts quand ses boucles de contour forment
    plusieurs composantes spatialement disjointes ET non imbriquées** — cas
    apparu concrètement seulement APRÈS le fix #38 (avant, la segmentation
    ne produisait quasiment jamais plus d'une région, donc ce cas de figure
    n'émergeait pas). Mesuré sur le modèle en T : ~20 couches sur 107
    contenaient 2 boucles de contour massivement disjointes (pied vs bout de
    bras, 16+ unités d'écart) partageant la même isovaleur — `Layer.local_frame`
    (`_compute_local_frame`) calculait un plan moyen UNIQUE pour les deux
    boucles combinées, géométriquement incohérent pour l'une comme pour
    l'autre (cause directe du faisceau de trajets convergents visible sur le
    rendu GLB).
    - **Invariant protégé, jamais régressé** : une boîte creuse (paroi
      extérieure + cavité interne) produit volontairement 2 boucles
      IMBRIQUÉES dans le MÊME `Layer` à mi-hauteur — ce cas ne doit JAMAIS
      être scindé (`tests/test_step5_layers.py::test_hollow_box_nested_contours_not_split_even_when_enabled`).
    - **Fix retenu** : nouveau paramètre `LayerParams.split_disjoint_clusters`
      (`Optional[float]`, défaut `None` = désactivé). Deux boucles restent
      dans le même cluster si IMBRIQUÉES (`point_in_polygon_2d`, priorité
      absolue) OU à une distance de bord-à-bord ≤
      `split_disjoint_clusters * layer_height_base` (seuil adaptatif à
      l'échelle du modèle, réutilise `layer_height_base` déjà présent
      plutôt qu'une nouvelle constante absolue). Regroupement par
      `nx.Graph` + `nx.connected_components` (`_cluster_contours`,
      `pipeline/step5_layers.py`).
    - **`region_ids` par cluster** : approximation explicitement documentée
      comme telle — la correspondance triangle→tétraèdre-source est perdue
      par `trimesh.outline()` (vérifié : les entités de `outline.entities`
      ne portent aucune référence aux indices de `faces` d'entrée). Recalcul
      par proximité spatiale (`_tets_for_cluster`) : tolérance dérivée de la
      taille CARACTÉRISTIQUE D'UN TÉTRAÈDRE (médiane des longueurs d'arête
      × 2.0), **PAS** une fraction de l'étendue globale du modèle — bug
      réel mesuré et corrigé pendant l'implémentation : une fraction de
      `scale` (~0.32 sur le modèle tripode de test) tombait EN DESSOUS de
      la distance réelle centroïde-tétraèdre → point de cluster le plus
      proche (mesurée à 0.328 sur ce même modèle), filtrant TOUS les
      tétraèdres et forçant un repli sur `unique_tets` entier
      (`region_ids` identiques pour tous les clusters, brisant
      silencieusement l'approximation).

---

## step5 : ordre des couches scindées non continu spatialement (alternance de branches)

40. **Le clustering spatial (limite #39) résout le mauvais plan local, mais
    PAS l'ordre final des `Layer` scindés — deux clusters issus de la MÊME
    région step3 (une région est une nappe de SURFACE connexe, pas
    nécessairement une seule composante spatiale du SOLIDE) restaient
    ordonnés uniquement par isovaleur croissante, les faisant alterner à
    quasi chaque isovaleur.** Mesuré sur le modèle en T (avec les fixes
    #38+#39 actifs) : la région "verticale" step3 a une bounding box
    `x∈[-2,18], y∈[-20,0], z∈[6,22]` — elle englobe À LA FOIS le pied ET la
    pointe de bras (surface d'un seul tenant, solide en Y/T). `_sort_key =
    (rank, isovalue)` ne trie que par rang de région (identique pour les
    deux clusters) puis isovaleur : `layers[9,10,11,12,13,14]` alternaient
    pied→pointe→pied→pointe à quasi chaque isovaleur, `step7_scheduling`
    insérant un trajet à vide (`step8`) à chaque bascule — le faisceau
    convergent résiduel signalé par l'utilisateur après le fix #38 seul.
    - **Fix retenu, généralisé (pas un cas "en T" codé en dur)** : graphe de
      continuité spatiale `_build_continuity_graph` (`nx.DiGraph`) sur
      `raw_layers` (avant tri final) — un nœud par cluster, une arête entre
      deux clusters d'isovaleurs DISTINCTES ADJACENTES (dans la séquence
      réellement présente, pas la grille théorique) s'ils sont connectés
      (même critère que le clustering intra-isovaleur, factorisé en
      `_loops_connected`, réutilisant le MÊME `gap_threshold` — pas une
      nouvelle constante). Le tri final devient un parcours DFS
      multi-racines (`_spatial_dfs_order`/`_dfs_order_within_rank`),
      **exécuté séparément par rang topologique de région** (jamais à
      travers un rang — partition STRICTE des nœuds avant tout parcours,
      pas une pondération : rend une violation du rang impossible par
      construction, pas seulement improbable).
    - **Ordre de visite à un branchement** : la branche de PLUS PETITE
      taille de sous-arbre en premier (`_subtree_sizes`, via
      `nx.descendants` — pas une somme naïve, pour gérer correctement les
      recombinaisons où deux parents convergent vers un même enfant).
      Justification : termine les branches courtes d'abord, jamais de
      retour en arrière sur une grande branche restée en suspens ; robuste
      au cas ≥3 enfants, contrairement à un critère de proximité au point
      de sortie du parent qui dégénère en quasi-égalité au moment même du
      split (par construction, les deux enfants sont chacun proches du même
      parent à un seuil comparable).
    - **Activation : PAS de nouveau paramètre séparé** — le tri DFS s'active
      automatiquement et exclusivement quand `split_disjoint_clusters is
      not None` : sans clustering actif, `raw_layers` a au plus un nœud par
      isovaleur, le graphe dégénère en chaîne simple, et le DFS produit
      EXACTEMENT le même ordre que `(rank, isovalue)`. Un flag séparé
      introduirait une combinatoire à cellules mortes.
    - **Cas testés au-delà du "T"** : `tests/test_step5_layers.py` ajoute une
      fixture `_portal_mesh`/`portal_data` (portique en U inversé, 2
      montants séparés + 1 traverse) pour prouver le cas MERGE
      (recombinaison, in-degree ≥ 2 dans le graphe), complémentaire au cas
      SPLIT déjà couvert par `tripod_data`. `test_spatial_order_is_a_permutation`
      teste `_spatial_dfs_order` unitairement (raw_layers synthétiques) pour
      garantir qu'aucun index n'est perdu ni dupliqué.
    - **Résultat mesuré sur le modèle en T** (`max_region_merge_angle_deg=30.0`
      + `split_disjoint_clusters=3.0`) : sur les MÊMES 157 couches
      (clustering identique, seul l'ordre change), grands sauts spatiaux
      entre clusters consécutifs réduits de **42 → 7** (-83%) par rapport à
      l'ordre `(rank, isovalue)` naïf reconstruit sur les mêmes couches.
    - **Effet de bord mesuré, non résolu** : le nombre d'orientations IK
      inatteignables (`step8_diagnostics`) reste à 528 sur ce modèle,
      inchangé par l'activation du clustering — ce fix ne touche que
      l'ORDRE des couches, pas les angles de tête résolus par step8 ; la
      cause de ces 528 diagnostics est un point ouvert, non investigué plus
      avant dans cette session.
    - **Piste "augmenter la résolution du maillage" testée et écartée pour
      les sauts résiduels (session ultérieure, après le fix limite #42)** :
      hypothèse initiale de l'utilisateur — le résidu de grands sauts
      (7 → 5-6 après le fix #42) pourrait venir d'un maillage encore trop
      grossier à certains endroits, et augmenter `resolution` globalement
      (ou étendre le raffinement local `add_region` du plateau aux zones de
      surplomb) pourrait le résorber. **Mesuré directement sur le modèle en
      T avant tout code : AUCUN effet.** `resolution` variant de 1.5 à 0.3
      (1913 → 57992 tétraèdres, ×30) : nombre de grands sauts constant
      (6, bruit à 7 sur le maillage le plus fin), distance totale de trajet
      à vide quasi inchangée (~96-102, sans tendance monotone). Diagnostic
      complémentaire (isovaleurs de part et d'autre de chaque saut) :
      plusieurs sauts sont des transitions à isovaleur DÉCROISSANTE
      (ex: 11.315→5.475, 11.419→7.769) — pas des continuations lisses au
      sein d'une branche, mais des points où le DFS achève une branche
      (pied ou un bras) et recommence sur une autre branche disjointe du
      même rang topologique. Un maillage plus fin ne change ni le nombre de
      branches ni leur disjonction spatiale réelle : structurellement, il
      ne peut pas faire disparaître ce type de saut. Confirme l'hypothèse
      déjà notée plus haut ("probablement structurellement inévitables") —
      désormais mesurée, pas seulement supposée. Écarté : raffiner
      localement les zones de surplomb nécessiterait de toute façon une
      tétraédrisation en deux passes (limite #36, non faite), un chantier
      disproportionné pour un résidu dont la cause n'est pas la résolution.

---

## TetGen.add_region ignore silencieusement le volume de la région "rest"

42. **`TetGen.add_region` avec deux régions de VOLUMES DIFFÉRENTS n'applique
    JAMAIS le volume de la région "rest" (le reste du modèle, censée suivre
    `resolution`) — c'est TOUJOURS le volume associé au SEED DU SLAB
    (`near_platform_resolution`, la fine tranche près du plateau) qui
    s'applique à TOUT le maillage combiné, quelle que soit sa valeur.**
    Bug réel du binding tetgen (0.8.4) découvert en investiguant un résidu
    de fragmentation en step5 sur le modèle en T
    (`res/test_models/l-triple-axis-overhanged-shape.stl`) : `resolution`
    n'avait **aucun effet mesurable** sur le maillage produit (1883
    tétraèdres identiques de `resolution=1.5` à `resolution=0.1`), y
    compris loin de toute zone plateau (bout de bras à z≈18-20). Diagnostic
    par élimination successive, 6 mesures isolées sur `_tetrahedralize`
    directement (pas supposé) :
    - Sans `add_region` (une seule région, un seul volume global) :
      `resolution` fonctionne parfaitement (419 → 57992 tétraèdres de 1.5 à
      0.3).
    - Avec `add_region`, les deux régions au MÊME volume : fonctionne aussi
      (177 → 6746 tétraèdres de 1.5 à 0.3) — le PLC combiné
      slab+rest (`_split_for_platform_refinement`) n'est donc PAS en cause.
    - Avec `add_region`, volumes DIFFÉRENTS : `slab=0.1/rest=1.5` → 150222
      tétraèdres (fin partout) ; `slab=1.5/rest=0.1` → 177 (grossier
      partout) ; `slab=2.0/rest=0.3` → 126 (grossier partout, contredit
      l'hypothèse "le plus petit volume domine") — dans TOUS les cas
      mesurés, c'est le volume du SLAB qui gagne, jamais celui du rest.
    - Ni l'ordre d'appel de `add_region` (rest appelé en premier ou en
      second) ni le numéro d'attribut de région utilisé (slab=1/rest=2 vs
      slab=2/rest=1) ne changent ce résultat — écarté comme cause possible
      par mesure directe.
    - **Conséquence pratique directe et significative** : depuis que
      `near_platform_resolution=0.5` est actif PAR DÉFAUT dans `main.py`
      (limite #36), le paramètre `resolution`/`FieldSolverParams.resolution`
      n'avait plus aucun effet réel sur la densité de maillage en dehors de
      la zone plateau, sur TOUT modèle utilisant les défauts du pipeline —
      pas seulement le modèle en T de ce diagnostic.
    - **Cause racine confirmée** (session suivante) : le seed de la coque
      "rest" (`_split_for_platform_refinement`, centroïde de la coque)
      **ne tombe pas à l'intérieur de sa propre coque** sur le modèle en T
      réel — `rest.contains([rest_seed])` → `[False]` (solide non convexe :
      pied + 2 bras à 90°). Quand un seed de région tombe hors de la coque
      qu'il représente, TetGen ne peut pas la localiser dans la
      triangulation Delaunay et le comportement dégénère silencieusement
      vers l'application du volume de l'AUTRE région à tout le maillage —
      confirmé par substitution manuelle d'un seed valide, qui restaure un
      comportement monotone correct (1913 → 13179 → 53807 tétraèdres pour
      `resolution` 1.5 → 0.5 → 0.3, contre ~1883 identique partout avant).
      Le docstring précédent de `_split_for_platform_refinement` affirmait
      à tort qu'un seed hors-coque "n'empêche pas TetGen de tétraédriser
      correctement... vérifié empiriquement" — cette affirmation était
      fausse, corrigée.
    - **Fix retenu** : nouvelle fonction `_find_guaranteed_interior_point(shell,
      fallback_point, context_label)` (`pipeline/step1_ingest.py`) — essaie
      d'abord le centroïde (chemin rapide, correct sur la grande majorité
      des coques), retombe sur `trimesh.sample.volume_mesh` (rejection
      sampling, **seed fixe** pour la reproductibilité) avec escalade
      déterministe du nombre d'essais si le tirage échoue, jusqu'à un
      plafond fini (`MeshIngestError` explicite au-delà). Retourne TOUJOURS
      un seul point individuellement valide, jamais une moyenne de
      plusieurs points acceptés — piège mesuré pendant la conception : la
      moyenne de plusieurs points valides peut elle-même retomber hors
      d'une coque non convexe. Garde-fou explicite pour une coque
      vide/dégénérée (évite un `TypeError` cryptique natif de
      `trimesh.sample.volume_mesh` sur ce cas, mesuré). Utilisé pour les
      DEUX seeds de `_split_for_platform_refinement` (le fixture
      synthétique de test, une forme en U, montre que `slab_seed` peut
      aussi échouer, pas seulement `rest_seed`), et réutilisé par
      `_find_internal_cavities` (même bug de fond atteignable sur une
      cavité interne non convexe, ex. en forme de tore — remplace
      l'échec immédiat précédent par ce même repli, pure amélioration de
      couverture).
    - **Résultat mesuré après fix** sur le modèle en T réel : `resolution`
      ∈ {1.5, 0.5, 0.3} avec `near_platform_resolution=0.5` produit
      1913/13179/53807 tétraèdres — progression strictement monotone
      restaurée, conforme à la mesure de validation faite pendant la
      conception du fix. 11/11 tests dans `tests/test_step1_ingest.py`
      (6 nouveaux, dont une fixture non convexe en U dédiée qui reproduit
      le mécanisme sans dépendre d'un fichier STL externe).

---

## TetGen refuse un mesh à auto-intersections exactes (death.stl)

43. **TetGen REFUSE tout maillage dont des triangles se traversent, même
    parfaitement watertight — et aucune réparation topologique existante ne
    corrige ça.** Découvert en tentant de trancher
    `res/test_models/death.stl` (12662 triangles, le premier modèle DENSE
    du dépôt) : `RuntimeError: The input surface mesh contain
    self-intersections`, alors que le maillage passe tous les critères
    classiques (`is_watertight=True`, `is_winding_consistent=True`,
    `body_count=1`, 0 face dégénérée, 0 face dupliquée).
    - **Ce qui NE corrige PAS le problème**, testé et écarté : les
      réparations déjà en place (`merge_vertices`, `nondegenerate_faces`,
      `unique_faces`, `fill_holes`, `fix_normals`, `fix_winding`) traitent
      la TOPOLOGIE, pas la géométrie qui s'interpénètre ; une union
      booléenne `manifold3d` rend le maillage à l'identique (même nombre de
      faces, même volume, mesuré) ; les commutateurs de tolérance de TetGen
      (`-T1e-6` à `-T1e-10`) échouent tous.
    - **Piège de diagnostic mesuré** : un test d'auto-intersection écrit
      avec la seule séparation par plan de face (condition NÉCESSAIRE mais
      pas SUFFISANTE) donne des faux positifs massifs — 41 paires sur le
      modèle en T, que TetGen accepte pourtant sans broncher, typiquement
      des faces qui se TOUCHENT le long d'une arête. Un vrai SAT (2
      normales + 9 produits vectoriels arête×arête) est nécessaire :
      `_triangles_cross`.
    - **Mais même le SAT correct rapporte 0 sur `death.stl`** : le défaut
      réel, révélé par `switches="d"` (mode détection de TetGen), est
      constitué de **3 triangles seulement** en contact EXACT (« Two
      segments exactly intersect », « A vertex lies on a facet ») — pas des
      traversées franches, donc invisibles pour tout test de croisement.
      **Conséquence actée : TetGen est le seul juge fiable de ce qui le
      bloque.**
    - **Fix retenu** : réparation **RÉACTIVE** dans `_tetrahedralize` — on
      tétraédrise, et on ne répare QUE si TetGen refuse (jamais
      préventivement : un maillage déjà accepté n'est jamais modifié).
      Trois tentatives en cascade, du moins au plus invasif : tel quel →
      `pymeshfix` → sans raffinement plateau. `_tetrahedralize_once`
      extraite pour rejouer la passe sans dupliquer la construction des
      régions/trous.
    - **Découverte au passage, non anticipée** :
      `_split_for_platform_refinement` **crée lui-même** des
      auto-intersections. Il recolle deux coques issues de
      `slice_plane(cap=True)` dont les faces de bouchage sont coplanaires
      et superposées le long de la frontière de coupe — inoffensif sur les
      modèles simples, bloquant sur une géométrie dense qui touche ce plan
      en de nombreux points. Mesuré sur `death.stl` réparé : **107037
      tétraèdres SANS raffinement plateau, refus net AVEC**. D'où le
      troisième repli (abandonner le raffinement local plutôt que de ne
      produire aucun maillage).
    - **Résultat mesuré** : `death.stl` passe de « échec total » à un
      pipeline step1→step8 complet en **180s** (107037 tétraèdres, 222
      régions, 693 couches, 181682 points de trajectoire), pour **0.27 %
      d'écart de volume** et 1014 faces retirées par la réparation. Les
      deux modèles existants sont **strictement inchangés** (1913 et 1755
      tétraèdres, identiques avant/après) : la réparation ne se déclenche
      jamais sur eux.
    - `IngestReport` porte `self_intersecting_faces_before` (compteur
      géométrique informatif, jamais un déclencheur) — donnée explicite
      plutôt qu'un silence, même discipline que
      `CollisionReport.chuck_checked`.

---

## Triangulation dupliquée du plan de coupe plateau (PLC platform-refinement)

44. **`_split_for_platform_refinement` (limite #36) trianguleait le MÊME
    contour de coupe DEUX FOIS indépendamment, et les deux triangulations
    pouvaient diverger sur un contour non convexe riche — bug réel, PAS
    propre à un seul modèle, découvert en diagnostiquant un défaut visuel
    signalé par l'utilisateur sur `res/test_models/death.stl` (benchy
    retourné avec disque de support rond) : le disque de support, qui
    devrait être une pile de couches PLANES horizontales (posé à plat sur
    le plateau, aucune raison géométrique d'être incliné), apparaissait
    fragmenté en un grand nombre de petites couches inclinées façon
    « courbes de niveau ».**
    Diagnostic par élimination, confirmé par mesure directe (scripts
    jetables, mesh nettoyé/réparé) : `_tetrahedralize` (limite #43) tentait
    3 stratégies en cascade et retombait SILENCIEUSEMENT sur la 3ᵉ (« sans
    raffinement plateau ») parce que TetGen refusait le PLC combiné produit
    par `_split_for_platform_refinement` — réintroduisant l'erreur de
    discrétisation FEM de la limite #34 exactement dans la zone qui
    devrait rester plate. Cause racine : `slice_plane(plane_origin=...,
    plane_normal=..., cap=True)` était appelée indépendamment sur `slab`
    et `rest`, partageant exactement le même contour de coupe (mêmes
    sommets XY à `z=platform_z+margin`, écart mesuré = 0.0) mais le
    TRIANGULANT séparément — `trimesh.creation.triangulate_polygon` choisit
    ses diagonales indépendamment à chaque appel dès que le polygone est
    non convexe. Mesuré sur les deux modèles de référence du projet :

    | | sommets au contour de coupe | triangles divergents (slab vs rest) | TetGen |
    |---|---:|---:|---|
    | `l_shape_overhang.stl` | 8 | 2 / 6 | acceptait (juste sous le seuil) |
    | `death.stl` | 47 | 12 / 47 | **refusait** |

    Plus le contour à `z=platform_z+margin` est riche géométriquement (ici :
    disque de support + jambes du benchy = 47 sommets non convexes), plus
    les deux triangulations indépendantes divergent, et les triangles
    superposés mais non identiques sont vus par TetGen comme de VRAIES
    auto-intersections (diagonales croisées du même polygone,
    « Two segments exactly intersect » / « A vertex lies on a facet » /
    « A duplicated triangle... is ignored »). **Pas un défaut propre à
    `death.stl`** : `l_shape_overhang.stl` était déjà à la limite (2
    triangles divergents sur 6) avant ce fix, sans jamais avoir échoué —
    un contour un peu plus riche aurait suffi à basculer.
    **Fix retenu** (`_shared_cut_cap`, nouvelle fonction, et réécriture de
    `_split_for_platform_refinement`, `pipeline/step1_ingest.py`) :
    trianguler le contour de coupe UNE SEULE FOIS (reproduit la logique
    interne de `cap=True` — dédup des sommets, extraction des arêtes de
    bord sur le plan via `trimesh.grouping.group_rows(require_count=1)`,
    `triangulate_polygon` — mais appliquée une fois, pas deux), puis
    l'ajouter comme facette **interne partagée** entre `slab` et `rest`
    (une seule copie dans le PLC combiné final, pas une par coque). Les
    seeds (`slab_seed`/`rest_seed`) restent calculés via
    `_find_guaranteed_interior_point` (limite #42, inchangée), mais sur des
    coques watertight temporaires (coque ouverte + cap orienté vers
    l'intérieur de CETTE coque) construites uniquement pour le test
    d'inclusion — jamais utilisées pour le PLC final. Signature de
    `_split_for_platform_refinement` strictement inchangée, aucun appelant
    (`_tetrahedralize_once`) n'a eu besoin de changer.
    **Résultat mesuré après fix** : `death.stl` — TetGen accepte le PLC
    combiné avec raffinement plateau réel, **0 warning, 0
    auto-intersection**, 112698 tétraèdres (contre 107037 en mode dégradé
    « sans raffinement plateau » avant ce fix). Confirmé de bout en bout via
    `scripts/export_slicing_preview.py` : le pipeline complet ne retombe
    plus jamais sur la 3ᵉ tentative de `_tetrahedralize` pour ce modèle — la
    2ᵉ tentative (mesh réparé, AVEC raffinement plateau) réussit
    directement. `l_shape_overhang.stl` : 1677 tétraèdres, **identique** au
    chiffre déjà documenté pour le raffinement plateau actif — non-régression
    confirmée. Nouveaux tests (`tests/test_step1_ingest.py`) : fixture
    synthétique tripode (plateforme circulaire + 3 jambes cylindriques
    placées de façon irrégulière, seed déterministe) reproduisant le
    mécanisme de divergence sans dépendre du fichier externe `death.stl` —
    mesuré : 13 triangles divergents sur 66 avec l'ANCIEN mécanisme,
    suffisant pour faire échouer TetGen (« Internal TetGen error within
    `recoversubfaces` », une variante différente mais également fatale de
    l'échec observé sur `death.stl`) ; le NOUVEAU mécanisme réussit sur
    cette même fixture (24530 tétraèdres). Piège rencontré en écrivant le
    test de bout en bout : une comparaison naïve de volume moyen
    « tétraèdres près du plateau vs loin du plateau » (même critère que
    `test_tetrahedralize_with_near_platform_resolution_refines_locally`
    sur `l_shape_overhang.stl`, un solide continu) donnait le résultat
    INVERSE de celui attendu sur cette fixture multi-branches — les jambes
    fines isolées loin du plateau ont un volume tétraédrique intrinsèquement
    élevé (peu de matière à mailler), sans rapport avec
    `near_platform_resolution`. Corrigé en comparant deux runs indépendants
    (avec/sans raffinement) DANS LA MÊME ZONE PHYSIQUE (plateforme, rayon
    XY > 6, hors jambes), plutôt qu'une comparaison near/far à l'intérieur
    d'un seul run — valide uniquement sur un solide continu unique comme
    `l_shape_overhang.stl`, pas sur une géométrie multi-branches.

---

## platform_z figé à 0.0 au lieu d'une détection automatique

45. **`platform_z` restait figé à son défaut `0.0` partout dans le pipeline
    (jamais surchargé par `main.py` ni aucun script), cassant en cascade
    tout modèle dont le point le plus bas n'est pas exactement `z=0` — bug
    réel découvert APRÈS le fix de la limite #44 : le disque de support de
    `death.stl` restait fragmenté malgré un raffinement plateau désormais
    réellement appliqué.**
    Diagnostic par élimination successive, chaque étape mesurée avant de
    passer à la suivante (pas supposée) :
    - Hypothèse « espacement d'isovaleurs global mal calibré » (limite
      #12) : mesurée et **écartée** — `‖∇φ‖` sur le disque (médiane 0.634)
      n'était pas anormalement faible, plutôt 18% AU-DESSUS de la médiane
      globale (0.536), composante verticale dominante comme attendu sur un
      disque plat (0.565 vs 0.235 horizontale).
    - Hypothèse « mélange de directions » (`direction_blend_radius=10.0`,
      limite #32) : testée et **écartée** par désactivation directe
      (`direction_blend_radius=0.0`) — le chaos mesuré sur `φ` restait
      identique (range 7.52 contre 7.52), aucun effet.
    - Cause réelle, confirmée par mesure directe : à `z=1.245` (hauteur
      strictement identique sur 101 nœuds voisins), `φ` variait de -13.35 à
      -8.42 SANS corrélation avec `z` (`corrcoef=0.15`) ni avec le rayon
      (`corrcoef=-0.10`) — pas un gradient lisse, du bruit angulaire. Tracé
      jusqu'à `res/test_models/death.stl` lui-même : son point le plus bas
      mesuré est `z=1.245`, PAS `z=0`
      (`trimesh.load(...).vertices[:,2].min() == 1.2450000047683716`).
    - Conséquence en cascade sur `_compute_platform_contact`
      (`step3_segmentation.py`, seuil `max_z_per_face < platform_z +
      tolerance`, tolérance par défaut `0.5`) : avec `platform_z=0.0` figé,
      AUCUNE face n'était jamais reconnue comme contact plateau sur ce
      modèle (`1.245 > 0.5`) — vérifié : les 325 faces du dessous du
      disque, pourtant les plus plates et les plus basses du modèle entier,
      retournaient `is_platform_contact=False`.
    - Deux effets en cascade, mesurés séparément :
      1. `_select_region_directions`/`_region_growing` (step3) perdaient la
         garantie de direction verticale prioritaire pour ces faces
         (`residual_risk` nul pour toute direction candidate quand une face
         EST reconnue contact plateau, `argmin` retombant alors sur
         `directions[0]`=verticale par construction) — sans cette
         reconnaissance, une région du disque récupérait une direction à
         **49.2°** au lieu de la verticale, mesurée directement.
      2. Plus grave : `_compute_platform_level` (`step4_field_solver.py`)
         dérive le niveau du plateau UNIQUEMENT depuis
         `RegionDependencyGraph.platform_contact_faces` — vide sur ce
         modèle, elle retournait `None`, et `_identify_platform_nodes`
         renvoyait alors un tableau de nœuds ancrés **vide**. Avec
         `enable_platform_anchor=True` (défaut) mais 0 nœud ancré,
         l'ancrage plateau du solveur FEM devenait un no-op **silencieux**
         — tout le champ scalaire de `death.stl` était résolu sans aucun
         point de référence z=0 stable, cause directe du chaos mesuré.
    - Root cause finale : **un problème de configuration** (le modèle
      source a un offset Z non nul), **pas un bug structurel** de
      step3/step4 — mais rien dans le pipeline ne détectait ni ne
      signalait ce décalage, il cassait entièrement silencieusement.
    **Fix retenu** : `SegmentationParams.platform_z` et
    `step1_ingest.run(platform_z=...)` deviennent `Optional[float] = None`
    — même discipline opt-in déjà en place ailleurs
    (`near_platform_resolution=None`, `max_region_merge_angle_deg=None`).
    `None` déclenche une détection automatique
    (`mesh.vertices[:, 2].min()`), suivant le même précédent de conception
    déjà acté dans le projet (`_find_guaranteed_interior_point`/
    `_find_internal_cavities` : dériver depuis la géométrie plutôt
    qu'exiger un paramètre externe supposé correct). La valeur RÉELLEMENT
    utilisée par step1 (fournie ou déduite) est remontée explicitement via
    `IngestReport.platform_z_used` — jamais un silence sur une valeur
    dérivée, même discipline que `self_intersecting_faces_before`
    (limite #43) et `CollisionReport.chuck_checked` (étape 9).
    `main.py::run_pipeline` et
    `scripts/export_slicing_preview.py::run_slicing_only` propagent cette
    MÊME valeur résolue par step1 vers step3 (`dataclasses.replace`,
    uniquement si `params.segmentation.platform_z is None` — un override
    explicite de l'appelant n'est jamais écrasé), pour que les deux étapes
    utilisent exactement le même nombre plutôt que de le re-dériver
    indépendamment deux fois (source de divergence silencieuse potentielle
    sinon).
    **Résultat mesuré après fix, sur `death.stl`** (mêmes params calibrés
    que la limite #44) : `platform_contact_faces` passe de **0 à 326
    faces marquées** (couvre bien le dessous du disque). Le nombre de
    `Layer` distincts touchant la zone du disque (r_xy>8, z<3) passe de
    **62 à 2** — exactement l'ordre de grandeur physiquement attendu (2-3
    couches pour ~0.5mm d'épaisseur à `layer_height_base=0.2mm`), les deux
    couches restantes proprement séparées en hauteur
    (`z=[1.348,1.351]` puis `z=[1.554,1.561]`) au lieu du chaos mesuré
    avant ce fix. Nouveaux tests : 4 dans `tests/test_step1_ingest.py`
    (détection sur mesh translaté, override explicite respecté,
    non-régression stricte sur mesh déjà à z=0, raffinement plateau
    ciblant la vraie hauteur détectée plutôt que z=0 fige), 1 dans
    `tests/test_step3_segmentation.py` (détection automatique produit
    EXACTEMENT le même `platform_contact_faces` qu'un `platform_z`
    explicite correct, sur `flat_box_data` dont la base est à z=-5 — pas
    z=0), et un nouveau fichier `tests/test_main_pipeline.py` (2 tests :
    propagation step1→step3 telle qu'orchestrée par `main.py`, et respect
    d'un override explicite jamais écrasé par la valeur détectée). Suite
    complète : 187/187 tests passants (180 avant ce fix + 7 nouveaux),
    zéro régression sur tous les modèles déjà posés à `z=0` (détection
    automatique confirmée être un no-op strict sur ce cas).

---

## Trajets à vide sans routage ni vérification de collision

46. **🔴 PROBLÈME OUVERT, PAS ENCORE CORRIGÉ — les trajets à vide insérés par
    `step8_machine_toolpath._iter_schedule_points` sont des lignes droites
    SANS AUCUN ROUTAGE NI VÉRIFICATION DE COLLISION, et peuvent traverser
    tout droit le vide en plein milieu de la trajectoire, à des distances
    largement supérieures à ce qu'un mouvement à vide légitime devrait
    parcourir.** Découvert en diagnostiquant un artefact visuel signalé par
    l'utilisateur sur `out/death_preview.glb` (pics verts/cyan qui dépassent
    nettement de la silhouette du mesh, près d'un renflement/surplomb
    latéral) — l'utilisateur a lui-même identifié que peu importe si le
    segment est extrudé ou à vide, **la tête n'a structurellement aucune
    raison de se déplacer en ligne droite au travers du vide au milieu d'une
    trajectoire**, ce qui est le vrai problème, indépendant de la couleur de
    rendu GLB.
    Diagnostic par élimination, chaque étape MESURÉE (scripts
    `scripts/diagnose_contour_overflow.py`,
    `scripts/diagnose_toolpath_jumps.py`,
    `scripts/diagnose_toolpath_color.py`, conservés dans le dépôt) :
    - Hypothèse initiale (débordement géométrique de `Layer.contour_3d`,
      sortie step5) : **écartée par mesure directe** — 0/64698 points de
      contour de couche hors du solide source (tolérance 0.05mm,
      `trimesh.proximity.ProximityQuery`, jamais `mesh.contains()` qui a
      explosé en RAM — >4.6 Go et >15 min sans terminer sur ~165k points,
      process tué — ni `mesh_data.kdtree`, limite #15). Les couches
      elles-mêmes sont géométriquement correctes.
    - Cause réelle localisée dans `ToolpathMachine` (sortie step8) : sur
      165278 segments, longueur médiane 0.48mm, mais les 827 segments les
      plus longs (>13mm, percentile 99.5) sont à 790/827 des trajets à vide
      (`extrusion=0.0` sur au moins une extrémité), jusqu'à **44.5mm** de
      long en ligne droite. Vérifié point par point sur les 10 plus longs
      segments : la couleur rendue par `_toolpath_path3d` est bien CYAN
      (`TRAVEL_COLOR`) dans tous les cas mesurés — le code de coloration
      n'est pas en cause pour ces segments précis, mais ça ne change rien à
      la gravité du problème sous-jacent.
    - Mécanisme en cause : `_iter_schedule_points` (step8) insère un point
      `extrusion=0.0` à CHAQUE frontière de `(layer_id, segment_index)`,
      "y compris entre deux passes non contiguës" (docstring déjà présente
      dans le code, jamais traité comme un défaut jusqu'à ce diagnostic) —
      **aucune vérification que le chemin en ligne droite reste dans
      l'espace libre**, aucun lever en Z avant déplacement (contrairement à
      un slicer 3-axes classique), aucune tentative de longer la surface du
      solide.
    - Sur une machine réelle, ce trajet traverserait potentiellement la
      pièce déjà imprimée ou le plateau — step9 (collision) est justement
      l'étape censée détecter ce genre de trajet dangereux, mais elle est
      actuellement court-circuitée par `scripts/export_slicing_preview.py`
      (vitesse d'itération) — donc ce problème de routage n'est aujourd'hui
      détecté par AUCUNE étape du pipeline dans ce script, ni signalé par
      aucun rapport structuré (`Step8Diagnostic` ne couvre que l'IK, pas la
      distance parcourue à vide).
    - **Pas encore corrigé.** Pistes à trancher avec l'utilisateur avant
      implémentation (pas tranchées à la date de cette note) : (a) lever la
      buse en Z avant tout trajet à vide long, palliatif simple inspiré des
      slicers 3-axes classiques, ne résout pas un vrai évitement de
      collision ; (b) un routage réel de trajet à vide (le long de la
      surface ou via un chemin validé sans collision), chantier plus lourd,
      recoupe step9/step10 ; (c) au minimum, un signal explicite
      (nouveau diagnostic structuré, même discipline que `Step8Diagnostic`)
      quand un trajet à vide dépasse un seuil de longueur, pour que
      l'absence de routage ne soit plus un silence.

**Note ajoutée lors du tri CLAUDE.md/FIXED.md (2026-08-23)** : cette entrée
affirmait « Pas encore corrigé » au moment de sa rédaction. Ce n'est plus le
cas — la limite #49 (toujours dans `CLAUDE.md`) documente le fix retenu
avec l'utilisateur : un palliatif lift-and-traverse
(`pipeline/step8b_travel_routing.py::route_long_travels`), orchestré entre
step8 et step9 dans `main.py` et `scripts/export_slicing_preview.py`,
validé par les mêmes primitives de collision que step9. Le mécanisme est
actif dès que `PipelineParams.travel_routing` est renseigné. Un résidu
honnête (60% de collisions encore signalées sur géométrie difficile) reste
documenté dans la limite #49 elle-même — ce n'est donc pas un problème
totalement clos, mais ce n'est plus un problème "pas encore corrigé, pistes
non tranchées" comme l'affirmait ce texte : la piste a été tranchée et
implémentée.

Le paragraphe narratif dupliquant ce même diagnostic
(« 🔴 PROBLÈME OUVERT À CORRIGER DÈS QUE POSSIBLE — trajets à vide sans
routage (limite #46) ») a été retiré de `CLAUDE.md` en même temps que cette
entrée, pour la même raison — il était devenu stale de façon identique.

---

## Champ scalaire non monotone en Z sur zone très détaillée (death.stl)

47. **Sur une zone très détaillée de `death.stl` (motif festonné, z≈62-67mm,
    couches ~180-227), le champ scalaire `φ` n'est PAS monotone en Z —
    plusieurs isovaleurs très différentes correspondent à des positions Z
    quasi identiques, produisant des dizaines de `Layer` quasi-superposées
    (écart Z mesuré jusqu'à 0.0024mm, soit ~80× plus fin que
    `layer_height_base=0.2mm`) qui créent un effet visuel de "faisceau
    dense" dépassant du solide — PAS un vrai bug géométrique.** Découvert
    en creusant un artefact visuel signalé par l'utilisateur (captures
    d'écran Open3DViewer ET Blender, deux moteurs de rendu différents
    montrant le même phénomène — écarte d'emblée toute hypothèse de bug de
    rendu spécifique à un viewer).
    Diagnostic par élimination, chaque étape MESURÉE (scripts
    `scripts/diagnose_contour_overflow.py`, `diagnose_loop_extent.py`,
    `diagnose_contour_vs_rendered_mesh.py`, `diagnose_layer_spacing.py`,
    `diagnose_gradient_variance.py`, `diagnose_isovalue_vs_z.py`, conservés
    dans le dépôt) :
    - Débordement géométrique de `Layer.contour_3d` hors du solide :
      **écarté**, max 0.0083mm mesuré (négligeable), y compris en comparant
      contre `result.mesh_data` (le mesh RÉELLEMENT affiché dans le GLB,
      pas le STL brut — les deux se sont avérés identiques sur ce modèle,
      pymeshfix jamais déclenché).
    - Boucles en "aller-retour" (zigzag anormal, ratio périmètre/bbox
      élevé) : **écarté**, ratio max 5.3 (normal pour un contour crénelé
      sur un motif festonné dense), corrélé à aucun débordement.
    - Divergence du gradient local `‖∇φ‖` par rapport à la médiane globale
      utilisée pour calibrer `delta` (limite #12) : **écarté** — médiane
      locale mesurée dans la zone (0.5306) quasi identique à la médiane
      globale (0.5063, ratio 1.05×), espacement Z théorique par isovaleur
      adjacente correctement calibré (~0.19mm, cohérent avec
      `layer_height_base=0.2mm`).
    - **Cause confirmée** : les paires de couches quasi-superposées en Z
      ont des isovaleurs très éloignées (`gap_isovalue` jusqu'à 3.14, soit
      ~31× `delta≈0.10`) et souvent des `region_ids` différents (ex: couche
      211 `regions=[9,11,17]` isovalue=-18.93 contre couche 180
      `regions=[9,11]` isovalue=-22.07, à seulement 0.013mm l'une de
      l'autre en Z). Le champ résolu par step4 "ondule" plusieurs fois en
      Z dans cette zone très détaillée — chaque couche individuelle reste
      géométriquement correcte (collée à la surface à <0.01mm), mais le
      tri DFS spatial de continuité (limite #40) place ces couches
      d'isovaleurs très différentes côte à côte dans la séquence finale
      parce qu'elles sont proches SPATIALEMENT, produisant l'illusion
      visuelle de dizaines de couches empilées au même endroit.
    - Recoupe probablement les limites #32/#34 (erreur de discrétisation
      FEM sur maillage grossier relatif à la complexité géométrique
      locale) et #36 dernier paragraphe (raffinement local `add_region`
      limité au cas plateau, aucun raffinement dédié aux jonctions
      multi-région/zones très détaillées) — non confirmé par mesure
      directe dans cette session, piste la plus probable pour une future
      investigation si ce point devient prioritaire.
    - **Fix tenté en session ultérieure, voir limite #48** — deux
      approches testées, l'une validée sur jonction localisée mais pas sur
      `death.stl` (zone trop étendue), l'autre sans effet mesuré.

---

## Diagnostic de jonction + raffinement local aux jonctions multi-région

48. **Livrable double pour la limite #47 : (1) un diagnostic structuré
    permanent, fonctionnel et validé sur tout modèle ; (2) un raffinement
    local de maillage aux jonctions, validé sur jonction LOCALISÉE mais
    confirmé NON APPLICABLE sur `death.stl` (la zone dense y couvre la
    quasi-totalité du volume) — deux approches de repli testées et sans
    effet mesuré sur ce modèle précis. Session longue, plusieurs bugs réels
    trouvés et corrigés en cours de route, aucun connu au préalable.**

    ### Livrable 1 — diagnostic structuré (`Layer.junction_variance_score`)
    `LayerParams.detect_junction_variance` (opt-in, `False` par défaut —
    même discipline que `split_disjoint_clusters=None`) formalise en
    diagnostic permanent la métrique calculée manuellement pendant le
    diagnostic de la limite #47 : pour chaque paire de couches consécutives
    (triées par Z), marque une anomalie si `gap_z < junction_variance_z_ratio
    * layer_height_base` **et** `gap_isovalue > junction_variance_isovalue_ratio
    * delta`. Champ pur sur `Layer`, ne modifie jamais `contour_3d`/
    `toolpath_2d`. **Validé et fonctionnel** : détecte 27/228 couches
    anormales sur `death.stl` (cohérent avec le diagnostic manuel de la
    limite #47), et reproduit l'anomalie sur la fixture `tripod_data` sans
    mélange de directions (`tests/test_step5_layers.py::
    test_junction_variance_detects_anomaly_at_tripod_junction`).
    `scripts/diagnose_isovalue_vs_z.py` réécrit pour l'utiliser (option
    `--assert-max-anomalies N`, non-régression exploitable en script).

    ### Livrable 2 — raffinement local aux jonctions (généralisation d'`add_region`)
    Nouveau module `pipeline/step1b_junction_detection.py::detect_junction_zones` —
    fonction pure `(MeshData, RegionDependencyGraph, angle_threshold_deg,
    margin) -> list[JunctionZone]`, identifie les faces frontière entre deux
    régions dont les `target_direction` divergent au-delà d'un seuil (même
    rôle que `SegmentationParams.max_region_merge_angle_deg`, valeur
    séparée), regroupe en clusters spatiaux. Nouvelle dataclass
    `core.mesh_types.JunctionZone` (centroïde + rayon).

    `main.run_pipeline`/`scripts/export_slicing_preview.py::run_slicing_only`
    orchestrent une **double passe explicite** (jamais cachée dans un
    step) : une passe légère (mesh + step2 + step3 SEULEMENT, jamais la
    tétraédrisation fine ni step4-10) localise les zones AVANT la vraie
    tétraédrisation. Nouveaux paramètres `PipelineParams.
    junction_refinement_resolution`/`_angle_deg`/`_margin` (`None` par
    défaut = désactivé, même discipline que `near_platform_resolution`).
    **Bug réel trouvé en écrivant ce fix** (pas anticipé) :
    `scripts/export_slicing_preview.py::run_slicing_only` avait sa PROPRE
    liste d'arguments `ingest_run(...)` distincte de `main.run_pipeline`,
    et ne transmettait initialement JAMAIS les nouveaux paramètres — deux
    runs "avant/après" via ce script produisaient un résultat identique,
    pas parce que le mécanisme ne fonctionnait pas, mais parce qu'il
    n'était jamais exercé. Corrigé (le script réplique maintenant la même
    double passe que `main.run_pipeline`).

    `pipeline/step1_ingest.py::_tetrahedralize_once` généralise
    `TetGen.add_region` : jusqu'à 3 régions actives simultanément (slab,
    junction, rest). Contrairement au plateau, une jonction n'a pas de plan
    de coupe unique — pas de `_shared_cut_cap` nécessaire, juste un seed
    interne supplémentaire.

    **Trois bugs réels de TetGen/logique trouvés et corrigés durant
    l'implémentation, aucun anticipé par le plan initial :**
    1. **Seed "rest" pouvant tomber DANS le rayon d'une zone de jonction**
       (le centroïde naïf du mesh n'a aucune raison d'être hors de la
       zone) — mesuré sur le tripode synthétique : centroïde à 3.42mm du
       centre d'une zone de rayon 6.0mm, causant une explosion de densité
       sur TOUT le maillage (905 → 44985 tétraèdres, loin de toute
       jonction : 34 → 2302 tétraèdres locaux) — même mécanisme de fond
       que la limite #42 (seed hors de sa région), pas un nouveau bug
       indépendant.
    2. **`TetGen.add_region` est sensible à l'ORDRE d'appel, pas
       seulement aux seeds** — découverte non documentée nulle part
       ailleurs dans le projet ni dans la doc TetGen : la région "rest"
       (grossière) DOIT être `add_region(1, ...)`, sinon la densité
       explose sur tout le maillage même avec des seeds individuellement
       corrects. Mesuré sur le tripode : `rest=1,junction=2` → 3080
       tétraèdres corrects ; `junction=1,rest=2` → 46000+, `far<3mm` passe
       de 34 à 1150-2277 selon la variante testée. **Confirmé aussi sur le
       cas plateau-seul déjà en production** (`slab=1,rest=2` vs
       `rest=1,slab=2` sur un modèle synthétique : `far<3mm` 847 vs 29) —
       mais **volontairement PAS corrigé pour ce cas historique**, qui
       reste validé sur `l_shape_overhang.stl` (1677 tétraèdres, identique
       au chiffre déjà documenté, non-régression confirmée) : inverser
       l'ordre partout aurait été un changement de comportement non
       mesuré sur un cas déjà en production. Seul le NOUVEAU cas à 3
       régions applique l'ordre correct (`rest` toujours en premier).
    3. **Fusion de zones proches non itérative** — sur `death.stl` (motif
       festonné, 24 régions), la détection brute produit 301 zones
       (distance médiane au plus proche voisin 1.18mm). Une première passe
       de fusion (`_merge_nearby_zones_once`) les réduisait à 2 zones dont
       la distance centre-à-centre (28.8mm) restait pourtant **en deçà**
       de la somme de leurs rayons (57.3mm) — chevauchement résiduel non
       détecté par une seule passe, parce que le rayon d'un groupe fusionné
       (englobant chaque membre à sa distance réelle) peut rapprocher deux
       groupes qui ne se touchaient pas via leurs zones D'ORIGINE. Corrigé
       par `_merge_nearby_zones` : itère `_merge_nearby_zones_once` jusqu'à
       convergence (nombre de zones stable) — 301 → 1 zone unique sur
       `death.stl` après ce fix.

    ### Confirmé non applicable sur `death.stl` — deux approches testées
    Après tous les fix ci-dessus, la zone unique fusionnée sur `death.stl`
    a un rayon de **47mm sur un modèle de 30×61×50mm** — la zone dense
    couvre la quasi-totalité du volume. `_find_rest_seed_outside_junction_zones`
    (réécrite pour chercher le sommet du mesh le plus éloigné de la zone,
    plutôt qu'un rejection-sampling aléatoire qui échouait après 4096
    essais malgré des seeds valides mais rares) lève une `MeshIngestError`
    **explicite** (pas un échec silencieux) : même le sommet le plus
    éloigné du mesh reste à l'intérieur de la zone (distance mesurée
    -12.73mm au bord de la zone). Le mécanisme point+rayon ne peut
    structurellement pas s'appliquer quand la zone dense n'est pas une
    petite jonction localisée mais une caractéristique dominante d'une
    grande partie du volume.
    **Approche de repli testée** : raffinement GLOBAL (`resolution` plus
    fine, sans aucun `add_region`) plutôt que local — **sans effet
    mesuré** : `resolution=0.7` (1.44× plus fin que le défaut 1.5,
    152889 tétraèdres contre 105880) donne 28/226 anomalies contre 27/228
    à `resolution=1.5` — quasiment identique, pas d'amélioration. Cohérent
    avec le précédent déjà documenté à la limite #40 (dernier point) :
    *"augmenter la résolution du maillage ×30 ne réduit pas du tout"* un
    problème structurellement similaire sur le modèle en T — confirme que
    la piste "juste raffiner plus" n'est structurellement pas la réponse à
    ce type de non-monotonie, ni en local ni en global, sur ce modèle
    précis.

    ### Validé et fonctionnel sur jonction LOCALISÉE (fixture tripode)
    `tests/test_step1_ingest.py::test_junction_refinement_increases_density_at_junction_not_elsewhere` :
    la densité de tétraèdres augmente nettement (>2×) près de la jonction
    physique et reste quasi inchangée (<2×) loin d'elle, sur la même
    fixture tripode que la limite #44 (pied + 3 bras non coplanaires,
    jonction physique compacte au sommet du pied). Le mécanisme fonctionne
    donc pour son cas d'usage visé — une jonction ponctuelle entre régions,
    pas une zone détaillée étendue sur une grande partie du modèle. La
    proportion de modèles utilisateur qui ressembleront à l'un ou l'autre
    cas n'est pas connue (un seul modèle dense/réaliste disponible dans le
    projet à ce stade, `death.stl` — pas d'échantillon pour estimer cette
    proportion).

    **Suite complète de tests : 197/197 passants** (187 avant cette
    session, +10 nouveaux : 2 diagnostic structuré
    `test_step5_layers.py`, 3 détection de jonction isolée
    `test_step1b_junction_detection.py`, 3 raffinement de maillage
    `test_step1_ingest.py`, 2 orchestration `test_main_pipeline.py`).
    Non-régression stricte confirmée sur `l_shape_overhang.stl` (1677
    tétraèdres, identique) — le nouveau chemin de code n'est jamais
    emprunté sans activation explicite de `junction_refinement_resolution`.

    Scripts de diagnostic conservés dans le dépôt (tous
    `scripts/diagnose_*`/`measure_*` cités ci-dessus), cohérent avec la
    discipline déjà en place : chaque hypothèse mesurée, jamais supposée,
    y compris (et surtout) les hypothèses qui se sont révélées fausses.

---

## step1 ne survivait pas à un crash natif de TetGen (isolation sous-processus)

50. **step1 doit survivre à N'IMPORTE QUEL défaut de maillage susceptible de
    faire échouer la tétraédrisation — y compris un CRASH NATIF de TetGen,
    pas seulement une exception Python — car un tel crash n'est pas
    rattrapable in-process : isolation de chaque tentative dans un
    sous-processus dédié.**

    ### Le vrai bug (pas un blocage, un crash)
    Sur `death.stl` modifié par l'utilisateur (2340 triangles, résolution
    réduite), le process se terminait par un crash natif (`EXIT=139` /
    `0xC0000005` selon la couche d'exécution) dès la 1ʳᵉ tentative de la
    cascade de réparation réactive déjà en place (limite #43 : tel quel →
    réparation `pymeshfix` → sans raffinement plateau), juste après les
    mêmes avertissements TetGen déjà documentés ("Two segments exactly
    intersect", "A vertex lies on a facet"). Le mécanisme existant
    intercepte `RuntimeError` — un mécanisme Python — mais un segfault/access
    violation tue tout l'interpréteur AVANT qu'aucun `except` ne puisse
    s'exécuter : la réparation `pymeshfix`, censée s'appliquer juste après,
    n'avait donc **jamais la moindre chance de tourner**, quel que soit son
    efficacité propre sur ce mesh.

    ### Fix retenu : isolation par sous-processus
    Nouveau module `pipeline/_tetrahedralize_worker.py` — exécute UNE
    tentative de `_tetrahedralize_once` dans un sous-processus dédié,
    lancé via `python -m pipeline._tetrahedralize_worker <input.pkl>
    <output.pkl>` (mesh sérialisé en vertices/faces bruts, pas un objet
    `Trimesh` complet — plus robuste au pickling inter-version). Nouvelle
    fonction `_tetrahedralize_once_isolated` (`pipeline/step1_ingest.py`) :
    lance le worker, détecte un code de sortie non nul (ou fichier de
    sortie absent/illisible) comme un crash natif et lève
    `_TetgenCrashError` (nouvelle exception, sous-classe de `RuntimeError`) ;
    une `RuntimeError`/`MeshIngestError` levée PROPREMENT par le worker est
    relayée telle quelle, comportement inchangé pour ce cas déjà géré.

    **CHAQUE tentative de la cascade** (y compris la première, "tel quel")
    tourne désormais dans son propre sous-processus — pas seulement en
    dernier recours. Nécessaire par construction : si seule la dernière
    tentative était isolée, un crash sur la 1ʳᵉ tentative (in-process)
    tuerait quand même tout le pipeline avant d'atteindre la logique de
    repli, exactement le bug observé sur `death.stl`. Décision actée
    explicitement avec l'utilisateur après clarification de cette
    contrainte technique (isoler seulement le dernier recours ne résout
    pas le cas réel qui a motivé le fix).

    ### Bug d'ordre des `except` trouvé ET corrigé pendant l'implémentation
    Premier essai sur le vrai cas : le crash était bien détecté
    (`_TetgenCrashError` levée avec le bon code de sortie), mais remontait
    quand même jusqu'à l'appelant au lieu de déclencher la tentative
    suivante. Cause : dans la boucle de cascade de `_tetrahedralize`,
    `except RuntimeError as exc: ... raise` était placé AVANT `except
    _TetgenCrashError`. Comme `_TetgenCrashError` est une sous-classe de
    `RuntimeError`, Python retient le PREMIER `except` applicable — toute
    instance de `_TetgenCrashError` tombait dans la branche générique, dont
    le message ("Le sous-processus TetGen s'est terminé anormalement...")
    ne contient jamais le texte `"self-intersection"`, déclenchant le
    `raise` immédiat prévu pour une VRAIE erreur non gérée. Corrigé en
    plaçant `except _TetgenCrashError` en premier — même discipline "mesuré
    pas supposé" que le reste du projet : le bug n'a été trouvé qu'en
    testant réellement le cas motivant, pas anticipé par relecture.

    ### Coût accepté, mesuré
    Chaque tentative relance un interpréteur Python complet (réimport de
    `tetgen`/`scipy`/`trimesh`) :

    | | avant isolation | après isolation |
    |---|---:|---:|
    | `tests/test_step1_ingest.py` (21 tests) | — | **16.84s** |
    | suite complète (206 tests) | **12.92s** | **77.16s** |

    Overhead ~×6 sur la suite complète — jugé acceptable au vu de la
    garantie de robustesse demandée explicitement par l'utilisateur
    ("step1 doit fixer TOUTES les erreurs... qui peuvent causer un échec
    de la tétraédralisation"), pas une dégradation de performance
    silencieusement acceptée sans mesure.

    ### Résultat sur le vrai cas (`death.stl` modifié par l'utilisateur)
    Après le fix : **step1 aboutit** (plus de crash), le pipeline progresse
    jusqu'à l'étape 6/8 (208 couches produites par step5) avant d'échouer
    sur un bug **préexistant et sans rapport**, découvert par cette
    progression nouvelle : `pipeline/step6_toolpath2d.py::_classify_contours`
    lève `KeyError: 1` — un trou topologiquement imbriqué dans un autre
    trou peut voir son "parent immédiat" (plus petit conteneur par aire)
    être lui-même classé comme trou, alors que `groups` n'a des clés que
    pour les contours NON-trou (`groups = {i: [] for i in range(...) if
    not is_hole[i]}`). Jamais rencontré avant sur aucun modèle de
    référence du projet — ce mesh est le premier à produire une couche
    avec cette configuration topologique précise.

    **CORRIGÉ dans la foulée** (même session, l'utilisateur ayant ensuite
    demandé de générer le GLB de `death.stl` — objectif bloqué par ce
    KeyError). Diagnostic complet, mesuré pas supposé :
    - Nouveau script permanent `scripts/extract_nested_hole_case.py` —
      instrumente `_classify_contours` (déjà corrigé, donc ne plante plus)
      pour DÉTECTER la condition pathologique sans planter, et exécute le
      pipeline complet sur `death.stl` pour capturer le vrai cas fautif.
      Résultat capturé sur la vraie couche fautive : 5 boucles, dont
      `containers = {0: [], 1: [0], 2: [1]}` / `is_hole = {0: False, 1:
      True, 2: True}` — la boucle 2 (trou) a pour **unique** conteneur la
      boucle 1 (elle-même un trou), **aucun contour de matière n'apparaît
      dans sa chaîne de conteneurs** — signature d'une isosurface
      légèrement incohérente topologiquement (mesh dense, motif déjà
      documenté ailleurs dans ce fichier comme source de ce genre
      d'irrégularité).
    - **Fix** (`pipeline/step6_toolpath2d.py::_classify_contours`) :
      `parent = min(matter_parents, key=...)` où `matter_parents =
      [p for p in containers[i] if not is_hole[p]]`, au lieu de
      `min(containers[i], key=...)` sans filtre. Si `matter_parents` est
      vide (comme le cas réel ci-dessus), le trou est ignoré — même
      discipline défensive déjà en place juste au-dessus dans la même
      fonction pour un trou sans AUCUN conteneur.
    - **Test de non-régression fidèle**, vérifié reproduire le vrai bug
      (pas juste un cas synthétique plausible) :
      `tests/test_step6_toolpath2d.py::test_classify_contours_hole_directly_nested_in_another_hole` —
      géométrie synthétique construite pour matcher EXACTEMENT le motif
      `containers`/`is_hole` réel (anneau non-convexe dont un lobe dépasse
      le contour extérieur, second trou placé dans ce lobe), confirmée par
      exécution manuelle de l'ancienne logique non filtrée AVANT d'écrire
      le test définitif (`KeyError: 1` reproduit à l'identique), puis
      confirmée passante avec le code corrigé.
    - **Résultat final sur `death.stl`** (`scripts/export_slicing_preview.py`,
      régénère `out/death_preview.glb`) : step1→step8 aboutit intégralement
      en **61.1s**, 208 couches, **206523 points de trajectoire**, **0
      orientation IK inatteignable**, GLB de **9 696 828 octets** (9.7 Mo).

    ### Tests
    `tests/test_step1_ingest.py` (21/21) + `tests/test_step6_toolpath2d.py`
    (21/21, dont le nouveau test ci-dessus) + suite complète
    (**207/207 passants**, 77.98s — 206 + 1 nouveau test). Non-régression :
    le comportement de la cascade sur un mesh déjà valide (jamais de
    crash, jamais de `RuntimeError`) est strictement inchangé — seul le
    CHEMIN d'exécution change (sous-processus au lieu d'in-process), jamais
    le résultat ; le comportement de `_classify_contours` sur toute couche
    déjà correctement imbriquée (le cas normal, largement majoritaire) est
    lui aussi strictement inchangé — seul le cas dégénéré (trou sans
    parent matière) change de comportement, passant d'un crash à un
    contournement silencieux du trou incohérent (même discipline que le
    cas "trou sans conteneur" déjà géré juste au-dessus dans le code).

---

## turkish.stl : trois défauts de maillage en cascade (coque non dominante, recoversubfaces)

51. **`res/test_models/turkish.stl` (10180 triangles, jamais traité avec
    succès avant cette session) a révélé TROIS défauts de maillage
    distincts en cascade — chacun corrigé, avec un premier fix qui s'est
    lui-même révélé dangereux et a dû être rejeté avant d'être remplacé.**
    Suite directe de la demande utilisateur "step1 doit fixer TOUTES les
    erreurs... qui peuvent causer un échec de la tétraédralisation"
    (limite #50) : chaque défaut a été traité au fur et à mesure qu'il
    bloquait la génération, jamais anticipé.

    ### Défaut 1 — mesh non-watertight, réparation standard insuffisante
    64 faces cassées sur 10180 : `trimesh.repair.fill_holes`/`fix_normals`/
    `fix_winding` ne suffisaient pas à rendre le mesh étanche.
    `_tetrahedralize` refusait immédiatement (`MeshIngestError`) SANS
    jamais essayer TetGen ni la cascade réactive existante (limite #43),
    qui suppose déjà un mesh étanche en entrée.

    **Premier fix tenté et REJETÉ** : ajout de `_repair_with_pymeshfix`
    (MeshFix, déjà présent dans le projet pour un autre usage — d'ailleurs
    renommée depuis `_repair_self_intersections` pour refléter sa portée
    réelle) comme repli. Accepté à tort comme un succès lors d'une première
    vérification superficielle (`is_watertight=True`) — mesure plus
    poussée : `MeshFix.repair()` réduisait le mesh de **10180 à seulement
    10 faces** (extentes ~2×1.8×1.1mm au lieu de 100×100×60.5mm), une perte
    quasi totale de la géométrie. Root cause réelle : le mesh combiné
    contient 16 fragments minuscules (76 faces, 0.31% du volume, volumes
    **négatifs** — normales inversées, probablement des artefacts d'export
    CSG) accolés au corps principal (10104 faces, **déjà individuellement
    watertight** avec les bonnes dimensions) par des jonctions
    non-manifold — ce n'étaient pas des "trous à combler", MeshFix
    échouait à les gérer proprement et détruisait tout.

    **Fix retenu** : nouvelle fonction `_extract_dominant_shell`
    (`pipeline/step1_ingest.py`) — isole via `mesh.split(only_watertight=False)`
    la plus grosse coque INDIVIDUELLEMENT watertight, si elle représente au
    moins `_DOMINANT_SHELL_VOLUME_RATIO=0.95` (95%) du **volume** total
    (mesuré sur turkish.stl : coque dominante à 99.69% du volume). Le seuil
    est basé sur le VOLUME, pas le nombre de faces — découvert pendant la
    conception des tests que le ratio de faces n'est PAS un critère fiable
    (un fragment 1000× plus petit en volume peut porter le même nombre de
    triangles qu'un cube primitif de référence). Tentée **dès le début**
    de `_clean_and_repair`, avant même le test `is_watertight` — couvre
    aussi le cas "techniquement watertight au sens strict de trimesh mais
    structurellement étrange" (rencontré en écrivant les tests : une
    coïncidence de sommet unique entre coques ne casse pas toujours le
    comptage d'arêtes manifold de trimesh, contrairement à l'intuition).
    `_repair_with_pymeshfix` reste en dernier recours, avec un **nouveau
    garde-fou anti-collapse** : son résultat n'est accepté que si
    `len(faces) >= 0.5 * len(faces_avant)` — empêche la répétition exacte
    de l'incident ci-dessus sur un futur mesh différent. Nouveau champ
    `IngestReport.dominant_shell_fragments_dropped` — donnée explicite,
    jamais une perte silencieuse de géométrie (même discipline que
    `self_intersecting_faces_before`).

    ### Défaut 2 — erreur TetGen "recoversubfaces" non reconnue par la cascade
    Une fois le mesh nettoyé (défaut 1 corrigé), TetGen levait `"Internal
    TetGen error within recoversubfaces"` — un message **déjà documenté**
    dans ce fichier (limite #44) comme "une variante différente mais
    également fatale" du refus "self-intersection" d'origine, mais
    **jamais relié au code de détection** de la cascade réactive : celle-ci
    ne cherchait que la sous-chaîne `"self-intersection"` dans le message
    d'erreur (`_tetrahedralize`), donc ne se déclenchait pas sur cette
    variante — l'erreur remontait brute, sans jamais essayer la réparation
    pymeshfix pourtant déjà fonctionnelle et disponible.

    **Fix** : nouvelle fonction `_is_recoverable_tetgen_error` généralisant
    la détection à un tuple de sous-chaînes reconnues
    (`_RECOVERABLE_TETGEN_ERROR_SUBSTRINGS = ("self-intersection",
    "recoversubfaces")`), utilisée par la boucle de cascade de
    `_tetrahedralize` à la place du test de sous-chaîne unique précédent —
    la cascade (réparation pymeshfix, puis sans raffinement plateau) a
    ainsi la même chance de résoudre les deux variantes, cohérent avec le
    fait qu'elles étaient déjà documentées comme la même classe de défaut.

    ### Résultat final sur `turkish.stl`
    step1→step8 aboutit intégralement en **63.6s**, 117 couches, **138759
    points de trajectoire**, **0 orientation IK inatteignable**, GLB de
    **6 172 648 octets** (6,2 Mo).

    ### Tests
    9 nouveaux tests dans `tests/test_step1_ingest.py` :
    `test_repair_with_pymeshfix_closes_holes_standard_repair_misses`,
    `test_clean_and_repair_falls_back_to_pymeshfix_when_standard_repair_fails`,
    `test_ingest_succeeds_on_mesh_pymeshfix_can_repair_but_standard_repair_cannot`
    (défaut 1, repli pymeshfix générique) ;
    `test_extract_dominant_shell_keeps_main_body_drops_negligible_fragments`,
    `test_extract_dominant_shell_returns_none_when_no_clear_dominance`,
    `test_clean_and_repair_reports_dropped_fragments` (défaut 1, coque
    dominante — garde-fou de non-dominance vérifié explicitement) ;
    `test_is_recoverable_tetgen_error_matches_self_intersection`,
    `test_is_recoverable_tetgen_error_matches_recoversubfaces`,
    `test_is_recoverable_tetgen_error_rejects_unrelated_error` (défaut 2).
    **216/216 tests passants** (207 avant cette investigation + 9
    nouveaux). Non-régression : `l_shape_overhang.stl` et
    `l-triple-axis-overhanged-shape.stl` (déjà propres, `split()` n'y
    révèle jamais qu'une seule composante) ne déclenchent ni
    `_extract_dominant_shell` ni la cascade `_is_recoverable_tetgen_error`
    — chemins de code inchangés sur tout mesh déjà valide.

---

## Garde-fou anti-effondrement pymeshfix manquant dans la cascade de tétraédrisation

52. **🔴 Un deuxième appel à `_repair_with_pymeshfix`, propre à la cascade de
    tétraédrisation (`_tetrahedralize`, limite #43) et INDÉPENDANT de celui
    déjà gardé dans `_clean_and_repair` (limite #51), n'avait AUCUN
    garde-fou anti-effondrement — il a laissé passer un résultat MeshFix
    qui tronquait silencieusement les DEUX TIERS SUPÉRIEURS de
    `turkish.stl`, produisant un GLB visuellement absurde sans jamais lever
    d'exception ni émettre le moindre avertissement. `turkish.stl` reste
    NON SLICEABLE à ce jour (`MeshIngestError` honnête, aucune réparation
    non-destructive trouvée pour son défaut réel).**

    Découvert sur signalement utilisateur direct (capture d'écran de
    `out/turkish_preview.glb` accompagnée de « C'est une CATASTROPHE !!! ») :
    la trajectoire rendue montrait un enchevêtrement dense de lignes
    vert/cyan concentré dans un coin du modèle, le reste (une rampe
    graduée en degrés, visuellement proéminente) totalement vierge de tout
    toolpath.

    **Diagnostic par élimination, chaque étape mesurée** (scripts
    `scripts/diagnose_turkish_manifold_repair.py`,
    `diagnose_turkish_manifold_tetgen.py`, `diagnose_turkish_pymeshfix_tetgen.py`,
    conservés dans le dépôt) :
    - Le toolpath rendu (`ToolpathMachine.points`) s'arrêtait net à
      `z=63.26mm`, alors que le mesh source s'étend jusqu'à `z=90.78mm` —
      86,6 % des points de trajectoire concentrés dans les 12 premiers mm
      du modèle (histogramme Z mesuré : 68020+52094 points sur 138759 dans
      les 2 premiers déciles de hauteur sur 10).
    - Remontée jusqu'à la source : `mesh_data.vertices` (sortie de
      `_clean_and_repair`, donc APRÈS le fix de la limite #51) couvre bien
      la pleine hauteur `z=[30.26, 90.78]` — le défaut n'est PAS dans le
      nettoyage de surface, déjà corrigé. Mais `tet_mesh.tet_vertices`
      (sortie de `_tetrahedralize`) ne couvre que `z=[30.26, 63.26]` — la
      troncature apparaît strictement ENTRE ces deux étapes, à l'intérieur
      de `_tetrahedralize` elle-même.
    - Reproduction directe de chaque branche de la cascade de repli
      (limite #43) :
      - **« tel quel »** : `TetGen` refuse ("51 input triangles are
        skipped due to self-intersections") puis crashe nativement (code
        de sortie détecté comme crash par `_TetgenCrashError`, limite #50).
      - **« sans raffinement plateau »** (mesh brut, sans
        `_split_for_platform_refinement`) : crashe identiquement — le
        défaut n'est pas propre au PLC combiné plateau, il préexiste sur
        le mesh nettoyé lui-même.
      - **« après réparation pymeshfix »** : `_repair_with_pymeshfix(mesh)`
        (l'appel de la ligne ~947, à l'intérieur de `_tetrahedralize`,
        PAS celui déjà gardé dans `_clean_and_repair`) réduit le mesh de
        **10104 à 1488 faces** — `is_watertight=True`, donc silencieusement
        accepté par la cascade — mais sa hauteur Z s'effondre de
        `[30.26, 90.78]` à `[30.26, 63.26]` : MeshFix a supprimé toute la
        moitié supérieure du modèle (la rampe graduée) plutôt que de
        réparer localement les triangles fautifs. TetGen tétraédrise ce
        mesh mutilé **sans une seule erreur ni un seul warning** — c'est un
        maillage parfaitement valide, juste géométriquement faux. Le
        pipeline entier (step2 à step8) traite ensuite fidèlement ce
        modèle amputé, produisant 117 couches et 138759 points de
        trajectoire tous confinés au tiers inférieur — un résultat
        « réussi » à tous les niveaux de vérification automatique
        existants (0 exception, 0 warning TetGen, 0 orientation IK
        inatteignable, 216/216 tests toujours verts), détecté uniquement
        par inspection visuelle humaine du GLB final.
    - **Cause structurelle** : la garde anti-effondrement
      (`len(faces) >= 0.5 * len(mesh.faces)`) ajoutée à la limite #51 ne
      protégeait QUE l'appel de `_clean_and_repair` (nettoyage de surface,
      avant tétraédrisation) — l'appel totalement indépendant de
      `_tetrahedralize` (retenter TetGen après un premier refus) avait été
      ajouté séparément (limite #43, avant l'incident de la limite #51) et
      n'avait jamais reçu cette même garde. Deux call sites du même
      mécanisme de repli, une seule protégée — angle mort classique de
      fix localisé à un seul endroit d'un mécanisme dupliqué.

    **Fix retenu** : le même garde-fou (`len(repaired.faces) < 0.5 *
    len(mesh.faces)` ⇒ résultat rejeté, `repaired = None`) appliqué à
    l'appel de `_tetrahedralize` — et propagé à TOUTE la cascade, pas
    seulement à la tentative dédiée « après réparation pymeshfix » : le
    repli suivant (« sans raffinement plateau », ligne ~976) utilisait
    aussi `repaired` comme source si non `None`, ce qui aurait réutilisé le
    même mesh mutilé même après un rejet explicite. `repaired = None`
    explicite en cas d'effondrement propage donc le rejet à toute la
    cascade en aval, un seul point de vérité.

    **Résultat mesuré après fix** : `step1_ingest.run()` sur `turkish.stl`
    lève désormais `MeshIngestError` explicite au lieu de réussir
    silencieusement sur une géométrie tronquée — comportement conforme à
    la discipline du projet (échec honnête plutôt que succès silencieux
    sur donnée corrompue, même principe que `CollisionReport.chuck_checked`,
    `IngestReport.self_intersecting_faces_before`).

    **Recherche d'une vraie réparation non-destructive, deux pistes
    testées et écartées, mesurées pas supposées** :
    - `pymeshfix` (déjà en place) : rend le mesh watertight mais en
      supprimant 87 % de sa géométrie (10104→1488 faces) — inutilisable.
    - `manifold3d` (`mesh.union(mesh, engine="manifold")`, déjà utilisé
      ailleurs dans le projet, limite #43) : préserve fidèlement la
      géométrie (10104→10372 faces, hauteur Z intacte `[30.26, 90.78]`) —
      contrairement à pymeshfix, AUCUNE perte de matière. Mais ne résout
      PAS le défaut sous-jacent : TetGen crashe toujours sur ce résultat
      avec exactement la même erreur (`Internal TetGen error within
      \`recoversubfaces\``, "Two segments exactly intersect" à répétition)
      — l'union booléenne recompose une géométrie topologiquement propre
      mais ne corrige pas les micro-contacts exacts (arêtes coïncidentes,
      sommets dupliqués à distance nulle) qui bloquent la reconstruction
      de frontière de TetGen. Confirme et étend la limite #43 : "TetGen
      est le seul juge fiable de ce qui le bloque" — un mesh topologiquement
      sain au sens de trimesh/manifold3d peut rester un mesh que TetGen
      refuse.
    - **Aucune réparation non-destructive trouvée à ce jour.** `turkish.stl`
      reste donc **NON SLICEABLE** — limitation honnête et documentée,
      pas contournée par un fix qui masquerait le problème (le
      comportement d'avant ce fix — silencieusement produire un résultat
      tronqué — aurait été objectivement pire que refuser explicitement).

    **`out/turkish_preview.glb` supprimé** (il montrait exactement le
    résultat tronqué diagnostiqué ci-dessus — laisser un artefact
    silencieusement faux dans le dépôt serait trompeur).

    ### Tests
    `tests/test_step1_ingest.py::test_tetrahedralize_cascade_rejects_collapsed_pymeshfix_repair` —
    monkeypatch `_repair_with_pymeshfix` pour simuler un effondrement
    (icosphere source à 320 faces, faux résultat une boîte à 12 faces,
    même ratio de perte que l'incident réel) sur un mesh par ailleurs
    valide, vérifie que `_tetrahedralize` tétraédrise correctement la VRAIE
    géométrie (icosphere) plutôt que d'accepter le faux résultat tronqué —
    reproduit le mécanisme du bug, pas seulement son symptôme sur
    `turkish.stl` (fichier externe, non utilisable dans un test unitaire
    rapide). **217/217 tests passants** (216 avant ce fix + 1 nouveau).
    Non-régression : sur tout mesh dont la cascade ne déclenche jamais le
    repli pymeshfix (le cas de loin le plus courant, tous les modèles de
    référence du projet), ce chemin de code reste inchangé.

**Note ajoutée lors du tri CLAUDE.md/FIXED.md (2026-08-23)** : cette entrée
affirmait à sa clôture « `turkish.stl` reste NON SLICEABLE à ce jour... 
Aucune réparation non-destructive trouvée à ce jour ». Ce n'est plus vrai —
une session ultérieure a résolu ce cas via un mécanisme distinct, déjà en
place dans `pipeline/step1_ingest.py` et `pipeline/_self_intersection_resolve.py` :
`resolve_self_intersections` (résolution GÉOMÉTRIQUE des auto-intersections,
préservant la topologie, contrairement à pymeshfix) suivi de
`_finish_geometric_resolution`/`_probe_tetgen_conflicts` (sondage de TetGen
lui-même en mode détection `switches='d'` dans un sous-processus isolé,
`pipeline/_tetgen_detect_worker.py`, pour cibler et réparer localement les
derniers segments nommés en conflit). La garde anti-effondrement décrite
dans cette entrée (le vrai sujet du fix ci-dessus) reste, elle, toujours
active et correcte — c'est uniquement l'affirmation finale sur le statut de
`turkish.stl` qui était devenue obsolète, sans jamais avoir été corrigée
dans `CLAUDE.md`.

---

## IK grille structurellement inatteignable sur modèle dense (56965/181682, 31.4%) — IK analytique

### Performance mesurée sur modèle DENSE (première mesure réelle)

Jusqu'ici tous les chiffres de performance du projet portaient sur des mesh
synthétiques de ~25 triangles. `res/test_models/death.stl` (12662 triangles,
débloqué par la limite #43) donne la première mesure sur une géométrie
réaliste, et elle **invalide une extrapolation naturelle** :

| | T-shape | death.stl | facteur |
|---|---:|---:|---:|
| triangles | 30 | 12662 | ×422 |
| tétraèdres | 1913 | 107037 | ×56 |
| points trajectoire | 15178 | 181682 | ×12 |
| **temps total** | **11.7s** | **180s** | **×15** |

Le temps ne suit PAS le nombre de triangles (×422 pour ×15). Mais il ne
suit pas non plus une seule étape : **le régime change complètement**.
step8 passe de 90 % à 70 % du total, et step5 explose de ×36 (0.77s →
27.96s) pour devenir le deuxième poste. step1 (11.8s) et step4 (9.4s),
négligeables sur les petits modèles, pèsent désormais 12 % à eux deux.

**Point ouvert et préoccupant** : **56965 orientations IK inatteignables
sur 181682 (31.4 %)**, contre 2.8 % sur le T-shape. Cause probable
identifiée par mesure, non encore corrigée : la grille IK de step8
(`ik_samples_per_axis=10`, 100 directions) a une **erreur médiane de
9.26°** alors que `max_unreachable_angular_error_deg=5.0` — la méthode de
résolution échoue structurellement son propre critère de rejet.

**IK analytique implémentée (session ultérieure) — résout le point ouvert
ci-dessus, mais PAS le gain de vitesse initialement espéré.** La piste
chiffrée envisageait une fermeture analytique pure (`B=acos(dz)`,
`C=atan2(dx,-dy)`, ~15000×) — valable seulement pour EXACTEMENT 2 axes
orthogonaux, donc rejetée : `_resolve_axis_angles` est documenté
génériquement pour N axes ("jamais 2 codés en dur", limite #18), un
fast-path 2-axes aurait contredit ce principe. Conçu à la place :
`core.geometry_utils.resolve_axis_angles_analytic`, une réduction par
« paire externe » — les axes internes (`axes[2:]`) sont repliés à une
valeur fixe (héritée du point précédent, `previous_angles`), puis la
paire externe (`axes[0]`, `axes[1]`) est résolue par bissection (angle
interne, borné par balayage vectorisé + bissection jusqu'à 1e-5°) suivie
d'un `atan2` dans le plan perpendiculaire à l'axe externe (angle
externe). Pour N=2 le repliement est vide : la méthode retombe
**identiquement** (pas approximativement) sur la fermeture 2-axes, d'où
la précision mesurée quasi identique à celle espérée.

**Risque confirmé réel avant l'implémentation** (pas supposé) : pour
N≥3, l'invariant géométrique utilisé (angle entre la direction et
`axes[0]` indépendant de l'angle d'`axes[0]`) ne suffit PAS à déterminer
une solution unique — comptage de degrés de liberté (cible = 2
contraintes sur S², N inconnues ⇒ famille de solutions de dimension N-2
dès N≥3). Sur des axes aléatoires non orthogonaux, un cas concret a été
produit où la paire externe seule est géométriquement incapable
d'atteindre la cible (résiduelle 55° avant repli) : l'angle atteignable
par rapport à `axes[0]` en balayant l'axe interne reste borné à un
ANNEAU `[|α-β|, α+β]` (α = angle entre la référence effective et l'axe
interne, β = angle entre les deux axes externes), pas `[0°,180°]` en
général. D'où un **repli explicite vers la méthode grille+raffinement
historique** (conservée intacte, jamais supprimée — toujours utilisée
indépendamment par `step10_adjustment._try_local_reorientation`),
déclenché uniquement quand le bracketing de la paire externe ne trouve
aucune racine. Jamais rencontré sur la machine à 2 axes orthogonaux
actuelle (le repliement est vide, la paire externe EST le système
complet).

**Résultat mesuré sur les deux modèles réels** (comparaison directe même
code, avant/après, `git stash`) :

| | l-triple-axis | death.stl |
|---|---:|---:|
| step8 avant (grille) | 13.02s | 142.12s |
| step8 après (analytique) | 6.38s | 103.94s |
| facteur | 2.0× | 1.4× |
| orientations inatteignables avant | 426 (2.8%) | 56965 (31.4%) |
| orientations inatteignables après | **0 (0.0%)** | **0 (0.0%)** |

Le gain de vitesse (1.4-2×) est très en-deçà des ~15000× de la piste
initiale — coût structurel de la bissection généralisée à N axes (~20-30
itérations par point, chacune une rotation de Rodrigues scalaire),
confirmé par profilage et non un défaut d'implémentation résiduel (deux
passages de profilage ont éliminé les surcoûts `np.cross`/`np.dot`/
`np.linalg.norm` génériques sur des vecteurs 3D triviaux, même correction
que la limite #37 déjà appliquée ailleurs dans le fichier). **Mais le
vrai objectif — l'échec structurel du critère de rejet — est résolu** :
0% d'orientations inatteignables sur les deux modèles, contre 31.4% sur
`death.stl` avant.

Le docstring de `_resolve_axis_angles` affirmant qu'une IK analytique
est « impossible à écrire pour toute configuration d'axes arbitraire »
a été corrigé — cette affirmation était fausse, comme identifié dans la
piste chiffrée initiale.

**Correction (session ultérieure)** : le paragraphe ci-dessus annonçait
« aucun nouveau test dédié » — corrigé dans la foulée, 8 nouveaux tests
ajoutés (`tests/test_geometry_utils.py`) : round-trip N=2 machine réelle,
N=1 et axe verrouillé, N=3 sur axes non orthogonaux (robustesse et
reconstruction exacte par repliement), axes parallèles, déclenchement du
repli grille. 174 tests au total après cet ajout.

---

## step5 : marching tetrahedra vectorisé (boucle Python par tétraèdre, np.cross générique)

### step5 vectorisé — marching tetrahedra sans boucle Python par tétraèdre

Objectif posé par l'utilisateur, session ultérieure : le pipeline complet
sur `death.stl` ≤ 4s (contre ~157s dans la config calibrée ci-dessus).
Mesuré et présenté avant tout code : même en vectorisant step5 entièrement
(le plus gros poste isolable après step8, déjà traité), le total resterait
très au-dessus de la cible — décision actée d'avancer sur ce gain mesuré
sans promettre 4s, puis remesurer et redécider de la suite.

**Cause** : `_march_tetrahedron` était appelée ~321 500 fois (une fois par
tétraèdre candidat, par isovaleur) sur `death.stl`, 89 % du temps de step5
en cumtime. `np.cross` appelé **421 545 fois sur des vecteurs 3×3
individuels** via `_orient_triangle` — même anti-pattern de dispatch numpy
générique déjà documenté et corrigé pour l'IK (limite #37) : le surcoût
de validation/permutation d'axes (`moveaxis`, `normalize_axis_tuple`)
domine une opération triviale appelée en boucle Python plutôt qu'en lot.

**Fix** : `_march_tetrahedra_batch` (`pipeline/step5_layers.py`) traite
**tous** les tétraèdres candidats d'une isovaleur en une seule passe
vectorisée — zéro boucle Python sur les tétraèdres ou les triangles.
Trois tables de correspondance statiques (construites une fois à l'import)
remplacent les indexations scalaires par ligne :
- **Cas "cap"** (`n_above ∈ {1,3}`, 1 triangle/tet) : `lone_idx` via
  `np.argmax` sur un masque booléen `(M,4)` (une seule ligne vraie par
  construction), `rest_idx` via `_CAP_EXCLUDE_TABLE` `(4,3)`.
- **Cas "quad"** (`n_above == 2`, 2 triangles/tet) : motif booléen encodé
  en entier 4 bits, deux tables `(16,2)` `_QUAD_TABLE_A`/`_QUAD_TABLE_B`.
- Les deux masques sont mutuellement exclusifs et exhaustifs (`n_above`
  n'a pas d'autre valeur possible pour un tétraèdre à 4 sommets), traités
  séparément puis concaténés avec leurs propres `local_tet_idx` pour
  préserver l'alignement `source_tets` (nécessaire à l'association
  couche→région en aval).

`_interpolate_edge_batch`/`_orient_triangle_batch` sont des batchings
directs des formules existantes, **aucun changement de logique** — la
justification géométrique de `reference_direction` (champ P1, `phi(hi) -
phi(lo) = grad(phi)·(p[hi]-p[lo])` exactement) reste strictement identique.
`_march_tetrahedron`/`_interpolate_edge`/`_orient_triangle` (scalaires)
sont **conservées intactes**, servent désormais d'oracle de test.

**Conception validée empiriquement AVANT l'implémentation** (pas seulement
tracée à la main) : un agent Plan dédié a transcrit son pseudo-code dans
un script jetable, exécuté contre `_march_tetrahedron` réel comme oracle
sur 30 essais / ~7129 tétraèdres synthétiques (mélange valeurs continues,
discrètes pour forcer des égalités, et dégénérées) — **0 écart** sur le
nombre de triangles, les coordonnées et l'ordre de winding. Le point de
risque signalé (ordre `reshape`/`repeat` du cas "cap", convention
tet-majeur/arête-mineur) a été revérifié directement dans le fichier
source (pas seulement accepté de l'agent) avant intégration.

Nouveau test permanent `test_march_tetrahedra_batch_matches_scalar_loop`
(`tests/test_step5_layers.py`), même méthodologie, plus deux tests de cas
limites (tous dégénérés, entrée vide). Aucun test d'intégration existant
modifié — le contrat public de `_extract_isosurface_triangles`/`run()` est
inchangé.

**Résultat mesuré** (config calibrée, comparaison directe avant/après) :

| | l-triple-axis | death.stl |
|---|---:|---:|
| step5 avant | ~0.78s | ~28s |
| step5 après | **0.22s** | **7.81s** |
| facteur | 3.5× | **3.6×** |

**Effet de bord non expliqué, mesuré reproductible (2 exécutions
identiques à la seconde près)** : step8 est passé de ~104s à ~51s sur
`death.stl` entre la mesure de la session précédente et celle-ci, sans
qu'aucun code de step8 n'ait été touché entre les deux. Cause non
investiguée (variance machine/thermique probable sur la mesure
précédente, pas un artefact de step5) — à garder en tête pour ne pas
sur-interpréter un futur écart de mesure sur step8 comme un effet réel
d'un changement de code.

**Nouveau total mesuré, pipeline complet step1→step8** :

| | l-triple-axis | death.stl |
|---|---:|---:|
| Avant cette session (IK seule) | ~4.8s* | ~157s |
| Après step5 vectorisé | **4.83s** | **80.38s** |

*(le chiffre "avant" pour l-triple-axis n'a jamais été mesuré isolément
avec l'IK analytique seule sans step5 vectorisé — la comparaison directe
n'existe que sur death.stl, où le gain total est net : 157s → 80s,
**×1.96**.)*

**Objectif de 4s toujours hors d'atteinte sur `death.stl`** (80s restants,
facteur ~20× encore requis) — pas une mission accomplie. step8 (~51s) et
step1+step4 (~18s combinés) restent les postes dominants. Prochaine
question à trancher avec l'utilisateur, une fois ce résultat communiqué :
réduire la résolution du maillage (coût qualité non mesuré à cette
échelle) vs un autre chantier algorithmique.

**Mesure sur une AUTRE machine (session ultérieure, `.venv` recréé à neuf,
Python 3.14)** : `scripts/export_slicing_preview.py res/test_models/death.stl`
— step1→step8 en **46.1s**, total script (avec chargement Python/TetGen/export
GLB) **51.9s**. Non comparable aux chiffres ci-dessus (tous mesurés sur la
même machine entre eux) — matériel différent, à ne PAS lire comme un gain
ou une régression de performance. 111978 tétraèdres, 0 warning TetGen
(confirme les fixes #43/#44 en usage réel), 176768 points de trajectoire,
**0 orientation IK inatteignable** (confirme le fix IK analytique).

---


---

## Ordre d'impression par hauteur de centroïde au lieu de la structure réelle (pièce pendante)

`_build_dependency_graph` (step3) décidait l'ordre entre régions adjacentes
en comparant la **hauteur de leurs centroïdes globaux** projetés sur
`up_vector` : la plus basse imprimée en premier. La spécification d'origine
de l'étape 3 prévoyait pourtant « détermine par ray-cast/test de contact qui
doit être imprimée avant qui » — l'heuristique de hauteur était une
simplification posée là, jamais la cible.

**Défaut réel** : la règle est fausse pour toute pièce **pendante**. Cas
rapporté par l'utilisateur — un T retourné avec une excroissance suspendue
sous une barre : l'excroissance est géométriquement plus basse, donc
imprimée *en premier*, alors qu'elle est accrochée à la barre et ne peut
exister qu'après elle. Constat structurel associé : deux régions voisines se
touchent par une **arête**, côte à côte et jamais empilées — aucune
comparaison de hauteur entre elles, globale ou locale, ne peut exprimer
« qui supporte qui ».

**Tentative écartée avant le fix retenu** (implémentée, mesurée, auditée,
puis annulée — conservée ici car elle documente pourquoi la voie évidente ne
marche pas) : test de contact par ray-cast vertical à la frontière partagée.
Deux défauts **mesurés** la rendaient arbitraire :
- les impacts de rayon n'étaient jamais attribués à une région — sur
  `l-triple-axis-overhanged-shape.stl`, l'impact décisif entre les régions 3
  et 0 tombait sur la face 9, appartenant à la région **2**, une troisième
  région sans rapport ;
- le seuil de contact, dérivé des longueurs d'arêtes, valait 20,49 / 34,44 /
  22,00 sur un modèle haut de 22 unités et 42 à 66 sur un modèle haut de 44 —
  un seuil supérieur à la pièce entière transforme le test en « y a-t-il
  quoi que ce soit en dessous ».
Les tests écrits pour cette tentative passaient, mais pour de mauvaises
raisons : sur la fixture pendante, dans les 4 paires de frontière, l'impact
décisif était le pendant touchant **son propre fond** (région 1 → région 1).
Leçon actée : valider le MÉCANISME, pas seulement l'assertion finale.

**Fix retenu — arborescence de croissance depuis le plateau.** Abandon total
du critère de hauteur. Parcours en largeur multi-sources
(`_compute_print_depths`, `nx.single_source_shortest_path_length` sur un
nœud source virtuel relié à toutes les racines) depuis les régions au
contact du plateau, à travers l'adjacence **non orientée** entre régions ;
`profondeur[région]` = nombre de sauts de région depuis le plateau ; chaque
arête d'adjacence est orientée de la profondeur faible vers la forte, et les
arêtes entre régions de même profondeur sont supprimées (régions sœurs
poussant en parallèle, aucune contrainte physique entre elles). Le cas
pendant est résolu **par construction** : le pendant n'atteint le plateau
qu'à travers son point d'accroche, donc profondeur strictement supérieure.

**Acyclicité désormais garantie structurellement** (la profondeur croît
strictement le long de chaque arête conservée) — garantie plus forte que
celle de la référence globale `up_vector` historique, qui existait
précisément pour éviter les cycles artificiels sur un anneau de régions.
Corollaire : la suppression forcée des arêtes entrantes des régions plateau
n'a plus lieu d'être (une racine est à profondeur 0, le minimum). Le
contrôle de cycle + `SegmentationError` sont conservés comme garde-fou de
contrat, désormais structurellement inatteignables.

**Deux cas limites traités explicitement, tous deux réels** :
- *aucune racine plateau* — le critère de contact exige que les **3 sommets**
  d'une face tiennent sous `platform_z + tolerance`, donc un objet en contact
  **tangent** (sphère, arête, pointe) en produit légitimement zéro ; déjà
  rencontré en production sur `death.stl`. Repli : la région dont le point le
  plus bas le long de `up_vector` est le plus bas sert de racine (minimum par
  région, pas la moyenne — ce qui compte pour démarrer une croissance est le
  point le plus proche du plateau).
- *régions inatteignables* — une cavité interne forme une nappe de surface ne
  partageant aucune arête avec la coque extérieure (cas nominal, cf.
  `_find_internal_cavities`) ; un mesh multi-corps aussi. Chaque composante
  connexe sans racine est traitée indépendamment avec sa propre racine de
  repli ; les profondeurs ne sont comparables qu'à l'intérieur d'une
  composante.

Nouveau champ `RegionDependencyGraph.region_print_depth: Optional[dict[int, int]]`
(`core/mesh_types.py`, défaut `None`) — rend l'arborescence inspectable
plutôt que seulement déductible des arêtes. Rétro-compatible : toutes les
constructions manuelles existantes (tests step1b/step4/step5) l'omettent.

**Validation du mécanisme AVANT les tests** (fixture « table sur pied +
excroissance pendante s'arrêtant au-dessus du plateau ») — arborescence
mesurée : `plateau(0,2) → pied(3) → dessous de table(4) → {dessus table(1),
parois pendant(6)} → fond pendant(5)`, profondeurs `0,0,1,2,3,3,4`. Les
parois du pendant sont à z≈22, plus bas que la table (z=30), et reçoivent
pourtant la profondeur 3 ; le fond du pendant (z=14, géométrie non-plateau
la plus basse du modèle) est imprimé **en dernier**.

**Preuve que les nouveaux tests attrapent l'ancien défaut** : l'ancienne
heuristique rejouée à l'identique sur cette fixture produit `(6→4)` et
`(5→6)` — pendant imprimé avant la table dont il pend, chaîne entière
inversée ; le nouveau mécanisme donne `4→6→5`.

**Mesures sur modèles réels** :

| modèle | arêtes (nouveau/ancien) | inversées | supprimées | profondeur | rangs step5 |
|---|---:|---:|---:|---:|---:|
| `l-triple-axis-overhanged-shape.stl` (5 régions) | 4 / 4 | **0** | 0 | 2 | 3 / 3 |
| `l_shape_overhang.stl` (1 région) | 0 / 0 | 0 | 0 | 0 | 1 / 1 |
| `death.stl` (22 régions) | 25 / 32 | **6** | 7 | 5 | 6 / 11 |

Aucun changement gratuit : le modèle en T et `l_shape_overhang` sont
strictement inchangés. Sur `death.stl`, les 6 inversions sont **le cas du
pendant sur un vrai modèle** — la région 3 est le corps principal (717 faces,
z=[61.5, 91.2], centroïde 79.3) et les régions 2/4/5/7/8/9 sont de minuscules
détails (2 à 6 faces) posés dessus, chacun dans une fine bande z *à
l'intérieur* du corps ; leur centroïde étant plus bas que celui du corps,
l'ancienne règle les imprimait avant le corps qui les porte.

**Impact aval mesuré** (`death.stl`, step1→step8) : trajet à vide total
23388,1 → **23354,4** (−0,14 %), grands sauts (>13 mm) 719 → **718**, plus
long saut inchangé (58,53 mm), 200 couches et 202416 points identiques. Le
fix corrige une erreur d'ordre réelle mais **n'améliore pas mesurablement les
trajets** — attendu, les 6 régions concernées faisant 2 à 6 faces chacune.
Fix de justesse, pas de performance : dit tel quel plutôt que survendu.

**Tests** : `tests/test_step3_segmentation.py` — 4 nouveaux
(`test_pendant_is_printed_after_what_it_hangs_from`,
`test_platform_regions_are_roots_of_the_tree`,
`test_depths_are_defined_for_every_region_even_without_platform_contact`,
`test_internal_cavity_regions_are_reachable_in_the_tree`) + fixture
`pendant_table_data`. Les 14 tests step3 existants passent sans
modification, dont `test_sphere_dependency_graph_has_no_cycle` (0 racine,
anneau de régions) et les deux tests exigeant qu'une région plateau n'ait
aucune arête entrante — désormais garanti par construction. **226/226** au
total (222 avant, zéro régression).

---

## Champ scalaire à croissance descendante — matière déposée dans le vide

**Symptôme observé par l'utilisateur** : sur le GLB de `death.stl` coloré par
ordre d'impression, une partie haute du modèle apparaissait imprimée **avant**
les parties censées la soutenir, elles-mêmes colorées comme imprimées en
dernier.

### Diagnostic (mesuré, jamais supposé)

Test de support sur la trajectoire machine réelle : pour chaque point extrudé,
y a-t-il de la matière **déjà déposée** en dessous (cylindre de 1 mm de rayon,
1,5 mm de hauteur) ? Grille 3D incrémentale — seuls les points d'indice
inférieur sont dans la grille au moment du test, donc « déjà imprimé » est
exact par construction.

- **22 393 / 198 999 points extrudés (11,25 %) déposés sans rien en dessous.**
- Sur ces 22 393, **99 % ont de la matière en dessous, mais imprimée plus tard**
  (défaut d'ordre corrigible) ; 1 % n'a aucun support à aucun moment (surplomb
  géométrique réel).

Ce test en Z est valide sur cette machine : le plateau ne tourne qu'autour de
Z, or tourner autour de l'axe de la gravité ne change jamais quel côté est
« le bas » pour la pièce. Le bas dans le repère pièce est donc toujours −Z,
y compris pendant l'impression.

Symptôme le plus net, vérifiable sans instrument : **l'impression de
`death.stl` commençait à z = 59,1 et descendait jusqu'au plateau à z = 51,5.**

**Cause racine** : 67,6 % du VOLUME du TetMesh avait `grad(φ)·up < 0`. Un champ
dont le gradient descend empile ses couches du haut vers le bas — chaque couche
atterrit dans le vide. Deux contributions :
- `step3` choisissait parfois une direction de région à composante Z négative
  pour une région non-racine (mesuré : 4 régions sur 22, dont la plus grande
  du modèle à 717 faces, et une racine plateau).
- `_blend_tet_directions` (step4) mélange les directions cibles de régions
  voisines par proximité spatiale **sans jamais garantir que le mélange reste
  montant** — mesuré : 61,6 % du volume avait une direction cible descendante.

### Mécanisme ÉCARTÉ (implémenté, mesuré, puis retiré)

Première approche, issue du plan initial : pénalité quadratique **repondérée
itérativement** — majorer `tet_weights` sur les tétraèdres dont
`grad(φ)·up < 0`, re-résoudre, recommencer. Le système restait linéaire
symétrique défini positif (aucun solveur QP introduit, aucune nouvelle
dépendance).

**Mesure : ça AGGRAVE le problème.** 67,6 % → 70,3 % de volume descendant
(3 passes), 69,3 % (6 passes).

La raison est structurelle, pas un défaut de réglage :
`_assemble_gradient_direction_term` pousse `grad(φ)` **VERS** `tet_directions`.
Majorer le poids d'un tétraèdre descendant **renforce donc sa cible** — et
61,6 % du volume a une cible qui pointe elle-même vers le bas. La pénalité
amplifiait exactement ce qu'elle devait supprimer. **Le levier correct est la
CIBLE, pas le poids.**

### Mécanisme RETENU

`FieldSolverParams.enforce_upward_growth` (opt-in, défaut `False` ;
activé dans `main.PipelineParams`) : les directions cibles par tétraèdre dont
la composante selon `up_vector` est négative sont **retournées** (`d → -d`)
avant l'assemblage du terme de direction
(`step4_field_solver._enforce_upward_directions`).

**Sans perte géométrique** : `φ` et `-φ` ont exactement les **mêmes
isosurfaces**, donc les mêmes couches. Retourner une direction cible ne déforme
aucune couche — cela n'inverse que le **sens de parcours** de la pile,
c'est-à-dire précisément le défaut corrigé. Coût de calcul nul (transformation
de tableau, aucune résolution supplémentaire).

Une direction exactement perpendiculaire à `up` est laissée telle quelle : elle
ne descend pas, et la retourner rendrait le résultat sensible au bruit
numérique du mélange spatial.

**Diagnostic explicite** : `ScalarField.upward_growth_diagnostic`
(`UpwardGrowthDiagnostic`) porte le nombre/volume de tétraèdres redressés et le
**résidu descendant mesuré sur le champ RÉELLEMENT résolu** — pas sur les
cibles. Redresser les cibles ne garantit pas mathématiquement un gradient
montant partout (le lissage FEM peut encore produire une descente locale), donc
le résidu est mesuré, jamais supposé. `None` = mécanisme désactivé, à
distinguer d'un résidu nul. Fractions pondérées par volume, jamais un compte
brut de tétraèdres (un maillage tétraédrique n'est pas uniforme).

### Résultat mesuré de bout en bout (`death.stl`)

| | avant | après |
|---|---:|---:|
| volume TetMesh à `grad(φ)·up < 0` | 67,64 % | **0,06 %** |
| **matière extrudée déposée en l'air** | **11,25 %** (22 393 / 198 999) | **0,57 %** (1 381 / 240 558) |
| couches | 200 | 257 |
| trajet à vide total | 23 354 mm | 25 950 mm (+11 %) |
| temps de slicing step1→8 | 66 s | 77 s (+17 %) |

Le nombre de couches augmente parce que le champ redressé produit un empilement
cohérent sur toute la pièce, au lieu de replier des zones sur elles-mêmes. Le
trajet à vide augmente de 11 % — coût assumé et mesuré, pas masqué.

### Tests

`tests/test_step4_field_solver.py`, 7 nouveaux tests — dont un test unitaire
**pur** de l'opérateur creux `G_up` découplé du solveur (vérifie que
`(G_up @ φ)` égale exactement la projection sur `up` du gradient calculé par
`_compute_per_tet_gradient`), et un test comparant le diagnostic rapporté à une
mesure **recalculée indépendamment** depuis le champ retourné.

Piège rencontré et corrigé pendant l'écriture des tests : la fixture initiale
ne forçait la direction descendante que sur les régions hors plateau — or les
fixtures synthétiques simples de ce fichier ne produisent qu'**une seule**
région, elle-même au contact du plateau. Le test passait sans jamais exercer le
mécanisme. Un garde-fou explicite (`assert descending_before > 0.2`, « fixture
invalide … le test ne discriminerait rien ») a détecté le problème et reste en
place pour empêcher la régression silencieuse de la fixture.

**233/233 tests** (226 avant, zéro régression).

---

## Étape 9b — vérification du support de dépôt sur la trajectoire réelle

Nouveau module `pipeline/step9b_support_check.py`, intercalé entre step8/step8b
et step9. Répond à une question qu'aucune étape existante ne posait : **chaque
point extrudé a-t-il de la matière DÉJÀ IMPRIMÉE sous lui au moment où il est
déposé ?**

### Pourquoi step2 ne suffisait pas

`step2_overhang` détecte les surplombs à partir de la seule **géométrie** du
modèle — il ne connaît pas l'ordre d'impression. Or une face parfaitement
supportée dans la pièce finie peut être déposée dans le vide si la matière qui
la porte est imprimée **après** elle. C'est précisément le défaut mesuré sur
`death.stl` : 11,25 % de la matière extrudée déposée sans rien en dessous, dont
99 % avait pourtant un support — imprimé trop tard.

### La donnée qui décide quoi corriger

`UnsupportedDepositEvent.has_future_support` sépare deux natures de défaut qui
appellent des corrections **incomparables** :
- **`True` — défaut d'ORDRE** : la matière porteuse existe, déposée trop tard.
  Corrigible par réordonnancement, sans toucher à la géométrie.
- **`False` — surplomb RÉEL** : aucun support à aucun moment de la trajectoire.
  Aucun ordre d'impression ne le résout ; il faut du support matériel.

Sans ce distinguo, un taux global de « matière en l'air » ne dit pas s'il faut
réordonner ou ajouter du support.

`SupportCheckReport.checked` est un signal explicite d'absence de vérification
(même famille que `CollisionReport.nozzle_model_checked`) : une liste `events`
vide serait sinon ambiguë entre « vérifié, rien à signaler » et « jamais
vérifié ». Le dénominateur de `unsupported_fraction` est le nombre de points
**réellement testés**, jamais la longueur de trajectoire — sinon la fraction
serait diluée par les trajets à vide, qui ne déposent rien.

### Réutilisation, jamais duplication

`_PrintedGeometryIndex`, `_grow_printed_index_to` et `_point_to_segment_distance`
sont importés **tels quels** depuis step9 — même patron que `step10_adjustment`,
qui importe déjà ces primitives. Aucune seconde implémentation de l'état
incrémental de la matière imprimée n'existe dans le projet. L'index grossit
**avant** chaque test, donc il ne contient que la matière strictement
antérieure, et jamais le cordon qui se termine au point testé (confondu avec la
buse — le compter comme support serait un faux négatif structurel).

### Validité du test vertical

Le test cherche du support vers le bas. C'est correct ici parce que le plateau
ne tourne qu'autour de l'axe de la gravité : tourner autour de cet axe ne change
jamais quel côté est « le bas » pour la pièce. Sur une machine capable de
basculer réellement le plateau, `up_vector` devrait être exprimé par point de
trajectoire — signalé dans la docstring du module plutôt que supposé en silence.

### Trois problèmes de performance rencontrés et corrigés (mesurés à chaque étape)

Sur le pire cas (`death.stl` non corrigé, 199 000 points testés, 22 393
événements) :

| état | coût |
|---|---:|
| version initiale (boucle Python par capsule) | **191 s** |
| + test de support vectorisé en numpy | 58,9 s |
| + grille spatiale dédiée (`grid_cell_size_multiplier=1.0`) | **41,2 s** |
| classification ordre/surplomb, version initiale | **ne terminait pas en 10 min** |
| classification restructurée | **9,7 s** |

1. **Boucle Python par candidat.** La grille retourne 432 candidats en médiane
   (jusqu'à 2338) par requête, pour ~199 000 requêtes. Vectorisé en une passe
   numpy sur tous les candidats — même nature de correction que
   `_points_to_segment_distances` de step9, pour exactement la même raison.
2. **Grille dimensionnée pour la mauvaise requête.** Le défaut de step9 (4,0 ×
   `nozzle_diameter` = 1,6 mm) est calibré pour l'AABB d'une tête de 76 mm ;
   une requête de support porte sur ~2×2×1,5 mm. Mesuré : 360 candidats
   médians à 4,0 contre **84 à 1,0**, à coût de construction inchangé. Même
   nature d'ajustement dédié que `routing_grid_cell_size_multiplier` de step8b,
   pour la raison symétrique (requêtes bien plus larges là-bas).
3. **Classification en O(événements × segments).** La version initiale
   re-testait tous les événements en attente à chaque segment. Restructurée en
   deux passes de coûts très différents : la **classification** ordre/surplomb
   se fait sur un index de l'état FINAL de la pièce, construit une fois — coût
   O(événements). La **localisation** du point porteur exact
   (`supporting_point_index`) exige de rejouer la croissance incrémentale, donc
   O(défauts × segments) : elle est désormais **opt-in et désactivée par
   défaut** (`locate_supporting_point`), alors que la distinction
   ordre/surplomb — la donnée qui décide quoi corriger — reste toujours
   calculée.

### Résultat mesuré (`death.stl`)

| | champ non corrigé | champ corrigé (défaut) |
|---|---:|---:|
| matière extrudée en l'air | **11,28 %** (22 452 / 198 999) | **0,53 %** (1 267 / 240 558) |
| dont défaut d'ORDRE | 22 243 | 1 073 |
| dont surplomb RÉEL | 209 | 194 |

Le module retrouve **11,28 %** sur la baseline, contre **11,25 %** mesurés à la
main lors du diagnostic initial par une méthode entièrement indépendante
(grille de hauteurs, sans capsules) — validation croisée du module lui-même.

Les 1 073 défauts d'ordre restants sont exactement la cible d'un
réordonnancement des couches ; les 194 surplombs réels ne s'en corrigent pas.

### Câblage et tests

`PipelineParams.support_check` (opt-in, `None` = désactivé), câblé
**simultanément** dans `main.run_pipeline` ET
`scripts/export_slicing_preview.run_slicing_only` — le piège des limites
#47/#48/#49 (un orchestrateur oublié) est ici vérifié **littéralement par un
test** sur la source des deux fonctions, plutôt que laissé à la convention.

`tests/test_step9b_support_check.py` (11 tests) + 3 tests d'orchestration dans
`tests/test_main_pipeline.py`. Dont un test comparant la projection sur segment
à `step9._point_to_segment_distance`, pour qu'une seconde géométrie ne puisse
pas diverger en silence.

**247/247 tests** (233 avant, zéro régression).
