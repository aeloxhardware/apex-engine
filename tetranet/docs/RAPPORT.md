# TetraNet — rapport de session (2026-10-26)

## Ce qui est fait (commité sur `claude/tender-hawking-wgt88y`)
- Mailleur BCC + « isosurface stuffing » avec découpes (Labelle–Shewchuk), nombre d'enroulement
  généralisé (dipôles/BVH), réparation d'orientation des faces, diagnostics, statut
  `InconsistentSurface` (> 5 % de points ambigus).
- Intégration step1 (`mesher="tetranet"`), code de sortie 3 + message R-016 dans `main.py`.
- Repli de grille de l'IK (step8) en natif.
- Vitesse : Benchy officiel (225 k faces) ≈ 2–4 s (TetGen : 411 s, échec en cascade).
- Tolérance mesurée sur le banc de cassures volontaires (`tools/eval_broken.py`).

## Comparaison TetGen / TetraNet, pipeline complet (résolution 1,5 sauf mention)
| modèle | mailleur | s | tétras | couches | φ_max | inatteignables % | collisions | sans appui % | défauts d'ordre |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| l_shape | TetGen | 0,6 | 1 638 | 177 | 35,38 | 25,2 | 26 | 3,87 | 0 |
| l_shape | TetraNet | 1,0 | 8 576 | 179 | 35,89 | 9,9 | 20 | 6,37 | 0 |
| l-triple | TetGen | 0,9 | 1 777 | 179 | 35,69 | 26,3 | 46 | 5,62 | 0 |
| l-triple | TetraNet | 2,1 | 12 160 | 181 | 36,16 | 10,8 | 37 | 6,10 | 0 |
| death | TetGen | 36,6 | 322 464 | 291 | 54,5 | 33,4 | 5 603 | 1,64 | 4 014 |
| death | TetraNet | 33,4 | 215 073 | 343 | 63,83 | 43,0 | 14 098 | 2,91 | 10 133 |
| death (1,0) | TetraNet | 54,0 | 688 732 | 339 | 63,55 | 41,1 | 13 716 | 2,76 | 14 171 |
| turkish | TetGen | 72,8 | 327 112 | 379 | 75,58 | 15,4 | 2 709 | 0,19 | 37 |
| turkish | TetraNet | 72,5 | 672 545 | 380 | 76,07 | 9,9 | 2 219 | 0,20 | 295 |

## Lecture honnête
- **Critère « pas de résultat significativement différent » : NON atteint sur death.stl**
  (collisions ×2,5 ; défauts d'ordre ×2,5 ; φ_max 63,8 contre 54,5, écart dont la référence géodésique
  ci-dessous suggère qu'il vient de TetGen), et partiellement sur turkish
  (défauts d'ordre 295 contre 37). Les petits modèles : écart de la matière sans appui
  (3,9 → 6,4 %) probablement lié au nombre de tétraèdres (1,6 k contre 8,6 k), pas à un défaut.
- Le raffinement en résolution (1,5 → 1,0) ne rapproche pas death de TetGen : l'écart n'est pas
  un problème de densité.
- Quel φ est le plus proche de la vérité sur death ? Référence géodésique indépendante
  (`/tmp/tn/geodesic_check.py` : Dijkstra sur grille de voxels intérieurs, 26 voisins, depuis le plateau).
  **death : 66,4 (h = 1,5) et 66,5 (h = 1,2)**, stable ; la grille surestime un peu la vraie géodésique
  (de l'ordre de quelques %), donc la vérité est probablement vers 62–66. TetraNet (φ_max 63,8) est
  **plus proche de cette référence que TetGen (54,5)** : l'écart de φ_max ne prouve PAS que TetraNet se
  trompe, TetGen sous-estime plutôt. À h = 0,8 le calcul est tué faute de mémoire. **l_shape** : la
  référence varie de 33,4 à 35,1 selon h (bruit de voxelisation de ±1,5), du même ordre que l'écart
  TetGen 35,38 / TetraNet 35,89 : rien à conclure. Réserve : φ du pipeline vient de la méthode de la
  chaleur, pas d'un Dijkstra ; et ce résultat n'explique PAS les collisions et défauts d'ordre plus
  nombreux de TetraNet sur death (14 098 contre 5 603 ; 10 133 contre 4 014), qui restent à comprendre.
- Écart de peau aux arêtes vives : boîte 0,5 mm à h = 1,5 ; Benchy p99 0,28 mm — pas 0,05 mm.
- Non portés : raffinement près du plateau, tolérance automatique (phase 5) ; zones de jonction
  retombent sur TetGen. Lecture STL/GLB toujours via trimesh.

## Décisions à prendre par l'utilisateur
1. Seuils de `corpus/CRITERES.md` (non fixés ; je ne les ai pas inventés).
2. Accepter ou non l'écart sur death avant d'activer TetraNet par défaut (aujourd'hui : opt-in).
3. Priorité : qualité de peau aux arêtes vives, ou raffinement plateau.

## Non fait / non touché
Pas de fusion dans `main`, pas de PR, aucun défaut ni seuil modifié, copies de `apex` et UI intactes.

## Fidélité de la coque : pistes C et D (2026-11-07, R-030 point 8)

Mesure : écart de chaque point de contour (sommets + 3 points par segment) à la surface STL d'origine, maille 1,5 mm,
`scripts/export_slicing_preview.py` (qui n'envoyait pas `mesher` : corrigé le même jour).

| modèle | variante | p99 (mm) | max (mm) | points > 0,1 mm | durée step1→8 | points de contour |
|---|---|---:|---:|---:|---:|---:|
| death | TetraNet seul | 0,294 | 0,541 | 7,07 % | 4,3 s | 507 461 |
| death | + C | 0,318 | 0,637 | 7,72 % | 4,7 s | 507 309 |
| death | + D | **0,017** | 0,466 | **0,21 %** | 5,3 s | 567 161 |
| death | + C + D | 0,015 | 0,442 | 0,13 % | 5,4 s | 552 501 |
| turkish | TetraNet seul | 0,232 | 0,557 | 8,09 % | 15,3 s | 980 369 |
| turkish | + C | 0,244 | 0,598 | 7,86 % | 15,6 s | 980 588 |
| turkish | + D | 0,136 | 0,555 | 1,88 % | 18,1 s | 1 134 581 |
| turkish | + C + D | 0,138 | 0,563 | 1,93 % | 20,0 s | 1 126 068 |

(TetGen : death p99 0,032 / max 0,069 ; turkish p99 0,007 / max 0,315.)

- **D** est la seule piste efficace sur les modèles réels. Plus de tours de raffinement (6, 9) ne changent presque rien
  sur turkish (1,7 %, 1,5 %).
- **Reliquat de D** : sur turkish, les segments restants (2,7 % des segments) relient deux points de faces aux normales
  OPPOSÉES (angle médian 180°, aucun sommet commun) et leur milieu est à 0,25 mm de la surface : des parois plus fines
  que la maille que le contour traverse. Ni D ni C ne peuvent les corriger ; il faut une maille locale plus fine ou une
  détection des parois minces.
- **C** : sur une boîte, p99 de peau 0,13 → 0,015 mm et rapport de volume 0,995 → 0,9996 à maille 1,5 mm ; à maille 1,0 mm
  l'angle dièdre minimal passe de 23° à 5°. Aucun gain sur death ni turkish. Désactivée par défaut.

## Écart dans les deux sens, raffinement de la peau, parois minces (2026-11-17)

Mesures du 2026-11-17, résolution 1,5 (maille 1,0 mm), noyau natif, 12 fils. Modèles : `death`, `turkish`, Benchy
officiel (225 143 faces). Les petits modèles de référence (`l_shape`, `l-triple`) ont un écart nul dans tous les cas.

### Ce que TetraNet ne mesurait pas
Le statut ne regardait que l'écart du BORD DU MAILLAGE à la surface. Mesuré dans l'autre sens (surface → bord), sur la
seule peau extérieure du solide (faces internes de coques qui se chevauchent écartées par le nombre d'enroulement,
validé contre l'enroulement exact : 3 465 contre ~3 100 points internes sur death, 44 836 contre ~45 700 sur turkish) :

| modèle | maillage → surface p99 / max | surface → bord p99 / max |
|---|---:|---:|
| death | 0,381 / 0,499 mm | 0,379 / 1,22 mm |
| turkish | 0,336 / 0,580 mm | 0,94 / 3,06 mm |
| Benchy | 0,276 / 0,428 mm | 0,44 / 1,34 mm |

Une plaque de 0,3 mm absente du maillage (volume 93 %) sortait « Ok » : écart maillage → surface 1,6e-14 mm. Sur le
Benchy, parmi les faces de bord à plus de 0,05 mm, 88 % coupent une arête vive, 11 % suivent mal la courbure, 1,5 %
enjambent une paroi mince ; au-delà de 0,2 mm, plus de la moitié des faces de death (1 792 / 3 167) et 46 % de celles de
turkish ENJAMBENT UNE PAROI d'épaisseur proche de la maille : deux nappes de surface dans un même tétraèdre du réseau,
que les gabarits de Labelle-Shewchuk ne représentent pas (la matière manque au-delà de la face).
**Gardé** : `check_cover` (défaut : actif), contrôle `cover_error` de `src/checks.cpp`, test
`tests/test_tetranet_cover.py`. Coût après réduction (test d'enroulement seulement au-delà de 0,25 h, un par triangle) :
+0,3 s sur le Benchy, +0,2 à 0,5 s sur turkish, +0,05 s sur death.

### Raffinement de la peau après remplissage : essayé, NON retenu
Bissection des arêtes du bord dont la face s'écarte de la surface (dans un sens ou dans l'autre), tous les tétraèdres
autour de l'arête coupés (maillage conforme), nouveau sommet sur l'arête vive (intersection des plans des deux faces les
plus proches) ou au point de surface le plus proche ; position refusée si un tétraèdre se retourne ou passe sous 5°.

| modèle | tours | maillage → surface p99 | > 0,1 mm | surface → bord p99 | durée du maillage |
|---|---:|---:|---:|---:|---:|
| death | 0 / 2 / 8 | 0,381 / 0,381 / 0,381 | 15,4 / 12,8 / 10,2 % | 0,379 / 0,298 / 0,283 | 0,7 / 1,6 / 5,4 s |
| turkish | 0 / 2 / 8 | 0,336 / 0,340 / 0,332 | 11,7 / 13,6 / 15,0 % | 0,94 / 0,94 / 0,94 | 4,2 / 7,7 / 20,7 s |
| Benchy | 0 / 2 / 8 | 0,276 / 0,246 / 0,223 | 12,6 / 9,7 / 7,2 % | 0,44 / 0,38 / 0,37 | 3,0 / 4,9 / 10,4 s |

Pourquoi : les tétraèdres du bord sortis du remplissage sont plats ; déplacer un sommet de 0,3 mm les retourne (31 576
coupes refusées sur le Benchy en 8 tours) ; les faces qui enjambent une paroi ne se corrigent pas en coupant (6 703
arêtes de paroi sur death). Une première version qui posait le sommet à mi-chemin quand la projection était refusée
EMPIRAIT le maillage (sommets de bord hors de la surface : 4 827 des 5 804 faces fautives de death après 4 tours).
Code retiré ; conception ci-dessus.

### Capture des parois minces dans le réseau : essayé, NON retenu
Un nœud de la bande dont le point de peau le plus proche est sur une paroi plus mince que la maille (rayon tiré à
travers, face opposée à moins de h) et à moins de 0,25 / 0,355 / 0,45 h du plan médian y est déplacé et compté dedans.
Nœuds capturés : death 11 / 19 / 31, Benchy 0 / 11 / 25, turkish 986 / 1 535 / 2 040 ; surface → bord p99 de turkish
0,94 → 0,93 / 0,93 / 0,84 mm, rien sur death ni le Benchy. Les parois fautives font 0,7 à 1 maille : leur plan médian est
à plus de 0,35 h de tout nœud extérieur. Code retiré.

### Maille plus fine partout (pour situer le coût)
| modèle | résolution | tétraèdres | durée | maillage → surface p99 | surface → bord p99 / max |
|---|---:|---:|---:|---:|---:|
| death | 1,5 / 1,0 / 0,75 | 0,22 / 0,69 / 1,61 M | 0,7 / 2,0 / 4,8 s | 0,381 / 0,182 / 0,127 | 0,379 / 0,236 / 0,168 ; max 1,22 / 0,82 / 0,76 |
| turkish | 1,5 / 1,0 / 0,75 | 0,67 / 2,19 / 5,05 M | 4,2 / 12,6 / 27,4 s | 0,336 / 0,250 / 0,147 | 0,94 / 0,90 / 0,74 ; max 3,06 / 2,96 / 2,71 |
| Benchy | 1,5 / 1,0 / 0,75 | 0,20 / 0,69 / 1,57 M | 3,0 / 5,2 / 8,1 s | 0,276 / 0,161 / 0,125 | 0,44 / 0,19 / 0,14 ; max 1,34 / 1,20 / 0,34 |

Turkish garde ~6 % de sa peau extérieure à plus de 0,3 mm même à maille 0,5 mm : détails de 0,1 à 0,6 mm d'épaisseur
(feuilles ouvertes : 62 points sur 9 265, négligeables), en partie plus fins qu'un cordon de buse.

### Projection des contours (étape 5c, piste D) sur le pipeline complet
Écart des points de contour (sommets + 3 par segment) au STL, étapes 1 à 8 :

| modèle | mailleur | p99 sans → avec D | > 0,1 mm sans → avec D | durée sans → avec D |
|---|---|---:|---:|---:|
| death | TetraNet | 0,294 → 0,017 mm | 7,07 → 0,21 % | 14,0 → 15,0 s |
| turkish | TetraNet | 0,232 → 0,136 mm | 8,09 → 1,88 % | 35,3 → 35,0 s |
| Benchy | TetraNet | 0,187 → 0,014 mm | 4,65 → 0,05 % | 17,3 → 19,2 s |
| death | TetGen | 0,032 → 0,001 mm (max 0,069 → 0,021) | 0 → 0 % | 17,9 → 16,4 s |
| turkish | TetGen | 0,006 → 0,000 mm (max 0,300 → 0,085) | 0,003 → 0 % | 28,3 → 32,9 s |

Matière sans appui, orientations inatteignables, couches : inchangées à quelques pour mille près. D corrige les arêtes
vives au niveau des contours (en 1D, peu coûteux) ; elle ne recrée pas la matière absente du maillage.

### Option `Y` de TetGen (sans point ajouté sur la surface) : NON retenue
Elle supprime en pratique le raffinement : l_shape 1 677 → 40 tétraèdres, death 333 950 → 5 425, turkish 320 372 →
15 790 ; orientations inatteignables 23 → 63 % des points (l_shape), 39 → 58 % (death) ; matière sans appui de death
1,5 → 3,0 %. Le gain de temps (death 17,9 → 5,2 s) vient d'un maillage dégradé.
