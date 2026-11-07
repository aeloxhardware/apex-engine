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
