# Element Aeolis (paquet `nonplanar/`) : tranchage non planaire « surface en hauteur »

Second algorithme du moteur, **indépendant du pipeline 5 axes** (`pipeline/`) : pas de
tétraédrisation, pas de champ scalaire, buse verticale. Les couches sont des surfaces
`z = S_k(x,y)` qui passent progressivement du dessous plat à la forme du dessus.

    python scripts/nonplanar_slice.py res/test_models/non_planar/plate_5.stl sortie.gcode --center 155,155

## Domaine

Un maillage dont chaque rayon vertical rencontre exactement un dessous et un dessus (plaque,
relief, coque sans surplomb). Les autres sont refusés (`NotAHeightfieldError`, contrôle par
volume intégré contre volume du maillage), jamais tranchés faux.

## Méthode (v0)

1. `heightfield.from_mesh` : dessous B et dessus T rastérisés depuis les triangles (grille de
   0,25 mm), empreinte XY par union des faces tournées vers le haut.
2. `layers.build` : couche 1 plane à B + t1 ; couches k = 2..N interpolées linéairement de
   la couche 1 à T. **Plancher et plafond d'épaisseur** (`min_layer_height` 0,1 et
   `max_layer_height` 0,25 mm par défaut) : ils viennent du profil d'imprimante (donc de
   l'interface), pas du moteur. Là où la pièce est trop basse pour le plancher, la colonne
   atteint T avant la dernière couche et n'y reçoit plus rien ; le reste est absorbé par sa
   dernière couche, qui finit exactement sur T. **Inter-couches** (`intermediate_layers`,
   activé par défaut) : N est dimensionné sur la colonne médiane et l'épaisseur nominale
   (0,2 mm) ; là où l'espacement dépasse le plafond, chaque intervalle est subdivisé et les
   couches intermédiaires ne couvrent que ces colonnes. Désactivé, N est dimensionné sur la
   colonne la plus haute.
3. `toolpath` : parois et remplissage zigzag en 2D (shapely), relevés sur S_k, segments de
   0,5 mm au plus ; couches pleines en bas et en haut.
4. `slicer` : G-code cartésien X Y Z E ; débit = longueur 3D × section locale du cordon (modèle
   « stade », épaisseur locale de la couche) ; dégagement au plus haut point de la surface.

## Parois : seulement sur le contour de la pièce

Une couche n'existe pas partout : une colonne terminée (la pièce y est trop basse pour le plancher
d'épaisseur) ou une inter-couche s'arrête sur une frontière INTERNE à la pièce. Les parois ne
suivent que le contour réel de la pièce ; aucune paroi ni retrait le long d'une frontière
interne (v0 : anneaux concentriques visibles sur le dessus, signalés par l'utilisateur le
2026-10-29 sur la plaque 6, aux creux et le long du bord). Le remplissage va jusqu'à la
frontière ; la couche suivante recouvre. L'angle du remplissage est compté depuis le haut de chaque
colonne, donc le dessus a partout la même direction.

## Collision de la vraie tête (`headcheck.py`, activée par défaut)

La tête (GLB `aetheris-4`, buse verticale) est réduite à son enveloppe inférieure (point le plus
bas par cellule polaire). Pour chaque point de dépôt, on teste si un point de la tête passe sous
la surface déjà imprimée (couche précédente ; la couche en cours n'est pas comptée, comme à
l'étape 9). Rapport : `SliceReport.head` (points testés, en collision, pénétration max, pente
d'un plan montant tolérée par la tête : 28°) et une ligne `; head_collisions` dans le G-code.
Mesuré le 2026-10-29 sur les plaques : 0 collision sur 2, 3, 5 et 1_1 (50 000 à 60 000 points
chacune), malgré des pentes de 38 à 71° au-dessus de 28° : ce sont des bosses courtes, le corps
large de la tête étant à plus de 11 mm au-dessus de la pointe. Limite : pas de correction
automatique (une collision est signalée, pas évitée) et pas d'inclinaison de la tête.

## Étapes, progression et temps de calcul (2026-10-30)

`slice_mesh(mesh, params, progress=cb)` et `slice_path(chemin, params, progress=cb)` appellent `cb(slug, fraction)`
pour les étapes de `slicer.AEOLIS_STEPS` (`python main.py steps --algorithm element-aeolis`) :
`import`, `heightfield` (modèle de surface en hauteur), `layer_surfaces` (surfaces des couches), `toolpath2d` (parois et
remplissage), `collision` (tête), `export` (écriture du G-code). Chaque étape s'ouvre à 0, finit à 1, et les boucles
émettent des fractions mesurées sur le travail fait (`core/progress.py`, R-017) : triangles rastérisés, couches
remplies, points testés contre la tête, points écrits. Aucune durée n'est estimée.

Le tranchage tourne en trois phases (chemins de toutes les couches, puis collision, puis écriture), pour que chaque
étape ait sa propre barre. Mesure (5 plaques de `res/test_models/non_planar/`, réglages par défaut) :

| plaque | avant | après |
|---|---:|---:|
| plate_1_1 | 31,0 s | 8,9 s |
| plate_2 | 14,0 s | 2,9 s |
| plate_3 | 20,0 s | 3,2 s |
| plate_5 | 16,3 s | 3,3 s |
| plate_6 | 37,2 s | 7,5 s |

D'où vient le gain : hauteur, épaisseur et dessous/dessus interpolés une fois par couche au lieu de ~200 000 petits
appels (`LayerStack.sample` / `Column`) ; rastérisation vectorisée sur les couples (triangle, nœud), **identique au bit**
(×40 à ×80) ; union de couverture de GEOS pour l'empreinte (6 s → 0,1 s, repli sur l'union robuste si les triangles du
dessus se recouvrent) ; lignes du zigzag écartées par un test d'intersection préparé avant l'overlay. Le G-code est
identique octet pour octet à celui d'avant, SAUF l'effet de l'empreinte : elle diffère d'environ 1e-7 mm² (la
version précédente la rétrécissait de 1e-7 mm), ce qui change le dernier chiffre de certains `E` et quelques points
(plate_5 : 58 208 → 58 214 points testés).

## Limites connues (mesurées le 2026-10-29)

- (Corrigé le 2026-10-29) L'épaisseur était non bornée en v0 (0,018 mm sur la plaque 1_1) ;
  elle reste désormais dans [plancher, plafond] sur les 5 plaques, avec ou sans inter-couches.
- Une couche intermédiaire commence par une marche (au bord de sa région) de l'épaisseur
  d'une sous-couche, et les inter-couches de nombres m différents se raccordent par des
  sauts de z inférieurs à une épaisseur ; l'effet à l'impression n'est pas essayé.
- Le plancher suppose plafond ≥ 2 × plancher (le reste d'une colonne est absorbé par sa
  dernière couche) ; non vérifié par une garde.
- Pente du dessus : 70,8° sur la plaque 2, sans limite imposée (seule la collision de la tête
  est contrôlée, voir plus haut).
- **Ne reproduit pas la fonction de couches d'IceSL** (`scripts/compare_nonplanar_icesl.py`,
  plaque 5, `plate_5_icesl_reference.gcode`) : l'interpolation d'IceSL n'est pas linéaire en k
  (écart de z jusqu'à 0,4 mm en milieu de pile) ; sa dernière couche dépasse le dessus d'environ
  0,18 mm. Hypothèses essayées, non retenues : linéaire en k, et max(échelle plane, échelle
  descendante du dessus). Le résultat de ce paquet est valide par construction (emboîtement,
  épaisseur bornée, fin exacte sur T), pas par imitation.
- Pas de gestion du plateau polaire ni de la cinématique 4 axes : G-code cartésien seulement.

## Sources et licences

Idées : CurviSlicer (AGPL-3.0), QuickCurve et Double QuickCurve (articles, pas de code lu).
IceSL (propriétaire) : seule sa sortie G-code a été lue, jamais son code.
