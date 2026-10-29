# nonplanar : tranchage non planaire « surface en hauteur »

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

## Limites connues (mesurées le 2026-10-29)

- (Corrigé le 2026-10-29) L'épaisseur était non bornée en v0 (0,018 mm sur la plaque 1_1) ;
  elle reste désormais dans [plancher, plafond] sur les 5 plaques, avec ou sans inter-couches.
- Une couche intermédiaire commence par une marche (au bord de sa région) de l'épaisseur
  d'une sous-couche, et les inter-couches de nombres m différents se raccordent par des
  sauts de z inférieurs à une épaisseur ; l'effet à l'impression n'est pas essayé.
- Le plancher suppose plafond ≥ 2 × plancher (le reste d'une colonne est absorbé par sa
  dernière couche) ; non vérifié par une garde.
- Pente du dessus : 70,8° sur la plaque 2 ; aucune limite de pente ni contrôle de collision
  de la vraie tête n'est encore appliqué.
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
