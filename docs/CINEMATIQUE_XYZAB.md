# Nouvelle cinématique XYZAB (décision de l'utilisateur, 2026-11-27)

« Il est temps de refaire la cinématique. » La machine devient une imprimante **XYZAB** :

- **bedslinger** : le plateau, **carré, 220 × 220 mm**, se déplace en Y ; la tête se déplace en X et Z ;
- **A et B sont sur la tête** : A tourne autour de X (moteur sur le chariot X) et porte B ; B tourne autour de Y en
  position de repos (lecture de la capture CAO vue de dessus, à confirmer) ; buse vers −Z au repos ;
- plages : **A et B ± 180°** ;
- **décalage pivot-pointe : 40 mm par défaut**, valeur provisoire (rappel permanent n° 3 de `CLAUDE.md`).

Pourquoi : sur l'Aetheris-4 (R C Z B, C imposé par la position, limite #59), la tête ne penche que vers l'extérieur du
plateau ; un plafond tenu des deux côtés (cabine du Benchy) ne peut pas être refermé en éventail (prototype
`PlanarParams.lean_deg`, DT du 2026-11-27).

## Étapes

1. Profil `res/machines/apex-xyzab.json` (données) : axes de la tête, plateau carré, décalage pivot-pointe.
2. Étape 4 : inclinaisons atteignables = toutes (± 180°) ; plus de projection radiale pour cette machine.
3. Étape 8 : angles A, B par cinématique inverse de la chaîne A(X)·B(Y) ; aucun axe imposé par la position ;
   continuité des angles (deux solutions par direction).
4. Étape 9 : plateau carré mobile en Y ; pose de la tête par A et B.
5. Export G-code (Python puis C++, test d'équivalence) : `G1 X Y Z A B`, position du PIVOT = pointe + 40 mm le long
   de l'axe de la buse ; description `apex-kinematics/1` ; contrat incrémenté.
6. Interface : profil, aperçu (demande au dépôt `apex`).
7. Éventail sous les plafonds (`lean_deg`) mesuré sur les 4 modèles avec la nouvelle machine, puis défaut.

L'Aetheris-4 reste la machine par défaut jusqu'à la fin de l'étape 5.
