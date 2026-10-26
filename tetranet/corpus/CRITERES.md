# Critères d'acceptation (à fixer avec l'utilisateur)

TetraNet est « équivalent » à TetGen sur une pièce si **tous** les seuils tiennent. Valeurs à
décider avant d'écrire du code ; aucune n'est arrêtée.

| mesure | seuil |
|---|---|
| écart de peau maximal (maillage ↔ STL) | ≤ tolérance (min(fraction de la taille, fraction de la buse)) |
| rapport de volumes | à fixer |
| angle dièdre minimal | à fixer |
| matière déposée sans appui (step9b) | écart relatif à fixer, par rapport à TetGen |
| matière flottante | idem |
| collisions de la tête | idem |
| temps d'impression estimé | idem |
| nombre de couches, épaisseur réelle | idem |
| temps de maillage | objectif : < 10 s sur le Benchy dense |

Un échec est un résultat : statut + raison, et repli sur TetGen côté moteur.
