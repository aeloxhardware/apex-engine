# Idée réservée pour la toute fin

_Extrait de l'ancien CLAUDE.md._

## IDÉE RÉSERVÉE POUR LA TOUTE FIN — NE PAS IMPLÉMENTER MAINTENANT

**Sélection manuelle des parties, avec paramètres et ordre d'impression par
partie.** Idée de l'utilisateur, actée explicitement comme **le tout dernier
chantier**, à n'ouvrir **que lorsque l'algorithme automatique sera au point**.
Ce n'est pas une piste à explorer quand un mécanisme automatique résiste : la
tentation sera précisément d'y recourir comme échappatoire à un défaut non
résolu, et ce serait déplacer le problème sur l'utilisateur au lieu de le
résoudre. **Tant que le pipeline automatique n'est pas satisfaisant, ne pas
commencer ce travail, même partiellement.**

Contenu visé, tel qu'énoncé :
- un **algorithme de sélection** permettant à l'utilisateur de désigner
  lui-même les différentes parties du modèle (pilier, bras, zone suspendue…),
  plutôt que de subir la segmentation automatique ;
- des **paramètres propres à chaque partie** (direction, hauteur de couche,
  et tout autre réglage aujourd'hui global) ;
- un **ordre d'impression choisi par l'utilisateur** entre ces parties.

**Ce qu'il faut préserver dès maintenant pour que ce soit implémentable plus
tard, sans rien coder de la fonctionnalité elle-même** — ces points ne coûtent
rien s'ils sont respectés au fil de l'eau, et coûteraient une refonte s'ils
étaient perdus :

1. **Toute pondération spatiale doit rester définie PAR TÉTRAÈDRE, jamais
   comme un scalaire global.** Un champ `(T,)` calculé automatiquement
   aujourd'hui est exactement la même structure qu'un champ peint à la main
   demain — la sélection utilisateur devient une simple **source alternative**
   du même tableau. Un scalaire global, lui, n'a aucun chemin de migration.
   Vaut pour la modulation de lissage, les pondérations de direction, et tout
   terme ajouté au système de step4.
2. **Ne jamais supposer qu'une partition du modèle provient de step3.**
   Partout où du code consomme un découpage (régions, branches, clusters),
   il doit le recevoir en argument, jamais aller le chercher lui-même dans
   `RegionDependencyGraph`. Une partition fournie par l'utilisateur doit
   pouvoir se substituer à celle de step3 sans modifier l'appelé.
3. **L'ordre d'impression doit rester une donnée, pas une conséquence
   implicite d'un tri.** `_frontier_order` est déjà sur ce modèle : il produit
   une **permutation explicite** des couches. Conserver cette forme — un ordre
   imposé de l'extérieur remplace alors la permutation calculée, sans toucher
   à la géométrie des couches.
4. **`core/constraints.py` est le point d'extension prévu** pour tout apport
   utilisateur au champ scalaire (registre `CONSTRAINT_ASSEMBLERS`). Une
   sélection de partie assortie de paramètres est, formellement, une
   contrainte de plus — ne pas court-circuiter ce registre en câblant des
   termes utilisateur directement dans step4.

Aucune de ces quatre disciplines ne demande d'écrire du code aujourd'hui :
ce sont des **contraintes de forme** sur le code qui sera écrit de toute façon.

