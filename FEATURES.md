# Idées de fonctionnalités futures — pistes exotiques

Ce fichier liste des extensions envisagées mais **non implémentées**, au-delà
du cœur algorithmique actuel (step1 à step10). Rien ici n'est planifié à
court terme — c'est une réserve d'idées à reprendre une par une, chacune
comme une tâche à part entière avec son propre plan, une fois le cœur du
pipeline validé visuellement et stabilisé.

Cohérence avec l'architecture existante : chaque idée doit rester compatible
avec le principe déjà acté du projet — modularité, aucune fonctionnalité
codée en dur pour une seule machine, chaque étape reste une fonction pure
`run(input, params) -> output`.

---

## Remplissage adaptatif à la structure et au sens d'impression

- [stated] Le remplissage (infill) doit se déformer selon la structure de
  l'objet et son sens d'impression, pour rester toujours parallèle aux
  couches — et doit pouvoir être ajusté pour répondre à des contraintes
  structurelles (donc devenir volontairement moins parallèle si la
  résistance mécanique locale l'exige).

## Personnalisation de l'algorithme de tranchage

- Solveur alternatif pour le champ scalaire (step4) basé sur une équation
  eikonale (`‖∇φ‖=1`, type fast marching) pour des couches à épaisseur
  vraiment constante, pas seulement approximative comme actuellement.
- Segmentation alternative (step3) par graph-cut multi-labels plutôt que
  region-growing, pour une segmentation potentiellement plus optimale sur
  des géométries très complexes (déjà évoqué comme alternative possible
  dès le premier plan de step3).
- Patterns d'infill au-delà de lignes/grille : gyroïde, nid d'abeille,
  densité adaptative selon la contrainte mécanique locale.

## Extensions du système de contraintes

- Contrainte de résistance mécanique réelle (issue d'une simulation FEM de
  charge) influençant automatiquement l'orientation des couches dans les
  zones à effort maximal — s'appuierait sur `ConstraintType.STRUCTURAL`,
  déjà présent dans `core/mesh_types.py` mais jamais encore exploité.
- Contrainte esthétique par "veine de bois" : orienter les couches pour
  suivre un motif décoratif choisi par l'utilisateur, au-delà du simple
  point d'attache.
- Import de contraintes depuis un fichier externe (résultat d'une
  simulation FEM tierce), pas seulement des points d'attache saisis
  manuellement.

## Planification avancée (au-delà de step7)

- Optimisation de trajet inter-îlots (type voyageur de commerce) : step6
  suit actuellement l'ordre d'apparition entre îlots séparés d'une même
  couche (limite déjà documentée), sans optimiser les trajets à vide.
- Simulation thermique simplifiée pour affiner `hole_cost` (actuellement
  un modèle de pénalité fixe, pas une vraie simulation de refroidissement
  entre deux passages de couleur).
- Multi-matériaux, pas seulement multi-couleurs : généraliser `ColorMap`/
  le mécanisme de repassage à différents filaments (PLA + support
  soluble, PLA + TPU), avec des contraintes de compatibilité thermique en
  plus du simple coût de purge.

## Extensions liées à la cinématique générique

- Optimisation globale de trajectoire pour minimiser les mouvements
  d'axes (pas seulement résoudre point par point comme step8 le fait
  actuellement) — utile si un axe est mécaniquement plus lent ou plus
  fragile qu'un autre, pour privilégier son immobilité.
- Bibliothèque de profils machine sauvegardés et sélectionnables dans
  l'UI, puisque `MachineProfile` est déjà générique par conception.

## Autres pistes exploratoires

- Détection automatique de la meilleure orientation de départ du modèle
  (avant step2), pour minimiser le nombre de régions nécessaires —
  actuellement l'orientation initiale du STL est prise telle quelle.
- Aperçu basse résolution en temps réel pour l'édition de points
  d'attache (évoqué dès la conception initiale comme mode dégradé
  possible, jamais formalisé en tâche concrète).
