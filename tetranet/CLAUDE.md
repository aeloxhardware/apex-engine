# TetraNet — contexte du projet

Mailleur tétraédrique du moteur Apex. Lire d'abord `README.md`. Il vit pour l'instant dans
`apex-engine/tetranet/` ; il deviendra le dépôt `aeloxhardware/tetranet`.

## Règles

1. **Mesurer avant de conclure.** Un résultat est « bon » s'il passe les seuils de
   `corpus/CRITERES.md`, pas parce qu'il semble correct. Ne jamais l'activer par défaut dans
   Apex sans avoir mesuré les quatre modèles de référence du moteur (`res/test_models/`).
2. **TetGen reste la référence et le repli.** Pas de dépendance à TetGen dans ce code : il est
   appelé par le moteur, dans son exécutable séparé.
3. **Une session par dépôt.** Quand ce dossier devient un dépôt, les règles de
   `apex/integration/RULES.md` s'appliquent. Toute demande à un autre dépôt passe par
   `apex/integration/REQUESTS.md`.
4. **Le numéro de version est dans `aeloxhardware/apex`** (`VERSION`, `VERSIONS.md`) ; ne pas
   en créer un second ici.
5. **Pas de `-ffast-math`.** L'équivalence avec le moteur repose sur des flottants prévisibles.
6. **Le contrat de sortie est celui de `TetMesh`.** Le changer est un changement de contrat
   (C1/C2 d'Apex), à demander, pas à faire seul.
7. **Dire les échecs.** Un cas hors tolérance est un résultat explicite (statut + raison),
   jamais un maillage dégradé en silence.

## État (2026-10-26)

- Phase 0-1 : référence TetGen (`tools/baseline_tetgen.py`), export VTK (`tools/vtk_export.py`).
- Phases 2-4 : réseau BCC ramené à la peau (`src/tetranet.cpp`), CLI (`tools/tetranet_cli.cpp`).
- **Dedans/dehors par nombre d'enroulement généralisé** (dipôles par paquets de triangles) et
  **remise en sens des faces retournées** (propagation sur les arêtes manifold, volume positif par
  composante). Tests : `tests/test_defects.cpp` (cube intact, trou, face manquante, face retournée,
  cube à l'envers).
- Diagnostic chiffré de la surface (trous, faces retournées, arêtes partagées, doublons,
  points incohérents/ambigus) et statut `InconsistentSurface`. Texte d'erreur : `ROADMAP.md` d'`apex`
  (R-016), code modèle `tools/diagnosis_message.py`.
- Banc de cassures : `tools/gen_broken.py`, `tools/eval_broken.py`, `tools/sweep_holes.py`.
- **Mesures (Benchy officiel, maille 1 mm)** : un trou jusqu'à 4 mm de rayon (8 000 faces) ne fait
  rien perdre (avant : fuite à partir de 3 mm, un tiers du volume perdu) ; faces retournées :
  identique à la pièce saine ; fentes et arêtes à trois faces : 1 à 2 % ; auto-intersections : 7 %
  de matière manquante (la géométrie elle-même est repliée, rien à deviner) ; fichier cumulant
  tout : 5 % manquant, 19 % en trop (refusé, `InconsistentSurface`).
- **Vitesse** : ~3,5 s pour 225 000 triangles (1,9 s sans nombre d'enroulement). Pas encore
  optimisé : la surface (cartes d'arêtes) et la classification des sommets hors bande.
- **Qualité des tétraèdres** : toujours non atteinte (angle dièdre minimal 0,7° au bord, tétraèdres
  plats aux arêtes vives). Prochaine étape : découpes de Labelle-Shewchuk.
- Reste : phase 5 (plateau, tolérance automatique), phase 7 (parallélisme déjà partiel), comparaison
  de bout en bout avec TetGen, test des cassures dans le pipeline complet.
