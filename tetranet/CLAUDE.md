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

- Phase 0 : référence TetGen mesurée avec `tools/baseline_tetgen.py` (death : 47 s ici, 318 183
  tétraèdres, binding Python, 4 cœurs) ; les grosses pièces tournent encore.
- Phase 1 : `tools/vtk_export.py` (VTK binaire + qualité des tétraèdres).
- Phases 2-4 (premier jet) : `src/tetranet.cpp` — BVH, signe par pseudonormales, composantes
  hors bande, réseau BCC, ramené à la peau, tri des tétraèdres, contrôles. CLI :
  `tools/tetranet_cli.cpp`.
- **Vitesse** : confirmée (death_hi-res, 226 156 faces, h = 0,6 : 3,9 s, dont 2,6 s de contrôles).
- **Qualité : NON atteinte.** Angle dièdre minimal ~0,5-0,7°, centaines de tétraèdres
  dégénérés écartés, arêtes de bord non manifold, écart de peau jusqu'à 1 mm sur les pièces
  fines. Cause identifiée : on ramène les sommets à la peau sans les découpes du schéma de
  Labelle-Shewchuk. Prochaine étape : ajouter ces découpes (points de coupe sur les arêtes
  qui traversent la peau) au lieu d'écarter les tétraèdres.
