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

## État

Basepack posé, aucun algorithme écrit. Prochaine étape : le plan du prototype (donné dans la
conversation du 2026-10-26), puis le corpus et ses seuils, puis le code.
