# TetraNet — contexte du projet

Mailleur tétraédrique du moteur Apex. Lire d'abord `README.md`. Il vit pour l'instant dans
`apex-engine/tetranet/` ; il deviendra le dépôt `aeloxhardware/tetranet`.

## Règles

**R-029, R-030 et R-031 sont FERMÉES le 2026-11-07** (décision de l'utilisateur, `apex/integration/RULES.md` §11.4, écrite
et poussée). La priorité absolue qui interdisait toute autre tâche est levée ; les **limites dures ci-dessous restent en
vigueur** (cliquet, pas d'état caché, machine en données, etc.), contrôlées en CI. Écarts déclarés à la fermeture, nommés et
suivis : R-030 point 2 → jalon « moteur 100 % C++ » suivi par R-014 ; R-030 point 8 → limite connue #86 (`docs/limites/086.md`).

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

8. **Limites de qualité d'Apex (2026-11-01), applicables ici sans exception** : Apex doit être
   le sommet des slicers, donc sans défaut de conception. Fichier source ≤ 800 lignes, fonction
   ≤ 80 lignes, signature ≤ 8 paramètres (`src/tetranet.cpp` fait déjà 725 lignes : il se découpe
   avant de grossir, jamais après). Aucune constante de tolérance cachée : tout seuil est un champ
   de paramètres nommé, documenté et testé. Aucun état global. Toute correction a un test qui
   échoue sans elle, et les seuils de `corpus/CRITERES.md` sont rejoués. La qualité des
   tétraèdres « non atteinte » ci-dessous est une dette à fermer, pas un état à accepter.
   Détail : `apex/integration/RULES.md` §11 et `REQUESTS.md` (R-029).

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
- **Qualité des tétraèdres** : les découpes de Labelle-Shewchuk sont faites (commit 05c71508 : angle dièdre
  minimal 0,7° → 19° sur la sphère ; 23 à 45° mesurés sur une boîte de 10 mm à maille 1,5 à 0,7 mm). Mise à jour du
  2026-11-07 : l'ancienne mention « toujours non atteinte » était périmée.
- **Fidélité de la coque (R-030 point 8, mesurée le 2026-11-07)** : écart du contour extérieur à la surface d'origine,
  maille 1,5 mm. death : p99 0,29 mm, max 0,54, 7,1 % des points > 0,1 mm ; turkish : p99 0,23, max 0,56, 8,1 %.
  TetGen : p99 0,03 / 0,007. Deux pistes essayées, détail dans `docs/RAPPORT.md` :
  - **D (projection des contours sur la surface, `element_ardor/step5c_skin_projection.py`, opt-in)** : death p99
    0,017 mm, 0,2 % > 0,1 mm ; turkish p99 0,13, 1,9 %. Durée +18 à +23 %. Reste : parois plus fines que la maille
    (segments reliant deux faces de normales opposées), que ni D ni C ne corrigent.
  - **C (protection des arêtes vives, `Params::feature_angle_deg`, désactivée par défaut)** : bon sur une boîte à maille
    1,5 mm (p99 de peau 0,13 → 0,015 mm), mais l'angle dièdre minimal chute de 23 à 5° à maille 1,0 mm, et aucun gain
    sur death ni turkish. Non fiable, à reprendre (test `tests/test_tetranet_features.py`, échec attendu documenté).
- Reste : phase 5 (plateau, tolérance automatique), phase 7 (parallélisme déjà partiel), comparaison
  de bout en bout avec TetGen, test des cassures dans le pipeline complet.
