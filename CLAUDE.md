# Apex Slicer — Contexte projet — FICHIER MODIFIABLE PAR L'IA

Moteur de slicing géométrique (Python) pour une imprimante 3D **4 axes** (3 linéaires
+ 1 rotatif). L'interface utilisateur est en **TypeScript**, séparée de ce moteur.
Ce fichier centralise tout le contexte nécessaire pour reprendre le développement
sans redemander les spécifications déjà actées. L'IA peut et doit modifier ce fichier lorsqu'elle opère.

> Note : les trajets à vide sans routage sont désormais corrigés par un
> palliatif lift-and-traverse — voir limite #49 pour le mécanisme, ses
> limites honnêtes (taux de résidu mesuré sur géométrie difficile), et
> pourquoi ce n'est pas un vrai path-planning 3D.

---

## VOIE UNIQUE : LA MÉTHODE DE LA CHALEUR (décision de l'utilisateur, 2026-10-27)

Le champ scalaire EST la distance au plateau à travers la matière, par la méthode de la chaleur
(`element_ardor/step4_field_solver.py`) ; les couches sont ses isosurfaces, ordonnées par isovaleur
(limite #79). **Le champ FEM directionnel et la voie Bird sont supprimés**, avec tout ce qui ne
servait qu'eux (champ multibranche, parties planes, contraintes utilisateur `Constraint[]`,
options de step3 et step5 liées). Leur code et leurs mesures restent dans l'historique git et
dans la branche `bird_old` ; ils ne sont plus cités que comme inspiration (S4 Slicer, Joshua
Bird, voir README et `CITATION.cff`). Contrat d'interface passé à la version 2
(`CONTRACT_CHANGES.md`). Les mesures chiffrées des limites ci-dessous qui parlent de « FEM » ou
de mécanismes retirés sont des archives, pas des comportements du code actuel. À refaire sur la
voie chaleur : points d'attache (contraintes utilisateur).

## NUMÉRO DE VERSION ET DÉPÔT COMMUN (décision de l'utilisateur, 2026-10-26)

L'application Apex (moteur + interface) a **un seul numéro de version**, rangé dans le dépôt
commun **`aeloxhardware/apex`** : fichier `VERSION` (source de vérité) et `VERSIONS.md`
(journal). **À chaque changement de version du moteur**, la session met aussi à jour ces deux
fichiers dans `apex` (récupérer d'abord son dernier état : l'interface a pu changer le numéro
entre-temps), et commence ses messages de commit par le numéro complet. **Copie aussi le numéro dans `config/APEX_VERSION`** (lu par `file_io/apex_signature.py`, écrit dans `meta.apex_version` du bloc chiffré des G-code). Format du numéro,
niveaux à monter et règle des conflits : `apex/integration/RULES.md`, §9. Le détail des
changements du moteur reste dans le DT et les limites numérotées de ce fichier.

Le moteur et l'interface se développent séparément, une session par dépôt :
`apex/integration/RULES.md` (règles et contrat entre les deux) et
`apex/integration/REQUESTS.md` (demandes ouvertes, à lire en début de session).

## RAPPELS PERMANENTS À L'UTILISATEUR (demande du 2026-10-18)

**À la fin de CHAQUE réponse, rappeler en une ligne le point ci-dessous, tant que
`DT.md` (section « RAPPELS PERMANENTS ») ne le marque pas fait `[x]`.** Ne cesser un rappel
que lorsque l'utilisateur dit l'avoir fait et que la case est cochée dans `DT.md`.
1. ~~Modifier les paramètres du firmware~~ : **fait** (confirmé par l'utilisateur, 2026-10-26, `DT.md` coché). Rappel arrêté. Salon MFR 2026 passé avec une récompense spéciale pour l'originalité du projet ; **le Benchy n'a toujours pas été imprimé**.
2. ~~Essayer le champ scalaire multibranche~~ : **fait le 2026-10-19** (limite #78, résultat mitigé,
   opt-in). Rappel arrêté.

## LIMITES DE QUALITÉ NON NÉGOCIABLES (décision de l'utilisateur, 2026-11-01)

**R-029, R-030 et R-031 sont FERMÉES le 2026-11-07** (décision de l'utilisateur, `apex/integration/RULES.md` §11.4, écrite
et poussée). La priorité absolue qui interdisait toute autre tâche est levée ; les **limites dures ci-dessous restent en
vigueur** (cliquet, pas d'état caché, machine en données, etc.), contrôlées en CI. Écarts déclarés à la fermeture, nommés et
suivis : R-030 point 2 → jalon « moteur 100 % C++ » suivi par R-014 ; R-030 point 8 → limite connue #86 (`docs/limites/086.md`).

**Apex doit être le sommet des slicers : un logiciel sans défaut de conception, pas un logiciel
« qui marche ».** Un audit externe (lecture seule, 2026-11-01) a relevé les écarts ci-dessous. Aucun
n'est acceptable à long terme ; tous doivent être corrigés, le moteur comme l'interface. Ce qui
suit l'emporte sur toute commodité de session. Détail : `apex/integration/RULES.md` §11 et `apex/integration/REQUESTS.md`.

### Limites dures (aucune exception sans décision écrite de l'utilisateur)
1. **Aucun fichier source > 800 lignes, aucune fonction > 80 lignes, aucune signature > 8 paramètres**
   (au-delà : une dataclass). Un script `scripts/check_limits.py` vérifie ces seuils et échoue sur
   tout dépassement nouveau.
   **Fait le 2026-11-04** : `scripts/limits_baseline.json` liste la dette (14 fichiers, 49 fonctions,
   15 signatures à cette date) ; `python scripts/check_limits.py --update` la resserre, `--report` la
   publie. `out/` est ignoré par git et `out/oldmods/` supprimé (R-030 point 7).
   **Fait le 2026-11-07** : plus aucune dette (0 fichier, 0 fonction, 0 signature hors limite).
   Prototypes et anciens modules : archivés dans la branche `archives` (voir `docs/ARCHIVES.md`),
   retirés de `main` le 2026-11-04 (14 `prototype_local_retriangulation*.py`, `out/oldmods/`).
2. **Cliquet** : un fichier ou une fonction déjà hors limite ne peut que rétrécir. Toucher une
   fonction hors limite oblige à la découper dans la même passe. Dettes connues à résorber, par
   ordre : `step1_ingest.py` (2 701 l.), `main.run_pipeline` (277 l.), `_tetrahedralize` (269 l.),
   `export_gcode` (257 l.), `element_aeolis/slicer._slice` (239 l.), `_self_intersection_resolve`
   (999 l.), `step10_adjustment.run` (211 l.), `step6_toolpath2d.py` (1 604 l.),
   `step8_machine_toolpath.py` (1 417 l.), `step5_layers.py` (1 269 l.).
3. **Aucune variable d'environnement ne change le résultat** (`APEX_NATIVE_FILL`, `APEX_NATIVE_LAYER`,
   `APEX_NATIVE_GRID`, `APEX_TETGEN_WORKER`, `APEX_USE_CPP`…) : tout passe par un objet d'options
   explicite. **Aucun `global` mutable** (`_LAST_NATIVE_FALLBACK_COUNT`, `_FAST_CHECKER_CACHE`,
   état de `core/progress.py`) : état porté par un objet passé en paramètre.
   **Fait le 2026-11-07** : `core/runtime.py` (`RuntimeOptions` dans `PipelineParams.runtime`, `RuntimeState` par
   run), contrôle `scripts/check_hidden_state.py`. Seule reste `APEX_CPP_LIB_DIR` (emplacement du binaire natif).
   Écart assumé : les fonctions profondes lisent le run COURANT (variable de contexte) au lieu d'un paramètre
   explicite ; `runtime.bind` propage le run aux fils.
4. **Aucune constante de machine dans le code** (`toolhead_id="aetheris-4"`, `PLATFORM_RADIUS=150.0`,
   noms de GLB en dur) : profil de données, jamais un défaut caché. Ajouter une machine ou une
   tête = ajouter un profil, sans toucher au cœur.
   **Fait le 2026-11-07** pour la cinématique, les limites d'axes, le plateau et la tête : `res/machines/*.json`,
   `config/machine_profiles.py`. Reste en code : vitesses et longueur de pivot du G-code (`GcodeParams`,
   `config/pivot_to_tip.json`), à passer en profil.
5. **Dépendances déclarées et figées** : `pyproject.toml` + fichier de verrouillage, licences
   (AGPL de TetGen et PyMeshFix comprises) listées. Aucun `pip install` sans version.
6. **Un correctif natif n'existe pas sans son test d'équivalence** Python/C++ ; la bascule C++
   a une date et une règle de fin (Python supprimé ou gelé comme oracle de test), écrites ici.
   **Règle : `docs/NATIF_ET_PYTHON.md` (jalon « moteur 100 % C++ », RULES §11.3, aucune date ; suivi par R-014).**
7. **Tout défaut qui change le G-code incrémente `APEX_CONTRACT_VERSION` dans le même commit** :
   plus de « à confirmer par l'utilisateur » laissé ouvert dans `CONTRACT_CHANGES.md`.
8. **Un palliatif porte un identifiant, un test qui le caractérise et une date de revue.** Le
   lift-and-traverse (limite #49), le repli natif→Python et tout contournement futur y sont soumis.
   **Fait le 2026-11-07** : `docs/PALLIATIFS.md`, contrôlé par `scripts/check_palliatives.py` (CI).
9. **Zéro code mort versionné** (`out/oldmods/`), zéro copie complète du projet à côté des dépôts
   (`*_local_*`), `out/` ignoré par git. Le dossier `scripts/` (80 fichiers) est trié : outil
   maintenu et testé, ou supprimé. **Fait le 2026-11-07** : 16 scripts maintenus (`scripts/README.md`, test
   `tests/test_scripts_registry.py`), 54 archivés (branche `archives`).
   **TetraNet (R-030 point 8) : limite connue #86 (2026-11-07).** Fidélité de la coque mesurée : la projection des contours
   sur la surface (`element_ardor/step5c_skin_projection.py`, opt-in) ramène death à p99 0,017 mm (0,2 % des points
   > 0,1 mm) et turkish à p99 0,13 mm (1,9 %) ; reste des parois plus fines que la maille. La protection des arêtes
   vives (piste C) n'est pas fiable. Détail : `tetranet/docs/RAPPORT.md`.
10. **Une seule source documentaire.** **Fait le 2026-11-07** (ce fichier : 4 515 → ~215 lignes ; voir le tableau
    « Où trouver le reste »). Ce fichier était découpé : règles et contrat
    ici (< 400 lignes), limites numérotées dans un registre structuré (`docs/limites/`), historique
    dans `FIXED.md`. Les imports passent par de vrais paquets (`element_ardor`, `element_aeolis`),
    sans alias `exec`.

### Définition de « terminé » pour toute session
Une tâche n'est finie que si : les tests passent **dans les deux positions du switch natif** ; la
mesure est rapportée sur les quatre modèles de référence ; `check_limits.py` ne signale aucune
régression ; aucune nouvelle dette n'est créée sans entrée au registre ; le numéro de version et
`CONTRACT_CHANGES.md` sont à jour. **Quand on touche une zone, on la laisse meilleure qu'on ne l'a
trouvée** : la dette ne se déplace pas, elle diminue à chaque passage.

## COMMENT LIRE CE FICHIER (à lire AVANT le reste)

Ce document est long, dense et écrit avec assurance. **Cette assurance n'est
pas une preuve, et sa structure n'est pas une contrainte.** Il décrit ce qui a
été fait, pas les limites de ce qui peut l'être. Lu naïvement, il induit un
biais mesuré en session réelle : optimiser *à l'intérieur* du découpage
existant au lieu de se demander si ce découpage pose la bonne question.

### Autorisation explicite de l'utilisateur

**Tout est modifiable pour peu que ça marche et que ce soit mesuré** : le
découpage en 10 étapes, les signatures `run(input, params)`, les dataclasses
de `core/mesh_types.py`, l'ordre des étapes, l'existence même d'une étape.
Ajouter une étape, en fusionner deux, faire remonter une donnée d'une étape
tardive vers une étape précoce, remplacer un critère par un autre : rien de
tout cela n'exige de demander la permission. Une fonctionnalité qui impose de
restructurer la moitié du pipeline se restructure. Les seules obligations qui
demeurent : **mesurer avant de conclure**, et **ne jamais activer par défaut
un mécanisme non mesuré sur les quatre modèles** de `res/test_models/`.

### La question à poser en premier quand un indicateur ne bouge pas

Sept mécanismes ont été écrits et mesurés pour réduire la matière déposée sans
appui. Six ont échoué. Tous les six modifiaient une **entrée** du circuit
existant (choix de direction, ancrage, ordre des couches). Aucun n'avait
interrogé le circuit lui-même. Le septième — celui qui a enfin fait bouger la
métrique — est né d'une question de l'utilisateur, pas du raisonnement mené
dans ce cadre : *« l'algorithme ne doit-il pas simplement vérifier, pour
chaque tétraèdre, s'il est en l'air ? »*

La réponse était que ce test existait déjà **deux fois**, et qu'il était
branché au mauvais endroit : `step2` évalue le surplomb contre la verticale
MACHINE, une seule fois, **avant que le champ scalaire n'existe** ; `step9b`
le mesure correctement mais *après tout*, sans rien réinjecter. Le solveur
n'était jamais pénalisé pour les surplombs qu'il créait lui-même.

**Règle de méthode qui en découle, à appliquer avant d'écrire un mécanisme de
plus :**

1. **Après deux ou trois mécanismes mesurés sans effet, arrêter d'en ajouter
   et auditer le circuit.** Un quatrième mécanisme sur une boucle ouverte ne
   la fermera pas.
2. **Pour chaque critère, demander : contre quelle RÉFÉRENCE est-il évalué, et
   est-ce celle qui compte à cet endroit du pipeline ?** Le surplomb mesuré
   contre la verticale machine et le surplomb mesuré contre `∇φ` portent le
   même nom et ne sont pas la même quantité.
3. **Suspecter tout ce qui est calculé UNE FOIS, en amont de ce qu'il
   contraint.** Une valeur figée avant l'existence de l'objet qu'elle régit
   est une boucle ouverte jusqu'à preuve du contraire.
4. **Un indicateur intermédiaire ne prouve rien.** Voir le 4ᵉ enseignement de
   la section « Prochaine étape » : direction, ancrage et ordre ont tous été
   améliorés sans que la métrique de référence (step9b) suive.

### Ce que ce fichier ne dit pas

Une limite documentée n'est une preuve de rien tant que le code n'a pas été
relu — la limite #50 a décrit pendant des sessions entières un mécanisme qui
n'existait pas. Les « pistes écartées » l'ont été dans un contexte donné ; une
piste écartée pour une raison qui ne tient plus se rouvre. Les sections
« ce qui reste ouvert » formulent des questions **déjà cadrées** : ce cadrage
est lui-même à remettre en cause en priorité, puisque c'est exactement là que
le raisonnement s'est enlisé six fois.


## OÙ TROUVER LE RESTE (le contexte détaillé a quitté ce fichier, R-030 point 7)

| Besoin | Fichier |
|---|---|
| Limites numérotées (#1 à #85, texte intégral, une par fichier) | `docs/limites/INDEX.md` |
| Spécification du produit, des 10 étapes, décisions techniques et dépendances | `docs/SPEC_PRODUIT_ET_ETAPES.md` |
| Détail module par module (`mesh_types`, `geometry_utils`, `spatial`, steps, `main`) | `docs/DETAIL_MODULES.md` |
| Prochaine étape, feuille de route, Maker Faire, temps de calcul | `docs/ROADMAP.md` |
| Idée réservée (sélection manuelle des parties), à ne pas implémenter | `docs/IDEE_RESERVEE_SELECTION_MANUELLE.md` |
| Convention de rédaction du Dev Tracker (`DT.md`) | `docs/DEV_TRACKER_CONVENTION.md` |
| Historique des correctifs | `FIXED.md` ; journaux de session : `docs/journal/` |
| Options d'exécution (natif, GPU, TetGen, fils), état de run | `core/runtime.py` |
| Profils de machine et de tête (données) | `res/machines/`, `res/3d/toolheads/` |
| Règle du commutateur natif et date de fin du Python | `docs/NATIF_ET_PYTHON.md` |
| Palliatifs (identifiant, test, date de revue) | `docs/PALLIATIFS.md` |

Remarque d'architecture : chaque algorithme est un vrai paquet Python `element_xxxx/` (`element_ardor`,
`element_aeolis`) ; les anciens alias `pipeline/` et `nonplanar/` (avec `exec`) sont supprimés le 2026-11-07. Les
identifiants publics du contrat (`element-ardor`, `element-aeolis`) ne changent pas. Détails dans `docs/SPEC_PRODUIT_ET_ETAPES.md` (arborescence).

## Conventions de code à respecter

- Chaque `element_ardor/stepN_*.py` expose `run(input, params) -> output`, fonction pure,
  sans état caché — permet à l'orchestrateur de ré-exécuter une portion du pipeline
  sans tout recalculer (utile pour l'étape 10, qui relance sélectivement 3 et 9).
- Un fichier de paramètres dataclass par étape (`OverhangParams`, `SegmentationParams`, etc.),
  actuellement définis dans chaque fichier `stepN_*.py` lui-même — à terme, `config/params.py`
  devra les centraliser sans casser les imports existants.
- Tests pytest par étape dans `tests/test_stepN_*.py`, mesh synthétiques générés à la
  volée (jamais de fichiers externes committés) via `trimesh.creation.*` et opérations
  booléennes (`union`/`difference`, moteur `manifold3d`).
- Toujours passer par un fichier temporaire + `step1_ingest.run()` dans les fixtures de
  test pour obtenir un `MeshData` réaliste (avec bvh/kdtree construits), plutôt que de
  mocker `MeshData` à la main.

