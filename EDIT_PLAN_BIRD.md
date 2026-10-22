# EDIT_PLAN_BIRD — Accélérer la voie Bird

Écrit le 2026-10-19. Sources : `CLAUDE.md` #75 à #78, `DT_Bird_CPP.md`, mesures de la session
(`D:/apex_multibranch/bird_results.jsonl`). Les gains annoncés sont des **estimations non mesurées**,
sauf mention contraire.

## 1. Question posée

Peut-on modifier la méthode Bird pour qu'elle soit plus rapide ? Est-elle entièrement portée en C++ ?

## 2. Réponse courte

- **Non, Bird n'est pas entièrement portée en C++.**
- **Oui, on peut l'accélérer.** Règle du projet : profiler avant de porter. Le profil détaillé du dino
  à jour n'existe pas encore.

## 3. État du portage (vérifié dans `cpp/` et `pipeline/bird/`)

### Déjà en C++
| Poste | Fichier | Équivalence |
|---|---|---|
| Moindres carrés de la déformation et des rotations (parallèle) | `bird_lsq.cpp` | tolérance déclarée / au bit en parallèle |
| Viterbi de `map_polyline` | `apex_core.viterbi_choice` | au bit |
| Contrôle de collision du lissage (P5, parallèle) | `skin_check.cpp` | au bit |
| Barycentriques des candidats | `bird_locate.cpp` | tolérance déclarée |
| Export G-code, découpage `max_tip_deviation`, contrôle des dégagements | `gcode_export.cpp` | texte identique |
| Socle Vulkan | `cpp/src/gpu/` | aucun noyau métier dessus |

### Encore en Python (4 760 lignes dans `pipeline/bird/`)
| Poste | Fichier | Remarque |
|---|---|---|
| Creusement du corps (P3) | `skin/skin_carve.py` | shapely, mémoire : bloque sur turkish |
| Dernier moment accessible par face (P2) | `skin/zones.py` | faces indépendantes, bon candidat multi-cœur |
| Chemins de lissage sur surfaces (P4) | `skin/surface_paths.py` (1 209 lignes) | dépend de shapely |
| Étape 8 Bird (boucle d'objets) | `machine_toolpath.py` | à passer en tableaux |
| Dijkstra ×2 + topologie | `deformation.py` | |
| Champ initial des rotations | `deformation.py` | ~9 s (benchy), ~27 s (death) |
| Orchestration (`loop`, `fold_judge`, `chain`) | `loop.py`, etc. | à porter en dernier |

## 4. Temps de calcul de référence

| Modèle | Bird | FEM |
|---|---:|---:|
| petits modèles (l_shape, l-triple, real_overhang) | 18 à 29 s | 3 s |
| dino | **698 s** (mesuré 2026-10-19) | 16 s |
| benchy | 150 à 230 s (mesures antérieures) | 18 s |
| death | ~714 s (mesure antérieure, avant les gains parallèles) | 26 s |

À noter : ces durées ont été prises sur une machine partiellement occupée à d'autres moments, avec une
variance d'environ 10 %.

## 5. Pistes d'accélération, par ordre

1. **Profiler le dino** (`scripts/profile_bird.py`), une fois le lot Bird terminé, pour connaître la
   répartition réelle. Aucun portage avant cela.
2. **Algorithmique, sans porter** (gain possible important, risque faible) :
   - arrêter les tours de lissage dès que le nombre de points en collision ne baisse plus
     (7 à 8 tours sur deer, plusieurs à faible gain) ;
   - réduire les essais de la boucle de déformation.
   À mesurer sur la qualité (matière sans appui, collisions) avant toute activation par défaut.
3. **Porter P2 puis P3** : plusieurs jours, P3 dépend de shapely. P4 en dernier.
4. **Étape 8 Bird en tableaux**, puis les deux Dijkstra.
5. **Ne pas lisser** les modèles qui n'en ont pas besoin (option existante, à tester modèle par modèle).
6. **GPU** : seulement si le C++ multi-cœur laisse plus de 1 s par appel sur un calcul. Aucun candidat
   mesuré à ce jour.

## 6. Règles à respecter (CLAUDE.md)

- Python = oracle ; équivalence au bit quand c'est possible, sinon tolérance déclarée et mesurée.
- Mesurer avant de conclure. Ne rien activer par défaut sans mesure sur les modèles de test.
- Reprofiler après chaque changement.

## 6 bis. Comparaison Bird / FEM déjà mesurée (matière sans appui, step9b)

À compléter par le tableau final du lot Bird (`D:/apex_multibranch/bird_results.jsonl`).

| Modèle | FEM (réf.) | Bird |
|---|---:|---:|
| real_overhang | 0,58 % | 0,32 % |
| l_shape | 2,19 % | 1,01 % |
| l-triple | 2,54 % | 1,52 % |
| dino | 0,019 % | 0,010 % |
| benchy, deer, fig, death | — | en cours |

Les collisions de Bird sont comptées après les paliers de réorientation, ce que la voie FEM n'applique
pas par défaut : seules les matières sans appui sont strictement comparables.

## 7. Décision en attente

Lancer le profil du dino (point 5.1) quand le lot Bird sera terminé.
