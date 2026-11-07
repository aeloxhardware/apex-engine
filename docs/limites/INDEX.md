# Registre des limites numérotées

Une limite = un fichier `NNN.md` (texte intégral, historique compris). Ce registre remplace la section
« Limites de conception importantes » de l'ancien `CLAUDE.md` (R-030 point 7, limite dure n° 10).
Les numéros sont ceux cités partout dans le code (`limite #NN`). Mesures chiffrées : archives datées,
à relire avec la voie unique « chaleur » (CLAUDE.md).

| N° | Titre |
|---|---|
| [1](001.md) | Le plateau d'impression n'est connu qu'à partir de l'étape 3 |
| [2](002.md) | `is_bridge` et `is_overhang` ne sont PAS mutuellement exclusifs |
| [3](003.md) | `bridge_max_gap_distance` doit être calibré relativement à l'échelle du modèle |
| [4](004.md) | La compatibilité face/direction (step3) se teste comme un surplomb, pas un |
| [5](005.md) | Le graphe de dépendances entre régions doit se baser sur une référence globale |
| [6](006.md) | Directions candidates 4 axes vs 5 axes |
| [8](008.md) | `TetMesh.tet_adjacency` (adjacence nœud-nœud) ne doit PAS servir de poids au |
| [9](009.md) | Un système Laplacien pur est singulier |
| [10](010.md) | Marching tetrahedra (étape 5) doit orienter ses triangles à partir des valeurs de |
| [11](011.md) | Le sens de parcours (winding) d'un contour extérieur et d'un trou intérieur |
| [12](012.md) | L'espacement des isovaleurs (étape 5) ne peut être qu'approximativement |
| [13](013.md) | Une paroi plus fine que `nozzle_width * wall_count` ne peut pas recevoir |
| [16](016.md) | Un trou temporaire n'est scellé que par la matière déposée APRÈS lui |
| [17](017.md) | Une `ColorMap` par face ne peut pas représenter une frontière plus fine |
| [18](018.md) | Le système doit décrire une cinématique de rotation ARBITRAIRE (N axes |
| [19](019.md) | Singularité géométrique au centre du plateau (`R=0`) : `C` est |
| [23](023.md) | Une ré-segmentation ne peut PAS être locale, contrairement à ce que |
| [24](024.md) | Une réorientation ponctuelle peut casser le lissage temporel que |
| [25](025.md) | Corollaire de la limite #21 pour la conception des FIXTURES de test |
| [27](027.md) | La boucle de distance point-capsule de `_check_printed_part` était en |
| [35](035.md) | Le mécanisme `bgmesh`/`metric` de TetGen (maillage de fond portant un |
| [41](041.md) | Repère visuel demandé pour le GLB de visualisation |
| [49](049.md) | Fix de la limite #46 (trajets à vide sans routage) — palliatif |
| [56](056.md) | PORTAGE C++ — pas 0 (socle de build) et pas 1 (index spatial) FAITS |
| [58](058.md) | Étape 8 : la tête était RETOURNÉE sur toutes les couches (corrigé le |
| [59](059.md) | Cinématique : sur la machine ACTUELLE, C est IMPOSÉ par la position du |
| [61](061.md) | Export G-code (`file_io/gcode_export.py`, écrit le 2026-09-30) |
| [62](062.md) | step9 : contrôle rapide par arbre k-d (`CollisionSimParams.fast_kdtree` |
| [63](063.md) | Étape 8 : filtre médian des axes d'orientation libres |
| [64](064.md) | Nouvelle référence step9b, mécanismes opt-in tranchés, et les |
| [65](065.md) | Étape 6 : des trous dans TOUTES les parois, depuis l'activation du |
| [66](066.md) | Étape 10 : en mode exact, aucune réorientation ne pouvait réussir |
| [67](067.md) | Étape 5b : ORDRE DES COUCHES CONTRAINT PAR LA TÊTE |
| [69](069.md) | Normale par point ACTIVÉE par défaut (2026-10-06, décision de |
| [71](071.md) | Suite autonome du 2026-10-06 (hors UI) : sept points traités |
| [72](072.md) | Points 6 a 9 du DT (2026-10-06, `config.g` confirme actif par l'utilisateur) |
| [73](073.md) | Pipeline complet de death : 124 s -> 37 s (2026-10-06), resultats identiques |
| [74](074.md) | Session autonome du 2026-10-07 (sans essai machine) : supports de l'étape 10 |
| [76](076.md) | Pose du modèle sur le plateau à l'ingestion (`place_on_platform`, 2026-10-13 |
| [79](079.md) | Champ = DISTANCE AU PLATEAU par la méthode de la chaleur (`use_heat_geodesic_field`) |
| [80](080.md) | Étape 6 : un anneau était coupé à CHAQUE oscillation numérique du champ |
| [81](081.md) | `InfillPattern.OSCILLATING_SPIRAL` (2026-10-24, opt-in, jamais activé par |
| [82](082.md) | Portage natif, session autonome du 2026-10-26 (poste Linux, g++ 13, Python 3.13) |
| [83](083.md) | Protocole d'appel du tranchage (R-003, contrats C2 et C3) et numéro de contrat (R-002) |
| [84](084.md) | Second algorithme « Element Aeolis » (nom choisi par l'utilisateur le 2026-10-29) : tranchage non planaire « surface en hauteur » (dossier ` |
| [84b](084b.md) | Element Aeolis tranche toute géométrie (2026-10-31, demande de l'utilisateur : « Aeolis doit supporter tous les |
| [85](085.md) | Famille de noms des algorithmes de tranchage (décision de l'utilisateur, 2026-10-29) : toujours préfixés |
| [86](086.md) | Limite connue : fidélité de la coque avec TetraNet, parois plus fines que la maille |
