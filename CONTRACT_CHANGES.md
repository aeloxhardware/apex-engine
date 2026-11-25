# Changements du contrat moteur <-> interface

Numero courant : voir `APEX_CONTRACT_VERSION` dans `config/schema.py` (copie dans
`config/params.schema.json`, cle `x-apex-contract`). Regles : `apex/integration/RULES.md`, §4.1 et §9.
Un changement INCOMPATIBLE incremente le numero et s'ajoute ici, avec ce que le consommateur doit faire.

## Ajout compatible (2026-11-25, alpha-0.2.3.1.9.1.0.221) : le numero du contrat NE CHANGE PAS (reste 9)
- `skin.conform_clearance` (nouveau, 1,5 mm) : le lissage non planaire (facultatif, `skin`) evite les faces proches
  d'une grande surface deja conforme (`field_solver.planar.conform`). Sans effet si `skin` est `null` (defaut) : le G-code
  par defaut ne change pas (limite #92).
- Ce que l'interface peut faire : `pnpm schema:sync` (nouveau champ facultatif).

## Contrat 9 (2026-11-25, alpha-0.2.3.1.9.0.1.221) : INCOMPATIBLE, le numero passe de 8 a 9
- Etape 5e (`planar_sections.surface_levels`) : le redecoupage a 0,3 mm n'est plus fait que sur les boucles courbes dans
  le repere de leur couche (`curved_tol`, nouveau, 0,05 mm) ; une isovaleur qui tombe pile sur des sommets du STL donne
  une ligne fermee (avant : ligne trouee, boucle gardee de l'isosurface). Benchy : G-code 376 000 -> 358 000 lignes,
  longueur deposee sans appui 30,2 -> 22,5 mm ; quatre modeles inchanges (limite #92).
- Ce que l'interface doit faire : `APEX_CONTRACT_VERSION = 9`, `pnpm schema:sync`. Ancien G-code :
  `"planar_sections": {"surface_levels": {"curved_tol": 0}}` (a la ligne de niveau trouee pres).

## Contrat 8 (2026-11-25, alpha-0.2.3.1.9.0.0.221) : INCOMPATIBLE, le numero passe de 7 a 8
Element Ardor, vers la qualite de surface d'Aeolis (limite #92) :
- `field_solver.planar.conform` (`ConformParams`, nouveau, ACTIF) : couches conformes sous les grandes surfaces en pente
  (2 a 45 deg, largeur 2A/P >= 5 mm) ; le toit et le pont du Benchy deviennent une surface continue.
- `planar_sections.surface_levels` (`SurfaceLevelParams`, nouveau, ACTIF ; etape 5e) : contour des couches courbes pris
  sur le STL (ligne de niveau du champ), redecoupe a 0,3 mm.
- `toolpath2d.run_split_min_slope_deg` = 20 (nouveau ; 0 = historique) : seuls les anneaux de plus de 20 deg sont
  decoupes en troncons montants (death : retractions -29 %).
- Ce qui change pour un meme modele : longueur deposee sans appui egale ou plus faible sur les quatre modeles ; G-code du
  Benchy plus long (~ +88 % de lignes) ; tranchage plus long (Benchy 30 -> ~57 s).
- Ce que l'interface doit faire : `APEX_CONTRACT_VERSION = 8`, `pnpm schema:sync`. Ancien G-code :
  `"field_solver": {"planar": {"conform": null}}, "planar_sections": {"surface_levels": null},
  "toolpath2d": {"run_split_min_slope_deg": 0}`.

## Contrat 7 (2026-11-25, alpha-0.2.3.1.8.0.7.221) : INCOMPATIBLE, le numero passe de 6 a 7
- Etape 5d (`planar_sections`) : une boucle plane prend la section du STL la plus proche, acceptee sous `match_tol`
  (0,5 mm) comme avant, ou sous `loose_match_tol` (1,0 mm, nouveau) si la suivante est au moins `ambiguity_ratio`
  (2,5, nouveau) fois plus loin. Benchy : dessus de la cheminee en anneaux au lieu d'un labyrinthe, 27 boucles de plus
  recoupees exactement (limite #91).
- Ce que l'interface doit faire : `APEX_CONTRACT_VERSION = 7`, `pnpm schema:sync`. Ancien G-code :
  `"planar_sections": {"loose_match_tol": 0.5}` (a tres peu pres : la plus proche est prise au lieu de la premiere).

## Contrat 6 (2026-11-24) : INCOMPATIBLE, le numero passe de 5 a 6
- `toolpath2d.min_fill_length` = 0,8 mm par defaut dans le pipeline (`Step6Params`, 0 = aucun filtre) : les cordons de
  remplissage (clairseme, plein, pont) plus courts sont retires (pates sur les dessus etroits en pente ; limite #90).
- Lissage (facultatif, ne change rien par defaut) : `skin` creuse sans recalculer toute la couche
  (`CarveParams.keep_projected`, interne, hors schema).
- Ce que l'interface doit faire : `APEX_CONTRACT_VERSION = 6`, `pnpm schema:sync`. Ancien G-code :
  `"toolpath2d": {"min_fill_length": 0}`.

## Contrat 5 (2026-11-19, publie avec alpha-0.2.3.1.8.0.3.221) : INCOMPATIBLE, le numero passe de 4 a 5 dans le commit qui publie ces defauts
Element Ardor imprime comme un trancheur classique partout ou l'impression multi-axes n'est pas necessaire (limite #88) :
- `field_solver.planar` (`PlanarParams`) ACTIF par defaut : couches planes loin des surplombs non imprimables a plat,
  direction lissee (20 iterations) dans la zone inclinee ;
- `toolpath2d.solid_top_layers` = 4 et `solid_bottom_layers` = 3 par defaut : couches pleines sous les dessus et sur les
  dessous (avant : aucune) ;
- `skin_projection.in_layer` = vrai : les contours sont ramenes sur le STL dans le plan de leur couche.
- `planar_sections` (`PlanarSectionParams`, etape 5d, ajoute le 2026-11-19, limite #89) ACTIF par defaut : chaque
  boucle plane loin de toute matiere inclinee est recoupee exactement dans le STL (simplifiee a 0,01 mm) ;
- zone pleine etroite en anneau (2A/P < 4 cordons) : anneaux concentriques au lieu du zigzag ;
- `skin` (lissage non planaire) reste `null` par defaut ; nouveaux champs de `SkinParams` : `stair_min_deg`,
  `stair_max_deg`, `proxy_loop`, `latest_band`.
- **Ce qui change pour un meme modele** : Benchy droit entierement plan, premiere couche plane et pleine ; matiere sans
  appui death 3,23 -> 1,20 %, turkish 0,36 -> 0,08 %, l_shape 6,39 -> 6,58 %, l-triple 6,25 -> 6,41 % ; filament +5 % (Benchy) a +34 % (turkish) ; 0 collision restante
  sur les quatre modeles. Detail : `docs/limites/088.md`.
- **Ce que le consommateur doit faire** : passer `APEX_CONTRACT_VERSION` a 5, recopier `config/params.schema.json`
  (nouveaux champs `planar`, `solid_top_layers`, `solid_bottom_layers`, `in_layer`, `planar_sections`, champs de
  `skin`), envoyer `"contract": 5`. Pour retrouver l'ancien G-code : `"field_solver": {"planar": null}`,
  `"toolpath2d": {"solid_top_layers": 0, "solid_bottom_layers": 0}`, `"skin_projection": {"in_layer": false}`,
  `"planar_sections": null` (les anneaux concentriques n'ont pas d'interrupteur). Lissage : `"skin": {}`.

## Contrat 4 (2026-11-17) : INCOMPATIBLE, le numero passe de 3 a 4
`PipelineParams.skin_projection` (etape 5c : contours ramenes sur la surface d'origine) est ACTIF par defaut
(`SkinProjectionParams()` au lieu de `null`). C'etait annonce dans l'ajout du 2026-11-07 ci-dessous : un defaut qui
change le G-code incremente le contrat.
- **Ce qui change pour un meme modele** (mesures du 2026-11-17, `tetranet/docs/RAPPORT.md`) : les contours de couche
  sont ramenes sur le STL. TetraNet : ecart p99 death 0,294 -> 0,017 mm, Benchy 0,187 -> 0,014, turkish 0,232 -> 0,136.
  TetGen : death 0,032 -> 0,001 mm. Points de trajectoire +7 a +9 %, duree 0 a +11 %, matiere sans appui et
  orientations inatteignables inchangees a quelques pour mille pres.
- **Ce que le consommateur doit faire** : passer `APEX_CONTRACT_VERSION` a 4, recopier `config/params.schema.json`
  (defaut de `skin_projection`), envoyer `"contract": 4` dans la requete (une requete au contrat 3 est refusee, code 2).
  Pour retrouver l'ancien G-code : `"skin_projection": null`. Demande correspondante dans `apex/integration/REQUESTS.md`.

## Ajout du 2026-11-14 (le numero reste 3, A CONFIRMER par l'utilisateur) : relief local du moteur general d'Aeolis
Le moteur general d'Element Aeolis (`report.mode = "general"`) ne deforme plus ses couches par un relief unique cale sur
le sommet de la piece : chaque surface tournee vers le haut (pente <= 45 deg et <= `nozzle_angle_deg`) est epousee par ses
propres couches (`element_aeolis/relief.py`, limite #87). Revu le 2026-11-15 : une PEAU de 3 couches epouse chaque
surface, les couches du dessous s'ecrasent (non imprimees) au lieu de se tordre, raccords jusqu'a 70 deg ; l'epaisseur
deposee reste dans [plancher, plafond] du profil (`report.thickness_*`), `report.relief_mm` vaut desormais la longueur
de buse (plus le quart de la hauteur). **Le G-code du mode general change pour un meme modele et les
memes reglages** ; aucun champ du schema, aucun defaut, aucun libelle ne change, et le moteur « surface en hauteur » est
intact. **Ce que le consommateur doit faire** : rien. Question ouverte, a trancher par l'utilisateur (regle dure n° 7) : ce
changement de G-code vient d'un correctif d'algorithme, pas d'un defaut ; s'il doit compter comme incompatible, le numero
passe a 4 dans le commit qui le publie.

## Ajout du 2026-11-07 (le numero reste 3) : `params.runtime`, plus de variable d'environnement
`PipelineParams.runtime` (`RuntimeOptions` : `native`, `gpu`, `grid_mode`, `tetgen_worker`, `tetgen_log`,
`speculative_tetgen`, `step6_threads`) remplace `APEX_USE_CPP`, `APEX_USE_GPU`, `APEX_CPP_GRID_MODE`,
`APEX_TETGEN_WORKER` et `APEX_TETGEN_LOG` (R-030 point 3). Ajout COMPATIBLE du schema : valeurs par defaut = ancien
comportement sans variable (Python, pas de GPU). **Ce que le consommateur doit faire** : un lanceur qui positionnait
`APEX_USE_CPP=1` dans l'environnement du processus doit envoyer `"params": {"runtime": {"native": true}}` ; la variable
n'est plus lue. Le G-code n'est pas modifie. Le profil de machine par defaut vient desormais de `res/machines/`
(memes valeurs qu'avant : aucun changement de G-code mesure par les tests).

## Ajout du 2026-11-07 (le numero reste 3) : `skin_projection` et `tetranet_feature_angle_deg`
Deux champs OPTIONNELS de `PipelineParams`, desactives par defaut (`skin_projection: null`, `tetranet_feature_angle_deg: 0`) :
sans eux le G-code est inchange. `skin_projection` (etape 5c) ramene les contours sur la surface d'origine ;
`tetranet_feature_angle_deg` protege les aretes vives dans TetraNet. Ajout COMPATIBLE du schema ; le jour ou un de ces
reglages devient actif par defaut, le G-code change et le numero de contrat doit etre incremente dans le meme commit.

## Contrat 3 (2026-11-04) : INCOMPATIBLE, le numero passe de 2 a 3
Decision de l'utilisateur (2026-11-04) : les defauts de buse tres bas d'Element Aeolis sont GARDES, pour eliminer des
defauts (conduits de ventilation serres, collisions). Un defaut qui change le G-code produit incremente le contrat
(`apex/integration/RULES.md` §4.1 et §11.1) : le numero passe donc de 2 a 3.
- **Ce qui change pour un meme modele** : avec `nozzle_angle_deg = 10`, `nozzle_reach = 2 mm` et `enforce_nozzle_limits = vrai`,
  les surfaces en hauteur dont les couches depassent ces limites (les plaques de `res/test_models/non_planar/`) passent au
  moteur general d'Element Aeolis ; `enforce_nozzle_limits = false` (ou 80° et 50 mm) rend l'ancien G-code.
- **Libelle renomme** : « Angle de la buse » devient **« Angle maximum buse-paroi »** (le nom du champ `nozzle_angle_deg` et son
  sens ne changent pas : pente maximale, depuis l'horizontale, de la surface deja imprimee autour de la pointe).
- **Ce que le consommateur doit faire** : l'interface passe `APEX_CONTRACT_VERSION` a 3 (`src/utils/engineContract.ts`), relance
  `schema:sync` et affiche le nouveau libelle. Une requete portant `contract: 2` est refusee (code 2).
- Rien d'autre ne change dans le schema (les deux reglages et `mode` existent deja).

## Ajout du 2026-11-01 (le numero reste 2) : types de ligne du G-code d'Element Aeolis
- Chaque polyligne extrudee d'Aeolis est precedee d'une etiquette `;TYPE:<type>` (convention de Cura, ecrite a chaque
  changement de type, remise a zero a chaque couche) : `WALL-OUTER` (paroi exterieure), `WALL-INNER` (parois interieures),
  `FILL` (remplissage creux), `SKIN` (peau pleine : dessus et dessous). Les deux moteurs d'Aeolis (surface en hauteur et
  general) l'ecrivent. Une visionneuse qui ignore les commentaires n'a rien a changer ; l'interface colore l'apercu d'apres
  ces etiquettes (R-028). Le G-code polaire (Element Ardor) n'en porte pas encore.

## Ajouts du 2026-10-31 (le numero reste 2 ; CONFIRMES le 2026-11-04 par l'utilisateur, le defaut qui change le G-code fait passer le contrat a 3, voir plus haut)
- Element Aeolis tranche TOUTE geometrie (`element_aeolis/general.py`). `AeolisSliceParams` gagne `mode` (`auto` par defaut,
  `heightfield`, `general`), `nozzle_angle_deg` (10) et `nozzle_reach` (2 mm), les deux reglages de buse de l'imprimante, et
  `enforce_nozzle_limits` (vrai). Le `report` gagne `mode`, `fallback_reason` et `relief_mm`.
- **Defaut qui change le G-code (confirme le 2026-11-04 : contrat 3)** : avec les defauts de buse tres bas, les surfaces en hauteur dont les couches
  depassent ces limites (les plaques de `res/test_models/non_planar/`) passent au moteur general ; `enforce_nozzle_limits=false`
  rend l'ancien resultat. Sans cela, le statut `not_a_heightfield` n'est plus renvoye par defaut : seulement avec
  `mode: heightfield`. Nouveau refus type (code 3) : `not_sliceable` (maillage vide ou ouvert, aucune section).
- Le controle de la buse utilise l'enveloppe des deux reglages quand `head_model` est vide (avant : le GLB `aetheris-4`).
- Interface : exposer les deux reglages de buse dans le profil d'imprimante (R-027).

## Ajouts du 2026-10-30 (le numero reste 2)
- Etapes de progression par algorithme (C3) : `main.py steps --algorithm <element-ardor|element-aeolis>` ; la reponse porte
  desormais `algorithm`. Sans option : `element-ardor`, liste inchangee. Element Aeolis : `import`, `heightfield`,
  `layer_surfaces`, `toolpath2d`, `collision`, `export` (libelles dans `apex/integration/REQUESTS.md`, R-026). Les fractions
  intermediaires d'Aeolis sont reelles (travail fait), comme R-017.
- Requete de tranchage (C2) : champ facultatif `algorithm` (`element-ardor` par defaut, ou `element-aeolis`) ; une valeur
  inconnue est refusee (code 2). Avec `element-aeolis`, `params` est un `AeolisSliceParams` (nouvelle racine du schema C1,
  `config/params.schema.json`, `cpp/include/apex/params.hpp`) et `gcode_params` est REFUSE (code 2) : le G-code est cartesien,
  sans pivot. Un champ de `params` inconnu reste une erreur.
- Resultat (C3) avec `element-aeolis` : meme objet que pour Ardor (`status`, `contract`, `gcode`, `report`, `warnings`) plus
  `algorithm`, `layers` (principales + inter-couches), `collisions_detected` (points de la tete en collision),
  `placement_offset_z`. Le `report` est celui d'Aeolis (`layer_count`, `intermediate_layer_count`, epaisseurs, pente max,
  `head`...). Avertissement `head_collisions` (`count`) si la vraie tete touche la matiere. Pas de `toolpath_points` ni de
  `step8_unreachable_orientations` (champs propres a Ardor).
- Refus type (C3, code 3) : statut `not_a_heightfield` quand le maillage n'est pas une surface en hauteur (surplomb, cavite,
  trou traversant debout), message lisible dans `message`.

## Ajouts du 2026-10-29 (le numero reste 2, a confirmer par l'utilisateur)
- G-code (C4, C5), signature : la PREMIERE ligne de tout G-code du moteur est `; 85Ejj-sV6dI-NXTGl-JiodA-0b8aZ` (en clair,
  sans signification, `file_io/apex_signature.py`). Le bloc chiffre porte aussi `meta.apex_version` (numero de version de
  l'application, lu dans `config/APEX_VERSION`, copie de `apex/VERSION` a mettre a jour a chaque montee de version).
- G-code (C4, C5), bloc d'en-tete chiffre (R-010) : `meta.machine` (`cartesian` ou `apex-polar-4axis`) et
  `meta.kinematics_description` (le meme objet `apex-kinematics/1` que la ligne `; apex_kinematics_json`) sont ecrits par
  les DEUX transcripteurs (polaire et non planaire, ce dernier n'ecrivait pas de bloc avant). Interface : lire ce bloc
  pour choisir le type de fichier au lieu de deviner d'apres la presence d'un mot B ou C (R-025).
- G-code (C4, C5), tous les algorithmes de tranchage : **chaque mouvement ecrit les six axes `X Y Z A B C`**, dans cet
  ordre, meme si l'un d'eux reste a 0 (decision de l'utilisateur). Machine Apex actuelle :
  `G00 X.. Y.. Z.. A0 B.. C0 [E..] F..` (A et C valent 0). Tranchage non planaire (cartesien) :
  `G1 X.. Y.. Z.. A0.000 B0.000 C0.000 E.. F..`. Un lecteur qui lit par lettre n'a rien a changer ; un lecteur a expression
  reguliere stricte doit accepter `A` et `C` (le moteur : `parse_gcode`, deja tolerant).
- G-code : description de la cinematique en tete de fichier, format `apex-kinematics/1` (`file_io/kinematics_description.py`) :
  une ligne lisible `; kinematics <machine> (apex-kinematics/1) : X=... Y=... ...` puis une ligne
  `; apex_kinematics_json {...}` (JSON sur une ligne, ASCII). Chaine serie de liaisons (cote pieces / cote outil), un
  axe par lettre de XYZABC (`linear`, `rotary` ou `unused`), vecteur directeur, plage, sens positif, decalage outil,
  conventions (compensation du pivot, signes). La ligne `; apex_kinematics pivot_to_tip=...` (R-007) reste.
  Champs ajoutes le meme jour : `part_origin` (point de la piece au pivot, axes a 0) et `conventions.rotation` (regle
  unique : main droite autour de `axis`). Exemple d'une machine decrite hors du code :
  `res/machines/corexy-tilting-bed.example.json` ; `kinematics_description.validate` refuse une description incoherente.

## Ajouts compatibles du 2026-10-27 (le numero reste 2)
Rien a faire pour l'interface qui lit deja le contrat 2 ; a connaitre :
- `machine_profile.toolhead_id` (texte, defaut `"aetheris-4"`, R-020) : identifiant du modele de tete,
  resolu par le moteur dans `res/3d/toolheads/<id>/<id>.glb` (jamais un chemin). Inconnu : journal, et
  `nozzle_model_checked=false` a l'etape 9.
- `params.allow_multi_mesh` vaut `true` par defaut (R-013) : un fichier a plusieurs objets est accepte, les
  transformations de nœud d'une scene (3MF, GLB) sont appliquees, une piece a volume positif n'est plus
  jetee comme « fragment ». Un fichier a un seul objet donne le meme resultat qu'avant.
- Resultat `surface_inconsistent` : `diagnostic` porte aussi `self_intersections` et `volume_pct` (R-016).
- G-code (C4, C5) : ligne en clair `; apex_kinematics pivot_to_tip=… b_sign=… tilt_sign=… z_reference=…`
  (R-007) et bloc chiffre `apex_block_begin v1` (R-010) portant `sliced` et `meta` (`kinematics` comprise),
  place apres l'en-tete en clair ; ecrit seulement si le paquet Python `cryptography` est present. Le
  marqueur `; apex_sliced = 1` reste.
- Progression (C3, R-017) : des lignes `{"type":"step","slug","fraction"}` avec une fraction intermediaire
  ]0, 1[ sont emises DANS chaque etape (croissantes, au plus une par pas de 1 % et par 100 ms, jamais
  estimees). Etapes a boucle : couches (5, 5b), remplissage (6), collision (9), ajustement (10), export ;
  etapes a phases indivisibles (maillage volumique, champ, cinematique inverse, ordonnancement) : quelques
  jalons seulement. `overhang` et `segmentation` n'en ont pas (instantanees). L'interface n'a rien a changer.
- Etapes (C3, R-022, 2026-10-28) : l'etape `ingest` est remplacee par `import` (lecture, nettoyage et reparation de la
  surface) puis `tetrahedralization` (TetGen ou TetraNet, puis structures d'acceleration), chacune de 0 a 1. Ancien
  identifiant `ingest` : disparu de `main.py steps` et des lignes de progression. Ajout compatible pour une interface
  qui affiche la liste de `steps` ; une interface qui testait le slug `ingest` en dur doit le retirer. Le numero de
  contrat reste 2.
- `params.skin` (objet ou `null`, defaut `null`, 2026-10-28) : lissage non-planaire sur la voie chaleur, opt-in
  (`element_ardor/skin/`). Champs : `thickness` (0,2 mm), `pattern` (`rectilinear` par defaut ; motifs de
  `surface_paths.PATTERNS`), `head_tolerance_deg` (45), `trim_radius`, `max_iter`, `loop_sampling`,
  `exclusion_rounds`, `plate_radius`, `face_sample`. Non mesure sur les quatre modeles : ne pas l'activer par defaut.
  Le G-code porte alors des passes `; --- lissage ---` ; `GcodeParams.skin_thickness` s'applique (avant : par erreur
  a tout point peint). Ajout compatible.
- `tests/vectors/kinematics_vectors.json` (R-008) : vecteurs de la cinematique pivot/pointe, generes par
  `scripts/gen_kinematics_vectors.py --check`.
- Modele de tete : celui de l'interface (v4, 16 237 sommets) remplace la v5 du moteur (R-009), copie par
  `scripts/sync_from_ui.py --check`. Effet mesure sur l'etape 9 : voir `DT.md`.

## 2 (2026-10-27) - voie unique : la methode de la chaleur ; FEM et Bird retires
INCOMPATIBLE. Decision de l'utilisateur : la voie par defaut (champ = distance au plateau par la
methode de la chaleur) devient la seule ; le champ FEM directionnel et la voie Bird sont supprimes.
Ce que l'interface doit faire :
- ne plus envoyer ni lire les champs retires de `params` : `field_solver.*` sauf `resolution`,
  `max_solver_iterations`, `convergence_tolerance`, `direct_solve_node_threshold`,
  `use_ilu_preconditioner`, `platform_anchor_weight`, `up_vector`, `heat_method_time_factor`,
  `radial_tilt_origin`, `heat_project_reachable` ; `constraints` (points d'attache : a refaire sur la
  voie chaleur) ; `segmentation.anchor_growth_weight` et `segmentation.branch_cut_threshold` ;
  `layers.support_driven_order`, `geodesic_order`, `frontier_order`, `nozzle_width`,
  `support_*` ;
- `layers.isovalue_order` vaut `true` par defaut ;
- ajout compatible : `collision.printed_part_margin` (mm, facultatif, defaut -0,15 dans le pipeline) ;
- un champ retire envoye dans une requete est une erreur (code 2), comme tout champ inconnu.

## 1 (2026-10-26) - premiere version numerotee
- Schema des parametres (C1) tel que genere par `scripts/gen_params.py`.
- Protocole d'appel du tranchage (C2, C3) : `main.py slice --request <requete.json>` et
  `main.py steps` (voir `slice_cli.py`). La requete porte `contract` ; un numero different est refuse
  (code de sortie 2).
