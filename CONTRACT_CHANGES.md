# Changements du contrat moteur <-> interface

Numero courant : voir `APEX_CONTRACT_VERSION` dans `config/schema.py` (copie dans
`config/params.schema.json`, cle `x-apex-contract`). Regles : `apex/integration/RULES.md`, §4.1 et §9.
Un changement INCOMPATIBLE incremente le numero et s'ajoute ici, avec ce que le consommateur doit faire.

## Ajout du 2026-11-04 (le numero reste 3) : options natives explicites (R-030 point 3)
- `Step6Params` gagne `native_fill` et `native_layer`, `Step8Params` gagne `native_grid_fallback` (tous vrais par defaut).
  Ils remplacent les variables d'environnement `APEX_NATIVE_FILL`, `APEX_NATIVE_LAYER` et `APEX_NATIVE_GRID`, qui n'ont plus
  aucun effet. Les valeurs par defaut etant celles d'avant, le G-code produit ne change pas : pas de hausse du contrat.
  Le consommateur n'a rien a faire (champs optionnels) ; `schema:sync` les fait apparaitre.

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
- Element Aeolis tranche TOUTE geometrie (`element-aeolis/general.py`). `AeolisSliceParams` gagne `mode` (`auto` par defaut,
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
  (`pipeline/skin/`). Champs : `thickness` (0,2 mm), `pattern` (`rectilinear` par defaut ; motifs de
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
