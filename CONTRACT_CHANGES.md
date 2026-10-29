# Changements du contrat moteur <-> interface

Numero courant : voir `APEX_CONTRACT_VERSION` dans `config/schema.py` (copie dans
`config/params.schema.json`, cle `x-apex-contract`). Regles : `apex/integration/RULES.md`, §4.1 et §9.
Un changement INCOMPATIBLE incremente le numero et s'ajoute ici, avec ce que le consommateur doit faire.

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
