# Changements du contrat moteur <-> interface

Numero courant : voir `APEX_CONTRACT_VERSION` dans `config/schema.py` (copie dans
`config/params.schema.json`, cle `x-apex-contract`). Regles : `apex/integration/RULES.md`, §4.1 et §9.
Un changement INCOMPATIBLE incremente le numero et s'ajoute ici, avec ce que le consommateur doit faire.

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
