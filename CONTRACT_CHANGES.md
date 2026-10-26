# Changements du contrat moteur <-> interface

Numero courant : voir `APEX_CONTRACT_VERSION` dans `config/schema.py` (copie dans
`config/params.schema.json`, cle `x-apex-contract`). Regles : `apex/integration/RULES.md`, §4.1 et §9.
Un changement INCOMPATIBLE incremente le numero et s'ajoute ici, avec ce que le consommateur doit faire.

## 1 (2026-10-26) - premiere version numerotee
- Schema des parametres (C1) tel que genere par `scripts/gen_params.py`.
- Protocole d'appel du tranchage (C2, C3) : `main.py slice --request <requete.json>` et
  `main.py steps` (voir `slice_cli.py`). La requete porte `contract` ; un numero different est refuse
  (code de sortie 2).
