# Element Ardor : pipeline principal 5 axes (méthode de la chaleur)

Les dix étapes du moteur (`step1_ingest.py` … `step10_adjustment.py`, plus `step5b`, `step8b`, `step9b` et
`skin/`). Orchestration : `main.py` à la racine ; briques partagées : `core/`, `file_io/`, `config/`, `cpp/`.

**Import** : le dossier porte un tiret, donc il n'est pas importable tel quel. Le paquet `pipeline/` à la racine le
redirige (`pipeline.step4_field_solver` → `element-ardor/step4_field_solver.py`). Ne pas déplacer ce code sans
mettre à jour cet alias. Détail et limites : `CLAUDE.md` à la racine.
