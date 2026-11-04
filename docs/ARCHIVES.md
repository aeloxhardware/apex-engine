# Archives (pointeur)

Les prototypes et anciens modules retirés de `main` vivent dans la branche **`archives`** (index :
`ARCHIVES.md` à sa racine). Règle : rien n'est effacé de `main` sans avoir d'abord été archivé là-bas.

Retirés le 2026-11-04 (R-030 point 7) : `scripts/prototype_local_retriangulation*.py` (14 fichiers),
`out/oldmods/` (3 modules), `scripts/check_step8_arrays.py`, et les 242 fichiers de sortie qui étaient
suivis dans `out/` (maintenant ignoré).

Retrouver un fichier : `git show archives:<chemin>` ou `git checkout archives -- <chemin>`.
