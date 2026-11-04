# Archives

Cette branche (`archives`) garde **tout ce qui a été retiré de `main`** parce que c'était un prototype, une
copie d'ancien module ou un outil jetable : rien n'est perdu, rien n'encombre `main` (R-030 point 7).

- **Elle part du dernier état de `main` qui les contenait** (`29de86a`) : chaque fichier archivé s'y trouve
  à son chemin d'origine, avec tout son historique.
- **Règle** : on n'efface rien de `main` sans l'avoir d'abord archivé ici ; on ajoute une entrée à la table
  ci-dessous, dans la branche `archives`, puis on retire le fichier de `main`.
- **Retrouver un fichier** : `git show archives:<chemin>` (lire), `git checkout archives -- <chemin>` (le
  ramener dans la copie de travail).
- Cette branche n'est pas fusionnée dans `main` et ne doit pas l'être : son contenu ne fait pas partie du
  produit, aucune règle de limites (`scripts/check_limits.py`) ne s'y applique.

| chemin | ce que c'est | retiré de `main` le | pourquoi |
|---|---|---|---|
| `scripts/prototype_local_retriangulation.py`, `…2.py` à `…14.py` (14 fichiers) | itérations numérotées d'un prototype de re-triangulation locale ; chaque version reprend la précédente (`…14` exécute `…13` par `exec`) | 2026-11-04 | prototype jamais branché dans le pipeline, aucun test, aucun appelant |
| `out/oldmods/step1_old.py`, `step6_old.py`, `step8_before_arrays.py` | copies d'anciennes versions de modules, gardées pour comparer | 2026-11-04 | code mort versionné ; l'historique git remplit ce rôle |
| `scripts/check_step8_arrays.py` | comparait `step8` à `out/oldmods/step8_before_arrays.py` | 2026-11-04 | sa seule référence a été archivée |
| `out/` (242 fichiers suivis, 608 Mo : journaux, G-code, GIF, mesures) | sorties générées | 2026-11-04 | sorties, pas du code : `out/` est désormais ignoré par git |
