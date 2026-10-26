# TetraNet

Mailleur tétraédrique du moteur Apex (Aelox Hardware). Il remplit le volume d'une pièce d'un
**réseau régulier de tétraèdres** (réseau cubique centré, BCC), puis **habille ce réseau à la
surface** — au lieu de construire les tétraèdres à partir de chaque triangle d'entrée comme
TetGen. Le temps dépend du nombre de cellules près de la peau, pas du nombre de triangles du STL.

> **Statut : prototype, avant tout code.** Ce dossier vit dans `apex-engine/tetranet/` le
> temps d'être validé ; il deviendra un dépôt à part (`aeloxhardware/tetranet`). TetGen reste le
> mailleur par défaut d'Apex tant que les mesures de `corpus/` ne disent pas le contraire.

## Objectif (reformulé par l'utilisateur, 2026-10-26)

Pour tout ce qui est tranchable avec Apex et imprimable, en réglant automatiquement la
tolérance selon la taille de la pièce **et la buse**, TetraNet ne donne pas un résultat
significativement différent de TetGen. « Significativement » = des seuils chiffrés, fixés dans
`corpus/CRITERES.md`, jamais une impression.

## Contrat avec le moteur

- Entrée : un maillage de surface (sommets `(N,3)`, faces `(M,3)`), une taille de maille, une
  tolérance de peau, optionnellement `platform_z` et une taille de maille locale près du plateau.
- Sortie : les mêmes tableaux que `core.mesh_types.TetMesh` (`tet_vertices (P,3)`,
  `tetrahedra (T,4)`), sans dépendre de Python. Les étapes 4, 5 et 8 du moteur ne changent pas.
- Il contrôle lui-même son résultat (couverture de la peau, rapport de volumes, qualité) et
  **dit** quand il est hors tolérance, pour que le moteur retombe sur TetGen.
- Unités : millimètres.

## Contenu

| dossier | rôle |
|---|---|
| `include/tetranet/` | API publique C++ (en-têtes) |
| `src/` | implémentation |
| `tests/` | tests (C++ et d'équivalence avec TetGen) |
| `corpus/` | pièces de référence et seuils d'acceptation |
| `tools/` | export VTK, comparaison TetGen / TetraNet |
| `docs/` | notes de conception |

## Construire

```
cmake -S . -B build -G Ninja && cmake --build build && ctest --test-dir build
```

C++20, `/fp:precise` (MSVC) ou `-ffp-contract=off` (GCC/Clang) : pas de `-ffast-math`.

## Licence

GPL-3.0 (`LICENSE`), comme Apex. Voir `NOTICE.md` pour les références.
