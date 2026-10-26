# Corpus de référence

Pièces sur lesquelles TetraNet est comparé à TetGen. Les STL ne sont pas copiés ici : ils sont
lus dans `apex-engine/res/test_models/` (ou donnés par chemin). Liste visée :

| pièce | ce qu'elle exerce |
|---|---|
| `l_shape_overhang.stl` | cas simple, surplomb |
| `l-triple-axis-overhanged-shape.stl` | parois fines, surplombs sur 3 axes |
| `death.stl` | pièce de référence (Benchy à l'envers, 2 418 faces) |
| `death_hi-res.stl` | maillage dense (226 156 faces) : le cas de vitesse |
| `turkish.stl` | auto-intersections, cascade de réparation |
| Benchy officiel | **à fournir** (≈ 225 000 triangles ; TetGen : 411 s) |
| pièces à parois fines, cavités, plusieurs coques | **à ajouter** |

Les seuils d'acceptation sont dans `CRITERES.md` (à rédiger avec l'utilisateur avant le code).
