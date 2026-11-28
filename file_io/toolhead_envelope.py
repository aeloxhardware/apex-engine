"""Enveloppe parametrique d'une tete d'impression, quand son modele 3D n'existe pas encore (machine XYZAB, limite #97).

Repere du modele de tete (comme les GLB de `res/3d/toolheads/`) : origine = pointe de la buse, corps vers +Z. L'enveloppe
est une union de formes simples, decrites dans le profil de la machine (`toolhead_envelope`), toutes en mm :
  - buse : cone de la pointe (rayon nul) a `nozzle_top_radius`, haut de `nozzle_length` ;
  - bloc de chauffe : pave `block_size` = [largeur X, profondeur Y, hauteur Z], pose sur la buse ;
  - corps : cylindre de rayon `body_radius`, du haut du bloc jusqu'a `body_top` au-dessus de la pointe.
Les dimensions sont des donnees (R-030 point 4), pas des constantes du code. L'etape 9 en echantillonne la surface ;
les formes se recouvrent sans union booleenne."""
from __future__ import annotations

import numpy as np
import trimesh


def envelope_mesh(spec: dict) -> trimesh.Trimesh:
    nozzle_len = float(spec["nozzle_length"])
    nozzle = trimesh.creation.cone(radius=float(spec["nozzle_top_radius"]), height=nozzle_len, sections=24)
    nozzle.apply_transform(np.diag([1.0, 1.0, -1.0, 1.0]))          # pointe en bas, a l'origine
    nozzle.apply_translation([0.0, 0.0, nozzle_len])
    sx, sy, sz = (float(v) for v in spec["block_size"])
    block = trimesh.creation.box(extents=(sx, sy, sz))
    block.apply_translation([0.0, 0.0, nozzle_len + sz / 2])
    bottom, top = nozzle_len + sz, float(spec["body_top"])
    body = trimesh.creation.cylinder(radius=float(spec["body_radius"]), height=top - bottom, sections=32)
    body.apply_translation([0.0, 0.0, (bottom + top) / 2])
    mesh = trimesh.util.concatenate([nozzle, block, body])
    mesh.fix_normals()
    return mesh
