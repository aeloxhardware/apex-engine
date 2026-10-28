"""
Insere les passes de lissage non-planaire dans la trajectoire machine du corps (apres l'etape 8).

L'ordre est celui de P5 (`print_order.order_and_check`) : chaque zone de lissage vient juste apres la derniere
couche du corps dont le rang est <= son moment, c'est-a-dire quand toute la matiere qui la porte est deposee. Les
points du lissage portent les angles de tete calcules par P5 (cinematique reelle, C impose par la position), la
premiere pointe de chaque polyligne est un trajet a vide (extrusion 0) et `is_skin` est vrai : l'export leur donne
l'epaisseur `skin_thickness` et releve la pointe de t/2 le long de l'axe de la tete, au sommet du cordon.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Optional

import numpy as np

from core.mesh_types import Layer, ToolpathMachine, ToolpathPoint


def inject(toolpath: ToolpathMachine, layers: list[Layer], passes: list[dict], feedrate: float,
           ) -> tuple[ToolpathMachine, np.ndarray]:
    """(trajectoire avec lissage, `old_to_new`) : `old_to_new[i]` est le nouvel indice du point `i` du corps (les
    diagnostics de l'etape 8 et toute reference par indice doivent etre remappes)."""
    ids = toolpath.layer_ids
    if ids is None:
        raise ValueError("la trajectoire n'a pas de layer_ids : lissage a inserer juste apres l'etape 8")
    n = len(toolpath.points)
    # pour chaque passe de lissage : couche du corps apres laquelle elle vient (index dans `layers`, -1 = avant tout)
    attach: dict[int, list[dict]] = {}
    last_body = -1
    for q in passes:
        if q["kind"] == "body":
            last_body = int(q["layer"])
        else:
            attach.setdefault(last_body, []).append(q)

    # position d'insertion (dernier point du corps appartenant a la couche, ou a la plus proche precedente)
    last_point_of = {}
    for i, lid in enumerate(ids):
        last_point_of[int(lid)] = i
    insert_after: dict[int, int] = {}      # indice de point du corps -> passes a inserer apres lui
    for layer_index, group in attach.items():
        pos = -1
        k = layer_index
        while k >= 0:
            pos = last_point_of.get(int(layers[k].id), -1)
            if pos >= 0:
                break
            k -= 1
        insert_after.setdefault(pos, []).extend(group)

    points: list[ToolpathPoint] = []
    new_ids: list[int] = []
    old_to_new = np.zeros(n, dtype=np.int64)

    def add_skin(group):
        for q in group:
            for path, angles in zip(q["points"], q["angles"]):
                path = np.asarray(path, dtype=float)
                angles = np.asarray(angles, dtype=float).reshape(len(path), -1)
                for j in range(len(path)):
                    points.append(ToolpathPoint(position=path[j], axis_angles=[float(a) for a in angles[j]],
                                                extrusion=0.0 if j == 0 else 1.0, feedrate=feedrate, is_skin=True))
                    new_ids.append(-1)

    if -1 in insert_after:
        add_skin(insert_after[-1])
    for i, p in enumerate(toolpath.points):
        old_to_new[i] = len(points)
        points.append(p)
        new_ids.append(int(ids[i]))
        if i in insert_after:
            add_skin(insert_after[i])
    return ToolpathMachine(points=points, layer_ids=np.asarray(new_ids, dtype=np.int64)), old_to_new
