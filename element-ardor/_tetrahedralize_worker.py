"""
pipeline/_tetrahedralize_worker.py

Worker isole (sous-processus dedie) pour UNE tentative de
`step1_ingest._tetrahedralize_once`.

Pourquoi ce module existe : `except RuntimeError` (dans
`step1_ingest._tetrahedralize`) ne peut rattraper que des erreurs Python.
TetGen est un binaire natif (via le binding `tetgen`) qui peut, sur
certaines geometries, planter par SEGFAULT plutot que de lever une
exception - un segfault tue tout le process Python qui l'a appele, sans
qu'aucun `except` ne puisse s'executer. Mesure concretement sur une version
de `res/test_models/death.stl` modifiee par l'utilisateur (2340 triangles,
recemment reduits en resolution) : le process se terminait avec le code de
sortie 139 (SIGSEGV) juste apres les memes avertissements TetGen deja
documentes CLAUDE.md limite #43 ("Two segments exactly intersect" etc.),
AVANT que `_repair_self_intersections` (pymeshfix) n'ait la moindre chance
de s'executer - le mecanisme de reparation reactive existant etait donc
inoperant sur ce cas precis, pas parce qu'il etait mal concu, mais parce
qu'un crash de process ne peut structurellement pas etre intercepte
in-process.

Ce worker fait tourner EXACTEMENT une tentative de
`step1_ingest._tetrahedralize_once` dans un sous-processus separe,
lance via `python -m pipeline._tetrahedralize_worker <input.pkl> <output.pkl>`.
L'appelant (`step1_ingest._tetrahedralize`) detecte un code de sortie
different de 0 (segfault inclus) exactement comme il detectait deja une
`RuntimeError` porteuse du message "self-intersection" - meme cascade de
repli (tel quel -> reparation pymeshfix -> sans raffinement plateau),
etendue pour couvrir aussi les crashs natifs, pas seulement les exceptions
Python.

Format d'echange (pickle, jamais de donnees sensibles, usage strictement
local au meme process de calcul) :
  - input.pkl : dict avec les memes arguments que `_tetrahedralize_once`
    (mesh serialise en vertices/faces bruts, pas un objet Trimesh complet -
    plus robuste au pickling inter-version).
  - output.pkl : soit {"ok": True, "tet_vertices": ..., "tetrahedra": ...},
    soit {"ok": False, "error": str(exception)} si TetGen a leve une
    exception Python NORMALE (le cas RuntimeError deja gere) - dans ce cas
    le worker se termine proprement avec le code 0, et c'est le contenu du
    pickle qui porte l'echec, pas le code de sortie. Un code de sortie NON
    NUL (ou pickle absent/corrompu) signale un crash natif, distinct d'une
    erreur Python normale.
"""

from __future__ import annotations

import pickle
import sys

import numpy as np
import trimesh


def _run(input_path: str, output_path: str) -> None:
    with open(input_path, "rb") as f:
        payload = pickle.load(f)

    from pipeline.step1_ingest import _tetrahedralize_once

    mesh = trimesh.Trimesh(
        vertices=payload["vertices"], faces=payload["faces"], process=False
    )

    # La sortie TEXTUELLE de TetGen (avertissements nommes - "Two segments
    # exactly intersect", indices de sommets/facettes precis) n'est PAS
    # capturable ici via `contextlib.redirect_stdout`/`redirect_stderr` -
    # tentative faite et mesuree comme inoperante : le binding tetgen est
    # une extension C qui ecrit DIRECTEMENT sur le descripteur de fichier
    # de l'OS (fd 1/2), en contournant l'objet Python `sys.stdout` que
    # `redirect_stdout` se contente de reassigner. La capture correcte se
    # fait au niveau du PROCESSUS PARENT, via `subprocess.run(...,
    # capture_output=True)` (deja en place dans
    # `step1_ingest._tetrahedralize_once_isolated`) qui, lui, capture bien
    # les descripteurs OS herites par CE sous-processus - aucune
    # modification necessaire ici, seulement cote parent.
    try:
        tet_mesh = _tetrahedralize_once(
            mesh,
            payload["resolution"],
            payload["near_platform_resolution"],
            payload["near_platform_margin"],
            payload["platform_z"],
            payload["junction_zones"],
            payload["junction_refinement_resolution"],
        )
        result = {
            "ok": True,
            "tet_vertices": np.asarray(tet_mesh.tet_vertices),
            "tetrahedra": np.asarray(tet_mesh.tetrahedra),
        }
    except Exception as exc:  # noqa: BLE001 - relaye tel quel au process parent
        result = {
            "ok": False,
            "error": str(exc),
            "error_type": type(exc).__name__,
        }

    with open(output_path, "wb") as f:
        pickle.dump(result, f)


if __name__ == "__main__":
    _input_path, _output_path = sys.argv[1], sys.argv[2]
    _run(_input_path, _output_path)
