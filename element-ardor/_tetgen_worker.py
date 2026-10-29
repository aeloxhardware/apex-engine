"""
Tetraedrisation de l autodocteur dans un SOUS-PROCESSUS (meme raison que
pipeline/_tetrahedralize_worker.py : TetGen peut planter nativement, et il ne
doit jamais etre lie dans le meme processus que le noyau, CLAUDE.md limite #56).

Appel : python -m pipeline._tetgen_worker <entree.npz> <sortie.npz>
entree : vertices, faces, switches (chaine vide = defauts du binding, comme
`TetGen.tetrahedralize()` sans argument dans la recherche).
sortie : node, elem -- les tableaux dont le binding construit sa grille.
"""
from __future__ import annotations

import sys

import numpy as np


def main(src: str, dst: str) -> None:
    import tetgen

    with np.load(src, allow_pickle=False) as npz:
        data = {k: npz[k] for k in npz.files}
    tet = tetgen.TetGen(np.asarray(data["vertices"], dtype=float), np.asarray(data["faces"], dtype=np.int64))
    switches = str(data["switches"])
    if switches:
        tet.tetrahedralize(switches=switches)
    else:
        tet.tetrahedralize()
    np.savez(dst, node=np.asarray(tet.node, dtype=float), elem=np.asarray(tet.elem, dtype=np.int64))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
