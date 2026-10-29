"""
pipeline/_tetgen_detect_worker.py

Worker isole (sous-processus dedie) pour UN sondage TetGen en mode
DETECTION (`switches='d'`), utilise par
`step1_ingest._probe_tetgen_conflicts`.

Pourquoi un worker SEPARE de `_tetrahedralize_worker.py` : ce dernier lance
la VRAIE tetraedrisation (`switches="pq1.2a..."`, mode qualite), dont la
sortie textuelle s'est averee TRONQUEE meme rediriges vers de vrais
fichiers (mesure directe sur turkish.stl : le texte s'arrete net juste
apres "Recovering boundaries...", avant les avertissements nommes eux-memes
- signature d'un tampon de la libc du binding jamais vide avant que le
sous-processus ne se termine, sur CE chemin de code precis). Le mode
DETECTION (`switches='d'`), plus leger (pas de maillage de qualite, juste
la phase de recouvrement de frontiere), a ete mesure comme flushant son
tampon correctement dans TOUS les cas testes - d'ou un worker dedie
n'executant QUE cette passe legere, jamais la vraie tetraedrisation.

Format d'echange : input.pkl = {"vertices":..., "faces":...} ; ce worker
n'ecrit AUCUN pickle de sortie (rien a en tirer, l'exception attendue est
systematique en mode detection sur un mesh non trivial) - seuls stdout et
stderr (rediriges vers de vrais fichiers par l'appelant, jamais des PIPE -
meme lecon que `_tetrahedralize_worker.py`) portent l'information utile.
Code de sortie toujours 0 sauf crash natif reel (le mode detection est
nettement plus leger que la tetraedrisation complete, un crash y est
possible mais pas la voie principale visee par ce module).
"""

from __future__ import annotations

import pickle
import sys

import tetgen


def _run(input_path: str) -> None:
    with open(input_path, "rb") as f:
        payload = pickle.load(f)

    tet = tetgen.TetGen(payload["vertices"], payload["faces"])
    try:
        tet.tetrahedralize(switches="d")
    except Exception:  # noqa: BLE001 - purement diagnostique, jamais relaye
        pass


if __name__ == "__main__":
    _run(sys.argv[1])
