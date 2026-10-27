"""
Copie dans le moteur ce que l'interface produit (R-009 et R-010, RULES §4.2, contrats C6 et C5).

  - `src/assets/models/aetheris-4.glb` (modele de tete, decision de l'utilisateur du 2026-10-27)
        -> `res/3d/toolheads/aetheris-4/aetheris-4.glb` ;
  - `docs/apex_cipher_vectors.json` (vecteurs du chiffrement de l'en-tete)
        -> `tests/vectors/apex_cipher_vectors.json`.

Une note `<fichier>.source.txt` (ou `SOURCE.txt` a cote du GLB) dit l'origine et le SHA-256. Un fichier copie
ne se modifie jamais a la main : une erreur dedans se corrige dans `apex-ui`.

Usage : python scripts/sync_from_ui.py [--check]
  APEX_UI_DIR : dossier d'apex-ui (defaut : ../apex-ui ; dans `apex` : slicer-ux).
  --check     : n'ecrit rien, sort en erreur si une copie n'est plus a jour.
"""
from __future__ import annotations

import hashlib
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UI_DIR = os.environ.get("APEX_UI_DIR") or os.path.join(os.path.dirname(ROOT), "apex-ui")
# (chemin dans apex-ui, chemin dans le moteur, fichier de note)
FILES = [
    (("src", "assets", "models", "aetheris-4.glb"), ("res", "3d", "toolheads", "aetheris-4", "aetheris-4.glb"),
     ("res", "3d", "toolheads", "aetheris-4", "SOURCE.txt")),
    (("docs", "apex_cipher_vectors.json"), ("tests", "vectors", "apex_cipher_vectors.json"),
     ("tests", "vectors", "apex_cipher_vectors.source.txt")),
]


def sha256(path: str) -> str:
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def note_text(ui_rel: tuple, digest: str) -> str:
    return (f"Genere depuis apex-ui:{'/'.join(ui_rel)} par scripts/sync_from_ui.py -- ne pas editer.\n"
            f"sha256 {digest}\n")


def main() -> int:
    check = "--check" in sys.argv
    status = 0
    for ui_rel, dst_rel, note_rel in FILES:
        src, dst, note = (os.path.join(UI_DIR, *ui_rel), os.path.join(ROOT, *dst_rel), os.path.join(ROOT, *note_rel))
        if not os.path.isfile(src):
            print(f"introuvable : {src} (definir APEX_UI_DIR)")
            return 2
        digest = sha256(src)
        if check:
            ok = os.path.isfile(dst) and sha256(dst) == digest and os.path.isfile(note) \
                and open(note, encoding="utf-8").read() == note_text(ui_rel, digest)
            if not ok:
                print(f"{'/'.join(dst_rel)} n'est plus celui de l'interface : relancer scripts/sync_from_ui.py")
                status = 1
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(src, dst)
        with open(note, "w", encoding="utf-8", newline="\n") as f:
            f.write(note_text(ui_rel, digest))
    return status


if __name__ == "__main__":
    sys.exit(main())
