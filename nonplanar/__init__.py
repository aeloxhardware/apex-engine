"""Alias d'import : le code d'Element Aeolis (non planaire) vit dans `element-aeolis/`.

Un dossier au nom à tiret n'est pas importable : ce paquet redirige `nonplanar.*` vers lui et
exécute son `__init__.py` (mêmes noms publics). Ne rien ajouter ici.
"""
from pathlib import Path

_real = Path(__file__).resolve().parent.parent / "element-aeolis"
__path__ = [str(_real)]
exec(compile((_real / "__init__.py").read_text(encoding="utf-8"), str(_real / "__init__.py"), "exec"))
