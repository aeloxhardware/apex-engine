"""Alias d'import : le code d'Element Ardor (pipeline 5 axes) vit dans `element-ardor/`.

Un dossier au nom à tiret n'est pas importable ; ce paquet redirige donc `pipeline.*` vers lui.
Tous les `from pipeline.stepN_... import ...` et `python -m pipeline._tetrahedralize_worker`
fonctionnent comme avant. Ne rien ajouter ici.
"""
from pathlib import Path

__path__ = [str(Path(__file__).resolve().parent.parent / "element-ardor")]
