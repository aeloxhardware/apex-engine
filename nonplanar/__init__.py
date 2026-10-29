"""Element Aeolis : tranchage non planaire « surface en hauteur » (couches z = S_k(x,y), buse verticale).

Second algorithme du moteur, indépendant du pipeline 5 axes (`pipeline/`) : il ne partage
avec lui ni le champ scalaire ni la tétraédrisation. Voir `nonplanar/README.md`.
"""
from .heightfield import Heightfield, NotAHeightfieldError, from_mesh
from .layers import LayerStack, build
from .slicer import SliceParams, SliceReport, slice_mesh

__all__ = ["Heightfield", "NotAHeightfieldError", "from_mesh", "LayerStack", "build",
           "SliceParams", "SliceReport", "slice_mesh"]
