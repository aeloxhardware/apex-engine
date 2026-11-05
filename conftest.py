"""
conftest.py — point d'entree des tests : installe les options d'execution (R-030 point 3).

La bibliotheque ne lit jamais l'environnement ; la suite le fait ICI, pour pouvoir etre lancee
dans les deux positions du noyau natif (`APEX_USE_CPP=0` / `=1`) et avec un module construit
ailleurs (`APEX_CPP_LIB_DIR`).
"""

from core.runtime import configure_from_environment

configure_from_environment()
