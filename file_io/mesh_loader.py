"""
io/mesh_loader.py

Chargement de maillages 3D depuis des fichiers externes (frontière du
pipeline vers le systeme de fichiers) - actuellement uniquement le modele
de tete d'impression, voir CLAUDE.md section "Modele 3D de tete
d'impression" pour le detail complet des conventions verifiees sur le
fichier reel (res/3d/toolheads/aetheris-4/aetheris-4.glb).
"""

from __future__ import annotations

import os
import re

import numpy as np
import trimesh

METERS_TO_MILLIMETERS = 1000.0

# Le corps de la tete s'etend selon +Y DANS LE FICHIER (axe long mesure :
# 77.2 mm en Y, contre 76.0 en X et 35.0 en Z), la pointe de buse etant a
# l'origine. Dans le repere machine, ce corps doit monter depuis le point de
# depot : +Y_fichier doit donc devenir +Z_machine, soit une rotation de +90
# deg autour de X.
#
# Correction posee sur MESURE, pas sur la convention documentee : sans elle,
# les bornes du mesh charge vont de z=-18.9 a z=+16.1 autour de la pointe,
# c'est-a-dire que 18.9 mm de tete pendent SOUS le point de depot. Effet
# constate en executant le pipeline complet : a une buse posee a z=1 mm,
# 47% des points d'echantillonnage de la tete tombent sous le plateau
# (z<0), 20% encore a z=11 mm - d'ou 3588 collisions detectees sur 3482
# points de trajectoire echantillonnes (~100%), rendant l'etape 10
# ininterpretable et interminable. Apres rotation : bornes z de 0.1 a 77.3
# (la tete est entierement AU-DESSUS de sa buse) et 0% de points sous le
# plateau a toute hauteur d'impression.
#
# ATTENTION : ceci CONTREDIT la convention "Axe Z = sens de la tete"
# documentee dans CLAUDE.md - la geometrie du fichier dit que l'axe long est
# Y, et Z (35 mm) est a peu pres centre sur la pointe. La mesure prime sur la
# note, mais si un futur modele de tete respecte reellement la convention Z,
# cette rotation devra devenir une propriete PAR MODELE (a cote du chemin du
# fichier) plutot qu'une constante globale.
TOOLHEAD_BODY_TO_MACHINE_UP = trimesh.transformations.rotation_matrix(
    np.pi / 2.0, [1.0, 0.0, 0.0]
)


def load_toolhead_mesh(path: str) -> trimesh.Trimesh:
    """
    Charge un GLB de tete d'impression et retourne un Trimesh UNIQUE, en
    espace millimetres, pret a etre positionne le long d'une trajectoire.

    Deux pieges evites, tous deux verifies concretement sur le fichier reel
    (voir CLAUDE.md) :

    1. `trimesh.load(path)` sur un GLB retourne une `Scene`, pas un
       `Trimesh` direct. Lire `scene.geometry[...]` brut donnerait les
       sommets dans l'espace LOCAL du noeud, PAS dans l'espace de la scene
       - ce fichier porte une rotation + mise a l'echelle composees sur la
       hierarchie de noeuds (verifie : extentes locales brutes tres
       differentes des extentes de scene). `scene.dump(concatenate=True)`
       applique cette transformation de noeud avant de retourner un Trimesh
       unique, seule facon correcte d'obtenir le mesh dans le repere
       documente (origine = pointe de buse, axes X/Y/Z de la tete).
    2. Le fichier source est en METRES (unite de base imposee par la spec
       glTF, confirme sur ce fichier) alors que tout le reste du pipeline
       travaille en millimetres - conversion x1000 obligatoire, faite ici
       une seule fois a la source plutot que dispersee chez chaque appelant.
    3. Le CORPS de la tete s'etend selon +Y dans le fichier, pas selon +Z -
       il faut donc une rotation de +90 deg autour de X pour l'amener dans
       le repere machine, ou la tete doit monter depuis la pointe de buse
       (voir TOOLHEAD_BODY_TO_MACHINE_UP ci-dessous).
    """
    scene = trimesh.load(path)
    if isinstance(scene, trimesh.Scene):
        mesh = scene.to_geometry()
    else:
        mesh = scene

    mesh = mesh.copy()
    mesh.apply_scale(METERS_TO_MILLIMETERS)
    mesh.apply_transform(TOOLHEAD_BODY_TO_MACHINE_UP)
    return mesh


TOOLHEADS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "res", "3d", "toolheads")


def toolhead_path(toolhead_id: str) -> str:
    """Chemin du GLB d'un modele de tete designe par son identifiant (R-020).

    L'identifiant vient du profil d'imprimante, donc de l'exterieur : seuls lettres, chiffres, `-` et `_`
    sont acceptes, pour qu'il ne puisse jamais sortir de `res/3d/toolheads/`."""
    if not re.fullmatch(r"[A-Za-z0-9_-]+", toolhead_id or ""):
        raise ValueError(f"identifiant de tete invalide : {toolhead_id!r}")
    return os.path.join(TOOLHEADS_DIR, toolhead_id, f"{toolhead_id}.glb")


def load_toolhead_by_id(toolhead_id: str) -> trimesh.Trimesh:
    """`load_toolhead_mesh` pour le modele designe par son identifiant ; `FileNotFoundError` s'il manque."""
    path = toolhead_path(toolhead_id)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"modele de tete inconnu : {toolhead_id!r} ({path})")
    return load_toolhead_mesh(path)
