"""
io/mesh_loader.py

Chargement de maillages 3D depuis des fichiers externes (frontière du
pipeline vers le systeme de fichiers) - actuellement uniquement le modele
de tete d'impression, voir CLAUDE.md section "Modele 3D de tete
d'impression" pour le detail complet des conventions verifiees sur le
fichier reel (res/3d/toolheads/aetheris-4/aetheris-4.glb).
"""

from __future__ import annotations

import trimesh

METERS_TO_MILLIMETERS = 1000.0


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
    """
    scene = trimesh.load(path)
    if isinstance(scene, trimesh.Scene):
        mesh = scene.to_geometry()
    else:
        mesh = scene

    mesh = mesh.copy()
    mesh.apply_scale(METERS_TO_MILLIMETERS)
    return mesh
