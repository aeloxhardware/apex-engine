"""
scripts/generate_test_model.py

Genere res/test_models/l_shape_overhang.stl : meme geometrie que
tests/test_step3_segmentation.py::l_shape_data (base verticale + bras
horizontal en porte-a-faux), reutilisee telle quelle - deja validee par la
suite de tests comme produisant un vrai surplomb detecte (step2) et une
direction de region alternative (step3), pas une geometrie inventee pour
l'occasion.

Execution : python scripts/generate_test_model.py
"""

from __future__ import annotations

import os

import trimesh


def build_l_shape_overhang() -> trimesh.Trimesh:
    vertical = trimesh.creation.box(extents=[4, 4, 20])
    vertical.apply_translation([0, 0, 10])

    horizontal = trimesh.creation.box(extents=[20, 4, 4])
    horizontal.apply_translation([8, 0, 20])

    return vertical.union(horizontal, engine="manifold")


def main() -> None:
    output_dir = os.path.join(os.path.dirname(__file__), "..", "res", "test_models")
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "l_shape_overhang.stl")

    mesh = build_l_shape_overhang()
    mesh.export(output_path)
    print(f"Modele de test ecrit : {os.path.abspath(output_path)}")
    print(f"  watertight: {mesh.is_watertight}, faces: {len(mesh.faces)}, extents: {mesh.extents}")


if __name__ == "__main__":
    main()
