"""Écrit deux petits G-code d'exemple (polaire et non planaire) pour tester la reconnaissance d'Apex par
l'interface : signature en clair, bloc d'en-tête chiffré, quelques mouvements.

    python scripts/gen_header_samples.py <dossier>
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def main() -> None:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    out.mkdir(parents=True, exist_ok=True)

    import trimesh
    from nonplanar import SliceParams, slice_mesh
    gcode, _ = slice_mesh(trimesh.creation.box([20, 20, 3]), SliceParams(grid_cell=1.0, check_head=False))
    (out / "sample_apex_nonplanar_cartesian.gcode").write_text("\n".join(gcode.splitlines()[:140]) + "\n")

    from main import _default_machine_profile
    from file_io.gcode_export import GcodeParams, export_gcode
    from tests.test_gcode_export import _toolpath
    machine = _default_machine_profile()
    lines, _ = export_gcode(_toolpath(machine), machine, GcodeParams())
    (out / "sample_apex_polar_4axis.gcode").write_text("\n".join(lines) + "\n")
    print("ecrit dans", out)


if __name__ == "__main__":
    main()
