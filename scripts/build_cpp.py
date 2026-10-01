"""
scripts/build_cpp.py

Construit le module natif `apex_core` (cpp/) et le depose dans `cpp/lib/`.

Pourquoi un script plutot qu'un simple appel cmake : sur cette machine,
`cl.exe` n'est PAS dans le PATH (Visual Studio est installe hors des
emplacements par defaut). Ce script localise l'installation via `vswhere`,
capture l'environnement produit par `vcvars64.bat`, puis lance cmake/ninja
dedans - de sorte qu'une seule commande suffise, sans terminal developpeur
prealable.

Usage :
    .venv\\Scripts\\python.exe scripts/build_cpp.py [--clean] [--debug]
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CPP_DIR = REPO_ROOT / "cpp"
BUILD_DIR = CPP_DIR / "build"
LIB_DIR = CPP_DIR / "lib"


def _find_vcvars() -> Path:
    """Chemin de vcvars64.bat, via vswhere (l'outil officiel de localisation)."""
    program_files_x86 = os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")
    vswhere = Path(program_files_x86) / "Microsoft Visual Studio" / "Installer" / "vswhere.exe"
    if not vswhere.is_file():
        raise SystemExit(f"vswhere introuvable : {vswhere}\nVisual Studio Build Tools requis.")

    result = subprocess.run(
        [
            str(vswhere), "-products", "*", "-latest",
            "-requires", "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
            "-property", "installationPath",
        ],
        capture_output=True, text=True, check=True,
    )
    for line in result.stdout.splitlines():
        candidate = Path(line.strip()) / "VC" / "Auxiliary" / "Build" / "vcvars64.bat"
        if candidate.is_file():
            return candidate
    raise SystemExit(
        "Aucune installation Visual Studio avec le workload C++ (VC.Tools.x86.x64) trouvee."
    )


def _msvc_environment() -> dict[str, str]:
    """
    Environnement complet apres vcvars64.bat.

    `vcvars64.bat` positionne une trentaine de variables (PATH, INCLUDE, LIB,
    ...) ; on ne peut pas les deviner, il faut les lire apres execution. D'ou
    le `&& set` et le parsing : c'est la methode standard, pas un contournement.
    """
    vcvars = _find_vcvars()
    result = subprocess.run(
        f'"{vcvars}" >nul 2>&1 && set',
        shell=True, capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise SystemExit(f"vcvars64.bat a echoue :\n{result.stderr}")

    env: dict[str, str] = {}
    for line in result.stdout.splitlines():
        if "=" in line:
            key, _, value = line.partition("=")
            env[key] = value
    if "INCLUDE" not in env:
        raise SystemExit("vcvars64.bat n'a pas produit d'environnement MSVC exploitable.")
    return env


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clean", action="store_true", help="Supprime cpp/build avant de configurer")
    parser.add_argument("--debug", action="store_true", help="Compile en Debug plutot qu'en Release")
    parser.add_argument("--dev-name", default="build_dev",
                        help="Nom du repertoire de --dev sous cpp/ (un second permet de construire pendant que "
                             "le premier est verrouille)")
    parser.add_argument("--dev", action="store_true",
                        help="Construit dans cpp/build_dev/ (module dans cpp/build_dev/lib) sans toucher "
                             "cpp/lib, verrouille par un calcul en cours ; charger avec "
                             "APEX_CPP_LIB_DIR=cpp/build_dev/lib")
    args = parser.parse_args(argv)
    global BUILD_DIR, LIB_DIR
    if args.dev:
        BUILD_DIR = CPP_DIR / args.dev_name
        LIB_DIR = BUILD_DIR / "lib"

    if args.clean and BUILD_DIR.exists():
        print(f"[build_cpp] suppression de {BUILD_DIR}")
        shutil.rmtree(BUILD_DIR)

    env = _msvc_environment()
    build_type = "Debug" if args.debug else "Release"

    # cmake et ninja sont installes DANS le venv (paquets pip) : on les appelle
    # par leur chemin absolu plutot que de dependre du PATH.
    scripts_dir = Path(sys.executable).parent
    cmake = scripts_dir / "cmake.exe"
    # Le lanceur Scripts/cmake.exe du paquet pip a cesse de fonctionner
    # (code 1 sans aucune sortie, constate le 2026-09-30) alors que le binaire
    # qu'il enveloppe marche : on appelle ce binaire directement s'il existe.
    try:
        import cmake as _cmake_pkg  # noqa: PLC0415
        real = Path(_cmake_pkg.CMAKE_BIN_DIR) / "cmake.exe"
        if real.is_file():
            cmake = real
    except (ImportError, AttributeError):
        pass
    ninja = scripts_dir / "ninja.exe"
    if not cmake.is_file():
        raise SystemExit(f"cmake introuvable dans le venv : {cmake}\nFaire : pip install cmake ninja")

    import pybind11  # noqa: PLC0415 - dependance de build, pas d'execution

    BUILD_DIR.mkdir(parents=True, exist_ok=True)
    LIB_DIR.mkdir(parents=True, exist_ok=True)

    configure = [
        str(cmake), "-S", str(CPP_DIR), "-B", str(BUILD_DIR),
        "-G", "Ninja",
        f"-DCMAKE_MAKE_PROGRAM={ninja}",
        f"-DCMAKE_BUILD_TYPE={build_type}",
        f"-DPython_EXECUTABLE={sys.executable}",
        f"-Dpybind11_DIR={pybind11.get_cmake_dir()}",
        f"-DAPEX_LIB_DIR={LIB_DIR.as_posix()}",
    ]
    print("[build_cpp] configuration...")
    if subprocess.run(configure, env=env).returncode != 0:
        return 1

    print("[build_cpp] compilation...")
    if subprocess.run([str(cmake), "--build", str(BUILD_DIR)], env=env).returncode != 0:
        return 1

    produced = sorted(LIB_DIR.glob("apex_core*.pyd"))
    if not produced:
        print(f"[build_cpp] ECHEC : aucun .pyd produit dans {LIB_DIR}")
        return 1
    for path in produced:
        print(f"[build_cpp] OK -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
