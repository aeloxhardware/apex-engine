"""Etape 1 : exécution de TetGen (sous-processus isolé, exécutable natif) et adjacence des tétraèdres.

Découpé de `step1_ingest.py` (R-030 point 1, limite de 800 lignes par fichier) : déplacement pur, aucun
changement de comportement. Les noms restent importables depuis `step1_ingest`.
"""

from __future__ import annotations

import os
import pickle
import re
import subprocess
import sys
import tempfile
from typing import Optional
import numpy as np
import scipy.sparse as sp
import trimesh
from core import _cpp
from core.mesh_types import JunctionZone, TetMesh
from pipeline.step1_tetgen_seeds import (  # noqa: F401
    _find_internal_cavities,
    _find_rest_seed_outside_junction_zones,
    _merged_junction_seed,
    _split_for_platform_refinement,
)
from pipeline.step1_types import (  # noqa: F401
    MeshIngestError,
)


# Sous-chaines de messages RuntimeError leves par TetGen (relayes via
# `_tetrahedralize_once_isolated`) que la cascade de reparation reactive
# (limite #43) doit traiter comme RECUPERABLES - une reparation
# (pymeshfix) ou un repli (sans raffinement plateau) a une vraie chance de
# les resoudre, contrairement a une erreur structurelle non liee au mesh.
# "self-intersection" couvrait le refus explicite initial de TetGen ; le
# crash `recoversubfaces` est une VARIANTE differente mais tout aussi
# fatale du meme probleme sous-jacent (deja documentee dans CLAUDE.md
# limite #44 comme telle, jamais reliee au code jusqu'a ce fix) - mesure
# concretement sur res/test_models/turkish.stl : la cascade ne se
# declenchait PAS sur ce message, remontant l'erreur brute a l'appelant
# sans jamais essayer la reparation pymeshfix, alors qu'elle etait deja
# tentee (et fonctionnelle) pour "self-intersection".
_RECOVERABLE_TETGEN_ERROR_SUBSTRINGS = (
    "self-intersection",
    "recoversubfaces",
    # Meme classe d'echec interne TetGen que "recoversubfaces" (limite
    # #44), variante distincte rencontree en traitant turkish.stl avec
    # raffinement plateau actif (`near_platform_resolution` non None) : le
    # PLC combine plateau/reste (`_split_for_platform_refinement`) peut
    # faire echouer une AUTRE phase interne de recouvrement de frontiere
    # que la reconstruction de sous-facettes. Meme repli deja en place
    # (cascade de reparation, puis "sans raffinement plateau" en dernier
    # recours) s'applique sans changement.
    "split_segment",
)


def _is_recoverable_tetgen_error(exc: RuntimeError) -> bool:
    message = str(exc).lower()
    return any(substring in message for substring in _RECOVERABLE_TETGEN_ERROR_SUBSTRINGS)


_TETGEN_BRACKETED_INDICES = re.compile(r"\[([\d,\s]+)\]")


def _parse_tetgen_conflict_vertices(tetgen_output: str) -> set[int]:
    """
    Extrait les indices de sommets nommes explicitement par TetGen dans un
    message d'avertissement/erreur ("Two segments exactly intersect. 1st
    seg [869,836]...", "A vertex lies on a facet. vertex : [4340]...",
    "Two facets exactly intersect. 1st facet triangle: [...]...").

    Approche volontairement generique (pas un parseur par type de message) :
    TOUT groupe `[a,b,c]` d'entiers separes par virgule dans la sortie est
    un indice de sommet dans TOUS les formats de message observes (segment,
    facette, sommet isole) - `tag(-1)` voisin utilise des parentheses, pas
    des crochets, donc aucune confusion mesuree. Un parseur specifique par
    motif de phrase serait plus fragile a la moindre variation de
    formatage entre versions de TetGen.
    """
    indices: set[int] = set()
    for match in _TETGEN_BRACKETED_INDICES.finditer(tetgen_output):
        for token in match.group(1).split(","):
            token = token.strip()
            if token.isdigit():
                indices.add(int(token))
    return indices


_GEOMETRIC_RESOLUTION_MAX_TETGEN_PROBES = 5


def _probe_tetgen_conflicts(mesh: trimesh.Trimesh) -> str:
    """
    Sonde TetGen en mode DETECTION (`switches='d'`, leger - pas de
    maillage de qualite) pour recuperer sa sortie textuelle complete,
    porteuse des indices de sommets/segments/facettes en conflit.

    Mesure directe sur turkish.stl : la sortie de la VRAIE tetraedrisation
    (`_tetrahedralize_once_isolated`, mode qualite `pq1.2a...`) s'arrete
    net juste apres "Recovering boundaries...", AVANT les avertissements
    eux-memes, meme redirigee vers de vrais fichiers - tampon de la libc du
    binding jamais vide sur ce chemin de code precis. Le mode DETECTION,
    plus leger, a ete mesure comme flushant correctement dans tous les cas
    testes - d'ou ce sondage SEPARE, dans son propre sous-processus isole
    (`pipeline/_tetgen_detect_worker.py`), uniquement pour obtenir le
    texte, jamais pour produire un maillage tetraedrique reel.

    Best-effort, ne leve jamais : retourne une chaine vide si le
    sous-processus crashe ou si les fichiers de sortie sont illisibles -
    l'appelant traite deja ce cas comme "aucun sommet nommable".
    """
    try:
        with tempfile.TemporaryDirectory(prefix="apex_tetgen_probe_") as tmp_dir:
            input_path = os.path.join(tmp_dir, "input.pkl")
            stdout_path = os.path.join(tmp_dir, "stdout.txt")
            stderr_path = os.path.join(tmp_dir, "stderr.txt")

            with open(input_path, "wb") as f:
                pickle.dump(
                    {
                        "vertices": np.asarray(mesh.vertices, dtype=np.float64),
                        "faces": np.asarray(mesh.faces, dtype=np.int64),
                    },
                    f,
                )

            with open(stdout_path, "wb") as stdout_file, open(stderr_path, "wb") as stderr_file:
                subprocess.run(
                    [sys.executable, "-m", "pipeline._tetgen_detect_worker", input_path],
                    cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    stdout=stdout_file,
                    stderr=stderr_file,
                )

            with open(stdout_path, "rb") as f:
                text = f.read().decode("utf-8", errors="replace")
            with open(stderr_path, "rb") as f:
                text += f.read().decode("utf-8", errors="replace")
            return text
    except OSError:
        return ""


class _TetgenCrashError(RuntimeError):
    """
    Leve par `_tetrahedralize_once_isolated` quand le sous-processus TetGen
    se termine par un crash natif (code de sortie non nul, typiquement
    SIGSEGV) plutot qu'une exception Python normale - distinct d'une
    `RuntimeError` "self-intersection" habituelle UNIQUEMENT pour le
    message d'erreur ; traitee de facon identique par la cascade de repli
    de `_tetrahedralize` (meme `except` logique, capture separee pour
    documenter explicitement la distinction).
    """


class _TetgenNamedConflictError(RuntimeError):
    """
    `RuntimeError` normale (echec Python, pas un crash) MAIS accompagnee du
    texte complet imprime par TetGen (`.tetgen_output`), qui peut nommer
    precisement les sommets/segments/facettes en conflit ("Two segments
    exactly intersect", "A vertex lies on a facet", etc. - meme motif que
    la limite #43). Ces avertissements ne survivent normalement pas au
    passage par le sous-processus isole (`str(exc)` seul remonte, message
    generique) - `_tetrahedralize_worker.py` les capture desormais via
    `redirect_stdout`/`redirect_stderr` autour de l'appel reel. Utilisee
    par `_tetrahedralize` pour cibler une reparation locale
    (`_repair_local_self_intersections(..., extra_seed_faces=...)`) sur
    exactement les faces nommees, au lieu d'un repli generique.
    """

    def __init__(self, message: str, tetgen_output: str = ""):
        super().__init__(message)
        self.tetgen_output = tetgen_output


def _native_tetgen_worker_path() -> Optional[str]:
    """
    Executable C++ `apex_tetgen_worker` (cpp/third_party/tetgen, AGPL-3.0),
    quand il est construit. `APEX_TETGEN_WORKER=0` le desactive (retour au
    worker Python, identique). Absence = jamais une erreur.
    """
    if os.environ.get("APEX_TETGEN_WORKER", "1") == "0":
        return None
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dirs = [os.environ.get("APEX_CPP_LIB_DIR", ""), os.path.join(root, "cpp", "lib")]
    for d in dirs:
        for name in ("apex_tetgen_worker", "apex_tetgen_worker.exe"):
            path = os.path.join(d, name) if d else ""
            if path and os.path.isfile(path):
                return path
    return None


class _NativeTetgen:
    """
    Memes methodes que `tetgen.TetGen` (add_hole, add_region, tetrahedralize,
    node, elem) mais execute le PLC dans l'executable C++ separe, sous
    sous-processus : crash natif -> `_TetgenCrashError`, erreur TetGen ->
    `_TetgenNamedConflictError` avec le texte complet (stdout + stderr).
    """

    def __init__(self, vertices, faces):
        self._vertices = np.ascontiguousarray(vertices, dtype=np.float64)
        self._faces = np.ascontiguousarray(faces, dtype=np.int32)
        self._holes: list[list[float]] = []
        self._regions: list[list[float]] = []
        self.node = None
        self.elem = None

    def add_hole(self, point) -> None:
        self._holes.append([float(v) for v in point])

    def add_region(self, id: int, point, max_vol: float = 0.0) -> None:
        self._regions.append([float(v) for v in point] + [float(id), float(max_vol)])

    def tetrahedralize(self, switches: str) -> None:
        exe = _native_tetgen_worker_path()
        sw = switches.encode("ascii")
        with tempfile.TemporaryDirectory(prefix="apex_tetgen_") as tmp_dir:
            input_path = os.path.join(tmp_dir, "input.bin")
            output_path = os.path.join(tmp_dir, "output.bin")
            stdout_path = os.path.join(tmp_dir, "stdout.txt")
            stderr_path = os.path.join(tmp_dir, "stderr.txt")
            holes = np.asarray(self._holes, dtype=np.float64).reshape(-1, 3)
            regions = np.asarray(self._regions, dtype=np.float64).reshape(-1, 5)
            with open(input_path, "wb") as f:
                np.array(
                    [len(self._vertices), len(self._faces), len(holes), len(regions), len(sw)],
                    dtype=np.int32,
                ).tofile(f)
                f.write(sw)
                self._vertices.tofile(f)
                self._faces.tofile(f)
                holes.tofile(f)
                regions.tofile(f)
            # Vrais fichiers (pas des pipes) : meme raison que pour le worker
            # Python (bufferisation de la libc de TetGen).
            with open(stdout_path, "wb") as out_f, open(stderr_path, "wb") as err_f:
                completed = subprocess.run(
                    [exe, input_path, output_path], stdout=out_f, stderr=err_f
                )
            log_path = os.environ.get("APEX_TETGEN_LOG")
            if log_path:  # temps par phase de TetGen, pour la mesure
                with open(log_path, "ab") as log_f, open(stdout_path, "rb") as o, open(
                    stderr_path, "rb"
                ) as e:
                    log_f.write(f"=== {switches}\n".encode() + o.read() + e.read())
            if completed.returncode != 0 or not os.path.isfile(output_path):
                with open(stderr_path, "rb") as f:
                    tail = f.read().decode("utf-8", errors="replace")[-500:]
                raise _TetgenCrashError(
                    f"Le sous-processus TetGen s'est termine anormalement "
                    f"(code de sortie {completed.returncode}, probablement un "
                    f"crash natif) - stderr : {tail}"
                )
            with open(output_path, "rb") as f:
                status = int(np.fromfile(f, dtype=np.int32, count=1)[0])
                if status != 0:
                    n = int(np.fromfile(f, dtype=np.int32, count=1)[0])
                    message = f.read(n).decode("utf-8", errors="replace")
                    with open(stdout_path, "rb") as g:
                        output = g.read().decode("utf-8", errors="replace")
                    with open(stderr_path, "rb") as g:
                        output += g.read().decode("utf-8", errors="replace")
                    raise _TetgenNamedConflictError(message, output)
                nn, nt = (int(v) for v in np.fromfile(f, dtype=np.int32, count=2))
                self.node = np.fromfile(f, dtype=np.float64, count=nn * 3).reshape(nn, 3)
                self.elem = np.fromfile(f, dtype=np.int32, count=nt * 4).reshape(nt, 4)


def _tetrahedralize_once_isolated(
    mesh: trimesh.Trimesh,
    resolution: float,
    near_platform_resolution: Optional[float],
    near_platform_margin: float,
    platform_z: float,
    junction_zones: Optional[list[JunctionZone]],
    junction_refinement_resolution: Optional[float],
) -> TetMesh:
    """
    Fait tourner UNE tentative de `_tetrahedralize_once` dans un
    sous-processus isole (`pipeline/_tetrahedralize_worker.py`), pour
    survivre a un crash NATIF de TetGen (segfault) que ni `except
    RuntimeError` ni aucun autre mecanisme Python ne peut intercepter
    in-process - un tel crash tue tout le process qui l'a appele.

    Mesure concretement sur une version modifiee de `res/test_models/death.stl`
    (2340 triangles) : le process se terminait avec le code de sortie 139
    (SIGSEGV) des la 1ere tentative de la cascade de reparation existante,
    juste apres les memes avertissements TetGen deja documentes a la limite
    #43 ("Two segments exactly intersect" etc.) - AVANT que la reparation
    pymeshfix n'ait la moindre chance de s'executer, puisque le process qui
    aurait du l'appeler etait deja mort. L'isolation par sous-processus
    s'applique donc a CHAQUE tentative de la cascade (y compris la 1ere,
    "tel quel"), jamais seulement en dernier recours - c'est la seule facon
    de garantir qu'un crash sur le mesh original n'empeche pas d'essayer la
    reparation ensuite.

    Cout accepte : overhead de lancement d'un sous-processus Python par
    tentative (import de tetgen/scipy/trimesh a chaque fois) - largement
    inferieur au temps de la tetraedrisation elle-meme sur tout modele non
    trivial, et strictement necessaire pour la garantie de robustesse
    demandee (voir CLAUDE.md limite #46/#49, section "point ouvert
    death.stl").

    Serialisation deliberement minimale (vertices/faces bruts, pas un objet
    Trimesh complet) : plus robuste au pickling inter-version que de
    transporter des objets trimesh a travers la frontiere de process.

    Raises:
        RuntimeError: relayee telle quelle si le worker a lui-meme leve une
            RuntimeError Python normale (ex: message "self-intersection" de
            TetGen) - meme comportement qu'un appel in-process, la cascade
            de `_tetrahedralize` la traite exactement comme avant.
        _TetgenCrashError: si le sous-processus s'est termine par un crash
            natif (code de sortie non nul, ou fichier de sortie absent/
            illisible) plutot qu'une exception Python.
        MeshIngestError: relayee telle quelle si le worker a leve cette
            exception (entree structurellement invalide, pas un defaut de
            maillage a reparer).
    """
    if _native_tetgen_worker_path() is not None:
        return _tetrahedralize_once(
            mesh,
            resolution,
            near_platform_resolution,
            near_platform_margin,
            platform_z,
            junction_zones,
            junction_refinement_resolution,
            _tetgen_factory=_NativeTetgen,
        )

    with tempfile.TemporaryDirectory(prefix="apex_tetgen_") as tmp_dir:
        input_path = os.path.join(tmp_dir, "input.pkl")
        output_path = os.path.join(tmp_dir, "output.pkl")
        stdout_path = os.path.join(tmp_dir, "stdout.txt")
        stderr_path = os.path.join(tmp_dir, "stderr.txt")

        with open(input_path, "wb") as f:
            pickle.dump(
                {
                    "vertices": np.asarray(mesh.vertices, dtype=np.float64),
                    "faces": np.asarray(mesh.faces, dtype=np.int64),
                    "resolution": resolution,
                    "near_platform_resolution": near_platform_resolution,
                    "near_platform_margin": near_platform_margin,
                    "platform_z": platform_z,
                    "junction_zones": junction_zones,
                    "junction_refinement_resolution": junction_refinement_resolution,
                },
                f,
            )

        # Les avertissements nommes de TetGen ("Two segments exactly
        # intersect [869,836]...") sortent de l'extension C du binding par
        # ses descripteurs OS bruts, en contournant l'objet Python
        # `sys.stdout` - une redirection Python (`contextlib.redirect_*`,
        # essayee et rejetee) ne les voit jamais. `subprocess.run(...,
        # capture_output=True)` (PIPE) ne suffit PAS non plus : mesure
        # directe sur turkish.stl, le texte s'arretait net juste apres
        # "Recovering boundaries..." (avant les avertissements eux-memes),
        # signature d'un tampon de la libc du binding jamais vide vers un
        # PIPE avant la fin du sous-processus. Rediriger vers de VRAIS
        # FICHIERS resout le probleme (verifie : texte complet recupere,
        # y compris les indices de sommets) - la nature du descripteur
        # (fichier vs pipe) change le comportement de bufferisation cote
        # OS, meme si Python ne voit aucune difference de son point de vue.
        with open(stdout_path, "wb") as stdout_file, open(stderr_path, "wb") as stderr_file:
            completed = subprocess.run(
                [sys.executable, "-m", "pipeline._tetrahedralize_worker", input_path, output_path],
                cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                stdout=stdout_file,
                stderr=stderr_file,
            )

        if completed.returncode != 0 or not os.path.isfile(output_path):
            with open(stderr_path, "rb") as f:
                stderr_tail = f.read().decode("utf-8", errors="replace")[-500:]
            raise _TetgenCrashError(
                f"Le sous-processus TetGen s'est termine anormalement "
                f"(code de sortie {completed.returncode}, probablement un "
                f"crash natif) - stderr : {stderr_tail}"
            )

        with open(output_path, "rb") as f:
            result = pickle.load(f)

        if not result["ok"]:
            error_type = result.get("error_type", "RuntimeError")
            message = result["error"]
            if error_type == "MeshIngestError":
                raise MeshIngestError(message)
            with open(stdout_path, "rb") as f:
                tetgen_output = f.read().decode("utf-8", errors="replace")
            with open(stderr_path, "rb") as f:
                tetgen_output += f.read().decode("utf-8", errors="replace")
            raise _TetgenNamedConflictError(message, tetgen_output)

    tet_vertices = result["tet_vertices"]
    tetrahedra = result["tetrahedra"]
    tet_adjacency = _build_tet_adjacency(tetrahedra, len(tet_vertices))
    return TetMesh(
        tet_vertices=np.asarray(tet_vertices, dtype=np.float64),
        tetrahedra=np.asarray(tetrahedra, dtype=np.int64),
        tet_adjacency=tet_adjacency,
    )


def _tetrahedralize_once(
    mesh: trimesh.Trimesh,
    resolution: float,
    near_platform_resolution: Optional[float],
    near_platform_margin: float,
    platform_z: float,
    junction_zones: Optional[list[JunctionZone]] = None,
    junction_refinement_resolution: Optional[float] = None,
    _tetgen_factory=None,
) -> TetMesh:
    """
    Une passe de tetraedrisation, sans aucune tentative de reparation :
    laisse remonter telle quelle l'erreur de TetGen. Separee de
    `_tetrahedralize` pour que celle-ci puisse la rejouer sur un maillage
    repare sans dupliquer la construction des regions/trous.

    Jusqu'a 3 regions TetGen actives simultanement : "slab" (pres du
    plateau), "junction" (jonctions multi-region fusionnees en UN seul
    seed, cf. `_merged_junction_seed`), "rest" (le reste, a `resolution`).
    Contrairement au plateau, une jonction n'a pas de plan de coupe : le
    PLC combine (`combined_vertices`/`combined_faces`) reste celui du
    plateau seul (ou le mesh brut si le plateau n'est pas actif) - le seed
    de jonction est juste un point interne supplementaire, jamais une
    decoupe de surface.
    """
    if _tetgen_factory is None:
        import tetgen

        _tetgen_factory = tetgen.TetGen

    max_tet_volume = resolution ** 3
    hole_points = _find_internal_cavities(mesh)

    merged_junction_seed = None
    if junction_zones and junction_refinement_resolution is not None:
        merged_junction_seed = _merged_junction_seed(mesh, junction_zones)

    if near_platform_resolution is not None:
        combined_vertices, combined_faces, slab_seed, rest_seed = _split_for_platform_refinement(
            mesh, near_platform_margin, platform_z
        )
        if merged_junction_seed is not None:
            # Garde-fou EXPLICITE (pas un silence) : si slab_seed ou
            # rest_seed tombent a l'interieur du rayon d'une zone de
            # jonction, TetGen appliquerait le volume FIN de la jonction a
            # toute la region plateau/rest - meme bug mesure et corrige
            # pour le cas jonction seule (_find_rest_seed_outside_junction_zones).
            # Non re-derive automatiquement ici (deplacer slab_seed/rest_seed
            # casserait la garantie deja etablie par _split_for_platform_refinement
            # d'etre chacun a l'interieur de SA PROPRE coque) - signale plutot
            # qu'un resultat silencieusement corrompu.
            for seed_name, seed in [("slab_seed", slab_seed), ("rest_seed", rest_seed)]:
                for zone in junction_zones:
                    if np.linalg.norm(seed - zone.centroid) <= zone.radius:
                        raise MeshIngestError(
                            f"{seed_name} ({seed}) tombe a l'interieur d'une zone de "
                            f"jonction (centre={zone.centroid}, rayon={zone.radius}) - "
                            "raffinement plateau et raffinement jonction se chevauchent "
                            "sur ce modele, non supporte (near_platform_margin ou les "
                            "zones detectees devraient etre ajustes)."
                        )
        tet = _tetgen_factory(combined_vertices, combined_faces)
        for hole_point in hole_points:
            tet.add_hole(hole_point)
        if merged_junction_seed is not None:
            # ORDRE CRITIQUE quand une 3e region (jonction) est active,
            # mesure pas suppose (scripts/measure_add_region_*.py, session
            # limite #48) : la region la plus GROSSIERE ("rest") doit etre
            # add_region(1, ...), sinon la densite explose sur TOUT le
            # maillage (mesure sur le tripode synthetique : far<3mm passe de
            # 34 a 1150 tetraedres selon l'ordre, memes seeds/volumes). Ce
            # meme defaut existe deja avec seulement 2 regions
            # (slab=1,rest=2 vs rest=1,slab=2 : far<3mm 847 vs 29 sur un
            # modele synthetique) mais N'EST PAS corrige ici pour le cas
            # SANS jonction (slab_seed/rest_seed) : `l_shape_overhang.stl`
            # (modele de reference du projet, deja valide) produit
            # exactement 1677 tetraedres avec l'ordre actuel slab=1/rest=2 -
            # identique au chiffre deja documente CLAUDE.md limite #36,
            # aucune regression mesuree sur ce cas. Inverser l'ordre pour
            # TOUS les appels aurait ete un changement de comportement non
            # mesure sur le cas deja valide - reserve au nouveau cas (3
            # regions) ou le defaut est confirme.
            tet.add_region(1, rest_seed, max_tet_volume)
            tet.add_region(2, slab_seed, near_platform_resolution ** 3)
            tet.add_region(3, merged_junction_seed, junction_refinement_resolution ** 3)
        else:
            tet.add_region(1, slab_seed, near_platform_resolution ** 3)
            tet.add_region(2, rest_seed, max_tet_volume)
        # -A active les regions (chaque region utilise SON PROPRE max_vol
        # declare via add_region) ; -a sans valeur numerique explicite ne
        # fixe PAS de plafond global unique - verifie fonctionnel, aucun
        # crash, resultat mesure (CLAUDE.md).
        tet.tetrahedralize(switches="pq1.2Aa")
    elif merged_junction_seed is not None:
        # Pas de raffinement plateau actif : 2 regions seulement ("junction"
        # + "rest"), sur le mesh brut (pas de PLC combine necessaire, aucune
        # decoupe de surface pour une jonction - voir _merged_junction_seed).
        # Seed "rest" garanti HORS de toute zone de jonction (pas seulement
        # a l'interieur du mesh) - voir docstring de
        # _find_rest_seed_outside_junction_zones pour le bug mesure sinon.
        rest_seed = _find_rest_seed_outside_junction_zones(mesh, junction_zones)
        tet = _tetgen_factory(mesh.vertices, mesh.faces)
        for hole_point in hole_points:
            tet.add_hole(hole_point)
        # ORDRE CRITIQUE, mesure pas suppose (scripts/measure_add_region_*.py) :
        # "rest" (grossier) DOIT etre la region 1, "junction" (fin) la region 2
        # - inverser l'ordre (junction=1, rest=2) fait EXPLOSER la densite sur
        # tout le maillage (mesure : 3080 tetraedres corrects avec rest=1, vs
        # 46000+ avec junction=1, memes seeds, meme mesh, meme volumes cibles).
        # Comportement TetGen non documente, mais reproductible a l'identique
        # sur ce binding (0.8.4) - meme discipline que la limite #42 (TetGen
        # est le seul juge fiable de son propre comportement d'add_region).
        tet.add_region(1, rest_seed, max_tet_volume)
        tet.add_region(2, merged_junction_seed, junction_refinement_resolution ** 3)
        tet.tetrahedralize(switches="pq1.2Aa")
    else:
        tet = _tetgen_factory(mesh.vertices, mesh.faces)
        for hole_point in hole_points:
            tet.add_hole(hole_point)
        tet.tetrahedralize(switches=f"pq1.2a{max_tet_volume:.10f}")

    # Lire depuis les attributs de l'objet plutot que la valeur de retour :
    # le nombre de valeurs retournees par tetrahedralize() varie selon la
    # version de tetgen installee, tandis que .node/.elem sont stables.
    tet_vertices = tet.node
    tetrahedra = tet.elem

    tet_adjacency = _build_tet_adjacency(tetrahedra, len(tet_vertices))

    return TetMesh(
        tet_vertices=np.asarray(tet_vertices, dtype=np.float64),
        tetrahedra=np.asarray(tetrahedra, dtype=np.int64),
        tet_adjacency=tet_adjacency,
    )


def _build_tet_adjacency(tetrahedra: np.ndarray, node_count: int) -> sp.csr_matrix:
    """
    Construit la matrice d'adjacence creuse noeud-noeud du maillage
    tetraedrique (1 si deux noeuds partagent une arete d'un tetraedre).
    Reutilisee par l'assemblage du Laplacien discret a l'etape 4 du pipeline.
    """
    if _cpp.enabled():
        # L'ordre des aretes n'a aucune importance ici : `csr_matrix` trie ses
        # indices a la construction, donc la matrice finale est identique quel
        # que soit l'ordre fourni.
        edge_array = _cpp.module().tet_unique_edges(np.asarray(tetrahedra, dtype=np.int64))
        if len(edge_array) == 0:
            return sp.csr_matrix((node_count, node_count))
        rows = edge_array[:, 0]
        cols = edge_array[:, 1]
    else:
        edges = set()
        for tet in tetrahedra:
            for i in range(4):
                for j in range(i + 1, 4):
                    a, b = int(tet[i]), int(tet[j])
                    edges.add((min(a, b), max(a, b)))

        if len(edges) == 0:
            return sp.csr_matrix((node_count, node_count))

        rows, cols = zip(*edges)
        rows = np.array(rows)
        cols = np.array(cols)

    data = np.ones(len(rows) * 2)
    all_rows = np.concatenate([rows, cols])
    all_cols = np.concatenate([cols, rows])

    adjacency = sp.csr_matrix(
        (data, (all_rows, all_cols)), shape=(node_count, node_count)
    )
    return adjacency
