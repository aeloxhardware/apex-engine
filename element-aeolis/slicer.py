"""Orchestration : maillage -> couches -> chemins -> G-code cartésien (X Y Z E).

Le débit suit l'épaisseur locale de la couche : section du cordon t·(w−t) + π·t²/4 (modèle
« stade »), volume = longueur 3D × section.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
import trimesh

from . import heightfield as hf
from . import layers as ly
from . import toolpath as tp
from core import progress as _progress
from file_io import kinematics_description as kin
from file_io.apex_header import header_block_lines
from file_io.apex_signature import SIGNATURE, common_meta
from file_io.gcode_axes import axes_words
from . import general as gen
from . import headcheck as hc
from .regions import mask_to_geometry


@dataclass
class SliceParams:
    # Hauteurs de couche : plancher et plafond viennent du profil d'imprimante (donc de l'UX).
    min_layer_height: float = ly.DEFAULT_MIN_LAYER
    max_layer_height: float = ly.DEFAULT_MAX_LAYER
    layer_height: float | None = None         # épaisseur nominale ; None : 0,2 (0,25 sans inter-couches)
    first_layer_height: float | None = None   # None : 0,2 borné par plancher et plafond
    # Buse de l'imprimante (profil d'imprimante, donc interface). Défauts TRÈS BAS : ils ne cassent pas une
    # imprimante à conduits de ventilation serrés (Anycubic...). À relever pour une tête dégagée.
    nozzle_reach: float = 2.0                 # « longueur maximale atteignable par la buse » : dénivelé max pointe / matière voisine (mm)
    nozzle_angle_deg: float = 10.0            # « angle de la buse » : pente max de la surface voisine depuis l'horizontale (°)
    enforce_nozzle_limits: bool = True        # le moteur « surface en hauteur » n'est gardé que s'il respecte ces limites
    mode: str = "auto"                        # "auto" | "heightfield" (refuse le reste) | "general" (tranchage général)
    check_head: bool = True                   # contrôle de la buse (enveloppe des deux réglages ci-dessus)
    head_model: str | None = None             # GLB d'une vraie tête (remplace l'enveloppe des deux réglages) ; None : réglages
    head_contact_clearance: float = 3.0       # zone de contact buse/cordon exclue du test (mm)
    head_safety_margin: float = 0.0
    intermediate_layers: bool = True          # inter-couches locales là où l'espacement dépasse le plafond
    grid_cell: float = 0.25
    toolpath: tp.ToolpathParams = field(default_factory=tp.ToolpathParams)
    filament_diameter: float = 1.75
    flow: float = 1.0
    print_speed: float = 40.0         # mm/s
    travel_speed: float = 100.0
    first_layer_speed: float = 20.0
    retract_length: float = 1.0
    retract_speed: float = 30.0
    z_hop: float = 0.6
    nozzle_temp: float = 200.0
    bed_temp: float = 50.0
    place_on_bed: bool = True         # pose le dessous de la pièce à z = 0
    center_xy: tuple[float, float] | None = None   # centre de l'empreinte sur le plateau


# Types de ligne écrits dans le G-code (`;TYPE:...`, convention de Cura, lue par la plupart des visionneuses) pour
# colorer l'aperçu : parois extérieure et intérieures, remplissage creux, peau pleine (dessus et dessous).
LINE_TYPES = {"wall_outer": "WALL-OUTER", "wall_inner": "WALL-INNER", "infill": "FILL", "solid_infill": "SKIN"}

# Étapes affichées pendant le tranchage Element Aeolis (`main.py steps --algorithm element-aeolis`) ;
# les libellés accentués sont ceux de l'interface (src/utils/sliceStepLabels.ts), ceux-ci en sont
# la forme du journal du moteur, sans accents, comme `main.PIPELINE_STEPS`.
AEOLIS_STEPS: tuple[tuple[str, str], ...] = (
    ("import", "Import du modele"),
    ("heightfield", "Modele de surface en hauteur"),
    ("layer_surfaces", "Surfaces des couches"),
    ("toolpath2d", "Parois et remplissage"),
    ("collision", "Simulation de collision"),
    ("export", "Ecriture du G-code"),
)


@dataclass
class SliceReport:
    layer_count: int                  # couches principales
    intermediate_layer_count: int
    fraction_below_min: float         # part de la longueur extrudée sous le plancher
    thickness_min: float
    thickness_max: float
    slope_max_deg: float
    extruded_length: float
    filament_length: float
    offset: tuple[float, float, float]
    head: hc.HeadReport = field(default_factory=hc.HeadReport)
    mode: str = "heightfield"         # moteur utilisé : "heightfield" ou "general"
    fallback_reason: str = ""         # pourquoi le moteur « surface en hauteur » n'a pas été gardé (mode auto)
    relief_mm: float = 0.0            # amplitude du relief des couches (tranchage général)


def _bead_area(t: np.ndarray, w: float) -> np.ndarray:
    t = np.minimum(t, w)
    return t * (w - t) + math.pi * t * t / 4.0


def _header_meta(p: "SliceParams", mode: str = "heightfield") -> dict:
    """Metadonnees du bloc d'en-tete chiffre (R-010) : la machine et sa cinematique, pour que l'apercu
    ne les devine pas d'apres les mots presents."""
    desc = kin.cartesian()
    return {**common_meta("Apex engine (Element Aeolis)"), "algorithm": "nonplanar-heightfield/1",
            "algorithm_name": "Element Aeolis",
            "machine": desc["machine"], "kinematics_description": desc,
            "min_layer_height": p.min_layer_height, "max_layer_height": p.max_layer_height,
            "intermediate_layers": p.intermediate_layers, "aeolis_mode": mode,
            "nozzle_reach": p.nozzle_reach, "nozzle_angle_deg": p.nozzle_angle_deg}


class _Stages:
    """Ouvre les étapes une à une et relaie leur progression RÉELLE (core/progress.py, R-017).

    `open(slug)` termine l'étape en cours (1.0), ouvre la suivante (0.0) et branche les boucles sur sa
    barre : une boucle appelle `_progress.report(fait, total)`, jamais une durée estimée. Sans
    rappel, tout est inerte."""

    def __init__(self, progress: Optional[Callable[[str, float], None]]):
        self.progress = progress
        self.current: Optional[str] = None

    def open(self, slug: Optional[str]) -> None:
        _progress.end()
        if self.progress is None:
            return
        if self.current is not None:
            self.progress(self.current, 1.0)
        self.current = slug
        if slug is not None:
            self.progress(slug, 0.0)
            _progress.begin(lambda fraction, slug=slug: self.progress(slug, fraction))

    def close(self) -> None:
        self.open(None)


def slice_path(path: str, params: SliceParams | None = None,
               progress: Optional[Callable[[str, float], None]] = None) -> tuple[str, SliceReport]:
    """Comme `slice_mesh`, mais lit le fichier (étape « import », avec sa propre barre)."""
    stages = _Stages(progress)
    try:
        stages.open("import")
        mesh = trimesh.load(path, force="mesh")
        return slice_mesh(mesh, params, progress=progress, _stages=stages)
    finally:
        _progress.end()


def slice_mesh(mesh: trimesh.Trimesh, params: SliceParams | None = None,
               progress: Optional[Callable[[str, float], None]] = None,
               _stages: "_Stages | None" = None) -> tuple[str, SliceReport]:
    """Tranche `mesh`. `progress(slug, fraction)` reçoit la progression des étapes de `AEOLIS_STEPS` :
    0 à l'ouverture, des fractions croissantes mesurées sur le travail fait, 1 à la fin."""
    stages = _stages or _Stages(progress)
    try:
        return _slice(mesh, params, stages)
    finally:
        stages.close()
        _progress.end()


def _place_mesh(mesh: trimesh.Trimesh, p: SliceParams):
    """Copie le maillage, le pose sur le plateau / le centre selon `p`. Retourne `(mesh, offset)`."""
    mesh = mesh.copy()
    off = np.zeros(3)
    if p.place_on_bed:
        off[2] = -mesh.bounds[0, 2]
    if p.center_xy is not None:
        c = mesh.bounds.mean(axis=0)
        off[0] = p.center_xy[0] - c[0]
        off[1] = p.center_xy[1] - c[1]
    mesh.apply_translation(off)
    return mesh, off


def _choose_engine(mesh: trimesh.Trimesh, p: SliceParams, stages: "_Stages"):
    """Choix du moteur : « surface en hauteur » si le modèle l'est ET que la buse y passe, sinon général.
    Retourne `(field_, stack, fallback_reason)` ; `field_` est None quand le moteur général s'impose."""
    stages.open("heightfield")
    field_ = None
    fallback_reason = ""
    if p.mode != "general":
        try:
            field_ = hf.from_mesh(mesh, p.grid_cell)
        except hf.NotAHeightfieldError as exc:
            if p.mode == "heightfield":
                raise
            fallback_reason = str(exc)
    stack = None
    if field_ is not None:
        stages.open("layer_surfaces")
        stack = ly.build(field_, p.layer_height, p.first_layer_height, min_layer=p.min_layer_height,
                         max_layer=p.max_layer_height, intermediate=p.intermediate_layers)
        _progress.milestone(0.5)
        if p.mode == "auto" and p.enforce_nozzle_limits:
            excess = gen.surfaces_violation(stack, field_, p.nozzle_angle_deg, p.nozzle_reach)
            if excess > p.min_layer_height:
                fallback_reason = (f"les couches de la surface en hauteur dépassent les limites de la buse "
                                   f"(angle {p.nozzle_angle_deg:g}°, longueur {p.nozzle_reach:g} mm) de {excess:.2f} mm")
                field_ = stack = None
    return field_, stack, fallback_reason


def _general_layers(mesh: trimesh.Trimesh, p: SliceParams, stages: "_Stages"):
    """Moteur général. Retourne `(plan, sequence, n, layers_data, layer_kinds)`."""
    analysis = dict(nozzle_angle_deg=p.nozzle_angle_deg, nozzle_reach=p.nozzle_reach, layer_height=p.layer_height,
                    first_layer=p.first_layer_height, min_layer=p.min_layer_height, max_layer=p.max_layer_height)
    if stages.current == "layer_surfaces":
        # Refus APRÈS le calcul des surfaces : on reste dans cette étape, sur la seconde moitié de sa barre.
        with _progress.phase(0.5, 1.0):
            plan = gen.analyze(mesh, **analysis)
            gen.compute_regions(plan)
    else:
        plan = gen.analyze(mesh, **analysis)
        stages.open("layer_surfaces")
        gen.compute_regions(plan)
    sequence = [(k, 0) for k in range(1, plan.count + 1)]
    stages.open("toolpath2d")
    layers_data, layer_kinds = gen.layer_paths(plan, p, p.first_layer_speed, p.print_speed)
    return plan, sequence, plan.count, layers_data, layer_kinds


def _heightfield_layer(p: SliceParams, field_, stack, grid, k: int, j: int, bounds: tuple):
    """Parois et remplissage d'une couche. Retourne `(speed, polylignes 3D, types)` ou None si vide.
    `grid` = `(gx, gy, grid_col)`, `bounds` = `(last_grid, split_grid)`."""
    xs, ys = field_.xs, field_.ys
    gx, _gy, grid_col = grid
    last_grid, split_grid = bounds
    mask = stack.active_mask(k, j, last_grid, split_grid)
    region = field_.footprint.intersection(mask_to_geometry(xs, ys, mask))
    if region.is_empty or region.area < 1e-3:
        return None
    zg = grid_col.layer_z(k, j).reshape(gx.shape)
    solid_mask = ((field_.top - zg) < p.toolpath.solid_top_mm) | \
        ((zg - field_.bottom) < p.toolpath.solid_bottom_mm)
    solid_mask = np.where(np.isnan(field_.top), False, solid_mask) & mask
    solid = field_.footprint.intersection(mask_to_geometry(xs, ys, solid_mask))
    speed = p.first_layer_speed if (k == 1 and j == 0) else p.print_speed
    # L'angle du remplissage alterne d'une couche à l'autre ; compté depuis le HAUT de chaque colonne pour que
    # la couche du dessus ait partout la même direction (sinon le dessus se compose de plaques d'orientations
    # différentes là où les colonnes finissent à des rangs différents).
    zones = None
    if j == 0:
        parity = np.where(np.isnan(last_grid), 0, (last_grid - k) % 2) == 0
        even = field_.footprint.intersection(mask_to_geometry(xs, ys, parity & mask))
        odd = field_.footprint.intersection(mask_to_geometry(xs, ys, ~parity & mask))
        zones = [(even, 45.0), (odd, 135.0)]
    path = tp.layer_polylines(region, solid, k * 100 + j, p.toolpath, field_.footprint, zones)
    xy_list = [tp.densify_xy(pl.xy, p.toolpath.max_segment) for pl in path.polylines]
    xy_all = np.vstack(xy_list) if xy_list else np.zeros((0, 2))
    z_all = stack.sample(xy_all[:, 0], xy_all[:, 1]).layer_z(k, j) if len(xy_all) else np.zeros(0)
    polys3d, kinds, lo = [], [], 0
    for pl, xy in zip(path.polylines, xy_list):
        pts = np.column_stack([xy, z_all[lo:lo + len(xy)]])
        lo += len(xy)
        if len(pts) < 2 or np.isnan(pts[:, 2]).any():
            continue
        polys3d.append(pts)
        kinds.append(pl.kind)
    return speed, polys3d, kinds


def _heightfield_layers(p: SliceParams, field_, stack, stages: "_Stages"):
    """Moteur « surface en hauteur » : parois et remplissage, une couche après l'autre.
    Retourne `(sequence, n, layers_data, layer_kinds)`."""
    sequence = stack.sequence()
    bounds = (stack.last_grid(), stack.split_grid())
    stages.open("toolpath2d")
    gx, gy = np.meshgrid(field_.xs, field_.ys)
    grid = (gx, gy, stack.sample(gx.ravel(), gy.ravel()))   # colonnes aux nœuds, calculées une fois
    layers_data = []   # (k, j, vitesse, polylignes 3D)
    layer_kinds = []   # types des polylignes, une liste par entrée de layers_data
    for done, (k, j) in enumerate(sequence):
        _progress.report(done, len(sequence))
        layer = _heightfield_layer(p, field_, stack, grid, k, j, bounds)
        if layer is None:
            continue
        speed, polys3d, kinds = layer
        layers_data.append((k, j, speed, polys3d))
        layer_kinds.append(kinds)
    return sequence, stack.count, layers_data, layer_kinds


def _head_collision_check(p: SliceParams, plan, field_, stack, layers_data) -> "hc.HeadReport":
    """Contrôle de la buse. Retourne le rapport (non vérifié si `check_head` est faux)."""
    head_env = None
    if p.check_head:
        head_env = hc.load_envelope(p.head_model, p.head_contact_clearance) if p.head_model \
            else hc.parametric_envelope(p.nozzle_angle_deg, p.nozzle_reach)
    head_rep = hc.HeadReport(checked=head_env is not None,
                             clearance_slope_deg=head_env.clearance_slope_deg if head_env else 0.0)
    if head_env is not None and plan is not None:
        gen.head_check(plan, layers_data, head_env, p.head_safety_margin, head_rep)
    elif head_env is not None:
        total_samples = sum(len(pts[::3]) for _, _, _, polys in layers_data for pts in polys) or 1
        sampled = 0
        for k, j, _speed, polys in layers_data:
            if not polys:
                continue
            samp = np.vstack([pts[::3] for pts in polys])

            def printed(x, y, k=k, j=j):
                col = stack.sample(x, y)
                return col.layer_z(k, j) - col.layer_thickness(k, j)

            hit, pen = hc.check_points(head_env, field_, None, None, samp[:, :2], samp[:, 2],
                                       p.head_safety_margin, printed_fn=printed)
            head_rep.points += len(samp)
            head_rep.colliding_points += hit
            head_rep.max_penetration = max(head_rep.max_penetration, pen)
            sampled += len(samp)
            _progress.report(sampled, total_samples)
    return head_rep


class _GcodeOut:
    """Écriture du G-code d'Aeolis : position courante, extrusion cumulée et statistiques d'épaisseur."""

    def __init__(self, p: SliceParams, report_mode: str, n: int, stack):
        self.p, self.stack = p, stack
        self.fil_area = math.pi * (p.filament_diameter / 2) ** 2
        self.out: list[str] = [
            SIGNATURE,
            "; apex_nonplanar = 1",
            *header_block_lines(_header_meta(p, report_mode)),
            "G21", "G90", "M82",
            f"M140 S{p.bed_temp:g}", f"M104 S{p.nozzle_temp:g}", "G28",
            f"M190 S{p.bed_temp:g}", f"M109 S{p.nozzle_temp:g}", "G92 E0",
            f"; layers = {n}",
            *kin.header_lines(kin.cartesian()),
        ]
        self.e = 0.0
        self.cur = None            # position courante (x, y, z)
        self.retracted = False
        self.total_len = 0.0
        self.thick_min, self.thick_max = np.inf, 0.0
        self.below_len = 0.0
        self.zero_axes = " ".join(axes_words({}, lambda v: f"{v:.3f}")[3:])   # « A0.000 B0.000 C0.000 »

    @staticmethod
    def axes6(x: float, y: float, z: float) -> str:
        """Les six axes XYZABC (file_io/gcode_axes.py) ; A, B, C valent 0 sur une machine cartésienne."""
        return " ".join(axes_words({"X": x, "Y": y, "Z": z}, lambda v: f"{v:.3f}"))

    def emit_travel(self, dest: np.ndarray, zc: float) -> None:
        p, out, cur = self.p, self.out, self.cur
        far = cur is None or np.linalg.norm(dest[:2] - cur[:2]) > 2.0
        if far and not self.retracted and cur is not None:
            out.append(f"G1 E{self.e - p.retract_length:.5f} F{p.retract_speed * 60:g}")
            self.retracted = True
        if cur is not None and far:
            out.append(f"G0 {self.axes6(cur[0], cur[1], zc)} F{p.travel_speed * 60:g}")
            out.append(f"G0 {self.axes6(dest[0], dest[1], zc)} F{p.travel_speed * 60:g}")
        out.append(f"G0 {self.axes6(dest[0], dest[1], dest[2])} F{p.travel_speed * 60:g}")
        if self.retracted:
            out.append(f"G1 E{self.e:.5f} F{p.retract_speed * 60:g}")
            self.retracted = False
        self.cur = dest.copy()

    def _layer_clearances(self, polys: list, k: int, j: int):
        """Hauteurs de dégagement et épaisseurs de TOUTE la couche en deux appels (les points de départ et
        d'arrivée des polylignes sont connus d'avance) : le dégagement d'une polyligne est le plus haut point
        de la couche entre la fin de la précédente et son départ. Retourne `(z_lines, t_flat)`."""
        cur, stack = self.cur, self.stack
        ramp = np.arange(12)[None, :, None]
        starts = np.array([pts[0, :2] for pts in polys])
        prev_end = np.array([(cur if i == 0 and cur is not None else polys[i - 1][-1])[:2]
                             if (i > 0 or cur is not None) else pts[0, :2]
                             for i, pts in enumerate(polys)])
        step_ = (starts - prev_end) / 11
        lines = ramp * step_[:, None, :] + prev_end[:, None, :]
        lines[:, -1, :] = starts
        col = stack.sample(lines[:, :, 0].ravel(), lines[:, :, 1].ravel())
        z_lines = np.nanmax(col.layer_z(k, j).reshape(len(polys), 12), axis=1)
        mids_all = np.vstack([(pts[:-1, :2] + pts[1:, :2]) / 2 for pts in polys])
        mcol = stack.sample(mids_all[:, 0], mids_all[:, 1])
        return z_lines, mcol.layer_thickness(k, j)

    def _polyline(self, pts: np.ndarray, t_all: np.ndarray, f_str: str) -> None:
        """Extrusion d'une polyligne (le déplacement jusqu'à son départ est déjà émis)."""
        p = self.p
        seg_len = np.linalg.norm(np.diff(pts, axis=0), axis=1)
        real = t_all > 1e-6
        if real.any():
            self.thick_min = min(self.thick_min, float(t_all[real].min()))
            self.thick_max = max(self.thick_max, float(t_all[real].max()))
            self.below_len += float(seg_len[real & (t_all < p.min_layer_height - 1e-6)].sum())
        vols = seg_len * _bead_area(np.maximum(t_all, 0.0), p.toolpath.line_width)
        e_run = self.e + np.cumsum(p.flow * vols / self.fil_area)
        self.total_len += float(seg_len.sum())
        self.out.extend([f"G1 X{x:.3f} Y{y:.3f} Z{z:.3f} {self.zero_axes} E{e_:.5f} F{f_str}"
                         for x, y, z, e_ in zip(pts[1:, 0].tolist(), pts[1:, 1].tolist(), pts[1:, 2].tolist(),
                                                e_run.tolist())])
        self.e = float(e_run[-1])
        self.cur = pts[-1].copy()

    def write_layers(self, layers_data: list, layer_kinds: list) -> None:
        p, out = self.p, self.out
        total_pts = sum(len(pts) for _, _, _, polys in layers_data for pts in polys) or 1
        written = 0
        for (k, j, speed, polys), kinds in zip(layers_data, layer_kinds):
            out.append(f"; <layer> {k}" + (f" intermediate {j}" if j else ""))
            if not polys:
                continue
            last_type = None
            z_lines, t_flat = self._layer_clearances(polys, k, j)
            t_lo = 0
            f_str = f"{speed * 60:g}"
            for i, pts in enumerate(polys):
                if self.cur is not None:
                    zc = max(float(z_lines[i]), float(self.cur[2])) + p.z_hop
                else:
                    zc = float(pts[0, 2]) + p.z_hop
                line_type = LINE_TYPES.get(kinds[i], "CUSTOM")
                if line_type != last_type:                    # une étiquette à chaque changement de type
                    out.append(f";TYPE:{line_type}")
                    last_type = line_type
                self.emit_travel(pts[0], zc)
                n_seg = len(pts) - 1
                self._polyline(pts, t_flat[t_lo:t_lo + n_seg], f_str)
                t_lo += n_seg
                written += len(pts)
                _progress.report(written, total_pts)

    def finish(self, head_env_checked: bool, head_rep) -> None:
        p, out = self.p, self.out
        cx, cy, cz = (self.cur[0], self.cur[1], self.cur[2]) if self.cur is not None else (0.0, 0.0, 0.0)
        out.append(f"G0 {self.axes6(cx, cy, cz + 5)} F{p.travel_speed * 60:g}")
        out.append("M104 S0\nM140 S0")
        if head_env_checked:
            out.insert(1, f"; head_collisions = {head_rep.colliding_points} / {head_rep.points} points, "
                          f"penetration max {head_rep.max_penetration:.2f} mm")


def _slice(mesh: trimesh.Trimesh, params: SliceParams | None, stages: "_Stages") -> tuple[str, SliceReport]:
    p = params or SliceParams()
    mesh, off = _place_mesh(mesh, p)
    if p.mode not in ("auto", "heightfield", "general"):
        raise ValueError(f"mode : 'auto', 'heightfield' ou 'general' attendu, pas {p.mode!r}")

    field_, stack, fallback_reason = _choose_engine(mesh, p, stages)
    report_mode = "heightfield" if field_ is not None else "general"
    plan = None
    if field_ is None:
        plan, sequence, n, layers_data, layer_kinds = _general_layers(mesh, p, stages)
        stack, field_ = plan.stack, plan.field
    else:
        sequence, n, layers_data, layer_kinds = _heightfield_layers(p, field_, stack, stages)

    stages.open("collision")
    head_rep = _head_collision_check(p, plan, field_, stack, layers_data)

    stages.open("export")
    writer = _GcodeOut(p, report_mode, n, stack)
    writer.write_layers(layers_data, layer_kinds)
    writer.finish(p.check_head and head_rep.checked, head_rep)
    out = writer.out
    report = SliceReport(
        layer_count=n, intermediate_layer_count=len(sequence) - n,
        fraction_below_min=writer.below_len / writer.total_len if writer.total_len else 0.0,
        thickness_min=float(writer.thick_min), thickness_max=float(writer.thick_max),
        slope_max_deg=plan.slope_deg if plan is not None else ly.max_slope_deg(field_, stack, stack.count),
        head=head_rep, extruded_length=writer.total_len, filament_length=writer.e,
        offset=(float(off[0]), float(off[1]), float(off[2])),
        mode=report_mode, fallback_reason=fallback_reason, relief_mm=plan.relief.amp if plan is not None else 0.0,
    )
    return "\n".join(out) + "\n", report
