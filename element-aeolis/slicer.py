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
from . import headcheck as hc
from .regions import mask_to_geometry


@dataclass
class SliceParams:
    # Hauteurs de couche : plancher et plafond viennent du profil d'imprimante (donc de l'UX).
    min_layer_height: float = ly.DEFAULT_MIN_LAYER
    max_layer_height: float = ly.DEFAULT_MAX_LAYER
    layer_height: float | None = None         # épaisseur nominale ; None : 0,2 (0,25 sans inter-couches)
    first_layer_height: float | None = None   # None : 0,2 borné par plancher et plafond
    check_head: bool = True                   # collision de la vraie tête (buse verticale)
    head_model: str | None = None             # GLB de la tête ; None : aetheris-4
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


def _bead_area(t: np.ndarray, w: float) -> np.ndarray:
    t = np.minimum(t, w)
    return t * (w - t) + math.pi * t * t / 4.0


def _header_meta(p: "SliceParams") -> dict:
    """Metadonnees du bloc d'en-tete chiffre (R-010) : la machine et sa cinematique, pour que l'apercu
    ne les devine pas d'apres les mots presents."""
    desc = kin.cartesian()
    return {**common_meta("Apex engine (Element Aeolis)"), "algorithm": "nonplanar-heightfield/1",
            "algorithm_name": "Element Aeolis",
            "machine": desc["machine"], "kinematics_description": desc,
            "min_layer_height": p.min_layer_height, "max_layer_height": p.max_layer_height,
            "intermediate_layers": p.intermediate_layers}


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


def _slice(mesh: trimesh.Trimesh, params: SliceParams | None, stages: "_Stages") -> tuple[str, SliceReport]:
    p = params or SliceParams()
    mesh = mesh.copy()
    off = np.zeros(3)
    if p.place_on_bed:
        off[2] = -mesh.bounds[0, 2]
    if p.center_xy is not None:
        c = mesh.bounds.mean(axis=0)
        off[0] = p.center_xy[0] - c[0]
        off[1] = p.center_xy[1] - c[1]
    mesh.apply_translation(off)

    stages.open("heightfield")
    field_ = hf.from_mesh(mesh, p.grid_cell)
    stages.open("layer_surfaces")
    stack = ly.build(field_, p.layer_height, p.first_layer_height, min_layer=p.min_layer_height,
                     max_layer=p.max_layer_height, intermediate=p.intermediate_layers)
    _progress.milestone(0.5)
    sequence = stack.sequence()
    last_grid, split_grid = stack.last_grid(), stack.split_grid()
    xs, ys = field_.xs, field_.ys
    n = stack.count
    w = p.toolpath.line_width
    fil_area = math.pi * (p.filament_diameter / 2) ** 2

    # -- parois et remplissage : une couche après l'autre --------------------------------------
    stages.open("toolpath2d")
    gx, gy = np.meshgrid(xs, ys)
    gxr, gyr = gx.ravel(), gy.ravel()
    grid_col = stack.sample(gxr, gyr)          # colonnes aux nœuds, calculées une fois pour toutes les couches
    layers_data: list[tuple[int, int, float, list[np.ndarray]]] = []   # (k, j, vitesse, polylignes 3D)
    for done, (k, j) in enumerate(sequence):
        _progress.report(done, len(sequence))
        mask = stack.active_mask(k, j, last_grid, split_grid)
        region = field_.footprint.intersection(mask_to_geometry(xs, ys, mask))
        if region.is_empty or region.area < 1e-3:
            continue
        zg = grid_col.layer_z(k, j).reshape(gx.shape)
        solid_mask = ((field_.top - zg) < p.toolpath.solid_top_mm) | \
            ((zg - field_.bottom) < p.toolpath.solid_bottom_mm)
        solid_mask = np.where(np.isnan(field_.top), False, solid_mask) & mask
        solid = field_.footprint.intersection(mask_to_geometry(xs, ys, solid_mask))
        speed = p.first_layer_speed if (k == 1 and j == 0) else p.print_speed
        # L'angle du remplissage alterne d'une couche à l'autre ; compté depuis le HAUT de chaque
        # colonne pour que la couche du dessus ait partout la même direction (sinon le dessus se
        # compose de plaques d'orientations différentes là où les colonnes finissent à des rangs
        # différents).
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
        polys3d, lo = [], 0
        for xy in xy_list:
            pts = np.column_stack([xy, z_all[lo:lo + len(xy)]])
            lo += len(xy)
            if len(pts) < 2 or np.isnan(pts[:, 2]).any():
                continue
            polys3d.append(pts)
        layers_data.append((k, j, speed, polys3d))

    # -- collision de la vraie tête -------------------------------------------------------------
    head_env = hc.load_envelope(p.head_model, p.head_contact_clearance) if p.check_head else None
    head_rep = hc.HeadReport(checked=head_env is not None,
                             clearance_slope_deg=head_env.clearance_slope_deg if head_env else 0.0)
    stages.open("collision")
    if head_env is not None:
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

    # -- G-code ----------------------------------------------------------------------------------
    stages.open("export")
    out: list[str] = [
        SIGNATURE,
        "; apex_nonplanar = 1",
        *header_block_lines(_header_meta(p)),
        "G21", "G90", "M82",
        f"M140 S{p.bed_temp:g}", f"M104 S{p.nozzle_temp:g}", "G28",
        f"M190 S{p.bed_temp:g}", f"M109 S{p.nozzle_temp:g}", "G92 E0",
        f"; layers = {n}",
        *kin.header_lines(kin.cartesian()),
    ]
    e = 0.0
    cur = None            # position courante (x, y, z)
    retracted = False
    total_len = 0.0
    thick_min, thick_max = np.inf, 0.0
    below_len = 0.0

    def axes6(x: float, y: float, z: float) -> str:
        """Les six axes XYZABC (file_io/gcode_axes.py) ; A, B, C valent 0 sur une machine cartésienne."""
        return " ".join(axes_words({"X": x, "Y": y, "Z": z}, lambda v: f"{v:.3f}"))

    def emit_travel(dest: np.ndarray, zc: float) -> None:
        nonlocal cur, retracted
        far = cur is None or np.linalg.norm(dest[:2] - cur[:2]) > 2.0
        if far and not retracted and cur is not None:
            out.append(f"G1 E{e - p.retract_length:.5f} F{p.retract_speed * 60:g}")
            retracted = True
        if cur is not None and far:
            out.append(f"G0 {axes6(cur[0], cur[1], zc)} F{p.travel_speed * 60:g}")
            out.append(f"G0 {axes6(dest[0], dest[1], zc)} F{p.travel_speed * 60:g}")
        out.append(f"G0 {axes6(dest[0], dest[1], dest[2])} F{p.travel_speed * 60:g}")
        if retracted:
            out.append(f"G1 E{e:.5f} F{p.retract_speed * 60:g}")
            retracted = False
        cur = dest.copy()

    zero_axes = " ".join(axes_words({}, lambda v: f"{v:.3f}")[3:])   # « A0.000 B0.000 C0.000 »
    total_pts = sum(len(pts) for _, _, _, polys in layers_data for pts in polys) or 1
    written = 0
    ramp = np.arange(12)[None, :, None]
    for k, j, speed, polys in layers_data:
        out.append(f"; <layer> {k}" + (f" intermediate {j}" if j else ""))
        if not polys:
            continue
        # Hauteurs de dégagement et épaisseurs de TOUTE la couche en deux appels (les points de départ
        # et d'arrivée des polylignes sont connus d'avance) : le dégagement d'une polyligne est le plus
        # haut point de la couche entre la fin de la précédente et son départ.
        starts = np.array([pts[0, :2] for pts in polys])
        prev_end = np.array([(cur if i == 0 and cur is not None else polys[i - 1][-1])[:2]
                             if (i > 0 or cur is not None) else pts[0, :2]
                             for i, pts in enumerate(polys)])
        step_ = (starts - prev_end) / 11
        lines = ramp * step_[:, None, :] + prev_end[:, None, :]
        lines[:, -1, :] = starts
        col = stack.sample(lines[:, :, 0].ravel(), lines[:, :, 1].ravel())
        z_lines = np.nanmax(col.layer_z(k, j).reshape(len(polys), 12), axis=1)
        mids = [(pts[:-1, :2] + pts[1:, :2]) / 2 for pts in polys]
        mids_all = np.vstack(mids)
        mcol = stack.sample(mids_all[:, 0], mids_all[:, 1])
        t_flat = mcol.layer_thickness(k, j)
        t_lo = 0
        f_str = f"{speed * 60:g}"
        for i, pts in enumerate(polys):
            if cur is not None:
                zc = max(float(z_lines[i]), float(cur[2])) + p.z_hop
            else:
                zc = float(pts[0, 2]) + p.z_hop
            emit_travel(pts[0], zc)
            seg_v = np.diff(pts, axis=0)
            seg_len = np.linalg.norm(seg_v, axis=1)
            t_all = t_flat[t_lo:t_lo + len(seg_len)]
            t_lo += len(seg_len)
            real = t_all > 1e-6
            if real.any():
                thick_min = min(thick_min, float(t_all[real].min()))
                thick_max = max(thick_max, float(t_all[real].max()))
                below_len += float(seg_len[real & (t_all < p.min_layer_height - 1e-6)].sum())
            vols = seg_len * _bead_area(np.maximum(t_all, 0.0), w)
            e_run = e + np.cumsum(p.flow * vols / fil_area)
            total_len += float(seg_len.sum())
            out.extend([f"G1 X{x:.3f} Y{y:.3f} Z{z:.3f} {zero_axes} E{e_:.5f} F{f_str}"
                        for x, y, z, e_ in zip(pts[1:, 0].tolist(), pts[1:, 1].tolist(), pts[1:, 2].tolist(),
                                               e_run.tolist())])
            e = float(e_run[-1])
            cur = pts[-1].copy()
            written += len(pts)
            _progress.report(written, total_pts)
    cx, cy, cz = (cur[0], cur[1], cur[2]) if cur is not None else (0.0, 0.0, 0.0)
    out.append(f"G0 {axes6(cx, cy, cz + 5)} F{p.travel_speed * 60:g}")
    out.append("M104 S0\nM140 S0")
    if head_env is not None:
        out.insert(1, f"; head_collisions = {head_rep.colliding_points} / {head_rep.points} points, "
                      f"penetration max {head_rep.max_penetration:.2f} mm")
    report = SliceReport(
        layer_count=n, intermediate_layer_count=len(sequence) - n,
        fraction_below_min=below_len / total_len if total_len else 0.0, thickness_min=float(thick_min), thickness_max=float(thick_max),
        slope_max_deg=ly.max_slope_deg(field_, stack, stack.count), head=head_rep, extruded_length=total_len,
        filament_length=e, offset=(float(off[0]), float(off[1]), float(off[2])),
    )
    return "\n".join(out) + "\n", report
