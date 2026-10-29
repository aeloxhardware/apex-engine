"""Orchestration : maillage -> couches -> chemins -> G-code cartésien (X Y Z E).

Le débit suit l'épaisseur locale de la couche : section du cordon t·(w−t) + π·t²/4 (modèle
« stade »), volume = longueur 3D × section.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import trimesh

from . import heightfield as hf
from . import layers as ly
from . import toolpath as tp
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


def slice_mesh(mesh: trimesh.Trimesh, params: SliceParams | None = None) -> tuple[str, SliceReport]:
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

    field_ = hf.from_mesh(mesh, p.grid_cell)
    stack = ly.build(field_, p.layer_height, p.first_layer_height, min_layer=p.min_layer_height,
                     max_layer=p.max_layer_height, intermediate=p.intermediate_layers)
    sequence = stack.sequence()
    last_grid, split_grid = stack.last_grid(), stack.split_grid()
    xs, ys = field_.xs, field_.ys
    n = stack.count
    w = p.toolpath.line_width
    fil_area = math.pi * (p.filament_diameter / 2) ** 2

    out: list[str] = [
        "; apex_nonplanar = 1",
        "G21", "G90", "M82",
        f"M140 S{p.bed_temp:g}", f"M104 S{p.nozzle_temp:g}", "G28",
        f"M190 S{p.bed_temp:g}", f"M109 S{p.nozzle_temp:g}", "G92 E0",
        f"; layers = {n}",
    ]
    e = 0.0
    cur = None            # position courante (x, y, z)
    retracted = False
    total_len = 0.0
    thick_min, thick_max = np.inf, 0.0
    below_len = 0.0
    head_env = hc.load_envelope(p.head_model, p.head_contact_clearance) if p.check_head else None
    head_rep = hc.HeadReport(checked=head_env is not None,
                             clearance_slope_deg=head_env.clearance_slope_deg if head_env else 0.0)

    def emit_travel(dest: np.ndarray, zc: float) -> None:
        nonlocal cur, retracted
        far = cur is None or np.linalg.norm(dest[:2] - cur[:2]) > 2.0
        if far and not retracted and cur is not None:
            out.append(f"G1 E{e - p.retract_length:.5f} F{p.retract_speed * 60:g}")
            retracted = True
        if cur is not None and far:
            out.append(f"G1 Z{zc:.3f} F{p.travel_speed * 60:g}")
        out.append(f"G0 X{dest[0]:.3f} Y{dest[1]:.3f} F{p.travel_speed * 60:g}")
        out.append(f"G1 Z{dest[2]:.3f} F{p.travel_speed * 60:g}")
        if retracted:
            out.append(f"G1 E{e:.5f} F{p.retract_speed * 60:g}")
            retracted = False
        cur = dest.copy()

    for k, j in sequence:
        mask = stack.active_mask(k, j, last_grid, split_grid)
        region = field_.footprint.intersection(mask_to_geometry(xs, ys, mask))
        if region.is_empty or region.area < 1e-3:
            continue
        gx, gy = np.meshgrid(xs, ys)
        zg = (stack.z(k, gx.ravel(), gy.ravel()) if j == 0
              else stack.z_sub(k, j, gx.ravel(), gy.ravel())).reshape(gx.shape)
        solid_mask = ((field_.top - zg) < p.toolpath.solid_top_mm) | \
            ((zg - field_.bottom) < p.toolpath.solid_bottom_mm)
        solid_mask = np.where(np.isnan(field_.top), False, solid_mask) & mask
        solid = field_.footprint.intersection(mask_to_geometry(xs, ys, solid_mask))
        out.append(f"; <layer> {k}" + (f" intermediate {j}" if j else ""))
        speed = p.first_layer_speed if (k == 1 and j == 0) else p.print_speed
        zfun = (lambda x, y, k=k: stack.z(k, x, y)) if j == 0 else \
            (lambda x, y, k=k, j=j: stack.z_sub(k, j, x, y))
        tfun = (lambda x, y, k=k: stack.thickness(k, x, y)) if j == 0 else \
            (lambda x, y, k=k, j=j: stack.thickness_sub(k, j, x, y))
        path = tp.layer_polylines(region, solid, k * 100 + j, p.toolpath)
        for pl in path.polylines:
            pts = tp.densify(pl.xy, stack, k, p.toolpath.max_segment, j)
            if len(pts) < 2 or np.isnan(pts[:, 2]).any():
                continue
            # Dégagement : plus haut point de la couche entre départ et arrivée, plus z_hop.
            if cur is not None:
                line = np.linspace(cur[:2], pts[0, :2], 12)
                zc = max(float(np.nanmax(zfun(line[:, 0], line[:, 1]))), float(cur[2])) + p.z_hop
            else:
                zc = float(pts[0, 2]) + p.z_hop
            emit_travel(pts[0], zc)
            if head_env is not None:
                samp = pts[::3]
                hit, pen = hc.check_points(head_env, field_, zfun, tfun, samp[:, :2], samp[:, 2],
                                           p.head_safety_margin)
                head_rep.points += len(samp)
                head_rep.colliding_points += hit
                head_rep.max_penetration = max(head_rep.max_penetration, pen)
            seg_v = np.diff(pts, axis=0)
            seg_len = np.linalg.norm(seg_v, axis=1)
            mids = (pts[:-1, :2] + pts[1:, :2]) / 2
            t_all = tfun(mids[:, 0], mids[:, 1])
            real = t_all > 1e-6
            if real.any():
                thick_min = min(thick_min, float(t_all[real].min()))
                thick_max = max(thick_max, float(t_all[real].max()))
                below_len += float(seg_len[real & (t_all < p.min_layer_height - 1e-6)].sum())
            vols = seg_len * _bead_area(np.maximum(t_all, 0.0), w)
            e_run = e + np.cumsum(p.flow * vols / fil_area)
            total_len += float(seg_len.sum())
            f_str = f"{speed * 60:g}"
            for b_, e_ in zip(pts[1:], e_run):
                out.append(f"G1 X{b_[0]:.3f} Y{b_[1]:.3f} Z{b_[2]:.3f} E{e_:.5f} F{f_str}")
            e = float(e_run[-1])
            cur = pts[-1].copy()
    out.append(f"G1 Z{(cur[2] if cur is not None else 0) + 5:.3f} F{p.travel_speed * 60:g}")
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
