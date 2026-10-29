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


@dataclass
class SliceParams:
    layer_height: float = 0.2
    first_layer_height: float | None = None   # None : = layer_height
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
    layer_count: int
    thickness_min: float
    thickness_max: float
    slope_max_deg: float
    extruded_length: float
    filament_length: float
    offset: tuple[float, float, float]


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
    stack = ly.build(field_, p.layer_height, p.first_layer_height)
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

    for k in range(1, n + 1):
        out.append(f"; <layer> {k}")
        speed = p.first_layer_speed if k == 1 else p.print_speed
        path = tp.layer_polylines(field_.footprint, k, n, p.toolpath)
        for pl in path.polylines:
            pts = tp.densify(pl.xy, stack, k, p.toolpath.max_segment)
            if len(pts) < 2 or np.isnan(pts[:, 2]).any():
                continue
            # Dégagement : plus haut point de la couche entre départ et arrivée, plus z_hop.
            if cur is not None:
                line = np.linspace(cur[:2], pts[0, :2], 12)
                zc = float(np.nanmax(stack.z(k, line[:, 0], line[:, 1]))) + p.z_hop
            else:
                zc = float(pts[0, 2]) + p.z_hop
            emit_travel(pts[0], zc)
            for a, b in zip(pts[:-1], pts[1:]):
                seg = float(np.linalg.norm(b - a))
                if seg < 1e-9:
                    continue
                mid = (a[:2] + b[:2]) / 2
                t = float(stack.thickness(k, mid[0], mid[1]))
                thick_min, thick_max = min(thick_min, t), max(thick_max, t)
                vol = seg * float(_bead_area(np.array(t), w))
                e += p.flow * vol / fil_area
                total_len += seg
                out.append(f"G1 X{b[0]:.3f} Y{b[1]:.3f} Z{b[2]:.3f} E{e:.5f} F{speed * 60:g}")
            cur = pts[-1].copy()
    out.append(f"G1 Z{(cur[2] if cur is not None else 0) + 5:.3f} F{p.travel_speed * 60:g}")
    out.append("M104 S0\nM140 S0")
    report = SliceReport(
        layer_count=n, thickness_min=float(thick_min), thickness_max=float(thick_max),
        slope_max_deg=ly.max_slope_deg(field_, n, n, stack.first_layer), extruded_length=total_len,
        filament_length=e, offset=(float(off[0]), float(off[1]), float(off[2])),
    )
    return "\n".join(out) + "\n", report
