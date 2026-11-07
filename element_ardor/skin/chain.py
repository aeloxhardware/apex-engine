"""
Enchainement du lissage non-planaire sur la voie chaleur (P2 -> P4 -> P5 en boucle de rognage, puis P3 et P5 exhaustif).

Porte de `pipeline/bird/skin/chain.py` (branche `bird_old`, decision de l'utilisateur du 2026-09-28) :

  1. couches du corps NON creuse : etape 6 (remplissage) sur les couches de l'etape 5 ;
  2. P2 : accessibilite de chaque face de la surface (cout : une fois) ;
  3. boucle : partition en zones (P2), trajectoires de surface (P4), ordre et collisions de la VRAIE tete contre le
     corps NON creuse, un point tous les `loop_sampling` mm (P5) ; les points de contact sont cumules et P2 retire les
     faces a moins de `trim_radius`, jusqu'a zero collision ;
  4. UN creusement (P3) sous les zones restantes, puis P5 exhaustif (tous les points) sur le corps creuse ; les zones
     encore en contact sont retirees et leur creusement annule (au plus `exclusion_rounds` tours).

Derive de S4 Slicer (Joshua Bird, GPL-3.0, https://github.com/jyjblrd/S4_Slicer).
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field, replace
from typing import Optional

import numpy as np

from core.mesh_types import Layer, MachineProfile, ScalarField, TetMesh
from element_ardor import step6_toolpath2d as s6
from element_ardor.skin import carve, paths as skin_paths, print_order, zones


@dataclass
class SkinParams:
    thickness: float = 0.2                 # epaisseur de la couche de lissage (mm)
    pattern: str = "rectilinear"           # motif des trajectoires de surface (surface_paths.PATTERNS)
    head_tolerance_deg: float = 45.0       # ecart maximal tete / normale
    trim_radius: float = 2.0               # rayon de rognage autour d'un point de collision (mm)
    max_iter: int = 8                      # tours de la boucle de rognage
    loop_sampling: float = 2.0             # un point teste tous les N mm pendant la boucle (le P5 final teste tous)
    exclusion_rounds: int = 3
    plate_radius: Optional[float] = None   # defaut : rayon du plateau de la machine, sinon 150 mm
    face_sample: Optional[int] = None      # n'evaluer qu'un echantillon de faces (mesure) ; None = toutes
    head_sample_spacing: float = 5.0       # echantillonnage de la tete pour l'accessibilite P2 (mm) ; P5 verifie tout apres
    workers: Optional[int] = None          # processus du calcul d'accessibilite (None = coeurs, 8 au plus)


@dataclass
class SkinResult:
    layers: list[Layer]                    # couches du corps CREUSE (contour_3d recalcule), toolpath_2d rempli
    body: list[dict]                       # idem, au format de P5 : {"isovalue", "toolpath_3d"}
    skin: dict                             # partition de P2
    paths: dict                            # trajectoires de surface (P4)
    print_order: dict                      # P5 : passes (corps + lissage), rapport, zones retirees
    iterations: int
    clean: bool
    excluded_zones: list[int]
    log: list[str] = field(default_factory=list)
    timings: dict = field(default_factory=dict)


def body_from_layers(layers: list[Layer]) -> list[dict]:
    """Couches remplies (etape 6) au format du controle P5 : trajectoires 3D du corps, par couche."""
    out = []
    for layer in layers:
        lines = [s6._lift_to_3d(np.asarray(p, dtype=float), layer.local_frame) for p in layer.toolpath_2d if len(p)]
        out.append({"isovalue": float(layer.isovalue), "toolpath_3d": lines})
    return out


@dataclass(frozen=True)
class SkinInputs:
    """Donnees des etapes amont consommees par `run`."""
    tet_mesh: TetMesh
    scalar_field: ScalarField
    layers: list[Layer]
    mesh_data: object
    overhang_map: object
    toolpath2d_params: s6.Step6Params


class _SkinRun:
    """Contexte d'une execution de `run` : entrees, journal et chronometrage."""

    def __init__(self, tet_mesh, scalar_field, layers, mesh_data, overhang_map, toolpath2d_params,
                 params, log: bool):
        self.tet_mesh, self.scalar_field, self.layers = tet_mesh, scalar_field, layers
        self.mesh_data, self.overhang_map, self.toolpath2d_params = mesh_data, overhang_map, toolpath2d_params
        self.params, self.log = params, log
        self.messages: list[str] = []
        self.timings: dict = {}

    def timed(self, name, fn, *args, **kwargs):
        t = time.time()
        try:
            return fn(*args, **kwargs)
        finally:
            self.timings[name] = self.timings.get(name, 0.0) + time.time() - t

    def say(self, msg):
        self.messages.append(msg)
        if self.log:
            print(f"[skin] {msg}", flush=True)

    def carved_body(self, skin, excluded=()):
        params = self.params
        cp = carve.CarveParams(thickness=params.thickness)
        keep = (skin["zone"] >= 0) & ~np.isin(skin["zone"], list(excluded))
        zone_faces = skin["faces"][keep]
        carved, _rep = self.timed("P3 creusement", carve.carve_layers, self.layers, self.tet_mesh, self.scalar_field,
                                  carve.ZoneDistance(skin["vertices"], zone_faces, cp), cp)
        filled_c = self.timed("remplissage du corps creuse (etape 6)", s6.run, carved, self.mesh_data,
                              self.overhang_map, self.toolpath2d_params)
        return filled_c, body_from_layers(filled_c)


def run(
    inputs: SkinInputs,
    machine: MachineProfile,
    sim,
    params: SkinParams = SkinParams(),
    log: bool = False,
) -> SkinResult:
    ctx = _SkinRun(inputs.tet_mesh, inputs.scalar_field, inputs.layers, inputs.mesh_data, inputs.overhang_map,
                   inputs.toolpath2d_params, params, log)
    timed, say = ctx.timed, ctx.say
    tet_mesh, scalar_field, layers = inputs.tet_mesh, inputs.scalar_field, inputs.layers
    mesh_data, overhang_map, toolpath2d_params = inputs.mesh_data, inputs.overhang_map, inputs.toolpath2d_params
    plate_r = params.plate_radius or getattr(machine, "platform_radius", None)
    if plate_r is None:
        raise ValueError("lissage de surface : ni SkinParams.plate_radius ni MachineProfile.platform_radius")
    offset = np.zeros(3)
    t0 = time.time()
    filled = timed("remplissage du corps (etape 6)", s6.run, layers, mesh_data, overhang_map, toolpath2d_params)
    body_uncarved = body_from_layers(filled)
    isovalues = np.sort(np.unique([l.isovalue for l in layers]))
    z = timed("P2 accessibilite", zones.compute, tet_mesh, scalar_field, isovalues, machine, sim,
              zones.SkinZoneParams(head_tolerance_deg=params.head_tolerance_deg,
                                   head_sample_spacing=params.head_sample_spacing),
              zones.ZoneSampling(face_sample=params.face_sample, surface=(mesh_data.vertices, mesh_data.faces),
                                 workers=params.workers))
    loop_sim = replace(sim, path_sampling_step=params.loop_sampling)
    full_sim = replace(sim, path_sampling_step=0.0)
    trim: list = []
    trace_cache: dict = {}
    clean, it = False, 0
    skin = traced = p5 = None
    for it in range(1, params.max_iter + 1):
        skin = timed("P2 partition", zones.partition, z, trim_points=trim, trim_radius=params.trim_radius)
        traced = timed("P4 trajectoires", skin_paths.trace, skin,
                       skin_paths.TraceOptions(pattern=params.pattern, thickness=params.thickness),
                       cache=trace_cache)
        p5 = timed("P5 controle (boucle)", print_order.order_and_check, body_uncarved, skin, traced, machine, loop_sim,
                   offset,
                   print_order.CheckOptions(tolerance=params.head_tolerance_deg, plate_r=plate_r))
        new = p5["collision_points"]
        say(f"tour {it} : {len(new)} point(s) de collision, {len(p5['excluded'])} zone(s) en contact")
        if not new:
            clean = True
            break
        trim.extend(new)

    excluded: set[int] = set()
    if clean:
        body_layers, body = ctx.carved_body(skin)
        p5 = timed("P5 controle (exhaustif)", print_order.order_and_check, body, skin, traced, machine, full_sim,
                   offset,
                   print_order.CheckOptions(tolerance=params.head_tolerance_deg, plate_r=plate_r))
        if not p5["collision_points"]:
            say(f"zero collision (verifie sur le corps creuse) en {time.time() - t0:.0f} s")
            return SkinResult(body_layers, body, skin, traced, p5, it, True, [], ctx.messages, ctx.timings)
        say(f"verification sur le corps creuse : {len(p5['collision_points'])} point(s) -> zones retirees")
    else:
        say(f"{params.max_iter} tours sans zero collision : zones restantes retirees")
    excluded |= set(int(zz) for zz in p5["excluded"])
    body_layers, body = None, None
    for r in range(1, params.exclusion_rounds + 1):
        body_layers, body = ctx.carved_body(skin, sorted(excluded))
        p5 = timed("P5 controle (exhaustif)", print_order.order_and_check, body, skin, traced, machine, full_sim,
                   offset,
                   print_order.CheckOptions(tolerance=params.head_tolerance_deg, plate_r=plate_r,
                                            exclude=tuple(sorted(excluded))))
        say(f"exclusion {r} : zones retirees {sorted(excluded)}, points restants {len(p5['collision_points'])}")
        if not p5["collision_points"]:
            break
        excluded |= set(int(zz) for zz in p5["excluded"])
    return SkinResult(body_layers, body, skin, traced, p5, it, clean, sorted(excluded), ctx.messages, ctx.timings)
