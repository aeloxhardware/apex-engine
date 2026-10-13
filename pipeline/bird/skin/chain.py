"""
Enchainement du lissage non planaire (porte de research/S4_Slicer/skin_loop.py
et de la partie --skin / --first-layer de bird_slice_deformed.py), EN MEMOIRE.

Boucle de rognage local (decision de l'utilisateur du 2026-09-28) :
P2 (zones) -> P4 (trajectoires) -> P5 (ordre + collisions de la VRAIE tete,
contre le corps NON creuse, un point tous les 2 mm), repetee tant que la tete
touche ; les points de contact sont cumules et P2 retire les faces a moins de
`trim_radius`. Puis UN SEUL creusement P3, verifie par un P5 exhaustif sur le
corps creuse ; les zones encore en contact sont retirees et leur creusement
annule (exclusion cumulee, 3 tours au plus).

Derive de S4 Slicer (Joshua Bird, GPL-3.0, https://github.com/jyjblrd/S4_Slicer).
"""
from __future__ import annotations

import os
import tempfile
import time
from dataclasses import dataclass, field, replace
from typing import Optional

import numpy as np

from pipeline.bird import deformed_slicing
from pipeline.bird.skin import paths as skin_paths
from pipeline.bird.skin import print_order, skin_carve, zones
from pipeline.step6_toolpath2d import _project_to_local_2d


@dataclass
class SkinParams:
    first_layer_pattern: Optional[str] = "concentric_conform"   # out/run_bird_chain.sh
    pattern: str = "rectilinear"
    skin_thickness: float = 0.2
    head_tolerance_deg: float = 45.0
    trim_radius: float = 2.0
    max_iter: int = 8                  # 8 tours (2026-10-01) : avec 3, le dino perdait le sommet de sa tete
    loop_sampling: float = 2.0         # P5 pendant la boucle (le P5 final teste tous les points)
    plate_radius: float = 150.0
    exclusion_rounds: int = 3
    carve_workers: Optional[int] = None


@dataclass
class SkinResult:
    body: list[dict]                   # couches du corps CREUSE (deformed_slicing), ordre d'impression
    skin: dict                         # zones.compute
    paths: dict                        # paths.trace
    print_order: dict                  # print_order.order_and_check (passes corps + lissage)
    iterations: int
    clean: bool
    excluded_zones: list[int]
    log: list[str] = field(default_factory=list)
    timings: dict = field(default_factory=dict)     # secondes cumulees par phase


def carve_hook(skin: dict, cells, deformed, original, thickness: float, excluded=(), workers: Optional[int] = None):
    """Creusement P3 sous les zones de lissage (bird_slice_deformed --skin) :
    rend une fonction couches -> couches pour `deformed_slicing.run(carve=...)`."""
    skin = dict(skin)
    if excluded:
        skin["zone"] = np.where(np.isin(skin["zone"], list(excluded)), -1, skin["zone"])

    def carve(layers):
        psi = skin_carve.skin_psi(original, skin, thickness)
        isovalues = [L.isovalue for L in layers]
        if len(set(isovalues)) != len(isovalues):
            raise ValueError("plusieurs couches par isovaleur : creusement non gere dans ce cas")
        todo = [i for i, L in enumerate(layers) if skin_carve.layer_is_touched(cells, deformed, psi, L.isovalue)]
        n_workers = workers if workers is not None else max(1, min(10, (os.cpu_count() or 2) - 2))
        if n_workers <= 1 or len(todo) < 2:
            exact = skin_carve.ExactDistance(skin, cells, deformed, original, thickness)
            for i in todo:
                region = skin_carve.carved_region(cells, deformed, psi, layers[i].isovalue, exact=exact, step=thickness)
                layers[i].contour_3d = skin_carve.region_loops(region, layers[i].isovalue)
        else:
            from concurrent.futures import ProcessPoolExecutor
            skin_d = {key: np.asarray(value) for key, value in skin.items()}
            with ProcessPoolExecutor(n_workers, initializer=skin_carve.carve_worker_init,
                                     initargs=(skin_d, cells, deformed, original, psi, thickness)) as pool:
                iso = [layers[i].isovalue for i in todo]
                for i, loops in zip(todo, pool.map(skin_carve.carve_worker_layer, iso,
                                                   chunksize=max(1, len(todo) // (4 * n_workers)))):
                    layers[i].contour_3d = loops
        return layers
    return carve


def first_layer_hook(pattern: str, toolpath2d_params):
    """Premiere couche pleine avec un motif (bird_slice_deformed --first-layer)."""
    def first_layer(layers):
        first = min(layers, key=lambda L: L.isovalue)
        loops = [_project_to_local_2d(np.asarray(c), first.local_frame) for c in first.contour_3d if len(c) >= 4]
        paths_ = skin_carve.first_layer_toolpath(loops, pattern, toolpath2d_params.nozzle_width,
                                                 toolpath2d_params.wall_count, toolpath2d_params.degenerate_area_ratio)
        first.toolpath_2d = [np.column_stack([p, np.zeros(len(p))]) for p in paths_]
        return layers
    return first_layer


def run(cells, original, deformed, offset, machine, sim, slicing_params: deformed_slicing.DeformedSlicingParams,
        params: SkinParams = SkinParams(), log: bool = False) -> SkinResult:
    messages: list[str] = []
    timings: dict = {}

    def timed(name, fn, *args, **kwargs):
        t = time.time()
        try:
            return fn(*args, **kwargs)
        finally:
            timings[name] = timings.get(name, 0.0) + time.time() - t

    def say(msg):
        messages.append(msg)
        if log:
            print(f"[skin] {msg}", flush=True)

    first = (first_layer_hook(params.first_layer_pattern, slicing_params.toolpath2d)
             if params.first_layer_pattern else None)

    def slice_body(skin=None, excluded=()):
        carve = (carve_hook(skin, cells, deformed, original, params.skin_thickness, excluded, params.carve_workers)
                 if skin is not None else None)
        name = "tranchage creuse (P3)" if skin is not None else "tranchage non creuse"
        return timed(name, deformed_slicing.run, cells, original, deformed, offset, slicing_params, carve=carve,
                     first_layer=first).layers

    t0 = time.time()
    body_uncarved = slice_body()
    isovalues = [L["isovalue"] for L in body_uncarved]
    loop_sim = replace(sim, path_sampling_step=params.loop_sampling)
    full_sim = replace(sim, path_sampling_step=0.0)
    trim: list = []
    trace_cache: dict = {}
    clean, it = False, 0
    with tempfile.TemporaryDirectory() as tmp:
        cache = os.path.join(tmp, "access_cache.npz")
        for it in range(1, params.max_iter + 1):
            skin = timed("P2 zones", zones.compute, cells, original, deformed, isovalues, machine, sim,
                                 head_tolerance=params.head_tolerance_deg, trim_points=trim,
                                 trim_radius=params.trim_radius, cache_path=cache)
            traced = timed("P4 trajectoires", skin_paths.trace, skin, pattern=params.pattern, thickness=params.skin_thickness,
                                      cache=trace_cache)
            p5 = timed("P5 controle (boucle)", print_order.order_and_check, body_uncarved, skin, traced, machine,
                       loop_sim, offset,
                                             tolerance=params.head_tolerance_deg, plate_r=params.plate_radius)
            new = p5["collision_points"]
            say(f"tour {it} : {len(new)} point(s) de collision, {len(p5['excluded'])} zone(s) en contact ; "
                f"P4 cache {trace_cache.get('_hits', 0)} reprise(s) / {trace_cache.get('_misses', 0)} calcul(s)")
            if not new:
                clean = True
                break
            trim.extend(new)

    excluded: set[int] = set()
    if clean:
        body = slice_body(skin)
        p5 = timed("P5 controle (exhaustif)", print_order.order_and_check, body, skin, traced, machine, full_sim,
                   offset, tolerance=params.head_tolerance_deg, plate_r=params.plate_radius)
        if not p5["collision_points"]:
            say(f"zero collision (verifie sur le corps creuse) en {time.time() - t0:.0f} s")
            return SkinResult(body, skin, traced, p5, it, True, [], messages, timings)
        say(f"verification sur le corps creuse : {len(p5['collision_points'])} point(s) -> zones retirees")
    else:
        say(f"{params.max_iter} tours sans zero collision : zones restantes retirees")
    excluded |= set(int(z) for z in p5["excluded"])
    for r in range(1, params.exclusion_rounds + 1):
        body = slice_body(skin, sorted(excluded))
        p5 = timed("P5 controle (exhaustif)", print_order.order_and_check, body, skin, traced, machine, full_sim,
                   offset, tolerance=params.head_tolerance_deg, plate_r=params.plate_radius,
                   exclude=sorted(excluded))
        say(f"exclusion {r} : zones retirees {sorted(excluded)}, points restants {len(p5['collision_points'])}")
        if not p5["collision_points"]:
            break
        excluded |= set(int(z) for z in p5["excluded"])
    return SkinResult(body, skin, traced, p5, it, clean, sorted(excluded), messages, timings)
