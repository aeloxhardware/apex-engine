# Plan: "Color by" selector for the G-code Preview legend

## 0. Scope decision - v1 vs deferred (15 criteria)

Grouped by real parsing/rendering cost, given the exact current-state facts (no temp/fan
parsing, no per-move numeric fields survive parseGCode(), no gradient utility exists anywhere in
the repo):

**V1 - ship now** (data either already available or cheap to add, no cross-dependencies):
1. Type de ligne (existing, default) - categorical, zero new work.
2. Vitesse (feed rate F, mm/s) - feeds[] already computed per push(); trivial to expose.
9. Acceleration - accels[] already computed per push(); trivial to expose (reuses existing modal
   tracking, just stop discarding it).
11. Inclinaison de la tete / axe B - GCodeMove.b already stored, zero parsing work.
14. Longueur du segment - trivially derivable from consecutive move coordinates (already computed
    as len[] inside planCartesianTimes for cartesian; needs an equivalent per-move chord/arc-length
    for polar, both cheap).
8. Temps ecoule - times[] (cumulative) already returned; purely a legend/colormap consumer, zero
   parser change.
15. Pression d'extrusion relative (E/mm) - eDeltas[] and segment length already available; ratio
    is one division.
7. Debit volumetrique (mm3/s) - E-volume (from eDeltas + filament diameter, already parsed)
   divided by per-move dt (already available via times[] deltas); no new parsing, just a derived field.
12. Rapport vitesse / vitesse max de la couche - derived in a post-pass once per-move actual speed
    exists (needs a "current layer" grouping by Z, a new but cheap aggregation, no new parsing).
6. Largeur de ligne - needs new comment regex (;WIDTH:) parsing (cheap, same handleComment pattern
   as ;TYPE:) plus a fallback estimate formula (E_volume / (distance x layer_height)); layer_height
   itself must be computed as part of criterion 4, so this depends on 4 shipping alongside it.
4. Hauteur de couche - NOT currently computed; needs per-move Z-delta-at-layer-change tracking
   (new bookkeeping, moderate complexity but self-contained, no cross-file dependency). Ship in v1
   because #6 depends on it.

**V1, shipped last within the milestone (more parsing work, still self-contained):**
3. Temperature (buse ; plateau si elle varie) - requires wholly new M104/M109/M140/M190 parsing
   (new modal state nozzleTemp/bedTemp, updated like feed/accelPrint today).
5. Vitesse du ventilateur - requires wholly new M106/M107 parsing (new modal fanSpeed, S-word
   0-255 -> %, M107 = 0).

**Deferred to phase 2 (stubbed extension points only):**
10. Vitesse angulaire de l'axe C - Apex-specific; derivable as dC/dt from consecutive polar moves,
    but has no cartesian-file analogue at all (would show as "n/a" for the vast majority of files
    users load, which are cartesian slicer output) - low value for v1 effort, defer.
13. Score de risque de defaut - explicitly the most complex: needs speed (#2), overhang detection
    (LineTypeId 'overhangWall' already exists as a proxy) and cooling (#5, fan) all landed AND
    validated first, plus a heuristic weighting formula with no existing precedent in this codebase
    to anchor design decisions on. Per the user's own reasoning this must land last, after the
    others prove out the plumbing (numeric per-move data -> gradient legend). Ship as an explicit
    stub (menu entry visible but marked "Bientot disponible", or omitted entirely from the v1 menu
    - recommend omitting, since a half-working risk score is worse than an honest absence) in phase 2.

**Recommendation flagged for user confirmation:** the above grouping is my judgment call, not
something the user specified numerically - the user only said "some now, some deferred, plan
should propose which." Please confirm before implementation, especially whether temperature/fan
(#3/#5) should really be v1 (non-trivial new parsing) vs. phase 2.

## 1. Data model extension (src/components/gcodeLayer.ts)

### 1.1 New parallel arrays inside parseGCode(), alongside existing feeds[]/eDeltas[]/etc.

All arrays are indexed identically to moves[]/times[] (index 0 unused/zero, matching existing
convention), pushed once per call to push() so arc subdivision (push() called per subdivided arc
segment in the G2/G3 handler) automatically produces correct per-segment values - no separate
handling needed for arcs, exactly as feeds/eDeltas/accels already work today.

New modal state (top of parseGCode, alongside feed, accelPrint, etc.):

    let nozzleTemp = 0;       // degC, from M104/M109 S
    let bedTemp = 0;          // degC, from M140/M190 S
    let fanSpeed = 0;         // 0-100 %, from M106 S (0-255 scaled) / M107 (=0)

New M-code cases in the "if (cmd[1] === 'M')" block:

    else if (num === 104 || num === 109) { if (w.S >= 0) nozzleTemp = w.S; }
    else if (num === 140 || num === 190) { if (w.S >= 0) bedTemp = w.S; }
    else if (num === 106) { if (w.S >= 0) fanSpeed = Math.round((w.S / 255) * 100); }
    else if (num === 107) { fanSpeed = 0; }

(M109/M190 are "wait" variants; for a static-preview colorizer we only care about the commanded
setpoint, not the wait - so treating them identically to M104/M140 is correct and simpler than
tracking "reached" state.)

New parallel arrays populated inside push() (mirroring feeds.push(feed)):

    temps.push(nozzleTemp);
    bedTemps.push(bedTemp);
    fans.push(fanSpeed);

### 1.2 Per-move speed (actual, not just commanded F)

feeds[] already holds the commanded F (mm/min) per move - sufficient for "Vitesse" in v1 as a
first pass (simpler, avoids touching planCartesianTimes). Flag for user: the *actual*
trapezoid-limited speed (vp inside planCartesianTimes) would be more accurate but requires
planCartesianTimes to also emit a speeds: Float64Array (extra return value, touches a
performance-sensitive function). Recommendation: ship v1 with commanded-F coloring (cheap, correct
for the vast majority of moves which do reach cruise speed) and note actual-speed as a phase-2
refinement - flag this trade-off to the user.

### 1.3 Layer height (#4) and line width (#6)

New bookkeeping in parseGCode, cartesian path only (polar has no discrete "layers" in the same
sense - treat layer-height/width as "n/a" for polar in v1, another judgment call to flag):
- Track layerStartZ / prevLayerZ: on each extruding move, if mv.z differs from the last extruding
  Z seen by more than 1e-4, that's a layer change; layerHeight = mv.z - prevLayerZ. Store into a
  new layerHeights: number[] parallel array (repeats the same value for every move within a layer -
  simplifies downstream gradient code, no need to look up "current layer" separately).
- New comment regex WIDTH_RE = /^\s*WIDTH\s*[:=]\s*([\d.]+)/i parsed in handleComment, alongside
  FEATURE_RE, setting a modal explicitWidth (mm) applied to subsequent extruding moves; reset to
  NaN if not re-stated (some slicers emit per-line ;WIDTH: comments before each new width region -
  treat as modal, matching how ;TYPE: is handled today).
- Fallback estimate when explicitWidth is NaN: width = eVolumeMm3 / (segmentLengthMm *
  layerHeightMm), computed per move from already-available deposits[k] (mm of filament) converted
  to volume via the existing filamentDiameter, plus the new layerHeights[k] and segment length
  (see 1.4). Store into widths: number[].

### 1.4 Segment length (#14) and volumetric flow (#7) and E/mm (#15)

- segLengths: number[] - chord length for cartesian (same formula as len[] in planCartesianTimes,
  but that function is cartesian-only and internal; add a small shared helper
  moveLength(a: GCodeMove, b: GCodeMove, kind: GCodeKind): number used both by planCartesianTimes
  (replacing its local len[] computation, so there's exactly one implementation) and by the new
  per-move array population after time-planning completes. For polar, reuse the existing
  nozzle-tip-arc-length formula already inlined in the polar times[] loop (lines ~586-597) -
  extract it into the same shared helper to avoid duplicating the arc-length-in-RCZ-space math.
- flowMm3s and eRatio are NOT stored as separate arrays - compute lazily in the legend/colormap
  layer from deposits[]/eDeltas[], segLengths[], and per-move dt (times[k]-times[k-1]), since
  they're pure derived ratios with no independent meaning; storing them separately would be
  redundant memory for large files. (Judgment call - flagging: this keeps ParsedGCode smaller at
  the cost of a few extra multiplications per mode-switch rebuild, negligible next to mesh rebuild
  cost.)

### 1.5 Relative speed within layer (#12)

Needs a per-layer max-speed lookup. Two-pass: first pass (right after times[]/feeds[] exist)
builds layerMaxFeed: Map<number /* layer index */, number> keyed by the same layer-change
detection as 1.3; second pass computes speedRatio[k] = feeds[k] / layerMaxFeed.get(layerOf(k)).
Store the per-move layerIndex: Uint32Array (small, reusable for any future per-layer aggregation)
alongside layerHeights.

### 1.6 GCodeMove vs. parallel arrays - recommendation

Keep GCodeMove itself unchanged ({x,y,z,b,extruding,type}) and return all new per-move numeric
data as new top-level Float32Array fields on ParsedGCode, not on GCodeMove, mirroring the existing
times: Float64Array pattern rather than the discarded local-array pattern. Rationale: GCodeMove[]
is a regular JS array of objects (one allocation per move); adding numeric fields there would
bloat per-object memory for every mode even when unused. Typed arrays returned once on ParsedGCode
are cheap, contiguous, and match the existing times convention exactly:

    export interface ParsedGCode {
      ...
      // Nouveaux champs numeriques par mouvement (Float32Array, un par critere continu) :
      feedMmS: Float32Array;      // vitesse commandee F, mm/s
      accelMmS2: Float32Array;    // acceleration modale au mouvement
      nozzleTempC: Float32Array;
      bedTempC: Float32Array;
      fanPct: Float32Array;
      layerHeightMm: Float32Array;
      widthMm: Float32Array;
      segLengthMm: Float32Array;
      layerIndex: Uint32Array;
      layerMaxFeedMmS: Float32Array;   // indexed by layerIndex, small lookup table, not per-move
      // b (inclinaison tete) already on GCodeMove; angular C speed deferred to phase 2.
      ranges: ColorByRanges;  // see section 3
    }

### 1.7 Per-file min/max ranges for gradient normalization (new aggregation pattern)

There is currently no continuous-value aggregation anywhere (lineStats is the only aggregation,
and it's discrete-bucket). Add a single pass after all arrays are built:

    export interface ColorByRange { min: number; max: number; avg: number }
    export type ColorByRanges = Record<ContinuousColorByMode, ColorByRange>;

computed once in parseGCode (cheap - one linear scan over already-built arrays, negligible next to
the existing O(n) passes already in the function) and returned on ParsedGCode. This is the "new
histogram/min-max aggregation" the requirements call out as currently absent.

## 2. Colormap utility - new file src/components/gcodeColormap.ts

Small, dependency-free module:

    export function normalize(value: number, range: ColorByRange): number; // clamp((v-min)/(max-min), 0, 1)
    export function viridisLike(t: number, out: BABYLON.Color3): BABYLON.Color3; // 4-stop RGB lerp
    export function colorForMode(mode: ColorByMode, k: number, parsed: ParsedGCode, out: BABYLON.Color3): BABYLON.Color3;

Recommendation: a 4-stop piecewise-linear RGB lerp (dark blue -> cyan -> yellow -> red), i.e. a
simplified viridis/turbo-style ramp, computed with plain arithmetic (no HSV conversion needed,
avoids hue-wraparound bugs, cheap enough for per-vertex use at parse/build time). Reuse
BABYLON.Color3.Lerp (confirmed unused today but perfectly fit for this) between the two bracketing
stops of t. This one function serves every continuous mode (#2,#3,#4,#5,#6,#7,#8,#9,#12,#14,#15) -
categorical mode (#1) keeps using LINE_TYPES[].color unchanged, never routed through this module.

## 3. Rendering strategy - decision

Hybrid, keyed by categorical vs. continuous, as the requirements doc's own analysis suggests:

- Categorical mode (#1, "Type de ligne"): keep the existing architecture completely unchanged -
  one mesh per LineTypeId, one flat StandardMaterial, setEnabled for visibility. Zero risk, zero
  regression surface, since this is the default/first mode and must remain exactly as-is per
  requirements.

- Continuous modes (all others in v1): switch to per-vertex vertex colors on ONE merged bead mesh
  + ONE merged travel/wipe line mesh, replacing the per-LineTypeId mesh split for the duration that
  a continuous mode is active. Justification:
  - Bucketing into N discrete bands (the alternative) would still require picking a band count,
    causes visible banding contrary to the "fluidity"/visual-quality bar for a gradient legend, and
    doesn't actually save implementation complexity here - buildBeadMesh/buildPathMesh already
    build one mesh per include() predicate, so a bucketed version needs the same amount of new code
    as a vertex-color version but produces a visibly worse result.
  - Per-vertex color avoids the "one mesh per category" ceiling (13 modes x potentially many
    buckets each would multiply draw calls); a single merged mesh per continuous mode is actually a
    net draw-call reduction vs. today's ~11 extrusion-type meshes, which helps fluidity, not hurts it.
  - Concrete Babylon changes: extend buildBeadMesh/buildPathMesh with an optional
    vertexColor: (k: number) => BABYLON.Color3 callback (mutually exclusive with the existing
    material/color param - when provided, write VertexData.colors (Float32Array, 4 floats/vertex,
    alpha 1) computed once per vertex written in the existing ring/position-writing loop, no extra
    pass). Material becomes a single shared "vertex-color aware" StandardMaterial: recommend
    mat.disableLighting = true; mat.emissiveColor = Color3.White() so vertex colors multiply
    emissive directly, avoiding any lighting-interaction surprises (verify mat.diffuseColor stays
    white if disableLighting is left false instead, so it doesn't tint the vertex colors); keep
    backFaceCulling = false as today.
  - buildPathMesh (LinesMesh) needs per-vertex color support too - Babylon's LinesMesh supports
    per-vertex colors natively via VertexData.colors plus mesh.useVertexColor = true - set that
    flag and write colors alongside positions in the existing vertex-writing loop.

- Mode switching / rebuild cost: switching colorBy requires a geometry rebuild only when crossing
  the categorical<->continuous boundary or between two different continuous modes (new vertex
  colors need different values). To respect the fluidity requirement (<=0.1mm/zero stutter):
  - Rebuild happens once, synchronously, on mode change - NOT per frame. This matches today's
    load() cost model (already a full rebuild on file load); a mode switch is strictly cheaper (no
    re-parsing, geometry/topology identical, only the color buffer changes).
  - Key optimization: cache the merged-mesh positions/normals/indices per kind (extrusion merged
    bead geometry, travel+wipe merged line geometry) built once at load() time regardless of
    initial mode, and on a continuous-mode switch only recompute+upload the VertexData.colors
    buffer via mesh.updateVerticesData(BABYLON.VertexBuffer.ColorKind, ...) - no
    VertexData.applyToMesh (positions/indices unchanged), avoiding the expensive parts of a full
    rebuild. This is the single most important perf decision for the "no stutter" requirement.
  - Only the very first switch away from categorical (or into it) needs a topology change (merging
    N per-type meshes into 1, or vice versa); cache BOTH geometry layouts (per-type AND merged)
    built once at load(), so subsequent toggling between categorical and any continuous mode is
    also just a setEnabled/color-buffer swap, never a rebuild. Given typical file sizes this
    roughly doubles peak geometry memory for the extrusion mesh set - acceptable trade for
    zero-stutter switching; flag this memory/latency trade-off to the user if very large files are
    a concern.

- Live trailing segment (updateLive, liveBead/liveLine): currently picks
  beadMaterials.get(m.type) / colors[type]. Extend to: in continuous mode, compute
  colorForMode(mode, next, parsed, tmpColor) once per frame (cheap, O(1), no allocation if
  tmpColor is a hoisted reusable BABYLON.Color3) and set liveBead.material.emissiveColor = tmpColor
  directly (since the live bead is a single small mesh, not vertex-colored, a flat material color
  update is fine and cheapest). Recommendation: flat live-color per frame - the live segment is
  tiny (<=1 move worth of geometry) so any visible mismatch with a gradient is imperceptible, and
  this avoids an updateVerticesData call every frame (better for fluidity than vertex-coloring the
  live strip too).

## 4. Event/state wiring

New CustomEvent, following apex-gcode-line-visibility's exact pattern (grepped - no existing
apex-* event carries a color-mode payload):

    // use-gcode-preview.ts
    export type ColorByMode = 'lineType' | 'speed' | 'accel' | 'headAngle' | 'segmentLength'
      | 'elapsedTime' | 'extrusionRatio' | 'volumetricFlow' | 'layerSpeedRatio'
      | 'lineWidth' | 'layerHeight' | 'temperature' | 'fanSpeed'; // v1 set per section 0

    const [colorBy, setColorByState] = useState<ColorByMode>('lineType');
    useEffect(() => {
      window.dispatchEvent(new CustomEvent('apex-gcode-color-by', { detail: { mode: colorBy } }));
    }, [colorBy]);
    const setColorBy = useCallback((mode: ColorByMode) => setColorByState(mode), []);

Viewport3D.tsx adds onGCodeColorBy mirroring onGCodeLineVisibility (add/remove listener next to
the existing ones, call a new gcodeLayer.setColorBy(mode) method on the GCodeLayer interface).

GCodeFileInfo (in use-gcode-preview.ts) gains the per-mode range data needed for legend display:

    colorByRanges: ColorByRanges; // from ParsedGCode.ranges, section 1.7

threaded through in loadText()'s setFileInfo({...}) call alongside the existing fields.

## 5. UI

Selector placement: reuse previewUpperContent - put the mode selector as a new small header
control directly above GCodeLineTypeLegend's existing header row. Recommendation: extend
GCodeLineTypeLegend in place (rather than splitting into a separate wrapper component) - it
already owns the legend table; adding the mode dropdown as its own header row above the grid keeps
the "swap legend content per mode" logic co-located with the legend body it swaps, per requirement
#5's own framing ("how GCodeLineTypeLegend's table body/columns change per mode").

Widget shape for 13 (v1) options: per the user's confirmed constraint (too many for a flat
segmented control), use a shadcn Popover + grouped list (Radix primitives already installed per
AI_RULES.md), trigger button styled like existing glassmorphism controls
(GCodeLineTypeLegend.module.scss / GCodeSpeedSelector.module.scss visual language - reuse the
checkboxRow/segment button treatment for the trigger's rest state). Group sections inside the
popover:
- "Type de ligne" (single item, always first, default)
- "Cinematique" - Vitesse, Acceleration, Inclinaison tete (B), Longueur segment, Temps ecoule
- "Extrusion" - Debit volumetrique, Pression E relative, Largeur de ligne, Hauteur de couche
- "Vitesse relative" - Rapport vitesse/vitesse max de couche
- "Thermique / refroidissement" - Temperature, Vitesse ventilateur
(Angular C speed / risk score omitted from v1 menu per section 0, added as a disabled/greyed
"Bientot disponible" group only if the user wants visible signposting - otherwise fully absent
until phase 2.)

Each popover row: label + tiny icon/swatch preview (for continuous modes, a small gradient bar
snippet using the same viridisLike ramp, for visual consistency with the eventual full legend
gradient bar).

Legend body change per mode:
- Categorical (lineType): unchanged - current grid/rows/checkboxes exactly as today.
- Continuous mode: replace the per-row grid with:
  - A single horizontal gradient bar (CSS linear-gradient, stops sampled from viridisLike at 5-6
    points, computed once per mode - hardcode the same stop colors as literal hex in the
    .module.scss/inline style so the CSS bar and the Babylon coloring are visually identical).
  - Min / Avg / Max readout row beneath it, sourced from fileInfo.colorByRanges[mode] (section
    1.7), formatted with the mode's unit (mm/s, degC, %, mm, s, mm3/s) - a new small
    formatColorByValue(mode, value) helper co-located in GCodeLineTypeLegend.tsx.
  - No per-row visibility checkboxes in continuous mode (there's no discrete category to toggle) -
    but the existing categorical visibility state (lineVisibility) must remain independently
    respected by the 3D layer regardless of colorBy (a continuous color mode still only shows
    line-types the user has toggled on) - confirm this stays orthogonal: setLineTypeVisibility and
    setColorBy are separate GCodeLayer methods, so visibility filtering (include() predicates) and
    coloring are already architecturally decoupled by this plan (section 3), good.

## 6. Extension points (16th+ criterion)

Given the chosen architecture, adding a new continuous criterion requires only:
1. A new parallel array (or derived-on-the-fly computation) added to ParsedGCode in gcodeLayer.ts.
2. One new entry in the ColorByMode union + one case in colorForMode() in gcodeColormap.ts.
3. One new popover row (with its group) in the selector.
4. One new unit/label mapping in formatColorByValue.
No changes to buildBeadMesh/buildPathMesh/updateLive/the event wiring are needed - they're already
generic over "a function that returns a color for move k." Adding a new categorical criterion
(rather than continuous) is more work (needs its own discrete-bucket mesh-per-value scheme like
LineTypeId today) - flag this asymmetry to the user; all 13 v1 criteria except #1 are continuous
by design, which is why the hybrid strategy optimizes for that case.

## 7. Verification plan

1. pnpm dev, open the Preview space.
2. Load a real polar Apex .gcode (B-word moves) - confirm: default "Type de ligne" mode renders
   identically to before (no regression); switching to "Vitesse"/"Acceleration"/"Inclinaison tete"
   recolors the toolpath with a visible gradient; angular-C-speed and risk-score are absent from
   the menu (phase 2); temperature/fan show "n/a" or flat color if the polar file has no
   M104/M106 lines (verify graceful fallback, not a crash, when a mode's data is all-zero/no
   variation).
3. Load a real cartesian slicer .gcode that includes temperature/fan/;WIDTH: comments (e.g. an
   Orca/Bambu export) - confirm temperature and fan modes show real variation (if the file
   actually changes temps/fan mid-print; many single-temperature prints won't show a gradient -
   verify the min=max case renders a solid single color, not divide-by-zero/NaN).
4. For each of the 13 v1 modes: confirm legend swaps from table-with-checkboxes to
   gradient-bar-with-min/avg/max, values look plausible (e.g. speed in mm/s roughly matches known
   print speeds; layer height matches the file's known layer height).
5. Confirm visibility checkboxes still work correctly in categorical mode (toggle a line type off,
   confirm it disappears) and that switching to a continuous mode and back to categorical restores
   the exact same visibility state (no accidental reset).
6. Confirm live playback trail recolors correctly while a continuous mode is active and during
   playback (press play, watch the trailing bead segment match the color of the static prefix it
   extends).
7. Switch rapidly between all modes with a large file loaded (many thousands of moves) and confirm
   no visible stutter/frame hitch (per the <=0.1mm fluidity requirement) - this specifically tests
   the "cache both geometry layouts, only swap color buffer" optimization in section 3; if a hitch
   IS visible, that is the first place to profile (likely updateVerticesData size or an accidental
   full VertexData.applyToMesh being triggered instead of the targeted buffer update).
8. Confirm mode selection behaves predictably across loading a new file in the same session -
   recommend resetting colorBy to 'lineType' on new file load, consistent with the existing
   setConversion(IDLE_CONVERSION) reset already done in loadText().

## Open judgment calls to flag to the user before implementation

- Whether temperature (#3) and fan (#5) really belong in v1 given they need net-new M-code parsing
  (vs. everything else in v1 which reuses already-tracked modal state).
- Commanded-F speed vs. actual trapezoid-limited speed for "Vitesse" mode (recommend commanded F
  for v1 simplicity).
- Layer height / width / relative-layer-speed as "n/a" for polar files (no discrete layer concept
  in the same sense) - confirm acceptable.
- Angular C speed and risk score fully deferred to phase 2, not even shown as disabled stub in the
  v1 menu - confirm acceptable, or user may prefer a visible "Bientot disponible" placeholder.
- Doubling peak geometry memory (caching both per-type and merged mesh layouts) to guarantee
  zero-rebuild mode switching - confirm acceptable given "fluidity above all" preference, especially
  for very large G-code files.

### Critical Files for Implementation
- R:/APEX_NEW_VERSION/slicer_ux/src/components/gcodeLayer.ts
- R:/APEX_NEW_VERSION/slicer_ux/src/components/gcodeColormap.ts (new)
- R:/APEX_NEW_VERSION/slicer_ux/src/hooks/use-gcode-preview.ts
- R:/APEX_NEW_VERSION/slicer_ux/src/components/GCodeLineTypeLegend.tsx
- R:/APEX_NEW_VERSION/slicer_ux/src/components/Viewport3D.tsx
