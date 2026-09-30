# Custom Notification Stack -- Implementation Plan

## 0. Decision: sonner stays installed, but fully unmounted from the app

Remove the `<Toaster>` mount and the `sonner` import from `src/App.tsx`. Do NOT
touch `src/components/ui/sonner.tsx` (CLAUDE.md forbids editing generated `ui/*`
files) and do NOT run `npm uninstall sonner` -- leave the dependency in
`package.json` untouched. Rationale: `ui/sonner.tsx` is shadcn-generated and
off-limits either way; removing the npm dependency is a separate, riskier,
unnecessary step (could affect lockfile / future shadcn regen) for zero
runtime benefit once nothing imports it. Net effect: sonner code is simply dead
weight in `node_modules`, nothing imports `sonner` or `ui/sonner.tsx` anymore.

## 1. New files

- `src/utils/notificationStore.ts` -- plain module-level store (no state lib),
  mirrors the module-state + subscriber pattern already used in this repo
  (e.g. how local hooks like `use-gcode-preview.ts` manage state without a
  store lib). Kept as a direct import between `toast.ts` and
  `NotificationStack.tsx` since they are tightly coupled 1:1 -- no need to
  round-trip through `window` events for this internal wiring; `window`
  events stay reserved for the existing cross-component signal bus per
  CLAUDE.md convention.
- `src/components/NotificationStack.tsx` -- the stack UI, mounted once in
  `App.tsx` in place of `<Toaster>`.
- `src/components/NotificationStack.module.scss` -- glass styling, following
  the `.bar` / `.toolbar` recipe from GCodePreviewBar / SliceSidebar.

`src/utils/toast.ts` is edited (not replaced) to call into the new store
instead of `sonner`, keeping exact existing exported signatures.

## 2. Store shape (`src/utils/notificationStore.ts`)

```ts
export type NotificationVariant = 'success' | 'error' | 'warning' | 'loading';

export interface NotificationItem {
  id: string;
  variant: NotificationVariant;
  message: string;
  createdAt: number; // Date.now(), used only for stable ordering
}

type Listener = (items: NotificationItem[], frontId: string | null) => void;

let items: NotificationItem[] = [];
// id of the notification currently pinned to front by the wheel gesture (or
// by having just been pushed). null = no override -> selector falls back to
// most recent item (see NotificationStack's orderedIds derivation).
let frontId: string | null = null;
const listeners = new Set<Listener>();

function emit() {
  const snapshot = items.slice();
  listeners.forEach((l) => l(snapshot, frontId));
}

let uid = 0;
function genId() {
  uid += 1;
  return `ntf-${Date.now()}-${uid}`;
}

export function pushNotification(variant: NotificationVariant, message: string): string {
  const id = genId();
  items = [...items, { id, variant, message, createdAt: Date.now() }];
  frontId = id; // a new notification always becomes the front card
  emit();
  return id;
}

/** Used by showLoading's returned id + dismissToast, and by the (x) button. */
export function removeNotification(id: string): void {
  const idx = items.findIndex((n) => n.id === id);
  if (idx === -1) return;
  items = items.filter((n) => n.id !== id);
  if (frontId === id) frontId = null; // stack re-derives front (see selector)
  emit();
}

/** Called by NotificationStack's wheel handler. */
export function setFrontNotification(id: string): void {
  if (!items.some((n) => n.id === id)) return;
  frontId = id;
  emit();
}

export function subscribeToNotifications(listener: Listener): () => void {
  listeners.add(listener);
  listener(items.slice(), frontId);
  return () => listeners.delete(listener);
}
```

Notes:
- The `loading` variant is kept only because `showLoading` / `dismissToast`
  must keep working (they currently have no call sites in the reviewed files,
  but are part of the public `toast.ts` API and must not regress). It renders
  like the others -- no spinner infra needed since these permanent cards have
  no timer to race against; a static "loading" tint on the icon is enough.
- Ordering is push order (`items`, oldest first, newest last). The stack
  always visually derives "who's on top" from `frontId ?? items[items.length-1]?.id`,
  so a fresh `pushNotification` always naturally becomes front.

## 3. `src/utils/toast.ts` (edited, same exported signatures)

```ts
import { pushNotification, removeNotification } from './notificationStore';

export const showSuccess = (message: string) => {
  pushNotification('success', message);
};

export const showError = (message: string) => {
  pushNotification('error', message);
};

export const showWarning = (message: string) => {
  pushNotification('warning', message);
};

export const showLoading = (message: string) => {
  return pushNotification('loading', message);
};

export const dismissToast = (toastId: string) => {
  removeNotification(toastId);
};
```

All ~25 existing call sites (AxisSidebar.tsx, use-gcode-preview.ts,
CategoriesSidebar.tsx, Index.tsx, ModelingSidebar.tsx,
PrinterConfigSidebar.tsx, TextureEditorModal.tsx, Viewport3D.tsx,
ViewCube.tsx) need ZERO changes -- they only import
`showSuccess`/`showError`/`showWarning`/`showLoading`/`dismissToast` from
`@/utils/toast`, whose signatures are unchanged.

## 4. `src/components/NotificationStack.tsx`

### Responsibilities
- Subscribe to `notificationStore` for `(items, frontId)`.
- Render the stack container, anchored bottom-right identically to the old
  sonner mount: same `right: 24px` / `bottom: 84px` offset documented in
  `App.tsx`'s old comment, reimplemented as fixed positioning in the SCSS
  module.
- Compute per-card stack position (`stackIndex`, 0 = front) and render only a
  capped number of "peek" cards behind the front one.
- Attach an `onWheel` handler on the stack container.
- Render nothing (`return null`) when `items.length === 0`.

### State shape inside the component

```tsx
const [items, setItems] = useState<NotificationItem[]>([]);
const [frontId, setFrontIdState] = useState<string | null>(null);

useEffect(() => subscribeToNotifications((its, front) => {
  setItems(its);
  setFrontIdState(front);
}), []);
```

### Deriving render order (front -> back)

```ts
const MAX_PEEK = 3; // front + up to 3 peeking behind = 4 cards max rendered

const orderedIds = useMemo(() => {
  if (items.length === 0) return [];
  const ids = items.map((n) => n.id);
  const resolvedFront = frontId && ids.includes(frontId) ? frontId : ids[ids.length - 1];
  // Everything else, newest-first, forms the rest of the cycle order.
  const rest = ids.filter((id) => id !== resolvedFront).reverse();
  return [resolvedFront, ...rest];
}, [items, frontId]);

const visibleIds = orderedIds.slice(0, MAX_PEEK + 1); // front + up to 3 peeking
```

Each id in `visibleIds` gets `stackIndex` = its position in that array (0 =
front, fully drawn/interactive; 1..3 = peeking slivers). DECISION: peeking
cards are not interactive at all -- only the front card has a live (x) and
receives wheel; to close a specific non-front notification the user must
first scroll it to front. This keeps hit-testing simple and matches the
"stacked deck, wheel brings one to front" mental model.

Cards beyond `MAX_PEEK+1` are not rendered (DOM stays small) but still exist
in `items` and become visible again as they cycle to a top-4 slot via wheel.

### Wheel exit direction -- user correction (overrides the generic spring-reflow described below)

The user specified precisely how the front-swap must look: **the current front
card slides out toward the bottom-right (off past the corner it's anchored
to) while the card behind it slides up into the front slot in the same
continuous motion** -- not a static in-place scale/fade crossfade. So the
wheel-driven transition is a directional slide-out + slide-up-into-place, not
a generic "recompute stackIndex, let spring interpolate" swap. This changes
the `animate`/`exit` targets in section 4's JSX and the framing of section 5
below: the card leaving front position gets an explicit `x`/`y` push
(down-right, positive values, away from the anchor corner) layered on top of
its normal fade-out, rather than just reusing the peek-position math for
"one slot further back."

### Per-card transform (the "peek" visual)

```ts
const PEEK_OFFSET_X = 10; // px per layer back, rightward
const PEEK_OFFSET_Y = -10; // px per layer back, upward (matches reference: top/right edges peek out)
const SCALE_STEP = 0.05;
const OPACITY_STEP = 0.18;

// for stackIndex i in 0..3:
transform: `translate(${i * PEEK_OFFSET_X}px, ${i * PEEK_OFFSET_Y}px) scale(${1 - i * SCALE_STEP})`,
opacity: Math.max(0.15, 1 - i * OPACITY_STEP),
zIndex: 100 - i,
```

Card 0 = full size/opacity, frontmost. Cards 1-3 progressively smaller, more
transparent, offset up-and-right, and painted behind (lower zIndex) so the
front card naturally covers most of their area except the intended peek
sliver.

Animate with `framer-motion` (`motion.div`, already a project dependency)
rather than raw CSS transitions, because:
- Cards need `layout`-aware animation when the stack reflows (item removed
  from the middle, or wheel changes front) -- `AnimatePresence` + per-card
  `animate` targets handle this cleanly across React key-based reorders.
- Use `<AnimatePresence mode="popLayout">` wrapping the mapped cards, each
  `motion.div` keyed by `item.id`:
  - `initial`: `{ opacity: 0, y: 20, scale: 0.9 }` (new card slides up from
    below the stack)
  - `animate`: the computed `{ x, y, scale, opacity }` for its current
    `stackIndex` (framer interpolates automatically on index change -- this
    is what produces the fluid wheel-driven slide, see section 5)
  - `exit`: `{ opacity: 0, x: 40, scale: 0.85, transition: { duration: 0.18 } }`
    (close (x) slides the card out to the right and fades)
  - `transition`: `{ type: 'spring', stiffness: 380, damping: 34 }` for
    stack-position changes.

### JSX sketch

Each card's `animate` target depends on whether it *was* the front card and
is now leaving (wheel-triggered), vs. a normal peek-position card advancing
toward the front. This requires tracking the previous front id across the
wheel gesture so the outgoing card gets the special exit-slide target instead
of the generic peek math for one frame of the transition.

```ts
const EXIT_PUSH_X = 60; // px, slides past the right edge (anchor corner)
const EXIT_PUSH_Y = 40; // px, slides down past the bottom edge

const prevFrontIdRef = useRef<string | null>(null);
const [leavingId, setLeavingId] = useState<string | null>(null);
```

```tsx
<div className={styles.stack} onWheel={handleWheel}>
  <AnimatePresence mode="popLayout">
    {visibleIds.map((id, stackIndex) => {
      const item = items.find((n) => n.id === id)!;
      const isLeaving = id === leavingId;
      return (
        <motion.div
          key={id}
          className={styles.card}
          data-variant={item.variant}
          style={{ zIndex: isLeaving ? 101 : 100 - stackIndex }}
          initial={{ opacity: 0, y: 20, scale: 0.9 }}
          animate={
            isLeaving
              ? {
                  // The outgoing front card: continues past the anchor
                  // corner (down-right) in the SAME motion the next card
                  // uses to slide up into front -- both driven by one
                  // wheel tick, one spring, opposite directions, so they
                  // read as a single continuous swap rather than two
                  // separate animations.
                  x: EXIT_PUSH_X,
                  y: EXIT_PUSH_Y,
                  scale: 1 - SCALE_STEP, // keeps receding-into-distance feel, not a snap-to-bigger jump
                  opacity: 0,
                }
              : {
                  x: stackIndex * PEEK_OFFSET_X,
                  y: stackIndex * PEEK_OFFSET_Y,
                  scale: 1 - stackIndex * SCALE_STEP,
                  opacity: Math.max(0.15, 1 - stackIndex * OPACITY_STEP),
                }
          }
          exit={{ opacity: 0, x: EXIT_PUSH_X, y: EXIT_PUSH_Y, scale: 0.85 }}
          transition={{ type: 'spring', stiffness: 380, damping: 34 }}
          onAnimationComplete={() => {
            if (isLeaving) setLeavingId(null); // unmount not required, card is simply gone (stackIndex > MAX_PEEK) or reflowed on the next store emit
          }}
        >
          {stackIndex === 0 && !isLeaving && (
            <>
              <NotificationIcon variant={item.variant} />
              <span className={styles.message}>{item.message}</span>
              <button
                className={styles.close}
                aria-label="Fermer"
                onClick={() => removeNotification(id)}
              >
                <X size={14} />
              </button>
            </>
          )}
          {stackIndex > 0 && !isLeaving && <div className={styles.peekFill} />}
          {isLeaving && <div className={styles.peekFill} />}
        </motion.div>
      );
    })}
  </AnimatePresence>
</div>
```

The card taking over the front slot (previously `stackIndex === 1`) is NOT
special-cased -- it already animates from its old peek offset
(`PEEK_OFFSET_X, PEEK_OFFSET_Y`, up-and-left-ish per section geometry) to
`stackIndex === 0`'s target (`{0, 0, 1, 1}`) via the normal branch, using the
same spring/duration as the outgoing card's exit push. Because both
transitions share one `transition` config and start on the same store emit,
they run in lockstep: the outgoing card recedes down-right while the incoming
one glides up into the vacated slot "in the same lancée" (per the user's
phrasing), rather than the incoming card just popping to full scale/opacity
in place.

(`peekFill` = an empty opaque-tinted glass rectangle for non-front cards --
no need to render icon/text/close for cards the user can't read anyway since
they are mostly hidden and scaled/faded down; this also avoids any
accessibility noise from off-screen-ish duplicate text. The leaving card also
switches to `peekFill` immediately so its own text doesn't visibly slide off
screen -- only the glass shape recedes.)

## 5. Wheel-cycling logic (`handleWheel`)

Interpretation chosen: scrolling down (`deltaY > 0`) cycles to the next card
further back in the current visual stack (the current `stackIndex 1` card
becomes front); scrolling up (`deltaY < 0`) cycles the other direction,
bringing back the previously-front card. This matches the "scroll down = go
forward/deeper" convention used elsewhere (e.g. carousels) and reverses
cleanly on scroll-up.

```ts
const wheelLockRef = useRef(false);

const handleWheel = (e: React.WheelEvent) => {
  if (orderedIds.length < 2) return;
  e.preventDefault(); // fixed-position overlay: never let this scroll page content behind it
  if (wheelLockRef.current) return; // one wheel "tick" = one step, not several
  wheelLockRef.current = true;
  setTimeout(() => { wheelLockRef.current = false; }, 220); // ~ spring settle time
  const delta = e.deltaY > 0 ? 1 : -1;
  const nextIdx = ((delta % orderedIds.length) + orderedIds.length) % orderedIds.length;
  const currentFront = orderedIds[0];
  setLeavingId(currentFront); // outgoing card takes the down-right exit path (see JSX sketch)
  setFrontNotification(orderedIds[nextIdx]);
};
```

Because `orderedIds` is recomputed from `frontId` + recency every render,
calling `setFrontNotification(nextId)` triggers a store emit -> re-render ->
`orderedIds` reshuffles so `nextId` is now at index 0 -> every card's
`stackIndex` (and therefore its `animate` target) changes -> framer-motion's
spring transition animates each card smoothly from its old slot to its new
slot. `leavingId` (set synchronously in the same handler, one React batch
before the store's re-render lands) makes sure the ex-front card's `animate`
branch is the down-right exit push instead of just landing on "whatever
`stackIndex` it ends up at" -- without it, the outgoing card would just be
treated like any other peek card sliding to a deeper slot, which is the
"generic spring-reflow" behavior the user explicitly asked NOT to have; with
it, the exit reads as one continuous down-right slide distinct from the
peek-position cards merely shifting one slot. Combined, this produces "la
notif courante glisse en bas vers la droite et la suivante prend sa place
dans la même lancée" -- both card animations are one spring transition
started on the same tick, so they visually read as a single swap.

### Interaction with new arrivals while a mid-stack item is focused to front

Per section 2/3, `pushNotification` unconditionally sets `frontId = <new id>`.
A new notification always interrupts and takes over the front slot, even if
the user had manually scrolled an older one to front. This is the simplest,
most predictable rule (new events are the most recent/important thing
happening) and matches how the old sonner stack behaved (new toast always
appears at the top position). The previously-front (user-pinned) card
animates back into the peek stack via the same `stackIndex` recompute +
spring, exactly like any other reordering.

## 6. Close (x) button behavior

`onClick={() => removeNotification(id)}` on the front card only (section 4
decision). `removeNotification`:
- Splices the item out of `items`.
- If the removed id was `frontId`, resets `frontId = null` so the selector
  falls back to "most recent remaining item" as the new front.
- Triggers `AnimatePresence`'s `exit` animation (slide right + fade) for the
  removed card while every remaining card's `stackIndex` shifts down by one
  and re-animates into its new slot via the same spring -- a natural "stack
  reflow".

## 7. Blur-without-ghosting -- resolution of the open question

CONCLUSION: no special stacking-context/isolation trickery is required.
`backdrop-filter` only samples pixels already painted behind the element
within the same compositing context, at the moment of paint. Because:
- Every card (front and peeking) is rendered with a translucent glass
  background (`rgba(10, 15, 26, 0.25)`, matching the existing recipe -- not
  fully opaque, by design), and
- Cards behind the front one are painted first (lower `zIndex`) and are
  positioned so the front card covers their bounding box except for the
  small peek slivers,

...the front card's `backdrop-filter: blur(4px)` correctly samples: (a) the
3D viewport wherever no card is behind it, and (b) the peeking cards' own
glass fill wherever they stick out from behind the front card. That second
case is NOT "ghosting" -- it is the intended visual (the peek slivers are a
deliberate part of the stacked-deck look). The real risk the task is
flagging -- peeking cards' TEXT/icons showing blurred through the front card
-- is structurally impossible here because ONLY the front card renders
text/icon content (section 4, `stackIndex === 0 &&` guard); peeking cards
render an empty `.peekFill` div, so there is no content behind the front card
that could "ghost" through the blur beyond the flat glass-tinted rectangles,
which is the desired peek effect itself.

The one real thing to add to the SCSS: give `.stack` `isolation: isolate;`
-- not to fix ghosting (there is none, per above), but purely for
defensive consistency with the existing codebase convention
(`GCodePreviewBar.module.scss`'s `.bar` also sets `isolation: isolate;`),
so the stack's z-index layering never interacts unpredictably with sibling
fixed-position chrome (SliceSidebar toolbar, GCodePreviewBar) that also
creates stacking contexts.

## 8. SCSS (`src/components/NotificationStack.module.scss`)

Follow the exact glass recipe from `GCodePreviewBar.module.scss` /
`SliceSidebar.module.scss`:

```scss
.stack {
  position: fixed;
  right: 24px;
  bottom: 84px;
  width: 0;
  height: 0;
  z-index: 60; // above viewport chrome, consistent with other floating panels
  isolation: isolate;
  pointer-events: none; // re-enabled per-card below so empty space isn't a hit target
}

.card {
  position: absolute;
  right: 0;
  bottom: 0;
  min-width: 280px;
  max-width: 360px;
  display: flex;
  align-items: center;
  gap: 10px;
  padding: 12px 14px;
  border-radius: 14px;
  border: none;
  box-shadow: 0 8px 32px 0 rgba(0, 0, 0, 0.4), 0 2px 8px 0 rgba(0, 0, 0, 0.3);
  background: rgba(10, 15, 26, 0.25);
  backdrop-filter: blur(4px);
  -webkit-backdrop-filter: blur(4px);
  color: rgb(212, 212, 216);
  pointer-events: auto;

  &:not(:first-child) {
    pointer-events: none; // peeking cards are inert, see sections 4/5
  }
}

.peekFill {
  width: 100%;
  height: 100%;
  border-radius: inherit;
}

.message {
  flex: 1;
  font-size: 13px;
  line-height: 1.4;
}

.close {
  flex: 0 0 auto;
  display: flex;
  align-items: center;
  justify-content: center;
  width: 22px;
  height: 22px;
  border-radius: 6px;
  background: transparent;
  border: none;
  color: rgba(212, 212, 216, 0.7);
  cursor: pointer;
  &:hover { background: rgba(255, 255, 255, 0.08); color: rgb(212, 212, 216); }
}

// Variant accent (left icon color only -- keep the glass bg identical across
// variants per "matching the rest of the chrome", don't reintroduce sonner's
// colored-background style):
.card[data-variant='success'] svg { color: #4ade80; }
.card[data-variant='error'] svg { color: #f87171; }
.card[data-variant='warning'] svg { color: #fbbf24; }
.card[data-variant='loading'] svg { color: #a1a1aa; }
```

`.card` uses `position: absolute` with transform-driven offsets (not
flex-gap stacking) so framer-motion's `animate` targets remain the single
source of truth for position -- `.stack` itself collapses to a 0x0 anchor
point at `right: 24px; bottom: 84px;`, and each `.card` anchors to that point
via `right: 0; bottom: 0;`, growing left/up as needed via
`min-width`/`max-width`.

## 9. `src/App.tsx` changes

- Remove `import { Toaster } from '@/components/ui/sonner';`
- Add `import NotificationStack from '@/components/NotificationStack';`
- Replace the `<Toaster ... />` JSX block with `<NotificationStack />`, kept
  in the same position in the tree (sibling to `<Routes>`, inside `<Router>`).
- Trim/replace the large doc-comment above the old `<Toaster>` -- the styling
  knowledge now lives in `NotificationStack.tsx` / `.module.scss` comments.

## 10. New event hook points and exact message text to add

All new dispatch sites import `showSuccess` from `@/utils/toast` (already the
established pattern) at the EXISTING funnel points identified in prior
research -- no new `window` events needed, since `toast.ts` already provides
the decoupling the app's event-bus convention is for; these are direct
function calls colocated with the action, matching how `Index.tsx:171`
already calls `showSuccess` directly after dispatching the slice event.

### a. Space changes -- `src/components/Sidebar.tsx`, `playSpaceTransition` (~139)

Both `handleTopTabClick` and `enterPreview` funnel through
`playSpaceTransition(space, onCovered)`, and `handleLayerClick`'s exit branch
also calls `playSpaceTransition('box', ...)` -- so a single call site covers
all space-change paths. CORRECTION vs. the naive placement: `label` is
computed synchronously at the top of `playSpaceTransition`, but the actual
tab switch only happens ~900ms later, inside each caller's `onCovered`
callback, once the curtain (`SpaceTransitionOverlay`) has fully covered the
screen (see the doc-comment on `playSpaceTransition`, lines 131-138 -- "c'est
là, et pas avant, que le vrai changement d'onglet doit avoir lieu"). Firing
`showSuccess` immediately after `const label = ...` would announce the space
change before the curtain even starts rising, out of sync with what the user
sees. Instead, wrap `playSpaceTransition` itself so the notification fires
inside `onCovered`, alongside the real switch:

```ts
const playSpaceTransition = (space: keyof typeof SPACE_LABELS, onCovered: () => void) => {
  const label = SPACE_LABELS[space];
  // ...unchanged geometry/sound logic...
  window.dispatchEvent(new CustomEvent('apex-space-transition', {
    detail: {
      label,
      onCovered: () => {
        onCovered();
        showSuccess(`Espace ${label} activé.`);
      },
    },
  }));
};
```

This keeps a single edit point (all four callers are unaffected) and fires
the notification at the same moment the space visually becomes active.

Resulting messages per space: "Espace Préparation activé.", "Espace
Modélisation activé.", "Espace Dev activé.", "Espace Paramètres activé.",
"Espace Prévisualisation activé."

### b. Tool changes -- `src/components/Sidebar.tsx`, `handleToolClick` (~184) and the keyboard shortcut block (~272-299)

Add a label map and a small helper, called from both the toolbar click
handler and the keydown handler (both already compute a `newTool` value):

```ts
const TOOL_LABELS: Record<'select' | 'move' | 'rotate' | 'scale' | 'texture' | 'layers', string> = {
  select: 'Selection',
  move: 'Deplacement',
  rotate: 'Rotation',
  scale: 'Echelle',
  texture: 'Texture',
  layers: 'Previsualisation',
};

const announceTool = (tool: keyof typeof TOOL_LABELS) => {
  showSuccess(`Outil ${TOOL_LABELS[tool]} active.`);
};
```

(Again, use proper accents -- "Sélection", "Déplacement", "Échelle",
"Prévisualisation", "activé" -- in the actual source.)

Call `announceTool(newTool)` at the end of `handleToolClick` and at the end
of each `m`/`r`/`s`/`t` branch in the keydown handler. For `Escape`, guard
with `if (activeTool !== 'select') announceTool('select');` before the state
update, to avoid a no-op notification when Escape is pressed with nothing
active.

### c. G-code color-by mode change -- `src/hooks/use-gcode-preview.ts`, `setColorBy` (line 169)

```ts
const setColorBy = useCallback((mode: ColorByMode) => {
  setColorByState(mode);
  const label = mode === 'lineType' ? 'Type de ligne' : COLOR_BY_LABELS[mode];
  showSuccess(`Coloration : ${label}.`);
}, []);
```

CORRECTION: no new label map needs to be created -- `COLOR_BY_LABELS: Record<ContinuousColorByMode, string>`
already exists in `src/components/gcodeColormap.ts` (imported by
`GCodeLineTypeLegend.tsx` today for the exact same purpose, e.g.
`COLOR_BY_LABELS[mode]` at GCodeLineTypeLegend.tsx:92,134). It only covers
`ContinuousColorByMode` (not the `'lineType'` member of the full
`ColorByMode` union, per `gcodeLayer.ts:70-75`), hence the ternary above.
Implementer imports `COLOR_BY_LABELS` from `./gcodeColormap` into
`use-gcode-preview.ts`.

### d. Model import success/failure -- `src/components/Viewport3D.tsx`, `onImportModel` (~3249-3366)

Exact insertion point verified in source: right after the existing
`window.dispatchEvent(new CustomEvent('apex-mesh-selected', ...))` call
(~3334-3341), still inside the `if (r) { ... }` block:

```ts
showSuccess(`Modèle importé : ${uniqueName}.`);
```

CORRECTION/addition vs. original plan: the `catch (err)` block (~3363-3365)
currently only does `console.error("Failed to import 3D model:", err)` --
completely silent to the user. Since this task is about surfacing events
including failures (not just successes), add a matching error notification
there too, consistent with the `showError` pattern used elsewhere (e.g.
use-gcode-preview.ts's load-failure toasts):

```ts
} catch (err) {
  console.error("Failed to import 3D model:", err);
  showError(`Échec de l'import du modèle : ${file.name}.`);
}
```

### e. Model delete (single) -- `src/components/Viewport3D.tsx`, `deleteSelection` (~2408-2430)

After the `selection.forEach(...)` loop completes, before `selection.clear()`
(need the count before clearing):

```ts
showSuccess(selection.size > 1
  ? `${selection.size} modeles supprimes.`
  : 'Modele supprime.');
```

(accented: "modèles supprimés" / "Modèle supprimé.")

### f. Model delete (all) -- `src/components/Viewport3D.tsx`, `onDeleteAll` (~3402-3423)

After the `all.forEach(...)` loop, before `selectedMeshesRef.current.clear()`:

```ts
showSuccess(`${all.length} modele(s) supprime(s) du plateau.`);
```

(accented: "modèle(s) supprimé(s)")

### g. Slicing start -- `src/pages/Index.tsx`, `handleConfirmSlice` (~163-167)

```ts
const handleConfirmSlice = () => {
  window.dispatchEvent(new CustomEvent('apex-slice-all-models'));
  setIsSliceModalOpen(false);
  setIsSlicingProgressOpen(true);
  showSuccess('Decoupage du plateau demarre...');
};
```

(accented: "Découpage du plateau démarré...")

This is additive alongside the existing completion toast
(`showSuccess("Slicing plate completed successfully!")` at line 171, inside
`handleSlicingFinished`) -- both remain, giving a start AND a completion
notification, each permanent until dismissed.

## 11. Porting existing showSuccess/showError/showWarning/showLoading call sites

No code changes needed at the ~25 existing call sites (section 3) -- they
automatically ride the new system once `toast.ts` is repointed. Verification
step for implementer: confirm no call site imports `toast` directly from
`sonner` (bypassing `utils/toast.ts`):
`grep -rn "from 'sonner'" src/` (excluding `src/components/ui/sonner.tsx` and
`src/utils/toast.ts` itself) should return zero results; port any hits found
to `@/utils/toast`.

## 12. Testing/validation checklist for implementer

- Multiple notifications queued -> verify peek offsets match the reference
  drawing orientation (up-and-right).
- Wheel over the stack cycles the front card with a visible spring animation;
  wheel outside the stack area does nothing (handler lives on `.stack`, not
  `window`).
- A new `showX` call while a non-front card is user-focused correctly
  interrupts and takes front (per section 5's chosen rule).
- Close (x) always closes the visually-front card only; remaining cards
  reflow.
- No card ever auto-dismisses (no `setTimeout` anywhere removes an item).
- Blur shows the 3D viewport through the front card in an empty area of the
  stack (no other cards behind it at that point), and shows only the
  peek-glass tint (never text) through the sliver areas.
- `sonner` is not imported anywhere reachable from the running app (App.tsx
  clean), while `ui/sonner.tsx` on disk is left untouched.

### Critical Files for Implementation
- src/utils/notificationStore.ts (new)
- src/components/NotificationStack.tsx (new)
- src/components/NotificationStack.module.scss (new)
- src/utils/toast.ts (edit: repoint to notificationStore)
- src/App.tsx (edit: swap Toaster for NotificationStack)
- src/components/Sidebar.tsx (edit: add space/tool change notifications)
- src/hooks/use-gcode-preview.ts (edit: add color-by change notification)
- src/components/Viewport3D.tsx (edit: add import/delete notifications)
- src/pages/Index.tsx (edit: add slicing-start notification)