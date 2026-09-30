# Plan — Transition Préparation/Prévisualisation ↔ Modélisation : corriger la cause réelle du flash

## Contexte

Pendant la transition vers Modélisation, on voit pendant **une frame** soit un « carré plus large »
(avant l'ajout du filet CSS), soit un panneau « transparent » (après). Les correctifs précédents
(unifier les boucles RAF, retirer les `motion.div`, etc.) visaient de mauvaises causes.

**Cause réelle (vérifiée dans le code) :** l'app utilise `createRoot` (React 18 concurrent). Dans
[ModelingSidebar.tsx:198-207](src/components/ModelingSidebar.tsx#L198-L207), la fin du morph fait, dans
le callback RAF :

```ts
el.style.clipPath = '';   // mutation DOM SYNCHRONE → visible à la prochaine peinture
setMorphStage(null);      // mise à jour React PLANIFIÉE → commitée dans une tâche ULTÉRIEURE
morphRef.current?.onDone();
```

Un `setState` hors gestionnaire d'événement React est traité par le Scheduler dans une tâche
séparée, donc **après** la peinture de cette frame. Résultat : pendant une frame, le clip-path est
retiré alors que `.morphing` et le style « stage » (364px+) du conteneur sont toujours là :
- sans filet CSS → panneau non découpé à la taille du stage = le **carré plus large** ;
- avec le filet `clip-path: inset(0 0 100% 0)` → panneau masqué, seule l'ombre `LiquidGlassShadow`
  (dessinée autour de la forme finale) reste = le **rectangle arrondi « transparent »** de la capture.

Le même défaut d'ordonnancement existe en fin de morph de sortie, aggravé par
`AnimatePresence mode={modelingMorph ? 'sync' : 'wait'}` dans Index.tsx, qui retarde le démontage de
l'ancien côté et le montage du nouveau d'au moins une tâche.

Deux décisions des tours précédents sont aussi à corriger :
- la pilule entière est estompée au clic (fausse piste « plus large ») alors que c'est **elle** la
  forme qui se déforme : seules ses icônes doivent disparaître ;
- le temps mort `ELEMENTS_FADE_MS` au début du chemin de morph est devenu inutile (les éléments
  s'estompent désormais dès le clic, ~900ms avant le début du morph), et en sortie il se retrouve à
  la **fin** (timeline inversée) au lieu du début.

## Principe de la correction

**Une seule source de vérité : la prop `morph`.** Tout ce qui décrit « un morph est en cours »
(style stage du conteneur, classe `.morphing`, ombre) est **dérivé de `morph` pendant le rendu**,
plus jamais d'un état local mis à jour depuis le RAF. Le RAF ne fait que dessiner les frames puis
appeler `onDone()`. Le retrait du clip-path se fait dans le **cleanup du `useLayoutEffect`**, qui
s'exécute pendant le commit, avant la peinture. La bascule entre les deux côtés se fait dans un
**seul commit**, sans `AnimatePresence`.

Séquence cible (les deux sens sont symétriques) :

| Étape | Vers Modélisation | Hors de Modélisation |
|---|---|---|
| Clic | ModelStatusBox, panneau du bas Preview et icônes de la pilule s'estompent en cascade ; la forme (pilule / panneau haut) reste | le contenu du panneau Modélisation s'estompe en cascade ; la forme reste |
| `onCovered` | même commit : CategoriesSidebar démontée, ModelingSidebar montée à la forme de départ ; la déformation commence | la déformation commence |
| Fin du morph | même commit : `.morphing`/stage retirés, clip effacé (layout effect) ; le contenu apparaît en cascade | même commit : ModelingSidebar démontée, CategoriesSidebar montée déjà dans sa forme finale ; icônes / panneau du bas / ModelStatusBox apparaissent |

## Modifications

### 1. [ModelingSidebar.tsx](src/components/ModelingSidebar.tsx) — le cœur du correctif
- Supprimer l'état `morphStage` / `setMorphStage`. Le dériver dans le rendu :
  `const stage = morph ? modelingMorphStage(window.innerHeight) : null;`
  Il sert au style du conteneur, à la classe `.morphing` et au montage de `<LiquidGlassShadow>`,
  qui est ainsi monté dès le commit où `morph` apparaît (ce qui corrige aussi l'ombre absente
  sur la 1re frame d'une sortie).
- `useLayoutEffect` (deps `[morph?.direction, morph?.fromPreview]`) :
  - le calcul de `stage` reste local à l'effet ; `apply(0)` synchrone ; boucle RAF inchangée ;
  - fin (`progress === 1`) : **seulement** `morphRef.current?.onDone()`, sans `el.style.clipPath = ''`
    ni `setState` de géométrie ;
  - **cleanup** : `cancelAnimationFrame(raf); el.style.clipPath = '';` — exécuté dans le commit où
    `morph` redevient `undefined` (même commit que le retrait de `.morphing`/stage) ou au démontage.
- Contenu en cascade :
  - écouter un nouvel événement `apex-leaving-modeling` → `setContentLeaving(true)` ;
  - `animate={contentRevealed && !contentLeaving ? 'visible' : 'hidden'}` ;
  - `contentContainerVariants.hidden` reçoit aussi `staggerChildren` (même ordre, du haut vers le bas) ;
  - en sortie, `revealed` démarre à `false` (le contenu est déjà parti au clic) ; en entrée, il
    reste `false` jusqu'à `t === 1` (`CONTENT_REVEAL_THRESHOLD = 1`, déjà en place).
- Garder le filet CSS `.panel.morphing { clip-path: inset(0 0 100% 0) }` (protection de la
  1re frame), ainsi que `.container`/`.panel` en simples `<div>`.

### 2. [spaceMorphGeometry.ts](src/utils/spaceMorphGeometry.ts)
- Retirer le temps mort du chemin : `modelingPathAt` redevient
  `stretchPath(fromPreview ? upper : PILL, target, t)`.
- Durées : `MODELING_MORPH_MS = 550` (depuis la pilule), `MODELING_MORPH_FROM_PREVIEW_MS = 480`.
- Garder `ELEMENTS_FADE_MS` (et son alias `FROM_PREVIEW_FADE_MS`) **uniquement** comme durée des
  fondus déclenchés au clic.

### 3. [Index.tsx](src/pages/Index.tsx)
- Remplacer `<AnimatePresence mode={…}>` autour des deux branches par un rendu conditionnel
  simple : bascule en un seul commit, sans démontage retardé. Rien n'en dépend :
  SliceSidebar n'utilise pas framer-motion, PrinterConfigSidebar a son propre `AnimatePresence`
  interne, et les fondus de CategoriesSidebar passent déjà par `animate`.
- Ajouter l'état `arrivingFromModeling` : passé à `true` dans la branche `exit` d'`onTabChange`,
  lu une seule fois au montage par CategoriesSidebar, remis à `false` dans un `useEffect` quand
  `modelingMorph` redevient `null`.

### 4. [CategoriesSidebar.tsx](src/components/CategoriesSidebar.tsx)
- **Annuler** le fondu de la pilule entière : revenir à `animate={{ opacity: 1, x: 0 }}` et à la
  transition 0.9s d'origine. Seules les icônes s'estompent au clic (`leavingToModeling`, déjà en
  place).
- Nouvelle prop `arrivingFromModeling`, figée au montage (`useState`), comme `wasInitialPreviewActive` :
  - conteneur de la pilule : `initial={arriving ? false : { opacity: 0, x: -20 }}`, pour qu'elle soit
    opaque dès la frame de bascule ;
  - icônes : `initial={arriving ? 'hidden' : false}`, `animate='visible'` → cascade d'apparition.
- Panneau du bas (arrivée en Prévisualisation) : la logique `wasInitialPreviewActive` actuelle est
  conservée (fondu-glissement d'entrée).

### 5. [Sidebar.tsx](src/components/Sidebar.tsx)
- Au clic, si `isModeling` et que la cible n'est pas Modélisation (dans `handleTopTabClick` pour
  `'box'`/`'developer'`, et dans la branche `entering` de `handleLayerClick`), dispatcher
  `apex-leaving-modeling` **avant** `playSpaceTransition`. C'est le symétrique d'`apex-entering-modeling`.

### Inchangé (déjà correct)
- `ModelStatusBox` : disparaît sur `apex-entering-modeling`, réapparaît sur
  `apex-mode-changed: 'slicer'` (démontage de ModelingSidebar, donc après la déformation).
- `LiquidGlassShadow` passif piloté par `applyLiquidGlassShadow` dans la même boucle RAF.
- `.container`/`.panel` sans `motion.div` (pas de `transform` sur les ancêtres du `backdrop-filter`).

## Vérification

1. `npx tsc --noEmit -p .` puis `npx eslint` sur les fichiers modifiés.
2. **Sonde d'invariant Playwright** (script jetable dans le scratchpad, avec throttling CPU via
   CDP `Emulation.setCPUThrottlingRate` pour se rapprocher d'une vraie machine) : un
   `requestAnimationFrame` injecté relève à chaque frame l'état de `.panel` et de son conteneur.
   Échec si une frame montre :
   - `.morphing` présent **et** `style.clipPath` vide (= le flash transparent) ;
   - `.morphing` absent **et** style stage encore présent (= le carré plus large) ;
   - aucun panneau gauche visible, ni pilule/panneaux Preview ni ModelingSidebar (= trou de bascule).
   À exécuter pour les 5 transitions : Préparation→Modélisation, Prévisualisation→Modélisation,
   Modélisation→Préparation, Modélisation→Prévisualisation, Modélisation→Dev.
3. Non-régression : mitose Préparation↔Prévisualisation (pas de fondu parasite du panneau du bas).
4. Validation manuelle par toi, en enregistrant à 120fps sur http://localhost:8082 (serveur dev déjà lancé).
