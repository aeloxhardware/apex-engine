# Plan : sélecteur "colorer par" pour la légende Prévisualisation G-code

## Contexte

Aujourd'hui, la Prévisualisation G-code (`GCodeLineTypeLegend.tsx`, montée dans
`Index.tsx` au-dessus du panneau `CategoriesSidebar`) ne colore le tracé 3D que par
**type de ligne** (paroi, remplissage, pont…), avec une légende tableau (couleur, durée,
%, longueur, poids) et une case d'affichage par ligne. L'utilisateur veut pouvoir
recolorer ce même tracé selon d'autres critères utiles au diagnostic d'impression
(vitesse, température, hauteur de couche, ventilateur…), un seul mode actif à la fois,
avec un sélecteur groupé (trop d'options — confirmé ~13 en v1 — pour un simple bouton
segmenté). Objectif : donner une vraie visibilité sur ce qui se passe machine pendant
l'impression, au-delà du seul type de ligne, sans dégrader la fluidité de la
prévisualisation (contrainte utilisateur existante : [[feedback_fluidity]], aucun
à-coup visible, tolérance ≤ 0,1 mm).

Cette tâche ne fait QUE poser l'architecture ; l'implémentation attend le feu vert de
l'utilisateur (autre tâche en cours en parallèle).

## 0. Périmètre v1 vs différé (15 critères proposés)

**V1 — 13 critères**, groupés par coût réel de parsing (confirmé avec l'utilisateur) :

Déjà disponibles ou dérivables à coût quasi nul (données suivies puis jetées aujourd'hui) :
1. **Type de ligne** — existant, reste le mode par défaut, inchangé.
9. **Accélération** — `accels[]` déjà calculé par `push()`, juste à ne plus le jeter.
11. **Inclinaison de la tête (axe B)** — déjà stocké sur `GCodeMove.b`.
14. **Longueur du segment** — dérivable des coordonnées consécutives (formule déjà
    inline dans `planCartesianTimes` pour le cartésien ; à extraire en helper partagé
    pour couvrir aussi le polaire).
8. **Temps écoulé** — `times[]` cumulé déjà retourné, aucun changement parseur.
15. **Pression d'extrusion relative (E/mm)** — ratio simple depuis `eDeltas[]` +
    longueur de segment.
7. **Débit volumétrique (mm³/s)** — E converti en volume (diamètre filament déjà connu)
    divisé par le dt du mouvement (déjà disponible via les deltas de `times[]`).

**Vitesse (#2) — décision utilisateur : vitesse RÉELLE, pas F commandé.** Contrairement
à la recommandation initiale (F commandé, plus simple), l'utilisateur a choisi la
vitesse réellement atteinte après limitation par accélération/jerk aux coins — c'est
justement dans les petits détails/coins que ce mode est le plus utile. Nécessite que
`planCartesianTimes` retourne aussi un tableau `speeds` (la vitesse de croisière `vp`/de
jonction déjà calculée en interne, simplement remontée) — touche une fonction sensible
aux perfs mais sans en changer l'algorithme, seulement ce qu'elle expose.

Nécessitent une nouvelle logique de détection de couche (auto-contenue) :
4. **Hauteur de couche** — y compris en **polaire** : décision utilisateur, l'épaisseur
   locale existe même sans palier Z global net. Calculée comme le delta Z entre le
   mouvement courant et le précédent mouvement extrudant au **même angle C (modulo
   360°, à une tolérance angulaire à définir, ex. ±2°)** — reflète l'épaisseur réellement
   déposée à cet endroit précis, cohérent avec la nature hélicoïdale de l'impression
   polaire. En cartésien, garde la détection par palier Z classique (delta Z entre
   couches).
6. **Largeur de ligne** — commentaire générateur (`;WIDTH:`) si présent, sinon estimée
   par `largeur = volume_E / (longueur_segment × hauteur_de_couche)` ; dépend donc de #4
   et doit être livré avec.
12. **Rapport vitesse / vitesse max de la couche** — nécessite un regroupement par
    couche (même détection que #4) puis un lookup du F max de cette couche.

Nécessitent un parsing M-code entièrement nouveau (confirmé maintenu en v1 par
l'utilisateur, livrés en dernier dans le lot car plus de risque parseur) :
3. **Température** (buse ; plateau si elle varie) — nouveau : `M104`/`M109` → temp buse,
   `M140`/`M190` → temp plateau. Traité comme un état modal (consigne commandée, pas
   distinction attente/atteinte — suffisant pour une prévisualisation statique).
5. **Vitesse du ventilateur** — nouveau : `M106 S<0-255>` → %, `M107` → 0.

**Différé phase 2** (hors périmètre de cette implémentation, pas même en entrée de menu
grisée — absence honnête plutôt qu'une option à moitié fonctionnelle) :
10. **Vitesse angulaire de l'axe C** — spécifique polaire, aucun équivalent pour les
    fichiers cartésiens (la majorité des fichiers chargés) → faible valeur en v1.
13. **Score de risque de défaut** — dépend de #2 (vitesse), du type de ligne (surplomb)
    et de #5 (ventilateur) déjà éprouvés, plus une formule de pondération sans précédent
    dans ce code — à concevoir une fois les autres critères validés à l'usage.

## 1. Extension du modèle de données (`src/components/gcodeLayer.ts`)

### 1.1 Nouveaux états modaux et parsing M-code

Dans `parseGCode()`, à côté de `feed`/`accelPrint` existants :

```ts
let nozzleTemp = 0;   // °C, M104/M109 S
let bedTemp = 0;      // °C, M140/M190 S
let fanPct = 0;        // 0-100 %, M106 S(0-255)→%, M107→0
```

Nouveaux cas dans le switch M :
```ts
else if (num === 104 || num === 109) { if (w.S >= 0) nozzleTemp = w.S; }
else if (num === 140 || num === 190) { if (w.S >= 0) bedTemp = w.S; }
else if (num === 106) { if (w.S >= 0) fanPct = Math.round((w.S / 255) * 100); }
else if (num === 107) { fanPct = 0; }
```
Poussés dans des tableaux parallèles à chaque `push()` (comme `feeds.push(feed)`
aujourd'hui), donc automatiquement corrects après subdivision d'arc.

### 1.2 Vitesse réelle (post-planification)

`planCartesianTimes()` (lignes ~186-253) calcule déjà `vp` (vitesse de croisière après
limitation accel/jerk) en interne pour chaque mouvement mais ne la retourne pas.
Ajouter un tableau de sortie `speeds: Float64Array` rempli avec `vp` (ou la vitesse
moyenne du mouvement) à côté du calcul de `times[k]`. Pour le polaire, dérive la vitesse
réelle directement de `distance / (times[k] - times[k-1])` (déjà tout ce qu'il faut,
pas de planificateur trapézoïdal séparé côté polaire).

### 1.3 Détection de couche (cartésien et polaire) — pour #4, #6, #12

- **Cartésien** : détection par palier Z existant à ajouter — sur chaque mouvement
  extrudant, si `z` diffère du dernier Z extrudant vu de plus de `1e-4`, nouvelle couche ;
  `hauteur = z - prevLayerZ`. Stocker dans un tableau parallèle `layerHeightMm[]`
  (répète la même valeur pour tous les mouvements de la couche) + `layerIndex:
  Uint32Array`.
- **Polaire** : pour chaque mouvement extrudant à l'angle `c`, chercher en arrière le
  dernier mouvement extrudant dont l'angle `c' ≡ c (mod 360°)` est à moins d'une
  tolérance (ex. ±2°, à ajuster à l'usage) ; `hauteur_locale = z_courant - z_de_ce_point`.
  Implémentation pragmatique : garder un petit historique indexé par angle (bucket de
  ~1-2° sur 360°) du dernier Z extrudant vu à cet angle, mis à jour au fil du parsing —
  évite une recherche arrière coûteuse (O(1) par mouvement au lieu de O(n)).
  `layerIndex` n'a pas vraiment de sens global en polaire (pas de palier net) : utiliser
  directement `layerHeightMm[]` comme signal continu, et pour #12 (rapport vitesse/max
  de couche) soit adapter la notion de "couche" au même découpage par bucket angulaire,
  soit — si trop bruité à l'usage — accepter que #12 reste approximatif en polaire (à
  valider une fois implémenté, sur un vrai fichier Apex).
- **Largeur (#6)** : nouveau regex `WIDTH_RE = /^\s*WIDTH\s*[:=]\s*([\d.]+)/i` dans
  `handleComment()` (même schéma que `;TYPE:`), état modal `explicitWidth`, appliqué aux
  mouvements suivants ; sinon fallback `largeur = volume_E / (longueur_segment ×
  hauteur_de_couche)` avec les tableaux de 1.1/1.3/1.4.

### 1.4 Longueur de segment (#14)

Extraire un helper partagé `moveLength(a, b, kind): number` — réutilisé par
`planCartesianTimes` (remplace son calcul local `len[]`, une seule implémentation) et
par le calcul polaire de longueur d'arc déjà inline dans la boucle `times[]` polaire
(lignes ~586-597, à extraire dans le même helper). Stocké dans `segLengthMm:
Float32Array`.

Débit volumétrique (#7) et E/mm (#15) : calculés à la volée (pas stockés) depuis
`deposits[]`/`eDeltas[]`, `segLengthMm[]` et les deltas de `times[]` — ce sont de purs
ratios dérivés, inutile de dupliquer la mémoire pour un gros fichier.

### 1.5 `GCodeMove` inchangé, nouveaux champs sur `ParsedGCode`

`GCodeMove` reste `{x, y, z, b, extruding, type}` — pas de champs numériques ajoutés par
mouvement (éviterait un alourdissement mémoire par objet). Les nouvelles données
vivent en `Float32Array`/`Uint32Array` parallèles sur `ParsedGCode`, exactement comme
`times: Float64Array` aujourd'hui :

```ts
export interface ParsedGCode {
  // ... champs existants ...
  speedMmS: Float32Array;       // vitesse réelle (post-planification), §1.2
  accelMmS2: Float32Array;
  nozzleTempC: Float32Array;
  bedTempC: Float32Array;
  fanPct: Float32Array;
  layerHeightMm: Float32Array;  // §1.3, "n/a" représenté par NaN si indéterminable
  widthMm: Float32Array;
  segLengthMm: Float32Array;
  layerIndex: Uint32Array;      // cartésien uniquement ; 0 partout en polaire
  layerMaxSpeedMmS: Float32Array; // table indexée par layerIndex (cartésien)
  ranges: ColorByRanges;        // §1.6
}
```

### 1.6 Agrégation min/max/moyenne par fichier (nouveau pattern)

Aucune agrégation continue n'existe aujourd'hui (`lineStats` est discret, par type).
Ajouter une passe unique après construction des tableaux :

```ts
export interface ColorByRange { min: number; max: number; avg: number }
export type ColorByRanges = Record<ContinuousColorByMode, ColorByRange>;
```

Un balayage linéaire supplémentaire, négligeable par rapport aux passes O(n) déjà
existantes dans `parseGCode`. Cas `min === max` (fichier à valeur constante, ex. une
seule température) à gérer explicitement en aval (couleur unique, pas de division par
zéro dans la normalisation).

## 2. Nouvel utilitaire de couleur — `src/components/gcodeColormap.ts`

```ts
export function normalize(value: number, range: ColorByRange): number; // clamp((v-min)/(max-min||1), 0, 1)
export function viridisLike(t: number, out: BABYLON.Color3): BABYLON.Color3; // rampe 4 points
export function colorForMode(mode: ColorByMode, k: number, parsed: ParsedGCode, out: BABYLON.Color3): BABYLON.Color3;
```

Rampe RGB à 4 points (bleu foncé → cyan → jaune → rouge), interpolation linéaire via
`BABYLON.Color3.Lerp` (confirmé inutilisé dans tout le repo aujourd'hui, mais adapté).
Pas de conversion HSV (évite les problèmes de teinte qui boucle). Une seule fonction
sert tous les modes continus ; le mode catégoriel (type de ligne) continue d'utiliser
`LINE_TYPES[].color` sans passer par ce module.

## 3. Stratégie de rendu Babylon — décision

**Hybride, selon catégoriel vs continu :**

- **Type de ligne (catégoriel)** : architecture actuelle intégralement conservée — un
  maillage par `LineTypeId`, un `StandardMaterial` plat, `setEnabled` pour la
  visibilité. Zéro risque de régression sur le mode par défaut.

- **Modes continus** : un maillage fusionné pour les cordons extrudés + un maillage
  fusionné pour déplacements/essuyage, colorés par **couleur par sommet**
  (`VertexData.colors`), au lieu d'un maillage par catégorie. Justification : le
  bucketing en bandes discrètes donnerait un effet de bandes visibles (contraire à la
  qualité visuelle attendue d'un dégradé) sans réduire la complexité d'implémentation ;
  un maillage fusionné réduit même le nombre d'appels de dessin par rapport
  aux ~11 maillages actuels par type d'extrusion.
  - `buildBeadMesh`/`buildPathMesh` reçoivent un callback optionnel
    `vertexColor: (k) => Color3`, exclusif du paramètre `material`/`color` existant ;
    écrit dans `VertexData.colors` (4 floats/sommet, alpha 1) au fil de la boucle
    d'écriture de sommets existante — pas de passe séparée.
  - Matériau partagé unique en mode continu : `disableLighting = true`,
    `emissiveColor = Color3.White()` pour que les couleurs par sommet ne soient pas
    teintées par l'éclairage.
  - `buildPathMesh` (LinesMesh) : Babylon supporte les couleurs par sommet nativement
    via `mesh.useVertexColor = true` + `VertexData.colors`.

- **Coût de bascule de mode — décision utilisateur (mise à jour du plan initial)** :
  contrairement à la recommandation initiale ("doubler la mémoire en cachant les deux
  géométries pour zéro à-coup, même au tout premier switch"), l'utilisateur a choisi
  **d'accepter une petite reconstruction ponctuelle au premier changement de mode**
  catégoriel↔continu (pas de doublement mémoire par défaut). Rationale acceptée par
  l'utilisateur : la règle de fluidité vise la lecture/animation en cours, pas une
  action manuelle ponctuelle de changement de mode.
  - Comportement : la géométrie fusionnée (pour les modes continus) n'est construite
    qu'à la première activation d'un mode continu, puis mise en cache pour le reste de
    la session sur ce fichier — les bascules suivantes (y compris retour au mode
    catégoriel) redeviennent instantanées (juste `setEnabled` + éventuellement mise à
    jour du buffer de couleur via `updateVerticesData(ColorKind, ...)`, jamais un
    `VertexData.applyToMesh` complet).
  - **Nouveau réglage utilisateur** : ajouter dans l'espace **Paramètres**
    (`src/components/ParamsPanel.tsx`, actuellement un panneau vide/placeholder — ce
    sera son premier contenu réel) une option pour basculer vers le comportement
    "toujours zéro à-coup, même au premier switch" (double cache mémoire), au cas où
    l'utilisateur constaterait un à-coup gênant sur de très gros fichiers. Stockée
    comme préférence simple (ex. `localStorage`, à confirmer avec le pattern de
    persistance déjà utilisé ailleurs dans l'app le cas échéant), lue par
    `gcodeLayer.ts` au moment de la construction du calque.

- **Segment de tête en cours de lecture (`liveBead`/`liveLine`)** : en mode continu,
  calculer `colorForMode(...)` une fois par frame (O(1), pas d'allocation avec une
  `Color3` réutilisée) et l'appliquer à `liveBead.material.emissiveColor` directement —
  le segment live est trop petit pour justifier une mise à jour de couleurs par sommet
  à chaque frame (moins coûteux, meilleur pour la fluidité).

## 4. Câblage événements/état

Nouvel événement, sur le modèle exact de `apex-gcode-line-visibility` :

```ts
// use-gcode-preview.ts
export type ColorByMode =
  | 'lineType' | 'speed' | 'accel' | 'headAngle' | 'segmentLength' | 'elapsedTime'
  | 'extrusionRatio' | 'volumetricFlow' | 'layerSpeedRatio' | 'lineWidth'
  | 'layerHeight' | 'temperature' | 'fanSpeed';

const [colorBy, setColorByState] = useState<ColorByMode>('lineType');
useEffect(() => {
  window.dispatchEvent(new CustomEvent('apex-gcode-color-by', { detail: { mode: colorBy } }));
}, [colorBy]);
```

`Viewport3D.tsx` ajoute un écouteur `onGCodeColorBy` (même schéma que
`onGCodeLineVisibility`) qui appelle une nouvelle méthode `gcodeLayer.setColorBy(mode)`.

`GCodeFileInfo` gagne `colorByRanges: ColorByRanges` (issu de `ParsedGCode.ranges`),
propagé dans `loadText()`. **Réinitialiser `colorBy` à `'lineType'` à chaque nouveau
chargement de fichier**, cohérent avec les autres resets déjà faits dans `loadText()`.

## 5. UI

- **Emplacement** : dans `previewUpperContent`, directement au-dessus de la grille de
  légende actuelle — étendre `GCodeLineTypeLegend.tsx` en place plutôt que créer un
  wrapper séparé, puisqu'il possède déjà la logique d'affichage du tableau qui doit
  changer selon le mode.
- **Widget** : un `Popover` shadcn (Radix déjà installé) avec sections groupées
  (trop d'options pour un bouton segmenté à plat, confirmé) :
  - **Type de ligne** (toujours en premier, item unique)
  - **Cinématique** — Vitesse, Accélération, Inclinaison tête (B), Longueur de segment, Temps écoulé
  - **Extrusion** — Débit volumétrique, Pression E relative, Largeur de ligne, Hauteur de couche
  - **Vitesse relative** — Rapport vitesse/vitesse max de couche
  - **Thermique / refroidissement** — Température, Vitesse ventilateur
  - Style du déclencheur et des lignes : réutiliser le langage visuel glassmorphism
    existant (`GCodeLineTypeLegend.module.scss`, `GCodeSpeedSelector.module.scss`).
  - Chaque ligne continue affiche un mini-aperçu en dégradé (même rampe que §2) pour
    cohérence visuelle avec la légende complète.
- **Corps de la légende selon le mode** :
  - Catégoriel (type de ligne) : grille actuelle inchangée (lignes + cases à cocher).
  - Continu : remplacé par une barre de dégradé horizontale (CSS `linear-gradient`,
    points d'arrêt identiques à la rampe Babylon pour cohérence visuelle) + un résumé
    Min / Moyenne / Max sous la barre, avec l'unité du mode (mm/s, °C, %, mm, s, mm³/s),
    via un petit helper `formatColorByValue(mode, value)` co-localisé dans
    `GCodeLineTypeLegend.tsx`.
  - Les cases de visibilité par type de ligne (`lineVisibility`) restent actives et
    respectées indépendamment du mode de coloration (un mode continu masque toujours
    les types de ligne décochés) — visibilité et coloration restent deux mécanismes
    orthogonaux dans `gcodeLayer.ts` (`setLineTypeVisibility` vs `setColorBy`).

## 6. Points d'extension (16e critère et suivants)

Pour un nouveau critère **continu** : un nouveau tableau/calcul sur `ParsedGCode`, une
entrée dans `ColorByMode` + un cas dans `colorForMode()`, une ligne dans le popover, un
mapping d'unité dans `formatColorByValue`. Aucun changement requis dans
`buildBeadMesh`/`buildPathMesh`/`updateLive`/le câblage d'événements — déjà génériques
sur "une fonction qui retourne une couleur pour le mouvement k". Un nouveau critère
**catégoriel** (comme le type de ligne) demande plus de travail (son propre découpage en
maillages par valeur discrète) — asymétrie assumée, la quasi-totalité des 13 critères v1
étant continus.

## 7. Vérification

1. `pnpm dev`, ouvrir l'espace Prévisualisation.
2. Charger un vrai fichier polaire Apex (mots B) : le mode par défaut "Type de ligne"
   doit rendre à l'identique d'avant (non-régression) ; basculer sur
   Vitesse/Accélération/Inclinaison tête doit recolorer visiblement le tracé ; vitesse
   angulaire C et score de risque absents du menu (phase 2) ; hauteur de couche/largeur
   doivent produire un résultat plausible via la détection par bucket angulaire.
3. Charger un fichier cartésien de slicer réel avec commentaires température/ventilateur/
   `;WIDTH:` (ex. export Orca/Bambu) : vérifier que température et ventilateur montrent
   une vraie variation si le fichier en a une ; vérifier le cas `min === max` (couleur
   unique, pas de NaN/division par zéro).
4. Pour chacun des 13 modes v1 : la légende bascule bien de tableau+cases à
   barre-dégradé+min/moy/max, valeurs plausibles (vitesse en mm/s cohérente avec les
   vitesses d'impression connues du fichier, hauteur de couche correspondant à la valeur
   connue du fichier).
5. Cases de visibilité : toujours fonctionnelles en mode catégoriel ; le fait de
   basculer vers un mode continu puis de revenir au catégoriel restaure exactement le
   même état de visibilité (pas de reset accidentel).
6. Le segment de tête en lecture se recolore correctement en mode continu, y compris
   pendant la lecture (bouton play).
7. Basculer rapidement entre tous les modes sur un gros fichier (plusieurs milliers de
   mouvements) : confirmer qu'après le tout premier switch catégoriel↔continu (léger
   à-coup ponctuel accepté), les bascules suivantes sont instantanées, sans à-coup —
   c'est le point le plus sensible vis-à-vis de la contrainte de fluidité
   ([[feedback_fluidity]]).
8. Vérifier le nouveau réglage dans l'espace Paramètres (§3) : par défaut,
   comportement "reconstruction au premier switch" ; activé, comportement "double
   cache, zéro à-coup dès le premier switch" — confirmer que les deux fonctionnent.
9. Charger un nouveau fichier dans la même session : `colorBy` doit revenir à
   `'lineType'` automatiquement.

## Fichiers concernés

- `src/components/gcodeLayer.ts` — parseur (nouveaux M-codes, tableaux parallèles,
  détection de couche polaire/cartésienne, agrégation min/max) + calque Babylon
  (maillages fusionnés + couleurs par sommet, câblage `setColorBy`).
- `src/components/gcodeColormap.ts` (nouveau) — rampe de couleur et normalisation.
- `src/hooks/use-gcode-preview.ts` — état `colorBy`, événement `apex-gcode-color-by`,
  extension de `GCodeFileInfo`.
- `src/components/GCodeLineTypeLegend.tsx` (+ `.module.scss`) — sélecteur Popover groupé,
  corps de légende conditionnel (tableau vs dégradé+stats).
- `src/components/Viewport3D.tsx` — écouteur du nouvel événement.
- `src/components/ParamsPanel.tsx` (+ `.module.scss`) — nouveau réglage de comportement
  mémoire/fluidité pour la bascule de mode (premier contenu réel de ce panneau).

## Écart d'implémentation validé (post-approbation)

En codant §3, j'ai proposé et l'utilisateur a validé une déviation par rapport au texte du
plan : au lieu de fusionner tous les types de ligne en un maillage unique par mode continu,
le calque garde sa structure actuelle (un maillage par `LineTypeId` extrudé + un pour
travel/wipe) et lui applique des couleurs par sommet (`VertexData.colors`) à la place du
matériau plat quand un mode continu est actif. Même rendu visuel, même mécanisme de
préfixe/lecture progressive (`cumIndices`/`indexCount`) inchangé et donc beaucoup moins de
risque de régression ; seul renoncement : pas de réduction du nombre d'appels de dessin
(mineure, la scène n'a déjà qu'une dizaine de maillages).

## Décisions actées avec l'utilisateur (ne pas rouvrir sans raison)

- Température et ventilateur restent en v1 malgré le nouveau parsing M-code requis.
- Vitesse = vitesse réelle post-planification (accel/jerk), pas le F commandé brut.
- Hauteur de couche/largeur/vitesse relative disponibles aussi en polaire, via détection
  par angle C (pas de "n/a" générique).
- Pas de doublement mémoire par défaut ; reconstruction ponctuelle acceptée au premier
  switch catégoriel↔continu, avec un réglage optionnel dans Paramètres pour revenir au
  comportement zéro-à-coup si besoin sur de gros fichiers.
- Vitesse angulaire C et score de risque de défaut restent en phase 2, absents du menu
  v1 (pas de placeholder grisé).
