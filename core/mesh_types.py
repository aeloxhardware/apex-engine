"""
core/mesh_types.py

Structures de données partagées entre toutes les étapes du pipeline de slicing 4 axes.

Ce fichier ne contient AUCUNE logique de calcul lourd (pas de résolution de champ,
pas de détection de collision, pas de segmentation). Il définit uniquement les
"contrats de données" échangés entre les fichiers pipeline/step*.py.

Chaque étape du pipeline prend un ou plusieurs de ces types en entrée et produit
un ou plusieurs de ces types en sortie, ce qui permet de tester chaque étape
indépendamment avec des instances construites à la main.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Optional

import numpy as np

# ---------------------------------------------------------------------------
# Alias sémantiques
# ---------------------------------------------------------------------------

# Un point ou vecteur 3D. Alias purement documentaire (numpy ne différencie
# pas les shapes au niveau du typage statique), mais clarifie l'intention
# dans les signatures de fonctions à travers tout le projet.
Vec3 = np.ndarray  # shape (3,)


# ---------------------------------------------------------------------------
# Géométrie de base
# ---------------------------------------------------------------------------

@dataclass
class AABB:
    """Boîte englobante alignée sur les axes."""
    min_corner: Vec3
    max_corner: Vec3

    def size(self) -> Vec3:
        return self.max_corner - self.min_corner

    def center(self) -> Vec3:
        return (self.min_corner + self.max_corner) * 0.5


# ---------------------------------------------------------------------------
# Étape 1 — Mesh source ingéré
# ---------------------------------------------------------------------------

@dataclass
class MeshData:
    """
    Résultat de l'ingestion du mesh (pipeline/step1_ingest.py).

    bvh et kdtree sont typés `Any` car ils enveloppent des objets de
    bibliothèques externes (trimesh / scipy.spatial) OU leur équivalent natif
    selon `APEX_USE_CPP` - on ne redéfinit pas leur structure ici.
    """
    vertices: np.ndarray          # shape (N, 3)
    faces: np.ndarray             # shape (M, 3) - indices dans `vertices`
    face_normals: np.ndarray      # shape (M, 3)
    vertex_normals: np.ndarray    # shape (N, 3)
    is_watertight: bool
    bvh: Any = None                # core.spatial.RayIntersector (trimesh ou natif)
    kdtree: Any = None             # core.spatial.KDTree sur les centres de faces

    @property
    def face_count(self) -> int:
        return self.faces.shape[0]

    @property
    def vertex_count(self) -> int:
        return self.vertices.shape[0]


# ---------------------------------------------------------------------------
# Étape 1 — Maillage tétraédrique (domaine de résolution du champ scalaire)
# ---------------------------------------------------------------------------

@dataclass
class TetMesh:
    """
    Maillage tétraédrique volumique, domaine sur lequel le champ scalaire
    continu (étape 4) est assemblé et résolu.

    Conventions alignées sur le pipeline tétraédrique CuPy existant, pour
    rester interopérable avec ce code déjà en place.
    """
    tet_vertices: np.ndarray      # shape (P, 3)
    tetrahedra: np.ndarray        # shape (T, 4) - indices dans `tet_vertices`
    tet_adjacency: Any = None      # topologie de voisinage précalculée (scipy.sparse)

    @property
    def tet_count(self) -> int:
        return self.tetrahedra.shape[0]

    @property
    def node_count(self) -> int:
        return self.tet_vertices.shape[0]


# ---------------------------------------------------------------------------
# Étape 2 — Détection de surplombs et ponts internes
# ---------------------------------------------------------------------------

@dataclass
class OverhangMap:
    """Résultat de pipeline/step2_overhang.py, un score/flag par face du MeshData."""
    risk_score: np.ndarray             # shape (M,) float
    is_overhang: np.ndarray            # shape (M,) bool
    is_bridge: np.ndarray              # shape (M,) bool
    bridge_span_direction: np.ndarray  # shape (M, 3) - direction PCA locale (NaN si non applicable)


# ---------------------------------------------------------------------------
# Étape 3 — Segmentation en régions et graphe de dépendances
# ---------------------------------------------------------------------------

@dataclass
class Region:
    """Un sous-volume du mesh partageant une direction d'impression candidate."""
    id: int
    face_ids: np.ndarray          # indices dans MeshData.faces
    target_direction: Vec3        # direction d'impression candidate (normalisée)


@dataclass
class RegionDependencyGraph:
    """
    Graphe de dépendances physiques entre régions : une arête (a, b) signifie
    que la région `a` doit être imprimée avant la région `b`.
    """
    regions: dict[int, Region] = field(default_factory=dict)
    dependencies: list[tuple[int, int]] = field(default_factory=list)
    # shape (M,) bool, indexé comme MeshData.faces : True si la face repose
    # sur le plateau d'impression. C'est l'étape 3 qui introduit la
    # connaissance du plateau dans le pipeline (limite de conception #1) ;
    # cette information est exposée ici parce que l'étape 4 en a besoin pour
    # ancrer le champ scalaire au plateau. Ne PAS la redériver ailleurs à
    # partir d'un second paramètre `platform_z` : deux sources de vérité qui
    # divergent en silence.
    platform_contact_faces: Optional[np.ndarray] = None
    # Profondeur de chaque région dans l'arborescence de croissance depuis le
    # plateau : 0 pour une région au contact du plateau (racine), n pour une
    # région atteinte en n sauts de région à travers la matière. C'est cette
    # profondeur qui oriente `dependencies` (étape 3) — l'ordre d'impression
    # suit la façon dont la matière rejoint réellement le plateau, jamais une
    # comparaison de hauteur (une pièce pendante est plus basse que ce à quoi
    # elle est suspendue, mais doit être imprimée après).
    # `None` = non renseignée. Les profondeurs ne sont comparables qu'à
    # l'intérieur d'une même composante connexe de régions (une cavité interne
    # forme sa propre composante, sans lien de surface avec la coque
    # extérieure) : seule l'orientation des arêtes en dépend, pas un ordre
    # global entre composantes.
    region_print_depth: Optional[dict[int, int]] = None


@dataclass
class JunctionZone:
    """
    Zone à raffiner localement lors d'une seconde tétraédrisation (CLAUDE.md
    limite #48) : un cluster de faces frontière entre deux régions dont les
    `target_direction` divergent au-delà d'un seuil angulaire (voir
    `pipeline/step1b_junction_detection.py::detect_junction_zones`).
    Généralise le mécanisme `TetGen.add_region` déjà utilisé pour le
    plateau (limite #36) à des jonctions multi-région en plein volume —
    contrairement au plateau, une jonction n'a pas de plan de coupe unique,
    juste un centre et un rayon d'influence.
    """
    centroid: Vec3
    radius: float


# ---------------------------------------------------------------------------
# Étape 4 — Champ scalaire continu
# ---------------------------------------------------------------------------

@dataclass
class ScalarField:
    """
    Champ scalaire résolu sur le TetMesh. Ses isosurfaces définissent les
    couches courbes (étape 5).
    """
    values: np.ndarray            # shape (P,) - une valeur par nœud du TetMesh
    # shape (T, 3) - un gradient par TÉTRAÈDRE, pas par nœud. Le champ est
    # linéaire par élément (éléments finis P1), donc son gradient est
    # constant sur chaque tétraèdre : c'est la représentation exacte, pas
    # une approximation. C'est aussi directement ce dont l'étape 5 a besoin
    # (norme moyenne du gradient local par couche).
    gradient: Optional[np.ndarray] = None


# ---------------------------------------------------------------------------
# Étape 5/6 — Couches et trajectoire 2D locale
# ---------------------------------------------------------------------------

@dataclass
class Layer:
    """
    Une couche courbe extraite du champ scalaire (étape 5), avec sa
    trajectoire 2D remplie ultérieurement par l'étape 6.
    """
    id: int
    contour_3d: list[np.ndarray]              # liste de polylignes fermées, shape (K,3) chacune
    region_ids: list[int]
    local_frame: tuple[Vec3, Vec3, Vec3]       # (origine, axe_u, axe_v) du repère local approximatif
    toolpath_2d: Optional[list[np.ndarray]] = None  # rempli à l'étape 6 (repère local)
    color_id: Optional[int] = None             # couleur dominante de la couche (affichage UI), étape 7
    # Isovaleur du champ scalaire dont cette couche EST l'isosurface. Donnée
    # d'identité de la couche, pas un doublon dérivé : l'étape 5 la calcule
    # (`_generate_isovalues`) et la jetait jusqu'ici. L'étape 7 en a besoin
    # pour convertir une valeur de champ en indice de couche (calcul de la
    # deadline d'accessibilité d'un trou temporaire).
    isovalue: Optional[float] = None
    # Une couleur par polyligne de `toolpath_2d` (même indexation, même
    # longueur), remplie par l'étape 7. `None` sur une entrée = matière
    # INTERNE, invisible depuis l'extérieur : elle peut être imprimée avec
    # n'importe quelle couleur déjà chargée, ce qui est le principal levier
    # de réduction des purges. `color_id` (une seule couleur par couche) ne
    # suffit pas : une couche traverse en général plusieurs zones peintes.
    segment_colors: Optional[list[Optional[int]]] = None
    # Score de non-monotonie locale du champ scalaire (CLAUDE.md limite
    # #47) : `None` sauf si `LayerParams.detect_junction_variance=True`.
    # Une valeur non-nulle signale que cette couche forme, avec une couche
    # Z-adjacente, une paire à écart Z anormalement petit ET écart
    # d'isovaleur anormalement grand — signature d'une jonction
    # multi-région où φ n'est pas monotone en Z (voir
    # `_compute_junction_variance_scores`). Champ pur, jamais utilisé pour
    # modifier `contour_3d`/`toolpath_2d` : un diagnostic, pas un fix.
    junction_variance_score: Optional[float] = None
    # Le dépôt de cette couche se fait dans le sens OPPOSÉ à la normale de son
    # repère local, c'est-à-dire le long de `-∇φ`.
    #
    # Le sens de dépôt n'est pas une propriété du champ mais de l'ORDRE dans
    # lequel on parcourt ses isosurfaces. Tant que l'ordre était « φ
    # croissant » (tout l'historique du projet), déposer selon `+∇φ` était vrai
    # par construction, et l'étape 8 pouvait le supposer sans le dire. Un ordre
    # qui redescend le long d'une pièce suspendue — une croissance de front — rend cette supposition fausse : on colle
    # alors de la matière SOUS ce qui existe déjà, ce qui demande à la tête de
    # pointer vers le haut (B ≈ 180° sur cette machine, qui couvre 0-185°).
    #
    # Renseigné par l'étape 5, qui seule connaît l'ordre final ; consommé par
    # l'étape 8, qui inverse alors la direction désirée. Le champ, les
    # isosurfaces et le repère local restent STRICTEMENT inchangés : seule
    # l'orientation de la tête bascule.
    deposition_reversed: bool = False


# ---------------------------------------------------------------------------
# Étape 7 — Carte de couleur et planification d'ordonnancement
# ---------------------------------------------------------------------------

@dataclass
class ColorMap:
    """Carte de couleur peinte par l'utilisateur sur le mesh source."""
    per_face_color: dict[int, int]    # face_id -> color_id
    merge_tolerance: float = 0.0       # lissage/regroupement avant ordonnancement


@dataclass
class PrintPass:
    """
    Un passage d'impression dans la séquence globale (étape 7).

    `segments` est un POINTEUR vers la géométrie, jamais une copie : les
    coordonnées restent dans `Layer.toolpath_2d` (source unique de vérité),
    une passe ne décrit que "quoi imprimer quand". Un couple (layer_id,
    index) désigne la polyligne `layer.toolpath_2d[index]` de la couche
    d'identifiant `layer_id`.

    `layer_range` seul ne suffisait pas à décrire une passe : une passe de
    couleur ne couvre qu'une PARTIE de chaque couche traversée (les segments
    de cette couleur), pas des couches entières contiguës — sans `segments`,
    l'étape 8 ne saurait pas quelles polylignes émettre.
    """
    segments: list[tuple[int, int]]    # (layer_id, index dans Layer.toolpath_2d)
    color_id: Optional[int]            # None = matière invisible, couleur libre
    is_hole_fill: bool = False         # True = passe de comblement d'un trou temporaire
    layer_range: Optional[tuple[int, int]] = None  # dérivé de `segments`, pour l'affichage UI
    deferred_regions: list[int] = field(default_factory=list)


@dataclass
class PrintSchedule:
    """
    Séquence globale d'impression, sortie de pipeline/step7_scheduling.py.

    L'ordre d'exécution est EXACTEMENT l'ordre de la liste `passes` :
    l'étape 8 la consomme séquentiellement, sans jamais réordonner ni
    re-dériver la logique d'ordonnancement (topologique + earliest-deadline
    -first), entièrement résolue et gelée ici.
    """
    passes: list[PrintPass] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Profil machine — cinématique générique à N axes de rotation composés
# ---------------------------------------------------------------------------

@dataclass
class RotationAxisSpec:
    """
    Un DOF de rotation : axe + plage angulaire, exprimé dans le repère du
    DOF parent (le repère déjà tourné par tous les axes précédents de la
    liste — même convention que l'ordre d'application matricielle).

    `angle_range_deg` peut être `(0, 360)` pour un tour complet, une plage
    restreinte (ex: `(0, 180)` pour la tête de la machine actuelle),
    ou `(v, v)` pour un DOF figé à une seule valeur — dégénérescence utile
    qui permet de décrire un axe verrouillé sans le retirer de la liste.
    """
    axis: Vec3
    angle_range_deg: tuple[float, float]


@dataclass
class MachineProfile:
    """
    Cinématique de la machine, décrite comme une liste ORDONNÉE et
    PARAMÉTRABLE d'axes de rotation — PAS un modèle figé sur une seule
    machine. La machine actuelle (plateau C + tête B∈[0°,180°])
    n'est qu'UN cas d'usage à valider parmi d'autres cinématiques futures
    (un seul DOF, deux DOF différemment agencés, plus de deux DOF) ; voir
    CLAUDE.md limite de conception #18.

    - liste vide : orientation libre (5 axes), aucune contrainte de rotation.
    - 1 élément, plage `(0,360)` : machine 4 axes classique à 1 DOF.
    - 2 éléments comme la machine actuelle : un exemple parmi d'autres
      configurations possibles, pas la référence du type lui-même.

    Commande firmware de la machine actuelle (gabarit symbolique,
    pas de valeurs numériques fictives) : `G00 Y<C> X<R> Z<Z> B<B> F<feedrate>`.
    """
    rotation_axes: list["RotationAxisSpec"] = field(default_factory=list)
    platform_origin: Vec3 = None     # centre de rotation du premier axe (si applicable) —
                                      # aussi le point de singularité où l'angle du premier
                                      # axe est indéterminé (limite #19)
    nozzle_model: Any = None         # maillage 3D tête (buse + porte-buse + carriage proche)
    # Identifiant du modele de tete (dossier de `res/3d/toolheads/<id>/<id>.glb`), lie au profil
    # d'imprimante (R-020). Jamais un chemin : il traverse le contrat, le moteur le resout chez lui.
    toolhead_id: str = "aetheris-4"
    nozzle_diameter: float = 0.4
    max_angular_velocity: float = 0.0     # rad/s
    max_angular_acceleration: float = 0.0
    build_volume: Optional["AABB"] = None
    z_clearance_min: float = 0.0
    # Rayon du plateau (mm). Machine actuelle : 300 mm de diametre
    # (utilisateur, 2026-09-29). None = inconnu : step9 retombe sur son disque
    # par defaut et aucun refus de taille n'est applique.
    platform_radius: Optional[float] = None
    # Axes dont l'angle est IMPOSÉ PAR LA POSITION du point (limite #59) :
    # quand la tête ne se déplace que dans un demi-plan (machine actuelle : R et
    # Z seulement, pas d'axe Y), le plateau doit tourner pour amener le point
    # dans ce demi-plan — C n'est plus un degré d'orientation. Seul le premier
    # axe (le plus externe) peut l'être : `[]` = orientation libre (machine
    # après upgrade, avec un axe linéaire de plus), `[0]` = machine actuelle.
    position_coupled_axes: list[int] = field(default_factory=list)
    # Direction +R (vers laquelle la tête se déplace en s'écartant de l'axe)
    # quand tous les angles sont nuls, dans le repère de la pièce. Requise si
    # `position_coupled_axes` est non vide : l'angle couplé est celui qui
    # oriente cette direction vers le point.
    radial_direction: Optional[Vec3] = None


# ---------------------------------------------------------------------------
# Étape 8 — Trajectoire machine (R, C, Z, B)
# ---------------------------------------------------------------------------

@dataclass
class ToolpathPoint:
    """
    Un point de trajectoire dans la cinématique GÉNÉRIQUE de la machine —
    `axis_angles` porte un angle par élément de `MachineProfile.rotation_axes`,
    même ordre, quel que soit N (0, 1, 2, ...). Remplace les anciens champs
    `radius`/`platform_angle`/`head_angle`, qui figeaient la structure sur
    exactement 2 DOF nommés "plateau"/"tête" — une hypothèse invalidée par la
    généralisation à N axes composés (voir CLAUDE.md limite #18).

    La traduction finale vers une commande firmware précise (ordre X/Y/Z/B,
    noms d'axes) est délibérément HORS de ce type : c'est la responsabilité
    d'un futur module d'export (`io/gcode_export.py`), spécifique à une
    machine donnée, jamais de step8 ni de ce contrat de données générique.
    """
    position: Vec3                  # cartésien machine, pour vérification/debug (étape 9)
    axis_angles: list[float]        # un angle en degrés par MachineProfile.rotation_axes, même ordre
    extrusion: float
    feedrate: float
    color_id: Optional[int] = None
    is_skin: bool = False           # couche de lissage non-planaire (epaisseur propre a l'export, `skin_thickness`)


@dataclass
class Step8Diagnostic:
    """
    Signal "orientation désirée inatteignable dans les plages d'axes
    disponibles à cette position" — donnée pure remontée par l'étape 8 vers
    l'étape 10 (ajustement), qui n'a besoin d'aucune connaissance de la
    résolution IK interne à l'étape 8 pour la consommer.
    """
    point_index: int                # index dans ToolpathMachine.points
    layer_id: int
    segment_index: int              # index dans Layer.toolpath_2d de la couche layer_id
    residual_angle_error_deg: float
    position: Vec3


@dataclass
class ToolpathMachine:
    """Trajectoire machine finale, prête pour simulation (étape 9) puis export G-code."""
    points: list[ToolpathPoint] = field(default_factory=list)
    # Couche (`Layer.id`) de chaque point du CORPS, issue de l'étape 8 ; `None` hors étape 8 (tests, trajectoire
    # relue). Sert à insérer le lissage non-planaire après la couche qui le porte (`pipeline/skin/inject.py`).
    layer_ids: Optional[np.ndarray] = None


@dataclass
class TravelRoutingDiagnostic:
    """
    Signal "trajet à vide routé en lift-and-traverse, mais toujours en
    collision avec la matière déjà imprimée même à la hauteur de lift
    maximale testée" — donnée pure remontée par
    `pipeline/step8b_travel_routing.py` (CLAUDE.md limite #46).

    Le mécanisme de routage n'est PAS un vrai path-planning 3D (voir la
    docstring du module) : c'est un palliatif lever/traverser/redescendre,
    documenté comme une limite du slicer à améliorer plus tard, pas une
    solution définitive. Ce diagnostic rend visible les cas où même ce
    palliatif échoue, plutôt que de les laisser silencieusement retomber
    sur le trajet en ligne droite d'origine.
    """
    point_index: int                # index dans ToolpathMachine.points (point de destination du trajet)
    travel_length: float            # longueur du trajet à vide d'origine (ligne droite), mm
    max_lift_tested: float          # hauteur de lift la plus grande essayée avant abandon, mm
    still_colliding: bool           # True si le meilleur candidat teste reste en collision


# ---------------------------------------------------------------------------
# Étape 9/10 — Simulation de collision et rapport
# ---------------------------------------------------------------------------

@dataclass
class CollisionEvent:
    point_index: int                # index dans ToolpathMachine.points
    obstacle_type: str              # "printed_part" | "chuck" | "platform"
    penetration_depth: float


@dataclass
class CollisionReport:
    """
    `chuck_checked`/`nozzle_model_checked` sont des signaux EXPLICITES
    d'absence de verification, jamais un silence - meme famille de motif
    que `Step8Diagnostic` (donnee pure, pas d'exception, pas de skip
    silencieux). Sans eux, une liste `events` vide serait ambigue entre
    "verifie, aucune collision" et "jamais verifie du tout".
    """
    events: list[CollisionEvent] = field(default_factory=list)
    chuck_checked: bool = False          # False = aucune geometrie de mandrin disponible - evenements "chuck" non exhaustifs
    nozzle_model_checked: bool = False   # False = nozzle_model absent/non chargeable - AUCUN test n'a eu lieu

    @property
    def has_collisions(self) -> bool:
        return len(self.events) > 0


@dataclass
class UnsupportedDepositEvent:
    """
    Un point de dépôt qui n'avait AUCUNE matière déjà imprimée en dessous au
    moment de son extrusion (étape 9b).

    Distinct d'un surplomb géométrique au sens de l'étape 2 : celle-ci ne
    connaît que la forme du modèle, pas l'ORDRE réel d'impression. Une face
    parfaitement supportée dans la pièce finie peut très bien être déposée
    dans le vide si la matière qui la porte est imprimée après elle.

    `has_future_support` sépare les deux natures de défaut, et c'est la
    donnée qui décide quoi corriger :
      - `True`  : la matière porteuse existe, mais elle est déposée plus tard
                  dans la trajectoire — défaut d'ORDRE, corrigible par un
                  réordonnancement, sans rien changer à la géométrie.
      - `False` : aucun support VERTICAL à aucun moment de la trajectoire —
                  candidat surplomb, à départager par `supported_by_tilt`.

    `supported_by_tilt` ajoute une troisième catégorie, nécessaire pour une
    machine dont la tête/le plateau peuvent s'incliner (mesuré sur
    res/test_models/death.stl : 194/194 points classés "surplomb réel" en
    ne testant QUE la verticale avaient en fait un support dès qu'on teste
    une direction inclinée atteignable par la cinématique réelle) :
      - `True`  : `has_future_support` est `False` (aucun support vertical),
                  mais un support existe le long d'une direction atteignable
                  par `machine_profile.rotation_axes` — la machine peut
                  déposer ce point sans surplomb en orientant la tête,
                  quel que soit l'ordre d'impression. Ni défaut d'ordre, ni
                  surplomb réel.
      - `False` : aucun support, dans AUCUNE direction testée — surplomb
                  géométrique réel, qu'aucun ordre ni aucune orientation de
                  tête ne peut résoudre.
    """
    point_index: int                 # index dans ToolpathMachine.points
    position: Vec3
    has_future_support: bool
    # Index du premier point (dans l'ordre de trajectoire) dont le dépôt
    # aurait supporté celui-ci. None si has_future_support est False.
    supporting_point_index: Optional[int] = None
    supported_by_tilt: bool = False


@dataclass
class SupportCheckReport:
    """
    Sortie de l'étape 9b. `checked` est un signal EXPLICITE d'absence de
    vérification (même famille que `CollisionReport.nozzle_model_checked` et
    `AdjustmentReport.converged`) : sans lui, une liste `events` vide serait
    ambiguë entre « vérifié, tout est supporté » et « jamais vérifié ».
    """
    events: list[UnsupportedDepositEvent] = field(default_factory=list)
    checked: bool = False
    # Dénominateur de `unsupported_fraction` : nombre de points RÉELLEMENT
    # testés (points en extrusion retenus par l'échantillonnage), jamais le
    # nombre total de points de la trajectoire — sinon la fraction serait
    # diluée par les trajets à vide, qui ne déposent rien.
    tested_extruded_points: int = 0

    @property
    def unsupported_fraction(self) -> float:
        if self.tested_extruded_points == 0:
            return 0.0
        return len(self.events) / self.tested_extruded_points

    @property
    def ordering_defect_count(self) -> int:
        """Défauts corrigibles par réordonnancement (le support existe, trop tard)."""
        return sum(1 for event in self.events if event.has_future_support)

    @property
    def true_overhang_count(self) -> int:
        """
        Surplombs réels : aucun ordre d'impression ni aucune orientation de
        tête atteignable ne peut les supporter. Exclut explicitement
        `supported_by_tilt` — sans quoi ce compte resterait pollué par des
        points que la machine sait déjà déposer correctement (mesuré :
        194/194 des surplombs "réels" de death.stl avant ce correctif).
        """
        return sum(
            1 for event in self.events if not event.has_future_support and not event.supported_by_tilt
        )

    @property
    def tilt_resolved_count(self) -> int:
        """Points sans support vertical, mais résolus par une orientation de tête atteignable."""
        return sum(1 for event in self.events if event.supported_by_tilt)


@dataclass
class AdjustmentLogEntry:
    """Une entrée du journal d'ajustement post-collision (étape 10)."""
    point_index: int
    method: str                     # "local_reorientation" | "region_resegmentation" | "fallback_support" | "orientation_deviation"
    detail: str = ""


@dataclass
class ResegmentationRequest:
    """
    Demande de ré-segmentation ciblée émise par l'étape 10 vers son
    appelant (orchestrateur), quand la réorientation locale ne suffit pas.

    L'étape 10 décide QUOI (quel point, quelle plage d'angles exclure) ;
    l'appelant décide COMMENT re-jouer la chaîne step3 -> step8. Ce
    découpage est délibéré : une ré-segmentation change `target_direction`
    d'une région, donc le terme de gradient du système STATIQUE de step4,
    donc le champ scalaire (résolu globalement par FEM), donc les couches,
    les trajectoires, l'ordonnancement et les angles — il n'existe aucune
    ré-exécution réellement locale. Faire porter cette chaîne par l'étape
    10 en ferait l'orchestrateur, rôle qui revient à `main.py`.
    """
    point_index: int
    position: Vec3
    # Une entrée par axe de MachineProfile.rotation_axes, même ordre.
    # None = aucune exclusion pour cet axe.
    excluded_axis_ranges_deg: list[Optional[tuple[float, float]]] = field(default_factory=list)
    reason: str = ""


@dataclass
class AdjustmentReport:
    """
    Sortie de l'étape 10. `converged` est un signal EXPLICITE, jamais un
    silence (même discipline que `CollisionReport.chuck_checked` et
    `Step8Diagnostic`) : une liste `unresolved_collisions` vide avec
    `converged=False` serait ambiguë sans lui.
    """
    entries: list[AdjustmentLogEntry] = field(default_factory=list)
    converged: bool = True
    unresolved_collisions: list["CollisionEvent"] = field(default_factory=list)
    iterations_used: int = 0