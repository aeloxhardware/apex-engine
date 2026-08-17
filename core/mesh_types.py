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

    bvh et kdtree sont typés `Any` car ce sont des objets de bibliothèques
    externes (trimesh / scipy.spatial) - on ne redéfinit pas leur structure ici.
    """
    vertices: np.ndarray          # shape (N, 3)
    faces: np.ndarray             # shape (M, 3) - indices dans `vertices`
    face_normals: np.ndarray      # shape (M, 3)
    vertex_normals: np.ndarray    # shape (N, 3)
    is_watertight: bool
    bvh: Any = None                # structure d'accélération pour ray-casting
    kdtree: Any = None             # scipy.spatial.cKDTree sur les centres de faces

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
# Contraintes utilisateur (points d'attache, esthétique, structurel)
# ---------------------------------------------------------------------------

class ConstraintType(Enum):
    ATTACH_POINT = "attach_point"
    AESTHETIC = "aesthetic"
    STRUCTURAL = "structural"


class ConstraintFalloff(Enum):
    LINEAR = "linear"
    GAUSSIAN = "gaussian"
    CUBIC = "cubic"


@dataclass
class Constraint:
    """
    Contrainte injectée dans le solveur de champ scalaire (étape 4).
    Envoyée dynamiquement depuis l'UI TypeScript, une par point d'attache
    ou par zone esthétique/structurelle.
    """
    type: ConstraintType
    position: Vec3
    target_field_value: Optional[float] = None
    target_normal: Optional[Vec3] = None
    influence_radius: float = 1.0
    falloff: ConstraintFalloff = ConstraintFalloff.GAUSSIAN
    weight: float = 1.0


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
    nozzle_diameter: float = 0.4
    max_angular_velocity: float = 0.0     # rad/s
    max_angular_acceleration: float = 0.0
    build_volume: Optional["AABB"] = None
    z_clearance_min: float = 0.0


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
class AdjustmentLogEntry:
    """Une entrée du journal d'ajustement post-collision (étape 10)."""
    point_index: int
    method: str                     # "local_reorientation" | "region_resegmentation" | "fallback_support"
    detail: str = ""