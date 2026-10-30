"""
Point d'acces UNIQUE a toutes les dataclasses de parametres du moteur.

Les classes restent definies dans leur etape (`pipeline/stepN_*.py`) : les
deplacer casserait les imports existants, ce que CLAUDE.md interdit
(« centraliser sans casser les imports existants »). Ce module les reexporte,
et `ROOTS` designe les racines du contrat, parcourues pour le schema JSON et l'en-tete C++
(`scripts/gen_params_hpp.py`).
"""
from core.mesh_types import AABB, ColorMap, MachineProfile, RotationAxisSpec  # noqa: F401
from file_io.gcode_export import GcodeParams  # noqa: F401
from main import PipelineParams  # noqa: F401
from pipeline.step2_overhang import OverhangParams  # noqa: F401
from pipeline.step3_segmentation import SegmentationParams  # noqa: F401
from pipeline.step4_field_solver import FieldSolverParams  # noqa: F401
from pipeline.step5_layers import LayerParams  # noqa: F401
from pipeline.step6_toolpath2d import Step6Params  # noqa: F401
from pipeline.step7_scheduling import HoleCostModel, PurgeCostModel, Step7Params  # noqa: F401
from pipeline.step8_machine_toolpath import Step8Params  # noqa: F401
from pipeline.step8b_travel_routing import TravelRoutingParams  # noqa: F401
from pipeline.step9_collision_sim import CollisionSimParams  # noqa: F401
from pipeline.step9b_support_check import SupportCheckParams  # noqa: F401
from pipeline.step10_adjustment import Step10Params  # noqa: F401

import importlib as _importlib
AeolisSliceParams = _importlib.import_module("nonplanar.slicer").SliceParams  # Element Aeolis (R-024)

# Racines du contrat UI : ce que l'interface envoie au moteur.
ROOTS = {"PipelineParams": PipelineParams, "GcodeParams": GcodeParams, "AeolisSliceParams": AeolisSliceParams}
