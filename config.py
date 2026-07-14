import numpy as np
import os

# --- Paramètres de calcul ---

NEIGHBOUR_LOSS_WEIGHT = 50
MAX_OVERHANG = 5
ROTATION_MULTIPLIER = 3
SET_INITIAL_ROTATION_TO_ZERO = True
INITIAL_ROTATION_FIELD_SMOOTHING = 30
MAX_POS_ROTATION = np.deg2rad(3600)
MAX_NEG_ROTATION = np.deg2rad(-3600)
ITERATIONS = 30
DEFORM_ITERATIONS = 1000

# --- Paramètres de posttraitement ---
SEG_SIZE = 0.6 # mm
MAX_ROTATION = 30 # degrees
MIN_ROTATION = -130 # degrees
NOZZLE_OFFSET = 42 # mm actuallt 41.5

# --- Paramètres système et export ---

SAVE_GIF = False
STEEP_OVERHANG_COMPENSATION = True

# Définit le dossier 'src' comme référence absolue
BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Construit le chemin complet vers ton modèle
model_path = os.path.join(BASE_DIR, 'input_model', 'model.stl')

# On construit les chemins absolus pour TOUS tes dossiers
INPUT_MODEL_DIR = os.path.join(BASE_DIR, 'input_model')
OUTPUT_MODEL_DIR = os.path.join(BASE_DIR, 'output_models')
INPUT_GCODE_DIR = os.path.join(BASE_DIR, 'input_gcode')
OUTPUT_GCODE_DIR = os.path.join(BASE_DIR, 'output_gcode')
PKL_DIR = os.path.join(BASE_DIR, 'pickle_files')

WORKSPACE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ux", "plate")
MODELS_DIR = os.path.join(WORKSPACE_DIR, "models")
ARRANGE_JSON = os.path.join(WORKSPACE_DIR, "arrange.json")
 
# Assurer l'existence des répertoires de stockage
os.makedirs(MODELS_DIR, exist_ok=True)

# On s'assure que les dossiers existent (très utile pour la prod)
for folder in [OUTPUT_MODEL_DIR, OUTPUT_GCODE_DIR]:
    if not os.path.exists(folder):
        os.makedirs(folder)