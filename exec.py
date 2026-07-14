import sys
import os

# On force Python à regarder dans le dossier 'src' pour trouver tes fichiers
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

# Imports directs, sans points, sans package
import config
import core

print("IMPORTING_GRAPHIC_PILOTS")

import cupy as cp
# Test rapide pour vérifier que le GPU est bien détecté
try:
    print(f"STATUS:CUDA_VERSION_{cp.cuda.runtime.runtimeGetVersion()}")
except Exception as e:
    print(f"STATUS:CUDA_ERROR_{e}")

print("PLEASE_WAIT")

def main():
    print("STATUS:START")

    print("STATUS:SETUP_START")
    core.SETUP()
    print("STATUS:SETUP_DONE")

    print("STATUS:SETUP_END_START")
    core.SETUP_END()
    print("STATUS:SETUP_END_DONE")
    
    print("STATUS:SHOW_GRAPHS")
    # core.SHOW_GRAPHS
    print("STATUS:SHOW_GRAPHS_DONE")

    print("STATUS:TET_MESH_DEFORM_START")
    core.START_TET_MESH_DEFORM()
    print("STATUS:TET_MESH_DEFORM_DONE")

    print("STATUS:PATH_PROCESSING_SETUP_START")
    core.PATH_PROCESSING_SETUP()
    print("STATUS:PATH_PROCESSING_SETUP_DONE")

    print("STATUS:N_MATRIX_SETUP_START")
    core.N_MATRIX_SETUP()
    print("STATUS:N_MATRIX_SETUP_DONE")

    print("STATUS:EXEC_DEFORMATIONS_START")
    core.EXEC_DEFORMATIONS()
    print("STATUS:EXEC_DEFORMATIONS_DONE")

    print("STATUS:UT_COPY_START")
    core.UT_COPY()
    print("STATUS:UT_COPY_DONE")

    print("STATUS:SAVE_STL_START")
    core.SAVE_STL()
    print("STATUS:SAVE_STL_DONE")

    print("STATUS:SAVE_DEFORM_START")
    core.SAVE_DEFORM()
    print("STATUS:SAVE_DEFORM_DONE")

    # --- POINT D'ARRÊT INTERMÉDIAIRE ---
    print("STATUS:WAITING_FOR_USER_INPUT")
    while True:
        user_input = input("--- PAUSE : Tapez 'run' pour continuer ou 'abort' pour annuler : ").strip().lower()
        
        if user_input == 'run':
            print("STATUS:RESUMING_EXECUTION")
            break
        elif user_input == 'abort':
            print("STATUS:USER_ABORT_REQUESTED")
            print("--- Annulation du programme. Fermeture sécurisée. ---")
            exit() # C'est la commande magique pour tuer ton programme proprement
        else:
            print("Commande inconnue. Tapez 'run' pour continuer, sinon 'abort' pour vider la mémoire à accès aléatoire et annuler l'opération.")

    print("STATUS:OPEN_DEFORM_START")
    core.OPEN_DEFORM()
    # ------------------------------------

    print("STATUS:OPEN_DEFORM_DONE")

    print("STATUS:EXEC_CALC_TET_ATTR_START")
    core.EXEC_CALC_TET_ATTR()
    print("STATUS:EXEC_CALC_TET_ATTR_DONE")

    print("STATUS:PREPARE_ROTATE_START")
    core.PREPARE_ROTATE()
    print("STATUS:PREPARE_ROTATE_DONE")

    print("STATUS:CALCULATE_VECTORIZED_ZSQUISH_START")
    core.CALCULATE_VECTORIZED_ZSQUISH()
    print("STATUS:CALCULATE_VECTORIZED_ZSQUISH_DONE")

    print("STATUS:OPEN_READ_GCODE_START")
    core.OPEN_READ_GCODE()
    print("STATUS:OPEN_READ_GCODE_DONE")

    print("STATUS:BRC_INTERPOL_PREPROCESS_START")
    core.BRC_INTERPOL_PREPROCESS()
    print("STATUS:BRC_INTERPOL_PREPROCESS_DONE")

    print("STATUS:ASSEMBLY_LOOP_START")
    core.ASSEMBLY_LOOP()
    print("STATUS:ASSEMBLY_LOOP_DONE")

    print("STATUS:SAVE_GCODE_START")
    core.SAVE_GCODE()
    print("STATUS:SAVE_GCODE_DONE")

    print("STATUS:PREV_TRAJECTORY_DOTS_START")
    core.PREV_TRAJECTORY_DOTS()
    print("STATUS:PREV_TRAJECTORY_DOTS_DONE")

    print("STATUS:PREV_TRAJECTORY_PATH_START")
    core.PREV_TRAJECTORY_PATH()
    print("STATUS:PREV_TRAJECTORY_PATH_DONE")

    print("STATUS:SUCCESS")
    
if __name__ == "__main__":
    main()