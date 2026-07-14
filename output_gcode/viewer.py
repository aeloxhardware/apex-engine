import sys
import os
import re
import json
import numpy as np
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QPushButton, QFileDialog, QLabel)

# Intégration de Matplotlib dans PyQt6 (Alternative robuste et native à Plotly dans PyQt)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class GCodeVisualizerApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Prévisualiseur de G-Code Cylindrique")
        self.setGeometry(100, 100, 900, 800)

        # Widget principal et Layout
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        self.layout = QVBoxLayout(self.main_widget)

        # Label d'information / Flux JSON simulé pour la progression
        self.label_status = QLabel("Aucun fichier chargé. Veuillez sélectionner un fichier G-code.")
        self.label_status.setStyleSheet("font-weight: bold; font-size: 14px; margin: 10px;")
        self.layout.addWidget(self.label_status)

        # Bouton pour ouvrir la boîte de dialogue native
        self.btn_open = QPushButton("Choisir et analyser un fichier G-Code")
        self.btn_open.setStyleSheet("padding: 10px; font-size: 14px; background-color: #2b5c8f; color: white; font-weight: bold;")
        self.btn_open.clicked.connect(self.open_file_dialog)
        self.layout.addWidget(self.btn_open)

        # Création de la figure Matplotlib intégrée
        self.canvas = FigureCanvas(Figure(figsize=(7, 7)))
        self.layout.addWidget(self.canvas)
        
        # Initialisation de l'axe 3D
        self.ax = self.canvas.figure.add_subplot(111, projection='3d')
        self.init_empty_plot()

    def init_empty_plot(self):
        """Initialise un graphique vide avec des axes visibles."""
        self.ax.clear()
        self.ax.set_xlabel('X (mm)')
        self.ax.set_ylabel('Y (mm)')
        self.ax.set_zlabel('Z (mm)')
        self.ax.grid(True)
        self.canvas.draw()

    def open_file_dialog(self):
        """Ouvre la boîte de dialogue native pour choisir le fichier."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Sélectionner le fichier G-Code",
            "",
            "Fichiers G-Code (*.gcode *.nc *.txt *.tap);;Tous les fichiers (*)"
        )

        if file_path:
            filename = os.path.basename(file_path)
            self.label_status.setText(f"Lecture réelle de : {filename} ...")
            # Force l'interface graphique à se mettre à jour avant le calcul
            QApplication.processEvents() 
            
            # Lance l'analyse stricte et le tracé
            self.process_and_plot_gcode(file_path)

    def process_and_plot_gcode(self, file_path):
        """Analyse le VRAI fichier sélectionné, extrait uniquement le G01 modal, 

        convertit le cylindrique en cartésien et l'affiche sans fausses trajectoires.
        """
        self.ax.clear()

        # Émission JSON de progression initiale au stdout comme demandé dans votre modèle
        sys.stdout.write(json.dumps({"type": "progress", "task": "Initialisation du fichier", "current": 0, "total": 2}) + "\n")
        sys.stdout.flush()

        # Mémoire de l'état de la machine (axes modaux)
        current_r = 0.0      # stocké dans X
        current_theta = 0.0  # stocké dans Y (en degrés)
        current_z = 0.0      # stocké dans Z
        current_b = 0.0      # stocké dans B
        
        active_command = None # Commande G-code active (G00, G01, etc.)
        
        cartesian_points = []

        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()

            total_lines = len(lines)
            
            for idx, line in enumerate(lines):
                # Nettoyage de la ligne et exclusion des commentaires
                line = line.strip().split(';')[0].split('(')[0]
                if not line:
                    continue

                # Détection du changement de commande (G0, G1, G00, G01)
                # Si absent, la commande précédente reste active (comportement modal du G-code)
                match_cmd = re.search(r'\b(G00|G01|G0|G1)\b', line)
                if match_cmd:
                    cmd = match_cmd.group(1)
                    if cmd in ["G0", "G00"]:
                        active_command = "G00"
                    elif cmd in ["G1", "G01"]:
                        active_command = "G01"

                # Extraction des valeurs des axes si présentes sur la ligne
                match_x = re.search(r'X\s*([-+]?\d*\.\d+|\d+)', line)
                match_y = re.search(r'Y\s*([-+]?\d*\.\d+|\d+)', line)
                match_z = re.search(r'Z\s*([-+]?\d*\.\d+|\d+)', line)
                match_b = re.search(r'B\s*([-+]?\d*\.\d+|\d+)', line)

                # Mise à jour des coordonnées absolues (si non spécifié, conserve l'ancienne valeur)
                if match_x: current_r = float(match_x.group(1))
                if match_y: current_theta = float(match_y.group(1))
                if match_z: current_z = float(match_z.group(1))
                if match_b: current_b = float(match_b.group(1))

                # Traitement exclusif de G01 (Adieu les G00 !)
                # On ne prend un point que si un mouvement a eu lieu sous la commande active G01
                if active_command == "G01" and (match_x or match_y or match_z or match_b or (match_cmd and active_command == "G01")):
                    # Conversion Cylindrique (R, Theta) -> Cartésien pour affichage correct
                    theta_rad = np.radians(current_theta)
                    visu_x = current_r * np.cos(theta_rad)
                    visu_y = current_r * np.sin(theta_rad)
                    visu_z = current_z
                    
                    cartesian_points.append([visu_x, visu_y, visu_z])

            # Émission JSON de fin de tâche
            sys.stdout.write(json.dumps({"type": "progress", "task": "Extraction et injection", "current": 2, "total": 2}) + "\n")
            sys.stdout.flush()

            cartesian_points = np.array(cartesian_points)

            if len(cartesian_points) == 0:
                self.label_status.setText("Aucune trajectoire G01 (Travail) trouvée dans ce fichier.")
                self.init_empty_plot()
                return

            # Séparation des axes cartésiens calculés
            x = cartesian_points[:, 0]
            y = cartesian_points[:, 1]
            z = cartesian_points[:, 2]

            # Tracé unique équivalent à votre configuration plotly (Lignes + Marqueurs en bleu)
            self.ax.plot(x, y, z, 
                         color="blue", 
                         linestyle='-', 
                         linewidth=1, 
                         marker='o', 
                         markersize=2, 
                         alpha=0.5, 
                         label="G01 (Travail)")

            # --- AJUSTEMENT DES PROPORTIONS SÉCURISÉ EN CAS D'AXE PLAT (Ex: Z constant) ---
            max_range = np.array([np.ptp(x), np.ptp(y), np.ptp(z)])
            max_range[max_range == 0] = 1.0 
            if np.max(max_range) > 0:
                aspect_ratio = max_range / np.max(max_range)
                self.ax.set_box_aspect((aspect_ratio[0], aspect_ratio[1], aspect_ratio[2]))

            # Habillage final des axes
            self.ax.set_xlabel('X (Cartésien - mm)')
            self.ax.set_ylabel('Y (Cartésien - mm)')
            self.ax.set_zlabel('Z (mm)')
            self.ax.legend(loc="upper left")
            self.ax.grid(True)
            
            filename = os.path.basename(file_path)
            self.label_status.setText(f"Affichage réussi : {filename} ({len(cartesian_points)} points G01)")

        except Exception as e:
            self.label_status.setText(f"Erreur lors de la lecture : {str(e)}")
            self.init_empty_plot()

        # Rafraîchissement requis de la vue PyQt
        self.canvas.draw()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = GCodeVisualizerApp()
    window.show()
    sys.exit(app.exec())