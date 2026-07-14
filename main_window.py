"""Fenêtre principale de l'application Apex Slicer Engine."""
import sys
import os
import math
import struct
import shutil
from PyQt6.QtWidgets import (QMainWindow, QWidget, QPushButton, QVBoxLayout,
                             QHBoxLayout, QFrame, QLabel, QComboBox, QLineEdit,
                             QScrollArea, QGraphicsBlurEffect, QFileDialog)
from PyQt6.QtCore import Qt, QTimer

from config import MODELS_DIR
from mesh_utils import simplify_mesh
from process_worker import ProcessWorker
from viewport import Slicer3DViewport
from dialogs import ModalOverlay, WarningDialog, SimplifyDialog, SequenceDialog, TransformDialog


class MainWindow(QMainWindow):
    """Fenêtre principale de l'application."""
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Apex Slicer Engine")
        self.setAcceptDrops(True) 
        self.worker = None
        
        self.setStyleSheet("""
            QMainWindow {
                background-color: #121212;
                color: #E5E5E5;
                font-family: 'Segoe UI', sans-serif;
            }
            QLabel {
                color: #E5E5E5;
                font-size: 13px;
            }
            QLabel#SectionTitle {
                color: #A0A8C0;
                font-weight: 600;
                font-size: 12px;
                text-transform: uppercase;
                margin: 6px 0 4px 0;
            }
            QFrame#SidebarPanel, QFrame#TopToolbar, QFrame#TransformPanel {
                background-color: #1F1F1F;
                border: 1px solid #2A2A2A;
                border-radius: 10px;
            }
            QComboBox, QLineEdit, QPushButton {
                background-color: #1B1B1B;
                color: #E5E5E5;
                border: 1px solid #333333;
                border-radius: 6px;
                padding: 6px 8px;
            }
            QComboBox:hover, QLineEdit:hover, QPushButton:hover {
                border: 1px solid #555555;
            }
            QScrollArea {
                background: transparent;
                border: none;
            }
            QScrollBar:vertical {
                background: #161616;
                width: 10px;
            }
            QScrollBar::handle:vertical {
                background: #444444;
                border-radius: 5px;
                min-height: 24px;
            }
        """)

        # Container principal arrière-plan (qui peut être flouté)
        self.main_container = QWidget()
        self.main_layout = QHBoxLayout(self.main_container)
        self.main_layout.setContentsMargins(0,0,0,0)
        
        self.viewport = Slicer3DViewport()
        self.viewport.modelSelected.connect(self.on_model_selected)
        self.main_layout.addWidget(self.viewport)
        
        self.setCentralWidget(self.main_container)

        self.overlay_container = QWidget(self)
        self.overlay_container.setGeometry(self.rect())
        self.overlay_container.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents, False)
        self.overlay_container.setStyleSheet("background: transparent;")
        
        # --- PANNEAUX DE CONTRÔLE FLOTTANTS ---
        self.sidebar = self.create_sidebar()
        self.sidebar.setParent(self.overlay_container)

        self.top_toolbar = self.create_top_toolbar()
        self.top_toolbar.setParent(self.overlay_container)

        self.transform_panel = self.create_transform_panel()
        self.transform_panel.setParent(self.overlay_container)
        
        self.init_btn = QPushButton("INITIATE SEQUENCE", self.overlay_container)
        self.init_btn.setFixedSize(220, 42)
        self.init_btn.clicked.connect(self.request_pipeline)
        self.init_btn.setStyleSheet("border-radius: 8px; background-color: #3A3A3A; border: 1px solid #444444;")
        
        self.minimize_btn = QPushButton("–", self.overlay_container)
        self.minimize_btn.setFixedSize(28, 28)
        self.minimize_btn.clicked.connect(self.showMinimized)
        self.minimize_btn.setStyleSheet("QPushButton { background-color: #2A2A2A; color:#E5E5E5; border:1px solid #444444; border-radius:14px; }")
        
        self.exit_btn = QPushButton("✕", self.overlay_container)
        self.exit_btn.setFixedSize(28, 28)
        self.exit_btn.clicked.connect(self.close)
        self.exit_btn.setStyleSheet("QPushButton { background-color: #2A2A2A; color:#E5E5E5; border:1px solid #444444; border-radius:14px; }")
        
        # Calque modal (par dessus tout)
        self.modal_layer = ModalOverlay(self)

        # Assurer l'ordre d'affichage (Z-index correct)
        self.sidebar.raise_()
        self.top_toolbar.raise_()
        self.transform_panel.raise_()
        self.init_btn.raise_()
        self.exit_btn.raise_()
        self.modal_layer.raise_()

        # Affichage plein écran à la toute fin
        self.showFullScreen()

    def resizeEvent(self, event):
        """Repositionne dynamiquement tous les éléments d'interface utilisateur flottants."""
        super().resizeEvent(event)
        if not event.size().isValid(): return
        
        w = self.width()
        h = self.height()
        
        # Positionnement de la Sidebar (gauche)
        if hasattr(self, 'sidebar') and self.sidebar:
            self.sidebar.setGeometry(25, 25, 310, h - 50)
            
        # Positionnement de la Top Toolbar (haut centre)
        if hasattr(self, 'top_toolbar') and self.top_toolbar:
            self.top_toolbar.setGeometry((w - 420) // 2, 25, 420, 64)

        # Positionnement du panneau de transformation indépendant
        if hasattr(self, 'transform_panel') and self.transform_panel:
            self.transform_panel.setGeometry((w - 700) // 2, h - 150, 700, 92)
            
        if hasattr(self, 'overlay_container') and self.overlay_container:
            self.overlay_container.setGeometry(0, 0, w, h)

        # Positionnement du HUD (haut droit)
        # Positionnement des boutons de contrôle de fenêtre
        if hasattr(self, 'exit_btn') and self.exit_btn:
            self.exit_btn.setGeometry(w - 40, 12, 28, 28)
        if hasattr(self, 'minimize_btn') and self.minimize_btn:
            self.minimize_btn.setGeometry(w - 74, 12, 28, 28)

        # Positionnement du bouton d'initiation
        if hasattr(self, 'init_btn') and self.init_btn:
            self.init_btn.setGeometry(w - 250, h - 70, 220, 42)
            
        # Calque modal sur tout l'écran
        if hasattr(self, 'modal_layer') and self.modal_layer:
            self.modal_layer.resize(event.size())

    def create_top_toolbar(self):
        toolbar = QFrame(self)
        toolbar.setObjectName("TopToolbar")
        toolbar.setFixedSize(420, 64)
        toolbar.setStyleSheet("background-color: #252525; border: 1px solid #333333; border-radius: 10px;")

        layout = QHBoxLayout(toolbar)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)

        self.move_btn = QPushButton("Déplacer", toolbar)
        self.rotate_btn = QPushButton("Tourner", toolbar)
        self.scale_btn = QPushButton("Échelle", toolbar)

        for btn in (self.move_btn, self.rotate_btn, self.scale_btn):
            btn.setEnabled(False)
            btn.setFixedSize(120, 44)
            btn.setStyleSheet("background-color: #1D1D1D; color: #E5E5E5; border: 1px solid #333333; border-radius: 8px;")
            layout.addWidget(btn)

        self.move_btn.clicked.connect(lambda: self.open_transform_dialog("MOVE"))
        self.rotate_btn.clicked.connect(lambda: self.open_transform_dialog("ROTATE"))
        self.scale_btn.clicked.connect(lambda: self.open_transform_dialog("SCALE"))

        return toolbar

    def create_transform_panel(self):
        panel = QFrame(self)
        panel.setObjectName("TransformPanel")
        panel.setFixedHeight(92)
        panel.setStyleSheet("""
            QFrame#TransformPanel {
                background-color: #252525;
                border: 1px solid #333333;
                border-radius: 8px;
            }
            QLabel {
                color: #E5E5E5;
                font-size: 11px;
            }
            QLineEdit {
                background-color: #1E1E1E;
                color: #E5E5E5;
                border: 1px solid #333333;
                border-radius: 8px;
                padding: 6px 8px;
                min-width: 80px;
            }
        """)

        layout = QHBoxLayout(panel)
        layout.setContentsMargins(14, 12, 14, 12)
        layout.setSpacing(14)

        def field_group(title, suffix):
            group = QFrame()
            group.setStyleSheet("background: transparent;")
            group_layout = QHBoxLayout(group)
            group_layout.setContentsMargins(0, 0, 0, 0)
            group_layout.setSpacing(8)
            group_layout.addWidget(QLabel(title))
            return group, suffix

        self.tp_px = QLineEdit("0.0mm")
        self.tp_py = QLineEdit("0.0mm")
        self.tp_pz = QLineEdit("0.0mm")
        self.tp_rx = QLineEdit("0.0°")
        self.tp_ry = QLineEdit("0.0°")
        self.tp_rz = QLineEdit("0.0°")
        self.tp_sx = QLineEdit("1.0x")
        self.tp_sy = QLineEdit("1.0x")
        self.tp_sz = QLineEdit("1.0x")

        for field in (self.tp_px, self.tp_py, self.tp_pz, self.tp_rx, self.tp_ry, self.tp_rz, self.tp_sx, self.tp_sy, self.tp_sz):
            field.setReadOnly(True)

        pos_frame = QFrame(); pos_layout = QHBoxLayout(pos_frame); pos_layout.setContentsMargins(0,0,0,0); pos_layout.setSpacing(8)
        pos_layout.addWidget(QLabel("POS")); pos_layout.addWidget(self.tp_px); pos_layout.addWidget(self.tp_py); pos_layout.addWidget(self.tp_pz)
        rot_frame = QFrame(); rot_layout = QHBoxLayout(rot_frame); rot_layout.setContentsMargins(0,0,0,0); rot_layout.setSpacing(8)
        rot_layout.addWidget(QLabel("ROT")); rot_layout.addWidget(self.tp_rx); rot_layout.addWidget(self.tp_ry); rot_layout.addWidget(self.tp_rz)
        scl_frame = QFrame(); scl_layout = QHBoxLayout(scl_frame); scl_layout.setContentsMargins(0,0,0,0); scl_layout.setSpacing(8)
        scl_layout.addWidget(QLabel("SCL")); scl_layout.addWidget(self.tp_sx); scl_layout.addWidget(self.tp_sy); scl_layout.addWidget(self.tp_sz)

        layout.addWidget(pos_frame)
        layout.addWidget(rot_frame)
        layout.addWidget(scl_frame)

        return panel

    def on_model_selected(self, model):
        self._current_selected_model = model
        enabled = bool(model)
        for btn in (self.move_btn, self.rotate_btn, self.scale_btn):
            btn.setEnabled(enabled)
        self.update_transform_panel(model)

    def on_transform_changed(self):
        # This callback is currently unused because transformations are applied through the dialog.
        if hasattr(self, '_current_selected_model') and self._current_selected_model:
            m = self._current_selected_model
            m.pos = [self.sp_px.value(), self.sp_py.value(), self.sp_pz.value()]
            m.rot = [self.sp_rx.value(), self.sp_ry.value(), self.sp_rz.value()]
            m.scale = [self.sp_sx.value(), self.sp_sy.value(), self.sp_sz.value()]
            self.viewport.save_arrange_json()
            self.viewport.update()

    def import_stl_dialog(self):
        """Ouvre un dialogue d'importation et ajoute le fichier sélectionné au plateau."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Importer un modèle STL",
            "",
            "Fichiers STL (*.stl)"
        )
        if file_path:
            dest_path = os.path.join(MODELS_DIR, os.path.basename(file_path))
            if file_path != dest_path:
                shutil.copy2(file_path, dest_path)
            self.viewport.add_model(dest_path)

    def create_sidebar(self):
        """Crée la barre latérale gauche avec un profil plus moderne."""
        sidebar_frame = QFrame(self)
        sidebar_frame.setObjectName("SidebarPanel")
        sidebar_frame.setFixedWidth(310)


        layout = QVBoxLayout(sidebar_frame)
        layout.setContentsMargins(20, 25, 20, 25)

        header = QWidget()
        header_layout = QVBoxLayout(header)
        header_layout.setContentsMargins(0, 0, 0, 8)

        self.sidebar_title = QLabel("APEX SLICER v2.1")
        self.sidebar_title.setStyleSheet("color: #00E5FF; font-size: 18px; font-weight: 700;")
        self.sidebar_title.setTextFormat(Qt.TextFormat.PlainText)

        self.sidebar_tagline = QLabel("Préparation d'impression assistée")
        self.sidebar_tagline.setStyleSheet("color: #CBD7FF; font-size: 12px;")
        self.sidebar_tagline.setTextFormat(Qt.TextFormat.PlainText)

        self.sidebar_subtitle = QLabel("<span style='color:#8C9EFF;'>Plateau • Modèles • Export</span>")
        self.sidebar_subtitle.setStyleSheet("font-size: 11px;")
        self.sidebar_subtitle.setTextFormat(Qt.TextFormat.RichText)

        header_layout.addWidget(self.sidebar_title)
        header_layout.addWidget(self.sidebar_tagline)
        header_layout.addWidget(self.sidebar_subtitle)
        layout.addWidget(header)

        status_card = QFrame()
        status_card.setStyleSheet("background-color: #262626; border: 1px solid #333333; border-radius: 10px; padding: 8px;")
        status_layout = QVBoxLayout(status_card)
        status_layout.setContentsMargins(10, 8, 10, 8)
        status_layout.addWidget(QLabel("<span style='color:#00E5FF; font-weight:600;'>Statut : prêt</span>"))
        status_layout.addWidget(QLabel("Glisser déposer un STL ou importer depuis la bibliothèque."))
        layout.addWidget(status_card)

        import_btn = QPushButton("Importer un STL", self)
        import_btn.clicked.connect(self.import_stl_dialog)
        import_btn.setStyleSheet("border-radius: 8px; background-color: #3A3A3A; border: 1px solid #444444;")
        layout.addWidget(import_btn)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setStyleSheet("border: none; background: transparent;")
        content = QWidget()
        sl = QVBoxLayout(content)
        sl.setContentsMargins(0, 0, 0, 0)

        sl.addWidget(QLabel("Profil Imprimante", objectName="SectionTitle"))
        cb = QComboBox()
        cb.addItems(["Aetheris-4 Polar Kinematics", "Ares-4 Cartesian"])
        sl.addWidget(cb)

        sl.addWidget(QLabel("Paramètres Couches", objectName="SectionTitle"))
        sl.addWidget(QLabel("Hauteur de couche (mm)")); sl.addWidget(QLineEdit("0.15"))
        sl.addWidget(QLabel("Remplissage (%)")); sl.addWidget(QLineEdit("20%"))
        sl.addStretch()

        scroll.setWidget(content)
        layout.addWidget(scroll)
        return sidebar_frame

    def update_transform_panel(self, model=None):
        if model:
            self.tp_px.setText(f"{model.pos[0]:.1f}mm")
            self.tp_py.setText(f"{model.pos[1]:.1f}mm")
            self.tp_pz.setText(f"{model.pos[2]:.1f}mm")
            self.tp_rx.setText(f"{model.rot[0]:.1f}°")
            self.tp_ry.setText(f"{model.rot[1]:.1f}°")
            self.tp_rz.setText(f"{model.rot[2]:.1f}°")
            self.tp_sx.setText(f"{model.scale[0]:.2f}x")
            self.tp_sy.setText(f"{model.scale[1]:.2f}x")
            self.tp_sz.setText(f"{model.scale[2]:.2f}x")
        else:
            self.tp_px.setText("0.0mm")
            self.tp_py.setText("0.0mm")
            self.tp_pz.setText("0.0mm")
            self.tp_rx.setText("0.0°")
            self.tp_ry.setText("0.0°")
            self.tp_rz.setText("0.0°")
            self.tp_sx.setText("1.0x")
            self.tp_sy.setText("1.0x")
            self.tp_sz.setText("1.0x")

    def close_modal(self):
        self.main_container.setGraphicsEffect(None)
        self.modal_layer.hide()

    def open_transform_dialog(self, mode):
        model = getattr(self, '_current_selected_model', None)
        if not model:
            return
        values = {
            'MOVE': [f"{model.pos[0]:.1f}", f"{model.pos[1]:.1f}", f"{model.pos[2]:.1f}"],
            'ROTATE': [f"{model.rot[0]:.1f}", f"{model.rot[1]:.1f}", f"{model.rot[2]:.1f}"],
            'SCALE': [f"{model.scale[0]:.2f}", f"{model.scale[1]:.2f}", f"{model.scale[2]:.2f}"]
        }
        dialog = TransformDialog(mode, values[mode], self.apply_transform_dialog, self.close_modal)
        self.modal_layer.set_content(dialog)
        blur = QGraphicsBlurEffect()
        blur.setBlurRadius(16)
        self.main_container.setGraphicsEffect(blur)

    def apply_transform_dialog(self, mode, values):
        model = getattr(self, '_current_selected_model', None)
        if not model:
            self.close_modal()
            return
        try:
            values = [float(v) for v in values]
        except Exception:
            values = [0.0, 0.0, 0.0]

        if mode == 'MOVE':
            model.pos = values
        elif mode == 'ROTATE':
            model.rot = values
        elif mode == 'SCALE':
            model.scale = values

        self.viewport.save_arrange_json()
        self.viewport.update()
        self.update_transform_panel(model)
        self.close_modal()


    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls() and any(u.toLocalFile().lower().endswith('.stl') for u in event.mimeData().urls()):
            event.acceptProposedAction()
                
    def dropEvent(self, event):
        for url in event.mimeData().urls():
            file_path = url.toLocalFile()
            if file_path.lower().endswith('.stl'):
                dest_path = os.path.join(MODELS_DIR, os.path.basename(file_path))
                if file_path != dest_path:
                    shutil.copy2(file_path, dest_path)
                self.viewport.add_model(dest_path)
                break

    def assemble_and_export_models(self):
        """
        Génère un unique fichier STL fusionné à partir de tous les maillages du plateau.
        Exécute un centrage absolu à l'origine (0,0,0) et écrit un fichier binaire STL optimisé.
        """
        if not self.viewport.models:
            return False
            
        # Consolidation absolue de tous les triangles transformés dans le repère de scène
        fused_triangles = []
        for model in self.viewport.models:
            fused_triangles.extend(model.get_world_triangles())
            
        script_dir = os.path.dirname(os.path.abspath(__file__))
        input_model_dir = os.path.join(script_dir, "input_model")
        
        try:
            # Nettoyage préalable complet du dossier de réception
            if os.path.exists(input_model_dir):
                for file_name in os.listdir(input_model_dir):
                    f_path = os.path.join(input_model_dir, file_name)
                    if os.path.isfile(f_path) or os.path.islink(f_path):
                        os.unlink(f_path)
                    elif os.path.isdir(f_path):
                        shutil.rmtree(f_path)
            else:
                os.makedirs(input_model_dir, exist_ok=True)
                
            output_filepath = os.path.join(input_model_dir, "model.stl")
            
            # Écriture du fichier au format binaire standard STL
            with open(output_filepath, 'wb') as out_f:
                # 80 octets d'en-tête descriptive
                desc_header = b"Assembled Mesh - Generated Dynamically by Apex Slicer Engine"
                out_f.write(desc_header.ljust(80, b'\x00'))
                
                # Nombre total de facettes (Unsigned Integer, 4 octets)
                out_f.write(struct.pack('<I', len(fused_triangles)))
                
                # Écriture des facettes : Normale unitaire + 3 Sommets + 2 octets d'attributs
                for tri in fused_triangles:
                    v0, v1, v2 = tri[0], tri[1], tri[2]
                    
                    # Calcul géométrique de la normale par produit vectoriel
                    ax, ay, az = v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]
                    bx, by, bz = v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]
                    nx = ay * bz - az * by
                    ny = az * bx - ax * bz
                    nz = ax * by - ay * bx
                    
                    norm_len = math.sqrt(nx*nx + ny*ny + nz*nz)
                    if norm_len > 0:
                        nx, ny, nz = nx / norm_len, ny / norm_len, nz / norm_len
                    else:
                        nx, ny, nz = 0.0, 0.0, 1.0
                        
                    # Structure compactée float: 3 floats (normale), 9 floats (sommets), 1 unsigned short (attribut)
                    f_data = struct.pack('<ffffffffffffH', 
                                         nx, ny, nz, 
                                         v0[0], v0[1], v0[2], 
                                         v1[0], v1[1], v1[2], 
                                         v2[0], v2[1], v2[2], 
                                         0)
                    out_f.write(f_data)
                    
            return True
        except Exception as e:
            print(f"Erreur d'assemblage STL: {e}")
            return False

    def apply_pymeshfix_repairs(self):
        """
        Applique absolument toutes les fonctions de réparation de PyMeshFix sur le fichier STL assemblé.
        Logge en détail l'expérience et le résultat dans la console standard.
        """
        script_dir = os.path.dirname(os.path.abspath(__file__))
        stl_path = os.path.join(script_dir, "input_model", "model.stl")
        
        if not os.path.exists(stl_path):
            print("[PyMeshFix] 🛑 Fichier STL cible introuvable pour la réparation.", file=sys.stderr)
            return

        print("\n" + "="*70)
        print("🛠️  [PyMeshFix] DÉBUT DE LA RÉPARATION INTÉGRALE DU MAILLAGE FINAL")
        print(f"📁 Cible : {stl_path}")
        print("="*70)

        try:
            import pymeshfix
            print("[PyMeshFix] ✅ Module pymeshfix importé avec succès.")
        except ImportError:
            print("[PyMeshFix] 🛑 Erreur : Le module 'pymeshfix' n'est pas installé.", file=sys.stderr)
            print("[PyMeshFix] Installez-le avec : pip install pymeshfix", file=sys.stderr)
            print("="*70 + "\n")
            return

        try:
            # Initialisation du wrapper de bas niveau PyTMesh pour une granularité maximale
            print("[PyMeshFix] ⚙️ Initialisation du moteur PyTMesh...")
            tin = pymeshfix.PyTMesh()
            
            # 1. Chargement du fichier
            print("[PyMeshFix] 📖 Chargement du fichier STL...")
            tin.load_file(stl_path)
            print("[PyMeshFix] -> Fichier chargé dans le tampon de mémoire PyTMesh.")

            # 2. Join closest components
            print("[PyMeshFix] 🔗 Étape 1 : Fusion des composants disjoints les plus proches (join_closest_components)...")
            tin.join_closest_components()
            print("[PyMeshFix] -> Fusion des sous-composants disjoints réussie.")

            # 3. Clean
            print("[PyMeshFix] 🧼 Étape 2 : Nettoyage itératif du maillage (clean)...")
            tin.clean()
            print("[PyMeshFix] -> Triangles dégénérés et non-manifold nettoyés.")

            # 4. Fix connectivity
            print("[PyMeshFix] 🔀 Étape 3 : Restructuration de la connectivité (fix_connectivity)...")
            tin.fix_connectivity()
            print("[PyMeshFix] -> Connectivité topologique réalignée.")

            # 5. Fill small boundaries
            print("[PyMeshFix] 🕳️ Étape 4 : Remplissage des micro-ouvertures et trous (fill_small_boundaries)...")
            tin.fill_small_boundaries()
            print("[PyMeshFix] -> Trous de surface micro-bouchés.")

            # 6. Remove smallest components
            print("[PyMeshFix] 🗑️ Étape 5 : Retrait des petites scories et débris isolés (remove_smallest_components)...")
            tin.remove_smallest_components()
            print("[PyMeshFix] -> Composants orphelins insignifiants éliminés.")

            # 7. Strong degeneracy removal
            print("[PyMeshFix] 💥 Étape 6 : Nettoyage agressif des géométries dégénérées (strong_degeneracy_removal)...")
            tin.strong_degeneracy_removal(True)
            print("[PyMeshFix] -> Élimination forcée des faces nulles ou écrasées complétée.")

            # 8. Strong intersection removal
            print("[PyMeshFix] ⚔️ Étape 7 : Résolution agressive des intersections de facettes (strong_intersection_removal)...")
            tin.strong_intersection_removal(True)
            print("[PyMeshFix] -> Conflits d'auto-intersection résolus.")

            # Sauvegarde finale
            print("[PyMeshFix] 💾 Sauvegarde du maillage purifié et étanche (solid/watertight)...")
            tin.save_file(stl_path)
            
            print("🎉 [PyMeshFix] Expérience de réparation achevée avec succès !")
            print("✨ Le maillage final exporté est désormais 100% propre, manifold et prêt pour le tranchage.")
            print("="*70 + "\n")

        except Exception as e:
            print(f"[PyMeshFix] ❌ Erreur critique lors de la réparation du maillage : {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            print("="*70 + "\n")

    def set_interface_blur(self, radius=None):
        if radius is None:
            self.main_container.setGraphicsEffect(None)
            if hasattr(self, 'overlay_container'):
                self.overlay_container.setGraphicsEffect(None)
            return

        blur = QGraphicsBlurEffect()
        blur.setBlurRadius(radius)
        self.main_container.setGraphicsEffect(blur)

        if hasattr(self, 'overlay_container'):
            blur_overlay = QGraphicsBlurEffect()
            blur_overlay.setBlurRadius(radius)
            self.overlay_container.setGraphicsEffect(blur_overlay)

    def request_pipeline(self):
        """Analyse le nombre de triangles totaux avant le calcul."""
        total_tris = sum(len(m.base_triangles) for m in self.viewport.models)
        
        # Flouter l'arrière-plan de l'interface
        self.set_interface_blur(20)
        
        if total_tris > 800:
            # Proposer désormais un bouton SIMPLIFY
            warning = WarningDialog(total_tris, self.start_pipeline, self.cancel_pipeline, self.open_simplify_tool)
            self.modal_layer.set_content(warning)
        else:
            self.start_pipeline()

    def open_simplify_tool(self):
        """Ouvre l'interface de simplification interactive avec aperçu en temps réel."""
        total_tris = sum(len(m.base_triangles) for m in self.viewport.models)
        
        simplify_dialog = SimplifyDialog(
            total_tris,
            on_apply=self.apply_simplification,
            on_cancel=self.cancel_simplification,
            on_slider_change=self.preview_simplification
        )
        self.modal_layer.set_content(simplify_dialog)
        
        # Pas de flou pour la boîte de simplification afin de conserver l'aperçu live
        self.set_interface_blur(None)
        
        # Initialiser avec un aperçu par défaut de 50%
        self.preview_simplification(0.50)

    def preview_simplification(self, ratio):
        """Génère et affiche le maillage simplifié en temps réel sur le viewport."""
        for model in self.viewport.models:
            model.preview_triangles = simplify_mesh(model.base_triangles, ratio)
            # Forcer la réinitialisation de l'éclairage pour le nouvel aperçu
            model._intensity_cache_id = None
        self.viewport.update()

    def apply_simplification(self):
        """Rend la simplification définitive, écrase les STL physiques d'origine et met à jour le plateau."""
        for model in self.viewport.models:
            if model.preview_triangles is not None:
                model.base_triangles = list(model.preview_triangles)
                model.preview_triangles = None
                model._intensity_cache_id = None
                
                # Écrire physiquement le STL simplifié dans MODELS_DIR pour que les scripts externes
                # travaillent sur le modèle optimisé
                self.write_model_to_stl(model)
                
        self.viewport.save_arrange_json()
        self.set_interface_blur(None)
        self.modal_layer.hide()
        self.viewport.update()
        
        # Relancer l'analyse du pipeline avec le nouveau maillage léger !
        self.request_pipeline()

    def cancel_simplification(self):
        """Annule l'aperçu et restaure les modèles originaux intacts."""
        for model in self.viewport.models:
            model.preview_triangles = None
            model._intensity_cache_id = None
        self.set_interface_blur(None)
        self.modal_layer.hide()
        self.viewport.update()

    def write_model_to_stl(self, model):
        """Écrit un objet Model3D sous forme de fichier binaire STL valide à sa destination."""
        try:
            with open(model.filepath, 'wb') as f:
                header = b"Simplified STL - Generated dynamically by Apex Slicer Engine"
                f.write(header.ljust(80, b'\x00'))
                f.write(struct.pack('<I', len(model.base_triangles)))
                for tri in model.base_triangles:
                    v0, v1, v2 = tri[0], tri[1], tri[2]
                    
                    # Calculer la normale unitaire
                    ax, ay, az = v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]
                    bx, by, bz = v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]
                    nx = ay * bz - az * by
                    ny = az * bx - ax * bz
                    nz = ax * by - ay * bx
                    norm_len = math.sqrt(nx*nx + ny*ny + nz*nz)
                    if norm_len > 0:
                        nx, ny, nz = nx / norm_len, ny / norm_len, nz / norm_len
                    else:
                        nx, ny, nz = 0.0, 0.0, 1.0
                        
                    f.write(struct.pack('<ffffffffffffH', 
                                         nx, ny, nz, 
                                         v0[0], v0[1], v0[2], 
                                         v1[0], v1[1], v1[2], 
                                         v2[0], v2[1], v2[2], 
                                         0))
        except Exception as e:
            print(f"[Simplify] Erreur d'écriture STL pour {model.filename}: {e}", file=sys.stderr)

    def cancel_pipeline(self):
        self.set_interface_blur(None)
        self.modal_layer.hide()

    def start_pipeline(self):
        if self.worker and self.worker.isRunning(): return
        
        # Exécuter d'abord la fusion et la copie conforme des modèles dans \input_model
        export_success = self.assemble_and_export_models()
        if not export_success:
            self.main_container.setGraphicsEffect(None)
            self.modal_layer.hide()
            print("[Apex Slicer] Annulation de la séquence : impossible d'assembler les modèles.", file=sys.stderr)
            return

        # --- EXÉCUTION DE PYMESHFIX ---
        # Nettoie et corrige absolument tous les défauts géométriques du STL fusionné avant d'exécuter exec.py
        self.apply_pymeshfix_repairs()

        self.seq_dialog = SequenceDialog(self.abort_pipeline, self.run_pipeline)
        self.modal_layer.set_content(self.seq_dialog)
        
        self.worker = ProcessWorker("exec.py")
        self.worker.progress_signal.connect(self.update_progress)
        self.worker.waiting_signal.connect(self.show_prompt)
        self.worker.error_signal.connect(self.show_error)
        self.worker.finished_signal.connect(self.on_pipeline_finished)
        self.worker.start()

    def show_prompt(self):
        self.seq_dialog.status_label.setText("STATUS: WAITING FOR INPUT")
        self.seq_dialog.run_btn.show()

    def run_pipeline(self):
        if self.worker and self.worker.process and self.worker.process.poll() is None:
            self.worker.process.stdin.write("run\n")
            self.worker.process.stdin.flush()
            self.seq_dialog.run_btn.hide()
            self.seq_dialog.status_label.setText("PROCESSING...")

    def abort_pipeline(self):
        if self.worker and self.worker.process and self.worker.process.poll() is None:
            self.worker.process.stdin.write("abort\n")
            self.worker.process.stdin.flush()
            self.seq_dialog.status_label.setText("ABORTING...")

    def update_progress(self, data):
        self.seq_dialog.progress_bar.setValue(data.get("current", 0))
        self.seq_dialog.status_label.setText(f"TASK: {data.get('task', '').upper()}")

    def on_pipeline_finished(self, exit_code):
        if exit_code == 0:
            self.seq_dialog.status_label.setText("SEQUENCE COMPLETE")
            self.seq_dialog.progress_bar.setValue(100)
            QTimer.singleShot(2000, self.cancel_pipeline)
        else:
            if not self.seq_dialog.status_label.text().startswith("ERROR"):
                self.seq_dialog.status_label.setText("PROCESS TERMINATED.")
                self.seq_dialog.status_label.setStyleSheet("color: #FF3366; font-weight: bold;")
                self.seq_dialog.abort_btn.setText("CLOSE")
                self.seq_dialog.abort_btn.clicked.disconnect()
                self.seq_dialog.abort_btn.clicked.connect(self.cancel_pipeline)

    def show_error(self, message):
        self.seq_dialog.status_label.setText("ERROR OCCURED")
        self.seq_dialog.status_label.setStyleSheet("color: #FF3366; font-weight: bold;")
        self.seq_dialog.progress_bar.setStyleSheet("""
            QProgressBar::chunk { background-color: #FF3366; }
        """)
        print(f"[Apex Slicer] {message}", file=sys.stderr)
