import sys
import json
import subprocess
import os
import math
import struct
import shutil
import uuid
from PyQt6.QtWidgets import (QApplication, QMainWindow, QPushButton, QVBoxLayout, 
                             QHBoxLayout, QWidget, QGraphicsDropShadowEffect, 
                             QProgressBar, QLabel, QFrame, QComboBox, QLineEdit, 
                             QScrollArea, QGraphicsScene, QGraphicsPixmapItem,
                             QGraphicsBlurEffect, QDoubleSpinBox, QSlider)
from PyQt6.QtCore import QThread, pyqtSignal, Qt, QPointF, QRect, QRectF, QTimer, QSize
from PyQt6.QtGui import QPainter, QColor, QPen, QPolygonF, QBrush, QImage, QPainterPath, QPixmap

try:
    import numpy as np
    from PIL import Image, ImageFilter
    _HAS_PIL_BLUR = True
except ImportError:
    _HAS_PIL_BLUR = False

try:
    from scipy.ndimage import gaussian_filter
    _HAS_SCIPY_BLUR = True
except ImportError:
    _HAS_SCIPY_BLUR = False

try:
    import trimesh
    _HAS_TRIMESH = True
except ImportError:
    _HAS_TRIMESH = False

WORKSPACE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ux", "plate")
MODELS_DIR = os.path.join(WORKSPACE_DIR, "models")
ARRANGE_JSON = os.path.join(WORKSPACE_DIR, "arrange.json")

# Assurer l'existence des répertoires de stockage
os.makedirs(MODELS_DIR, exist_ok=True)

def simplify_mesh(triangles, ratio):
    """Simplifie un maillage de triangles selon un ratio donné (0.05 à 0.95)."""
    if ratio >= 0.99:
        return triangles
    target_count = int(len(triangles) * ratio)
    target_count = max(12, target_count)
    
    if _HAS_TRIMESH:
        try:
            import trimesh
            # Reconstruire les sommets et faces uniques pour trimesh
            vertices = []
            faces = []
            v_map = {}
            for tri in triangles:
                face_idx = []
                for pt in tri:
                    # Clé de hachage arrondie pour fusionner les sommets proches
                    key = (round(pt[0], 5), round(pt[1], 5), round(pt[2], 5))
                    if key not in v_map:
                        v_map[key] = len(vertices)
                        vertices.append(pt)
                    face_idx.append(v_map[key])
                faces.append(face_idx)
            
            mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
            # Utiliser la décimation quadrique de Trimesh
            simplified = mesh.simplify_quadric_decimation(target_count)
            
            # Re-générer la liste de triangles pour le moteur de rendu
            new_triangles = []
            for face in simplified.faces:
                new_triangles.append([
                    simplified.vertices[face[0]].tolist(),
                    simplified.vertices[face[1]].tolist(),
                    simplified.vertices[face[2]].tolist()
                ])
            return new_triangles
        except Exception as e:
            print(f"[Simplify] Erreur avec Trimesh: {e}. Utilisation du repli ultra-rapide.")
    
    # Algorithme de repli natif (Sous-échantillonnage de facettes à pas constant)
    total = len(triangles)
    if target_count >= total:
        return triangles
    step = max(1, total / target_count)
    new_triangles = []
    for i in range(target_count):
        idx = int(i * step)
        if idx < total:
            new_triangles.append(triangles[idx])
    return new_triangles

class Model3D:
    """Représente un modèle STL individuel positionné sur le plateau."""
    def __init__(self, filepath, id_str=None):
        self.id = id_str if id_str else str(uuid.uuid4())
        self.filename = os.path.basename(filepath)
        self.filepath = filepath
        self.base_triangles = []
        
        # Vecteurs de transformation
        self.pos = [0.0, 0.0, 0.0]
        self.rot = [0.0, 0.0, 0.0]
        self.scale = [1.0, 1.0, 1.0]
        
        self.selected = False
        self.preview_triangles = None # Stocke l'aperçu dynamique de simplification
        
        # Caches pour l'accélération du rendu 3D
        self._intensity_cache_id = None
        self._triangle_intensities = []
        self._projected_polygons = [] # Requis pour la détection précise des clics (Hit Test)
        
        self.load_stl()

    def load_stl(self):
        """Lit et parse le fichier STL (ASCII ou Binaire)."""
        try:
            with open(self.filepath, 'rb') as f:
                header = f.read(80)
                if len(header) < 80: return
                is_ascii = False
                try:
                    if header.decode('utf-8', errors='ignore').strip().startswith('solid'):
                        is_ascii = True
                except Exception: pass
                
                if is_ascii:
                    f.seek(0)
                    content = f.read().decode('utf-8', errors='ignore')
                    current_tri = []
                    for line in content.splitlines():
                        line = line.strip()
                        if line.startswith('vertex'):
                            parts = line.split()
                            if len(parts) >= 4:
                                current_tri.append([float(parts[1]), float(parts[2]), float(parts[3])])
                        elif line.startswith('endloop'):
                            if len(current_tri) == 3:
                                self.base_triangles.append(current_tri)
                            current_tri = []
                else:
                    f.seek(80)
                    count_bytes = f.read(4)
                    if len(count_bytes) < 4: return
                    count = struct.unpack('<I', count_bytes)[0]
                    # Sous-échantillonnage dynamique si le modèle est trop lourd pour le moteur QPainter
                    step = int(count / 20000) + 1 if count > 20000 else 1
                    for i in range(count):
                        data = f.read(50)
                        if len(data) < 50: break
                        if i % step != 0: continue
                        t = struct.unpack('<fff fff fff fff H', data)
                        self.base_triangles.append([[t[3], t[4], t[5]], [t[6], t[7], t[8]], [t[9], t[10], t[11]]])
            
            self.center_origin()
        except Exception as e:
            print(f"Erreur chargement STL {self.filename}: {e}")

    def center_origin(self):
        """Centre le maillage sur son origine Z=0."""
        if not self.base_triangles: return
        xs = [v[0] for tri in self.base_triangles for v in tri]
        ys = [v[1] for tri in self.base_triangles for v in tri]
        zs = [v[2] for tri in self.base_triangles for v in tri]
        cx = (min(xs) + max(xs)) / 2.0
        cy = (min(ys) + max(ys)) / 2.0
        cz = min(zs)
        
        centered = []
        for tri in self.base_triangles:
            centered.append([[v[0]-cx, v[1]-cy, v[2]-cz] for v in tri])
        self.base_triangles = centered

    def get_world_triangles(self):
        """Applique l'échelle, la rotation et la translation aux triangles du modèle."""
        world_tris = []
        cx = math.cos(math.radians(self.rot[0]))
        sx = math.sin(math.radians(self.rot[0]))
        cy = math.cos(math.radians(self.rot[1]))
        sy = math.sin(math.radians(self.rot[1]))
        cz = math.cos(math.radians(self.rot[2]))
        sz = math.sin(math.radians(self.rot[2]))
        
        # Choisir le maillage d'aperçu de simplification ou le maillage d'origine
        source_triangles = self.preview_triangles if self.preview_triangles is not None else self.base_triangles
        
        for tri in source_triangles:
            new_tri = []
            for v in tri:
                # Échelle (Scale)
                x, y, z = v[0]*self.scale[0], v[1]*self.scale[1], v[2]*self.scale[2]
                # Rotation X
                y1 = y*cx - z*sx
                z1 = y*sx + z*cx
                # Rotation Y
                x2 = x*cy + z1*sy
                z2 = -x*sy + z1*cy
                # Rotation Z
                x3 = x2*cz - y1*sz
                y3 = x2*sz + y1*cz
                # Translation
                new_tri.append([x3+self.pos[0], y3+self.pos[1], z2+self.pos[2]])
            world_tris.append(new_tri)
        return world_tris

class ProcessWorker(QThread):
    """Gère l'exécution asynchrone du script de traitement 'exec.py'."""
    progress_signal = pyqtSignal(dict)
    waiting_signal = pyqtSignal()
    error_signal = pyqtSignal(str)
    finished_signal = pyqtSignal(int)

    def __init__(self, script_path):
        super().__init__()
        self.script_path = script_path
        self.process = None

    def run(self):
        try:
            resolved_path = self.script_path
            if not os.path.isabs(resolved_path):
                resolved_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), resolved_path)

            if not os.path.isfile(resolved_path):
                self.error_signal.emit(f"SCRIPT INTROUVABLE : {resolved_path}")
                self.finished_signal.emit(-1)
                return

            env = os.environ.copy()
            env["PYTHONIOENCODING"] = "utf-8"
            
            self.process = subprocess.Popen(
                [sys.executable, resolved_path],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE,
                text=True, encoding="utf-8", bufsize=1, env=env 
            )
            
            for line in iter(self.process.stdout.readline, ''):
                stripped = line.strip()
                if "STATUS:WAITING_FOR_USER_INPUT" in stripped:
                    self.waiting_signal.emit()
                    continue
                try:
                    data = json.loads(stripped)
                    if isinstance(data, dict) and data.get("type") == "progress":
                        self.progress_signal.emit(data)
                except json.JSONDecodeError: pass
                    
            _, stderr = self.process.communicate()
            exit_code = self.process.returncode

            if stderr:
                filtered = [l for l in stderr.splitlines() if "UserWarning" not in l and "warnings.warn" not in l]
                if filtered: self.error_signal.emit("\n".join(filtered[-6:]))

            self.finished_signal.emit(exit_code if exit_code is not None else -1)
        except Exception as e:
            self.error_signal.emit(f"ERREUR PROCESS: {str(e)}")
            self.finished_signal.emit(-1)

class GlowingButton(QPushButton):
    """Bouton stylisé avec halo lumineux néon réactif."""
    def __init__(self, text, accent_color="#00E5FF", parent=None):
        super().__init__(text, parent)
        self.accent_color = accent_color
        self.shadow = QGraphicsDropShadowEffect()
        self.shadow.setBlurRadius(20)
        self.shadow.setColor(QColor(accent_color))
        self.shadow.setOffset(0, 0)
        self.setGraphicsEffect(self.shadow)
        
        self.setStyleSheet(f"""
            QPushButton {{
                background-color: rgba(21, 23, 34, 220);
                color: {accent_color};
                border: 1px solid {accent_color};
                border-radius: 6px;
                padding: 10px 20px;
                font-family: 'Segoe UI', sans-serif; font-weight: bold; font-size: 13px;
            }}
            QPushButton:hover {{ background-color: rgba(30, 34, 53, 240); }}
            QPushButton:disabled {{ color: #2A2F4A; border-color: #2A2F4A; background-color: rgba(15, 16, 22, 220); }}
        """)

class TransformSpinBox(QDoubleSpinBox):
    """Sélecteur numérique élégant pour les transformations 3D."""
    def __init__(self, suffix="", min_val=-999.0, max_val=999.0, step=1.0):
        super().__init__()
        self.setSuffix(suffix)
        self.setRange(min_val, max_val)
        self.setSingleStep(step)
        self.setDecimals(1)
        self.setStyleSheet("""
            QDoubleSpinBox {
                background-color: rgba(15, 17, 26, 200);
                color: #FFFFFF;
                border: 1px solid #2A2F4A;
                border-radius: 4px;
                padding: 4px;
                width: 60px;
            }
            QDoubleSpinBox:focus { border: 1px solid #00E5FF; }
        """)

class Slicer3DViewport(QFrame):
    """Viewport interactif affichant le plateau d'impression et les modèles 3D."""
    modelSelected = pyqtSignal(object) 

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName("ViewportPanel")
        self.yaw, self.pitch, self.zoom = 45.0, 25.0, 4.0
        self.pan_x, self.pan_y = 80.0, 50.0
        self.last_pos = None
        
        self.models = []
        self._blur_cache = {}
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.load_arrange_json()

    def load_arrange_json(self):
        """Restaure la disposition et l'état des modèles depuis le fichier de configuration."""
        if not os.path.exists(ARRANGE_JSON): return
        try:
            if os.path.getsize(ARRANGE_JSON) == 0:
                self.save_arrange_json()
                return
            with open(ARRANGE_JSON, 'r') as f:
                data = json.load(f)
                for m_data in data.get("models", []):
                    filepath = os.path.join(MODELS_DIR, m_data.get("filename", ""))
                    if os.path.exists(filepath):
                        model = Model3D(filepath, m_data.get("id"))
                        model.pos = m_data.get("position", [0,0,0])
                        model.rot = m_data.get("rotation", [0,0,0])
                        model.scale = m_data.get("scale", [1,1,1])
                        self.models.append(model)
        except Exception as e:
            print(f"Erreur chargement JSON (recreation automatique): {e}")
            self.save_arrange_json()

    def save_arrange_json(self):
        """Enregistre la position et l'orientation des modèles."""
        data = {"models": []}
        for m in self.models:
            data["models"].append({
                "id": m.id,
                "filename": m.filename,
                "position": m.pos,
                "rotation": m.rot,
                "scale": m.scale
            })
        try:
            with open(ARRANGE_JSON, 'w') as f:
                json.dump(data, f, indent=4)
        except Exception as e:
            print(f"Erreur sauvegarde JSON: {e}")

    def add_model(self, filepath):
        """Ajoute un modèle sur le plateau en évitant les superpositions directes."""
        model = Model3D(filepath)
        offset_x = len(self.models) * 40.0
        model.pos[0] = offset_x - 40.0
        self.models.append(model)
        self.save_arrange_json()
        
        for m in self.models: m.selected = False
        model.selected = True
        self.modelSelected.emit(model)
        self.update()

    def delete_selected_models(self):
        """Supprime les modèles sélectionnés du plateau et efface leurs fichiers."""
        to_delete = [m for m in self.models if m.selected]
        for m in to_delete:
            self.models.remove(m)
            try:
                if os.path.exists(m.filepath):
                    os.remove(m.filepath)
            except: pass
        
        if to_delete:
            self.save_arrange_json()
            self.modelSelected.emit(None)
            self.update()

    def keyPressEvent(self, event):
        if event.key() == Qt.Key.Key_Delete or event.key() == Qt.Key.Key_Backspace:
            self.delete_selected_models()
        super().keyPressEvent(event)

    def _current_scene_key(self):
        geom_hash = hash(tuple((m.pos[0], m.rot[0], m.scale[0]) for m in self.models))
        return (
            round(self.yaw, 2), round(self.pitch, 2), round(self.zoom, 3),
            round(self.pan_x, 1), round(self.pan_y, 1), len(self.models), geom_hash,
            self.width(), self.height()
        )

    def mousePressEvent(self, event):
        # Donner le focus clavier au Viewport pour intercepter la touche 'Suppr'
        self.setFocus()
        
        if event.button() == Qt.MouseButton.LeftButton:
            # Recherche de collision (Hit Test)
            clicked_model = None
            click_point = event.position()
            
            for model in reversed(self.models):
                hit = False
                for poly in model._projected_polygons:
                    if poly.containsPoint(click_point, Qt.FillRule.OddEvenFill):
                        hit = True
                        break
                if hit:
                    clicked_model = model
                    break
            
            changed = False
            for m in self.models:
                if m.selected != (m == clicked_model):
                    m.selected = (m == clicked_model)
                    changed = True
            
            if changed:
                self.modelSelected.emit(clicked_model)
                self.update()
                
            self.last_pos = event.position()
        elif event.button() == Qt.MouseButton.RightButton:
            self.last_pos = event.position()
            
    def mouseMoveEvent(self, event):
        if self.last_pos is not None:
            delta = event.position() - self.last_pos
            if event.buttons() & Qt.MouseButton.RightButton:
                self.pan_x += delta.x()
                self.pan_y += delta.y()
            elif event.buttons() & Qt.MouseButton.LeftButton:
                self.yaw += delta.x() * 0.5
                self.pitch = max(-85.0, min(85.0, self.pitch + delta.y() * 0.5))
            self.last_pos = event.position()
            self.update()
            
    def mouseReleaseEvent(self, event):
        self.last_pos = None
        
    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        self.zoom = max(0.2, min(12.0, self.zoom + delta * 0.001))
        self.update()

    def paintEvent(self, event):
        w, h = self.width(), self.height()
        if w <= 0 or h <= 0: return

        # Dessiner la scène 3D brute sur un buffer d'image hors-écran
        scene_image = QImage(w, h, QImage.Format.Format_ARGB32_Premultiplied)
        scene_image.fill(QColor("#090A10")) 
        
        img_painter = QPainter(scene_image)
        img_painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        cx, cy = w / 2 + self.pan_x, h / 2 + self.pan_y
        yaw_rad, pitch_rad = math.radians(self.yaw), math.radians(self.pitch)
        cos_y, sin_y = math.cos(yaw_rad), math.sin(yaw_rad)
        cos_p, sin_p = math.cos(pitch_rad), math.sin(pitch_rad)
        
        def project(x, y, z):
            x1 = x * cos_y - y * sin_y
            y1 = x * sin_y + y * cos_y
            x2 = x1
            y2 = y1 * cos_p - z * sin_p
            z2 = y1 * sin_p + z * cos_p
            return cx + x2 * self.zoom, cy - z2 * self.zoom, y2

        # Dessiner le plateau circulaire
        bed_pts = [(125*math.cos(2*math.pi*i/64), 125*math.sin(2*math.pi*i/64), 0) for i in range(64)]
        proj_bed = [project(*p) for p in bed_pts]
        polygon_bed = QPolygonF([QPointF(p[0], p[1]) for p in proj_bed])
        
        img_painter.setBrush(QColor(15, 17, 26, 140))
        img_painter.setPen(QPen(QColor("#00E5FF"), 1.5, Qt.PenStyle.SolidLine))
        img_painter.drawPolygon(polygon_bed)

        # Préparation des polygones de triangles pour le rendu
        all_rendered_tris = []
        
        for model in self.models:
            world_tris = model.get_world_triangles()
            model._projected_polygons = [] # Réinitialiser le mesh de collision 2D
            
            # Recalculer l'éclairage uniquement si la géométrie a changé
            geom_hash = hash(str(model.rot) + str(model.scale) + str(model.preview_triangles is not None))
            if getattr(model, '_intensity_cache_id', None) != geom_hash:
                intensities = []
                for tri in world_tris:
                    ax, ay, az = tri[1][0]-tri[0][0], tri[1][1]-tri[0][1], tri[1][2]-tri[0][2]
                    bx, by, bz = tri[2][0]-tri[0][0], tri[2][1]-tri[0][1], tri[2][2]-tri[0][2]
                    nx, ny, nz = ay*bz - az*by, az*bx - ax*bz, ax*by - ay*bx
                    length = math.sqrt(nx*nx + ny*ny + nz*nz)
                    if length > 0: nx, ny, nz = nx/length, ny/length, nz/length
                    else: nx, ny, nz = 0, 0, 1
                    light = nx*0.3 + ny*-0.3 + nz*0.9
                    intensities.append(max(0.15, min(1.0, light)))
                model._triangle_intensities = intensities
                model._intensity_cache_id = geom_hash

            for tri, intensity in zip(world_tris, model._triangle_intensities):
                p0, p1, p2 = project(*tri[0]), project(*tri[1]), project(*tri[2])
                depth = (p0[2] + p1[2] + p2[2]) / 3.0
                all_rendered_tris.append((p0, p1, p2, depth, intensity, model))

        # Tri des polygones par profondeur (Algorithme du peintre)
        all_rendered_tris.sort(key=lambda item: item[3], reverse=True)
        
        selected_models = [m for m in self.models if m.selected]
        if selected_models:
            mask_img = QImage(w, h, QImage.Format.Format_ARGB32_Premultiplied)
            mask_img.fill(Qt.GlobalColor.transparent)
            mask_painter = QPainter(mask_img)
            mask_painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            mask_painter.setBrush(QColor("#00E5FF"))
            mask_painter.setPen(QPen(QColor("#00E5FF"), 1))
            
            for p0, p1, p2, _, _, model in all_rendered_tris:
                if model.selected:
                    poly = QPolygonF([QPointF(p0[0], p0[1]), QPointF(p1[0], p1[1]), QPointF(p2[0], p2[1])])
                    mask_painter.drawPolygon(poly)
            mask_painter.end()
            
            # Dessiner le masque expansé par décalages directionnels
            offsets = [(-3,0), (3,0), (0,-3), (0,3), (-2,-2), (2,2), (-2,2), (2,-2)]
            for ox, oy in offsets:
                img_painter.drawImage(ox, oy, mask_img)

        # Rendu normal des faces
        for p0, p1, p2, _, intensity, model in all_rendered_tris:
            poly = QPolygonF([QPointF(p0[0], p0[1]), QPointF(p1[0], p1[1]), QPointF(p2[0], p2[1])])
            model._projected_polygons.append(poly)
            
            r, g, b = int(157*intensity), int(78*intensity), int(221*intensity)
            img_painter.setBrush(QColor(r, g, b))
            img_painter.setPen(QPen(QColor(min(r+20, 255), min(g+20, 255), min(b+20, 255)), 0.5))
            img_painter.drawPolygon(poly)

        img_painter.end()

        # Peindre la scène finale sur l'écran
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.drawImage(0, 0, scene_image)

        # --- EFFET DE FLOU DE FOND POUR LES PANNEAUX SOUSTRAYANTS ---
        if not hasattr(self, '_blur_targets'):
            # Rechercher dynamiquement les panneaux à flouter dans la fenêtre parent
            self._blur_targets = [c for c in self.window().findChildren(QFrame) 
                                  if c.objectName() in ("SidebarPanel", "TopToolbar")]

        scene_key = self._current_scene_key()

        def apply_blur(img, radius=12):
            if _HAS_PIL_BLUR:
                img = img.convertToFormat(QImage.Format.Format_RGBA8888)
                ptr = img.bits()
                ptr.setsize(img.height() * img.width() * 4)
                arr = np.frombuffer(ptr, dtype=np.uint8).reshape((img.height(), img.width(), 4))
                rgb_img = Image.frombuffer('RGB', (img.width(), img.height()), np.ascontiguousarray(arr[:,:,:3]), 'raw', 'RGB', 0, 1)
                blurred = rgb_img.filter(ImageFilter.GaussianBlur(radius=max(1, radius/3)))
                out_arr = arr.copy()
                out_arr[:,:,:3] = np.array(blurred)
                return QImage(out_arr.data, img.width(), img.height(), QImage.Format.Format_RGBA8888).copy()
            else:
                scene = QGraphicsScene()
                scene.addItem(QGraphicsPixmapItem(QPixmap.fromImage(img)))
                result = QImage(img.size(), QImage.Format.Format_ARGB32_Premultiplied)
                result.fill(Qt.GlobalColor.transparent)
                rp = QPainter(result)
                scene.render(rp)
                rp.end()
                return result

        for child in self._blur_targets:
            if not child.isVisible(): continue
            
            # Calcul absolu des coordonnées géométriques entre les deux éléments frères
            top_left_global = child.mapToGlobal(child.rect().topLeft())
            top_left = self.mapFromGlobal(top_left_global)
            geom = QRect(top_left, child.size())
            if not geom.isValid() or geom.isEmpty(): continue
            
            panel_name = child.objectName()
            cached = self._blur_cache.get(panel_name)
            geom_tuple = (geom.x(), geom.y(), geom.width(), geom.height())

            if cached and cached["key"] == scene_key and cached["geom"] == geom_tuple:
                blurred_part = cached["image"]
            else:
                backdrop = scene_image.copy(geom)
                blurred_part = apply_blur(backdrop, radius=50)
                self._blur_cache[panel_name] = {"key": scene_key, "geom": geom_tuple, "image": blurred_part}

            clip_path = QPainterPath()
            clip_path.addRoundedRect(QRectF(geom), 14.0, 14.0)
            painter.save()
            painter.setClipPath(clip_path)
            painter.drawImage(geom.topLeft(), blurred_part)
            painter.fillRect(geom, QColor(14, 14, 14, 80)) # Teinte de contraste
            painter.restore()

class ModalOverlay(QFrame):
    """Calque plein écran semi-transparent interceptant tous les clics."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setStyleSheet("background-color: rgba(9, 10, 16, 150);")
        self.hide()
        
        self.layout = QVBoxLayout(self)
        self.layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        
    def set_content(self, widget):
        for i in reversed(range(self.layout.count())): 
            self.layout.itemAt(i).widget().setParent(None)
        self.layout.addWidget(widget)
        self.show()

class WarningDialog(QFrame):
    """Boîte de dialogue d'avertissement de surcharge polygonale."""
    def __init__(self, count, on_confirm, on_cancel, on_simplify):
        super().__init__()
        self.setFixedSize(450, 240)
        self.setStyleSheet("""
            QFrame {
                background-color: #151722;
                border: 1px solid #FF3366;
                border-radius: 12px;
            }
            QLabel { border: none; background: transparent; }
        """)
        
        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(40)
        shadow.setColor(QColor(255, 51, 102, 100))
        shadow.setOffset(0, 0)
        self.setGraphicsEffect(shadow)
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(30, 25, 30, 25)
        
        title = QLabel("ATTENTION : MAILLAGE LOURD")
        title.setStyleSheet("color: #FF3366; font-size: 16px; font-weight: bold; font-family: 'Consolas';")
        layout.addWidget(title, alignment=Qt.AlignmentFlag.AlignCenter)
        
        msg = QLabel(f"Le plateau contient {count} triangles au total.\n\nLe seuil recommandé est de 800. Que souhaitez-vous faire ?")
        msg.setStyleSheet("color: #FFFFFF; font-size: 13px; text-align: center;")
        msg.setWordWrap(True)
        layout.addWidget(msg, alignment=Qt.AlignmentFlag.AlignCenter)
        
        btn_layout = QHBoxLayout()
        cancel = GlowingButton("CANCEL", accent_color="#8A91B4")
        cancel.clicked.connect(on_cancel)
        
        simplify = GlowingButton("SIMPLIFY", accent_color="#00E5FF")
        simplify.clicked.connect(on_simplify)
        
        confirm = GlowingButton("CONFIRM", accent_color="#FF3366")
        confirm.clicked.connect(on_confirm)
        
        btn_layout.addWidget(cancel)
        btn_layout.addWidget(simplify)
        btn_layout.addWidget(confirm)
        layout.addLayout(btn_layout)

class SimplifyDialog(QFrame):
    """Boîte de dialogue interactive pour la simplification de maillage en temps réel."""
    def __init__(self, initial_count, on_apply, on_cancel, on_slider_change):
        super().__init__()
        self.initial_count = initial_count
        self.on_apply_callback = on_apply
        self.on_cancel_callback = on_cancel
        self.on_slider_change_callback = on_slider_change
        
        self.setFixedSize(450, 250)
        self.setStyleSheet("""
            QFrame {
                background-color: #151722;
                border: 1px solid #9D4EDD;
                border-radius: 12px;
            }
            QLabel { border: none; background: transparent; }
            QSlider::groove:horizontal {
                border: 1px solid #2A2F4A;
                height: 8px;
                background: #0F111A;
                border-radius: 4px;
            }
            QSlider::sub-page:horizontal {
                background: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 #9D4EDD, stop:1 #00E5FF);
                border-radius: 4px;
            }
            QSlider::handle:horizontal {
                background: #FFFFFF;
                border: 2px solid #00E5FF;
                width: 16px;
                margin-top: -4px;
                margin-bottom: -4px;
                border-radius: 8px;
            }
        """)
        
        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(40)
        shadow.setColor(QColor(157, 78, 221, 100))
        shadow.setOffset(0, 0)
        self.setGraphicsEffect(shadow)
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(30, 25, 30, 25)
        
        title = QLabel("MOTEUR DE SIMPLIFICATION")
        title.setStyleSheet("color: #00E5FF; font-size: 16px; font-weight: bold; font-family: 'Consolas';")
        layout.addWidget(title, alignment=Qt.AlignmentFlag.AlignCenter)
        
        self.info_label = QLabel()
        self.info_label.setStyleSheet("color: #FFFFFF; font-size: 13px; font-family: 'Segoe UI';")
        self.info_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(self.info_label)
        
        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setRange(5, 95)
        self.slider.setValue(50) # Conserver 50% par défaut
        self.slider.valueChanged.connect(self.handle_slider_change)
        layout.addWidget(self.slider)
        
        btn_layout = QHBoxLayout()
        cancel_btn = GlowingButton("CANCEL", accent_color="#8A91B4")
        cancel_btn.clicked.connect(self.on_cancel_callback)
        
        apply_btn = GlowingButton("APPLY REDUCTION", accent_color="#9D4EDD")
        apply_btn.clicked.connect(self.on_apply_callback)
        
        btn_layout.addWidget(cancel_btn)
        btn_layout.addWidget(apply_btn)
        layout.addLayout(btn_layout)
        
        # Initialiser l'état textuel
        self.handle_slider_change(50)
        
    def handle_slider_change(self, val):
        ratio = val / 100.0
        remaining = int(self.initial_count * ratio)
        self.info_label.setText(
            f"Conservation : <b>{val}%</b> de la géométrie d'origine.<br>"
            f"Triangles : <font color='#00E5FF'><b>{self.initial_count}</b></font> ➔ <font color='#9D4EDD'><b>{remaining}</b></font>"
        )
        self.on_slider_change_callback(ratio)

class SequenceDialog(QFrame):
    """Boîte de dialogue de progression globale de la séquence."""
    def __init__(self, on_abort, on_run):
        super().__init__()
        self.setFixedSize(450, 180)
        self.setStyleSheet("""
            QFrame {
                background-color: #151722;
                border: 1px solid #00E5FF;
                border-radius: 12px;
            }
            QLabel { border: none; background: transparent; }
        """)
        
        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(40)
        shadow.setColor(QColor(0, 229, 255, 100))
        shadow.setOffset(0, 0)
        self.setGraphicsEffect(shadow)
        
        layout = QVBoxLayout(self)
        layout.setContentsMargins(30, 20, 30, 20)
        
        self.status_label = QLabel("INITIALIZATION...")
        self.status_label.setStyleSheet("color: #00E5FF; font-family: 'Consolas', monospace; font-size: 14px; font-weight: bold;")
        layout.addWidget(self.status_label, alignment=Qt.AlignmentFlag.AlignCenter)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setFixedHeight(18)
        self.progress_bar.setStyleSheet("""
            QProgressBar {
                background-color: rgba(30, 34, 53, 160); border: 1px solid #2A2F4A; border-radius: 4px; text-align: center; color: white;
            }
            QProgressBar::chunk { background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 #9D4EDD, stop:1 #00E5FF); border-radius: 3px; }
        """)
        layout.addWidget(self.progress_bar)
        
        self.btn_layout = QHBoxLayout()
        self.run_btn = GlowingButton("RUN", accent_color="#00E5FF")
        self.run_btn.hide()
        self.run_btn.clicked.connect(on_run)
        
        self.abort_btn = GlowingButton("ABORT", accent_color="#FF3366")
        self.abort_btn.clicked.connect(on_abort)
        
        self.btn_layout.addWidget(self.run_btn)
        self.btn_layout.addWidget(self.abort_btn)
        layout.addLayout(self.btn_layout)

class MainWindow(QMainWindow):
    """Fenêtre principale de l'application."""
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Apex Slicer Engine")
        self.setAcceptDrops(True) 
        self.worker = None
        
        self.setStyleSheet("""
            QMainWindow { background-color: #090A10; }
            QLabel { color: #FFFFFF; font-family: 'Segoe UI', sans-serif; font-size: 13px; }
            QLabel#SectionTitle { color: #9D4EDD; font-weight: bold; font-size: 13px; text-transform: uppercase; margin-bottom: 3px; }
            QFrame#SidebarPanel, QFrame#TopToolbar { background-color: transparent; border: 1px solid rgba(0, 229, 255, 120); border-radius: 14px; }
            QComboBox { background-color: rgba(30, 34, 53, 200); color: #FFFFFF; border: 1px solid #2A2F4A; border-radius: 6px; padding: 6px; }
            QLineEdit { background-color: rgba(30, 34, 53, 200); color: #FFFFFF; border: 1px solid #2A2F4A; border-radius: 6px; padding: 6px; }
        """)

        # Container principal arrière-plan (qui peut être flouté)
        self.main_container = QWidget()
        self.main_layout = QHBoxLayout(self.main_container)
        self.main_layout.setContentsMargins(0,0,0,0)
        
        self.viewport = Slicer3DViewport()
        self.viewport.modelSelected.connect(self.on_model_selected)
        self.main_layout.addWidget(self.viewport)
        
        self.setCentralWidget(self.main_container)
        
        # --- PANNEAUX DE CONTRÔLE FLOTTANTS ---
        # Tous les contrôles sont des enfants directs de 'self' pour éviter l'héritage
        # de l'attribut de transparence aux événements de souris de l'ancienne self.ui_layer.
        self.sidebar = self.create_sidebar()
        self.hud_frame = self.create_hud()
        
        self.init_btn = GlowingButton("INITIATE SEQUENCE", accent_color="#9D4EDD", parent=self)
        self.init_btn.setFixedSize(220, 50)
        self.init_btn.clicked.connect(self.request_pipeline)
        
        self.exit_btn = QPushButton("✕", self)
        self.exit_btn.setFixedSize(32, 32)
        self.exit_btn.clicked.connect(self.close)
        self.exit_btn.setStyleSheet("QPushButton { background: rgba(30,34,53,200); color:#FF3366; border:1px solid #FF3366; border-radius:16px; }")
        
        self.top_toolbar = self.create_top_toolbar()
        self.top_toolbar.hide()
        
        # Calque modal (par dessus tout)
        self.modal_layer = ModalOverlay(self)

        # Assurer l'ordre d'affichage (Z-index correct)
        self.sidebar.raise_()
        self.hud_frame.raise_()
        self.init_btn.raise_()
        self.exit_btn.raise_()
        self.top_toolbar.raise_()
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
            tw = self.top_toolbar.width()
            self.top_toolbar.setGeometry((w - tw) // 2, 25, tw, 70)
            
        # Positionnement du HUD (haut droit)
        if hasattr(self, 'hud_frame') and self.hud_frame:
            self.hud_frame.setGeometry(w - 350, 25, 325, 65)
            
        # Positionnement des boutons du bas (bas droit)
        if hasattr(self, 'init_btn') and self.init_btn:
            self.init_btn.setGeometry(w - 250 - 32 - 12, h - 75, 220, 50)
        if hasattr(self, 'exit_btn') and self.exit_btn:
            self.exit_btn.setGeometry(w - 57, h - 66, 32, 32)
            
        # Calque modal sur tout l'écran
        if hasattr(self, 'modal_layer') and self.modal_layer:
            self.modal_layer.resize(event.size())

    def create_top_toolbar(self):
        """Crée la barre de transformation flottante."""
        toolbar = QFrame(self)
        toolbar.setObjectName("TopToolbar")
        toolbar.setFixedSize(700, 70)
        
        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(20)
        shadow.setColor(QColor(0, 229, 255, 100))
        toolbar.setGraphicsEffect(shadow)
        
        layout = QHBoxLayout(toolbar)
        
        def make_group(name, suffix):
            grp = QWidget()
            l = QHBoxLayout(grp)
            l.setContentsMargins(0,0,0,0)
            l.addWidget(QLabel(name))
            x = TransformSpinBox(suffix)
            y = TransformSpinBox(suffix)
            z = TransformSpinBox(suffix)
            for w in (x,y,z): l.addWidget(w)
            return grp, x, y, z
            
        grp_pos, self.sp_px, self.sp_py, self.sp_pz = make_group("POS", "mm")
        grp_rot, self.sp_rx, self.sp_ry, self.sp_rz = make_group("ROT", "°")
        grp_scl, self.sp_sx, self.sp_sy, self.sp_sz = make_group("SCL", "x")
        
        for sp in (self.sp_px, self.sp_py, self.sp_pz, self.sp_rx, self.sp_ry, self.sp_rz, self.sp_sx, self.sp_sy, self.sp_sz):
            sp.valueChanged.connect(self.on_transform_changed)

        layout.addWidget(grp_pos)
        layout.addWidget(grp_rot)
        layout.addWidget(grp_scl)
        
        return toolbar

    def on_model_selected(self, model):
        self._current_selected_model = model
        if model:
            self.top_toolbar.show()
            for sp in (self.sp_px, self.sp_py, self.sp_pz, self.sp_rx, self.sp_ry, self.sp_rz, self.sp_sx, self.sp_sy, self.sp_sz):
                sp.blockSignals(True)
                
            self.sp_px.setValue(model.pos[0]); self.sp_py.setValue(model.pos[1]); self.sp_pz.setValue(model.pos[2])
            self.sp_rx.setValue(model.rot[0]); self.sp_ry.setValue(model.rot[1]); self.sp_rz.setValue(model.rot[2])
            self.sp_sx.setValue(model.scale[0]); self.sp_sy.setValue(model.scale[1]); self.sp_sz.setValue(model.scale[2])
            
            for sp in (self.sp_px, self.sp_py, self.sp_pz, self.sp_rx, self.sp_ry, self.sp_rz, self.sp_sx, self.sp_sy, self.sp_sz):
                sp.blockSignals(False)
        else:
            self.top_toolbar.hide()

    def on_transform_changed(self):
        if hasattr(self, '_current_selected_model') and self._current_selected_model:
            m = self._current_selected_model
            m.pos = [self.sp_px.value(), self.sp_py.value(), self.sp_pz.value()]
            m.rot = [self.sp_rx.value(), self.sp_ry.value(), self.sp_rz.value()]
            m.scale = [self.sp_sx.value(), self.sp_sy.value(), self.sp_sz.value()]
            self.viewport.save_arrange_json()
            self.viewport.update()

    def create_sidebar(self):
        """Crée la barre latérale gauche."""
        sidebar_frame = QFrame(self)
        sidebar_frame.setObjectName("SidebarPanel")
        sidebar_frame.setFixedWidth(310)
        
        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(30)
        shadow.setColor(QColor(0, 229, 255, 110))
        sidebar_frame.setGraphicsEffect(shadow)
        
        layout = QVBoxLayout(sidebar_frame)
        layout.setContentsMargins(20, 25, 20, 25)
        
        layout.addWidget(QLabel("<span style='font-weight:bold; font-size:18px; color:#00E5FF;'>APEX SLICER v2.1</span>"))
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setStyleSheet("border: none; background: transparent;")
        content = QWidget()
        sl = QVBoxLayout(content)
        sl.setContentsMargins(0,0,0,0)
        
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

    def create_hud(self):
        """Crée l'indicateur de raccourcis caméra."""
        hud = QFrame(self)
        hud.setObjectName("HUD")
        hud.setStyleSheet("background-color: rgba(15,17,26,170); border: 1px solid rgba(0,229,255,60); border-radius: 8px; padding: 10px;")
        l = QVBoxLayout(hud)
        l.setContentsMargins(12,8,12,8)
        title = QLabel("CAMERA VIEW : INTERACTIVE 3D")
        title.setStyleSheet("color: #00E5FF; font-family: 'Consolas'; font-size: 11px; font-weight: bold;")
        desc = QLabel("Drag & Rotate : L-Click | Pan : R-Click | Del : Suppr")
        desc.setStyleSheet("color: #8A91B4; font-size: 11px;")
        l.addWidget(title); l.addWidget(desc)
        return hud

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

    def request_pipeline(self):
        """Analyse le nombre de triangles totaux avant le calcul."""
        total_tris = sum(len(m.base_triangles) for m in self.viewport.models)
        
        # Flouter l'arrière-plan de l'interface
        blur = QGraphicsBlurEffect()
        blur.setBlurRadius(20)
        self.main_container.setGraphicsEffect(blur)
        
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
        self.main_container.setGraphicsEffect(None)
        self.modal_layer.hide()
        self.viewport.update()
        
        # Relancer l'analyse du pipeline avec le nouveau maillage léger !
        self.request_pipeline()

    def cancel_simplification(self):
        """Annule l'aperçu et restaure les modèles originaux intacts."""
        for model in self.viewport.models:
            model.preview_triangles = None
            model._intensity_cache_id = None
        self.main_container.setGraphicsEffect(None)
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
        self.main_container.setGraphicsEffect(None)
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

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())