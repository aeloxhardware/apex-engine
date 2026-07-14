"""Viewport 3D interactif : rendu du plateau d'impression, interactions souris/clavier, persistance JSON."""
import os
import json
import math
from PyQt6.QtWidgets import QFrame, QGraphicsScene, QGraphicsPixmapItem
from PyQt6.QtCore import Qt, QPointF, QRect, QRectF, pyqtSignal, QTimer
from PyQt6.QtGui import QPainter, QColor, QPen, QPolygonF, QImage, QPainterPath, QPixmap

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

from config import MODELS_DIR, ARRANGE_JSON
from model3d import Model3D


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
        self._gizmo_hit_targets = []
        self._gizmo_active_axis = None
        self._gizmo_drag_start = None
        self._gizmo_initial_pos = None
        self._gizmo_dirty = False
        self._viewport_intro_active = False
        self._viewport_intro_phase = 0
        self._viewport_intro_tick = 0
        self._viewport_intro_point_count = 0
        self._viewport_intro_line_progress = 0.0
        self._viewport_intro_box_progress = 0.0
        self._viewport_intro_text_index = 0
        self._viewport_intro_text = "FORMATION DES BOÎTES"
        self._viewport_intro_fade = 1.0
        self._viewport_intro_timer = QTimer(self)
        self._viewport_intro_timer.setInterval(36)
        self._viewport_intro_timer.timeout.connect(self._advance_intro)
        self._setup_intro_geometry()
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setStyleSheet("""
            QFrame#ViewportPanel {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                    stop:0 #0A0E17, stop:1 #06070C);
                border: 1px solid rgba(255,255,255,0.06);
                border-radius: 24px;
            }
        """)
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
            click_point = event.position()

            # Priorité aux zones de gizmo cliquables
            for axis, poly in self._gizmo_hit_targets:
                if poly.containsPoint(click_point, Qt.FillRule.OddEvenFill):
                    self._gizmo_active_axis = axis
                    self._gizmo_drag_start = event.position()
                    self._gizmo_dirty = False
                    selected = next((m for m in self.models if m.selected), None)
                    self._gizmo_initial_pos = selected.pos.copy() if selected else None
                    self.update()
                    return

            # Recherche de collision (Hit Test)
            clicked_model = None
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
        if self._gizmo_active_axis is not None and event.buttons() & Qt.MouseButton.LeftButton:
            selected = next((m for m in self.models if m.selected), None)
            if selected and self._gizmo_drag_start and self._gizmo_initial_pos is not None:
                delta = event.position() - self._gizmo_drag_start
                screen_axis = self._gizmo_screen_direction(self._gizmo_active_axis)
                if screen_axis is not None:
                    move_amount = delta.x() * screen_axis[0] + delta.y() * screen_axis[1]
                    move_amount *= 0.16 / max(0.2, self.zoom)
                    axis_idx = {'X': 0, 'Y': 1, 'Z': 2}[self._gizmo_active_axis]
                    selected.pos[axis_idx] = self._gizmo_initial_pos[axis_idx] + move_amount
                    self._gizmo_dirty = True
                    self.update()
            return

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
        if self._gizmo_dirty:
            self.save_arrange_json()
            self._gizmo_dirty = False
        self._gizmo_active_axis = None
        self._gizmo_drag_start = None
        self._gizmo_initial_pos = None
        self.update()
        
    def _setup_intro_geometry(self):
        self._intro_point_norms = [
            (-0.34, -0.18), (-0.22, -0.12), (-0.08, -0.28), (0.06, -0.14), (0.22, -0.32),
            (-0.30, 0.12), (-0.16, 0.22), (0.02, 0.18), (0.24, 0.04), (0.38, 0.20)
        ]
        self._intro_h_lines = [(-0.36, -0.18, 0.72), (-0.32, 0.0, 0.64), (-0.28, 0.18, 0.56)]
        self._intro_v_lines = [(-0.18, -0.36, 0.72), (0.06, -0.32, 0.64), (0.28, -0.28, 0.56)]
        self._intro_boxes = [
            (-0.34, -0.26, 0.32, 0.20),
            (0.04, -0.30, 0.28, 0.18),
            (-0.30, 0.14, 0.22, 0.20),
            (0.10, 0.08, 0.26, 0.24)
        ]

    def start_intro_sequence(self):
        self._viewport_intro_active = True
        self._viewport_intro_phase = 0
        self._viewport_intro_tick = 0
        self._viewport_intro_point_count = 0
        self._viewport_intro_line_progress = 0.0
        self._viewport_intro_box_progress = 0.0
        self._viewport_intro_text_index = 0
        self._viewport_intro_fade = 1.0
        self._viewport_intro_timer.start()
        self.update()

    def _advance_intro(self):
        self._viewport_intro_tick += 1
        if self._viewport_intro_phase == 0:
            self._viewport_intro_point_count = min(len(self._intro_point_norms), self._viewport_intro_tick)
            if self._viewport_intro_tick >= 10:
                self._viewport_intro_phase = 1
                self._viewport_intro_tick = 0
        elif self._viewport_intro_phase == 1:
            self._viewport_intro_line_progress = min(1.0, self._viewport_intro_tick / 16.0)
            if self._viewport_intro_tick >= 16:
                self._viewport_intro_phase = 2
                self._viewport_intro_tick = 0
        elif self._viewport_intro_phase == 2:
            self._viewport_intro_line_progress = 1.0
            self._viewport_intro_box_progress = min(1.0, self._viewport_intro_tick / 18.0)
            if self._viewport_intro_tick >= 18:
                self._viewport_intro_phase = 3
                self._viewport_intro_tick = 0
        elif self._viewport_intro_phase == 3:
            self._viewport_intro_text_index = min(len(self._viewport_intro_text), self._viewport_intro_tick)
            if self._viewport_intro_tick >= len(self._viewport_intro_text) + 8:
                self._viewport_intro_phase = 4
                self._viewport_intro_tick = 0
        elif self._viewport_intro_phase == 4:
            self._viewport_intro_fade = max(0.0, 1.0 - self._viewport_intro_tick / 10.0)
            if self._viewport_intro_tick >= 10:
                self._viewport_intro_active = False
                self._viewport_intro_timer.stop()
        self.update()
        selected = next((m for m in self.models if m.selected), None)
        if not selected:
            return None

        center = [sum(p[i] for tri in selected.get_world_triangles() for p in tri) /
                  len([p for tri in selected.get_world_triangles() for p in tri]) for i in range(3)]
        axis_offset = {'X': [1,0,0], 'Y': [0,1,0], 'Z': [0,0,1]}[axis]
        base = self._project_point(center)
        end = self._project_point([center[0] + axis_offset[0] * 100, center[1] + axis_offset[1] * 100, center[2] + axis_offset[2] * 100])
        if not base or not end:
            return None
        dx = end[0] - base[0]
        dy = end[1] - base[1]
        length = math.hypot(dx, dy)
        if length == 0:
            return None
        return (dx / length, dy / length)

    def _project_point(self, point):
        w, h = self.width(), self.height()
        if w <= 0 or h <= 0:
            return None
        cx, cy = w / 2 + self.pan_x, h / 2 + self.pan_y
        yaw_rad, pitch_rad = math.radians(self.yaw), math.radians(self.pitch)
        cos_y, sin_y = math.cos(yaw_rad), math.sin(yaw_rad)
        cos_p, sin_p = math.cos(pitch_rad), math.sin(pitch_rad)
        x1 = point[0] * cos_y - point[1] * sin_y
        y1 = point[0] * sin_y + point[1] * cos_y
        x2 = x1
        y2 = y1 * cos_p - point[2] * sin_p
        z2 = y1 * sin_p + point[2] * cos_p
        return (cx + x2 * self.zoom, cy - z2 * self.zoom)

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

        # Fond de scène et grille légère pour un rendu plus soigné
        img_painter.setPen(QPen(QColor(255, 255, 255, 18), 1))
        for x in range(0, w, 28):
            img_painter.drawLine(x, 0, x, h)
        for y in range(0, h, 28):
            img_painter.drawLine(0, y, w, y)
        img_painter.setPen(QPen(QColor(0, 229, 255, 35), 1.5))
        img_painter.drawEllipse(int(w/2 - 140), int(h/2 - 140), 280, 280)

        if self._viewport_intro_active:
            self._draw_intro_sequence(img_painter, w, h)
            img_painter.end()

            painter = QPainter(self)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.drawImage(0, 0, scene_image)
            return

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

        if self.models and any(m.selected for m in self.models):
            self._draw_transform_gizmo(img_painter, w, h)

        img_painter.end()

        # Peindre la scène finale sur l'écran
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.drawImage(0, 0, scene_image)

        # --- EFFET DE FLOU DE FOND POUR LES PANNEAUX SOUSTRAYANTS ---
        if not hasattr(self, '_blur_targets'):
            # Rechercher dynamiquement les panneaux à flouter dans la fenêtre parent
            self._blur_targets = [c for c in self.window().findChildren(QFrame) 
                                  if c.objectName() in ("SidebarPanel", "TopToolbar", "TransformPanel")]

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

    def _draw_intro_sequence(self, painter, width, height):
        center_x = width * 0.5 + self.pan_x * 0.15
        center_y = height * 0.52 + self.pan_y * 0.15
        center = QPointF(center_x, center_y)
        span = min(width, height) * 0.48

        # Points initiaux
        painter.setPen(Qt.PenStyle.NoPen)
        for idx in range(self._viewport_intro_point_count):
            nx, ny = self._intro_point_norms[idx]
            px = center_x + nx * span
            py = center_y + ny * span
            radius = 3 + (idx % 2)
            painter.setBrush(QColor(0, 229, 255, 230))
            painter.drawEllipse(QPointF(px, py), radius, radius)

        # Lignes horizontales en expansion
        if self._viewport_intro_phase >= 1 or self._viewport_intro_line_progress > 0:
            painter.setPen(QPen(QColor(0, 180, 255, 190), 2, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap))
            for x_norm, y_norm, length_norm in self._intro_h_lines:
                y = center_y + y_norm * span
                x_start = center_x + x_norm * span
                x_end = x_start + length_norm * span * self._viewport_intro_line_progress
                painter.drawLine(QPointF(x_start, y), QPointF(x_end, y))

        # Lignes verticales en expansion
        if self._viewport_intro_phase >= 2 or self._viewport_intro_line_progress > 0.5:
            painter.setPen(QPen(QColor(0, 229, 255, 160), 2, Qt.PenStyle.SolidLine, Qt.PenCapStyle.RoundCap))
            line_progress = max(0.0, min(1.0, (self._viewport_intro_line_progress - 0.4) / 0.6))
            for x_norm, y_norm, length_norm in self._intro_v_lines:
                x = center_x + x_norm * span
                y_start = center_y + y_norm * span
                y_end = y_start + length_norm * span * line_progress
                painter.drawLine(QPointF(x, y_start), QPointF(x, y_end))

        # Boîtes apparentes
        if self._viewport_intro_phase >= 2:
            painter.setPen(QPen(QColor(0, 229, 255, 210), 2, Qt.PenStyle.SolidLine, Qt.PenCapStyle.SquareCap))
            for box in self._intro_boxes:
                bx, by, bw, bh = box
                x = center_x + bx * span
                y = center_y + by * span
                w_box = bw * span
                h_box = bh * span
                progress = self._viewport_intro_box_progress
                self._draw_partial_rect(painter, QRectF(x, y, w_box, h_box), progress)

        # Texte de diagnostic futuriste
        if self._viewport_intro_phase >= 3:
            painter.setPen(QColor(115, 235, 255, int(220 * self._viewport_intro_fade)))
            font = painter.font()
            font.setPointSize(16)
            font.setFamily('Consolas')
            painter.setFont(font)
            text = self._viewport_intro_text[:self._viewport_intro_text_index]
            painter.drawText(QRectF(center_x - span*0.45, center_y + span*0.34, span*0.9, 40), Qt.AlignmentFlag.AlignCenter, text)

        if self._viewport_intro_phase == 4:
            painter.setBrush(QColor(5, 10, 18, int(180 * (1.0 - self._viewport_intro_fade))))
            painter.setPen(Qt.PenStyle.NoPen)
            painter.drawRect(0, 0, width, height)

    def _draw_partial_rect(self, painter, rect, progress):
        x, y, w, h = rect.x(), rect.y(), rect.width(), rect.height()
        if progress <= 0:
            return
        length = 2 * (w + h)
        target = length * progress
        # top
        top_len = min(target, w)
        painter.drawLine(QPointF(x, y), QPointF(x + top_len, y))
        target -= top_len
        if target <= 0: return
        # right
        right_len = min(target, h)
        painter.drawLine(QPointF(x + w, y), QPointF(x + w, y + right_len))
        target -= right_len
        if target <= 0: return
        # bottom
        bottom_len = min(target, w)
        painter.drawLine(QPointF(x + w, y + h), QPointF(x + w - bottom_len, y + h))
        target -= bottom_len
        if target <= 0: return
        # left
        left_len = min(target, h)
        painter.drawLine(QPointF(x, y + h), QPointF(x, y + h - left_len))

    def _draw_transform_gizmo(self, painter, width, height):
        selected = next((m for m in self.models if m.selected), None)
        if not selected:
            return

        # Calculer le point central du modèle pour placer le gizmo
        all_points = [pt for tri in selected.get_world_triangles() for pt in tri]
        if not all_points:
            return
        center = [sum(p[i] for p in all_points) / len(all_points) for i in range(3)]

        def project(x, y, z):
            cx, cy = width / 2 + self.pan_x, height / 2 + self.pan_y
            yaw_rad, pitch_rad = math.radians(self.yaw), math.radians(self.pitch)
            cos_y, sin_y = math.cos(yaw_rad), math.sin(yaw_rad)
            cos_p, sin_p = math.cos(pitch_rad), math.sin(pitch_rad)
            x1 = x * cos_y - y * sin_y
            y1 = x * sin_y + y * cos_y
            x2 = x1
            y2 = y1 * cos_p - z * sin_p
            z2 = y1 * sin_p + z * cos_p
            return cx + x2 * self.zoom, cy - z2 * self.zoom, y2

        base = project(*center)
        axis_len = 100 / max(1.0, self.zoom)
        axes = {
            'X': ([center[0] + axis_len, center[1], center[2]], QColor(255, 95, 95)),
            'Y': ([center[0], center[1] + axis_len, center[2]], QColor(95, 255, 120)),
            'Z': ([center[0], center[1], center[2] + axis_len], QColor(120, 170, 255)),
        }

        self._gizmo_hit_targets = []

        for label, (end_pt, color) in axes.items():
            end = project(*end_pt)
            painter.setPen(QPen(color, 3.5))
            painter.drawLine(QPointF(base[0], base[1]), QPointF(end[0], end[1]))

            arrow = QPolygonF([
                QPointF(end[0], end[1]),
                QPointF(end[0] - 12, end[1] - 6),
                QPointF(end[0] - 12, end[1] + 6)
            ])
            painter.setBrush(color)
            painter.setPen(QPen(color, 1.8))
            painter.drawPolygon(arrow)

            label_pos = QPointF(end[0] + 10, end[1] + 10)
            painter.setPen(QPen(QColor(255, 255, 255, 220)))
            painter.setFont(painter.font())
            painter.drawText(label_pos, label)

            # Build a larger hit area around the axis arrow
            hit_area = QPolygonF([
                QPointF(base[0], base[1]),
                QPointF(end[0] - 14, end[1] - 10),
                QPointF(end[0] + 18, end[1] + 10)
            ])
            self._gizmo_hit_targets.append((label, hit_area))
