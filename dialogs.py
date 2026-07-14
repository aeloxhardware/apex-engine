"""Boîtes de dialogue modales : avertissement de surcharge, simplification interactive, séquence de traitement."""
from PyQt6.QtWidgets import (QFrame, QVBoxLayout, QHBoxLayout, QLabel,
                             QProgressBar, QSlider, QGraphicsDropShadowEffect,
                             QLineEdit, QPushButton)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor

from widgets import GlowingButton


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
        self._drag_offset = None

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

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            self._drag_offset = event.position()
            event.accept()

    def mouseMoveEvent(self, event):
        if self._drag_offset is not None:
            diff = event.position() - self._drag_offset
            self.move(self.x() + int(diff.x()), self.y() + int(diff.y()))
            event.accept()

    def mouseReleaseEvent(self, event):
        self._drag_offset = None
        event.accept()

    def handle_slider_change(self, val):
        ratio = val / 100.0
        remaining = int(self.initial_count * ratio)
        self.info_label.setText(
            f"Conservation : <b>{val}%</b> de la géométrie d'origine.<br>"
            f"Triangles : <font color='#00E5FF'><b>{self.initial_count}</b></font> ➔ <font color='#9D4EDD'><b>{remaining}</b></font>"
        )
        self.on_slider_change_callback(ratio)

class TransformDialog(QFrame):
    """Boîte de dialogue pour modifier position, rotation ou échelle."""
    def __init__(self, mode, values, on_apply, on_cancel):
        super().__init__()
        self.setFixedSize(440, 220)
        self.mode = mode
        self.on_apply_callback = on_apply

        self.setStyleSheet("""
            QFrame {
                background-color: #151722;
                border: 1px solid #00E5FF;
                border-radius: 14px;
            }
            QLabel { color: #FFFFFF; border: none; background: transparent; }
            QLineEdit {
                background-color: rgba(20, 24, 38, 240);
                border: 1px solid rgba(255,255,255,0.08);
                border-radius: 10px;
                padding: 10px 12px;
                color: #FFFFFF;
                font-size: 13px;
            }
            QPushButton {
                background-color: rgba(0, 229, 255, 0.18);
                color: #E5F4FF;
                border: 1px solid rgba(0, 229, 255, 0.40);
                border-radius: 10px;
                padding: 10px 16px;
                font-weight: 600;
            }
            QPushButton:hover { background-color: rgba(0, 229, 255, 0.30); }
        """)

        shadow = QGraphicsDropShadowEffect()
        shadow.setBlurRadius(40)
        shadow.setColor(QColor(0, 229, 255, 100))
        shadow.setOffset(0, 0)
        self.setGraphicsEffect(shadow)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(24, 24, 24, 24)
        layout.setSpacing(14)

        title_label = QLabel(f"{mode} : Modifier les valeurs")
        title_label.setStyleSheet("font-size: 16px; font-weight: bold; color: #00E5FF;")
        layout.addWidget(title_label)

        axes_frame = QFrame()
        axes_layout = QHBoxLayout(axes_frame)
        axes_layout.setContentsMargins(0, 0, 0, 0)
        axes_layout.setSpacing(12)

        self.axis_inputs = []
        for axis, value in zip(("X", "Y", "Z"), values):
            sub_frame = QFrame()
            sub_frame.setStyleSheet("background: rgba(15, 18, 30, 220); border: 1px solid rgba(255,255,255,0.08); border-radius: 10px;")
            sub_layout = QVBoxLayout(sub_frame)
            sub_layout.setContentsMargins(12, 10, 12, 10)
            sub_layout.setSpacing(6)
            label = QLabel(axis)
            field = QLineEdit(value)
            self.axis_inputs.append(field)
            sub_layout.addWidget(label)
            sub_layout.addWidget(field)
            axes_layout.addWidget(sub_frame)

        layout.addWidget(axes_frame)

        button_row = QHBoxLayout()
        apply_btn = QPushButton("APPLIQUER")
        cancel_btn = QPushButton("ANNULER")
        apply_btn.clicked.connect(self.apply)
        cancel_btn.clicked.connect(on_cancel)
        button_row.addWidget(cancel_btn)
        button_row.addWidget(apply_btn)
        layout.addLayout(button_row)

    def apply(self):
        values = [field.text().strip() for field in self.axis_inputs]
        self.on_apply_callback(self.mode, values)


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
