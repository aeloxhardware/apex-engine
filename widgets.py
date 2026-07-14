"""Composants d'interface réutilisables : boutons et champs numériques stylisés."""
from PyQt6.QtWidgets import QPushButton, QGraphicsDropShadowEffect, QDoubleSpinBox, QLabel
from PyQt6.QtCore import QEasingCurve, QPropertyAnimation, Qt, QTimer
from PyQt6.QtGui import QColor


class GlowingButton(QPushButton):
    """Bouton stylisé avec halo lumineux néon réactif et animation douce."""
    def __init__(self, text, accent_color="#00E5FF", parent=None):
        super().__init__(text, parent)
        self.accent_color = accent_color
        self.shadow = QGraphicsDropShadowEffect(self)
        self.shadow.setBlurRadius(20)
        self.shadow.setColor(QColor(accent_color))
        self.shadow.setOffset(0, 0)
        self.setGraphicsEffect(self.shadow)

        self.setStyleSheet(f"""
            QPushButton {{
                background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 rgba(25, 28, 42, 240), stop:1 rgba(12, 15, 25, 240));
                color: {accent_color};
                border: 1px solid {accent_color};
                border-radius: 8px;
                padding: 10px 20px;
                font-family: 'Segoe UI', sans-serif; font-weight: 600; font-size: 13px;
            }}
            QPushButton:hover {{
                background: qlineargradient(x1:0, y1:0, x2:1, y2:1, stop:0 rgba(35, 40, 60, 250), stop:1 rgba(18, 22, 35, 250));
                border-color: #7DE2FF;
            }}
            QPushButton:pressed {{
                padding-top: 11px;
                padding-bottom: 9px;
            }}
            QPushButton:disabled {{ color: #2A2F4A; border-color: #2A2F4A; background-color: rgba(15, 16, 22, 220); }}
        """)

        self._hover_anim = QPropertyAnimation(self.shadow, b"blurRadius", self)
        self._hover_anim.setDuration(140)
        self._hover_anim.setEasingCurve(QEasingCurve.Type.OutCubic)

    def enterEvent(self, event):
        self._hover_anim.stop()
        self._hover_anim.setStartValue(self.shadow.blurRadius())
        self._hover_anim.setEndValue(32)
        self._hover_anim.start()
        super().enterEvent(event)

    def leaveEvent(self, event):
        self._hover_anim.stop()
        self._hover_anim.setStartValue(self.shadow.blurRadius())
        self._hover_anim.setEndValue(20)
        self._hover_anim.start()
        super().leaveEvent(event)


class AnimatedRevealLabel(QLabel):
    """Label de texte avec apparition progressive caractère par caractère."""
    def __init__(self, text="", parent=None):
        super().__init__(text, parent)
        self.full_text = text
        self._reveal_index = 0
        self._reveal_timer = QTimer(self)
        self._reveal_timer.setInterval(36)
        self._reveal_timer.timeout.connect(self._reveal_step)
        self.setText("")

    def setRevealText(self, text):
        self.full_text = text
        self._reveal_index = 0
        self.setText("")

    def start_reveal(self, delay_ms=0):
        self._reveal_timer.stop()
        self._reveal_index = 0
        self.setText("")
        if delay_ms > 0:
            QTimer.singleShot(delay_ms, self._reveal_timer.start)
        else:
            self._reveal_timer.start()

    def _reveal_step(self):
        if self._reveal_index >= len(self.full_text):
            self._reveal_timer.stop()
            return
        self._reveal_index += 1
        self.setText(self.full_text[:self._reveal_index])
        if self._reveal_index >= len(self.full_text):
            self._reveal_timer.stop()


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
                background-color: rgba(15, 17, 26, 220);
                color: #FFFFFF;
                border: 1px solid rgba(255,255,255,0.10);
                border-radius: 6px;
                padding: 5px 7px;
                min-width: 70px;
            }
            QDoubleSpinBox:focus { border: 1px solid #00E5FF; }
        """
)
