"""A matplotlib canvas with its toolbar and a 'stale' banner, themed to the application."""
from __future__ import annotations

import matplotlib
from matplotlib import rc_context
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget

from . import theme


class PlotPane(QWidget):
    def __init__(self, theme_name: str = "light", parent=None):
        super().__init__(parent)
        self.theme_name = theme_name
        self.fig = Figure(figsize=(6, 4))
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        self.banner = QLabel("Model, data or window changed since this result: run it again.")
        self.banner.setObjectName("stale")
        self.banner.hide()
        self.hint = QLabel("")
        self.hint.setObjectName("muted")
        self.hint.setWordWrap(True)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.banner)
        lay.addWidget(self.canvas, 1)
        lay.addWidget(self.hint)
        self._apply_face()
        self.has_result = False

    def _apply_face(self):
        self.fig.set_facecolor(theme.mpl_rc(self.theme_name)["figure.facecolor"])

    def set_theme(self, name: str):
        self.theme_name = name
        self._apply_face()
        self.canvas.draw_idle()

    def draw(self, fn, *args, **kwargs):
        """Call a ``fitstability.plots`` function with this pane's figure, in the pane's theme."""
        with rc_context(theme.mpl_rc(self.theme_name)):
            fn(*args, fig=self.fig, **kwargs)
        self._apply_face()
        self.canvas.draw_idle()
        self.has_result = True
        self.banner.hide()

    def show_message(self, text: str):
        """Replace the figure with a centred message."""
        with rc_context(theme.mpl_rc(self.theme_name)):
            self.fig.clear()
            self.fig.text(0.5, 0.5, text, ha="center", va="center", color=matplotlib.rcParams["text.color"])
        self._apply_face()
        self.canvas.draw_idle()
        self.has_result = False
        self.banner.hide()

    def set_stale(self, stale: bool):
        self.banner.setVisible(bool(stale) and self.has_result)
