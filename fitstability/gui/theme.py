"""Light and dark themes: Qt palette, a small style sheet, and matching matplotlib colours."""
from __future__ import annotations

from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import QApplication

THEMES = {
    "light": dict(bg="#f4f5f7", surface="#ffffff", border="#cfd4dc", text="#1c2430", muted="#667085",
                  accent="#1f6feb", accent_text="#ffffff", warn="#b54708", ok="#067647", danger="#b42318"),
    "dark": dict(bg="#161b22", surface="#1f2630", border="#3a4252", text="#e6e9ef", muted="#98a2b3",
                 accent="#4c8dff", accent_text="#0b1220", warn="#f79009", ok="#32d583", danger="#f97066"),
}


def mpl_rc(name: str) -> dict:
    t = THEMES[name]
    return {
        "figure.facecolor": t["surface"], "axes.facecolor": t["surface"], "savefig.facecolor": t["surface"],
        "axes.edgecolor": t["muted"], "axes.labelcolor": t["text"], "text.color": t["text"],
        "xtick.color": t["text"], "ytick.color": t["text"], "grid.color": t["border"],
        "legend.labelcolor": t["text"],
    }


def apply(app: QApplication, name: str) -> dict:
    """Apply the theme to the application; returns its colour tokens."""
    t = THEMES[name]
    app.setStyle("Fusion")
    pal = QPalette()
    for role, key in ((QPalette.Window, "bg"), (QPalette.Base, "surface"), (QPalette.AlternateBase, "bg"),
                      (QPalette.Button, "surface"), (QPalette.ToolTipBase, "surface"),
                      (QPalette.WindowText, "text"), (QPalette.Text, "text"), (QPalette.ButtonText, "text"),
                      (QPalette.ToolTipText, "text"), (QPalette.PlaceholderText, "muted"),
                      (QPalette.Highlight, "accent"), (QPalette.HighlightedText, "accent_text"),
                      (QPalette.Link, "accent")):
        pal.setColor(role, QColor(t[key]))
    for role in (QPalette.WindowText, QPalette.Text, QPalette.ButtonText):
        pal.setColor(QPalette.Disabled, role, QColor(t["muted"]))
    app.setPalette(pal)
    app.setStyleSheet(f"""
        QGroupBox {{ border: 1px solid {t['border']}; border-radius: 6px; margin-top: 10px; padding-top: 8px; }}
        QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 4px; font-weight: 600; }}
        QPushButton {{ padding: 5px 12px; border: 1px solid {t['border']}; border-radius: 5px; background: {t['surface']}; }}
        QPushButton:hover {{ border-color: {t['accent']}; }}
        QPushButton:disabled {{ color: {t['muted']}; }}
        QPushButton#primary {{ background: {t['accent']}; color: {t['accent_text']}; border-color: {t['accent']}; font-weight: 600; }}
        QPushButton#primary:disabled {{ background: {t['border']}; border-color: {t['border']}; color: {t['muted']}; }}
        QTabWidget::pane {{ border: 1px solid {t['border']}; }}
        QHeaderView::section {{ background: {t['bg']}; border: 0; border-bottom: 1px solid {t['border']}; padding: 4px; }}
        QLabel#muted {{ color: {t['muted']}; }}
        QLabel#warn {{ color: {t['warn']}; }}
        QLabel#danger {{ color: {t['danger']}; }}
        QLabel#stale {{ background: {t['warn']}; color: {t['accent_text']}; padding: 3px 8px; border-radius: 3px; }}
    """)
    return t
