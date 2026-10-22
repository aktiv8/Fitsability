"""Results of a batch run: one row per block, with a CSV export."""
from __future__ import annotations

from PySide6.QtWidgets import (QDialog, QDialogButtonBox, QFileDialog, QHeaderView, QPushButton, QTableWidget,
                               QTableWidgetItem, QVBoxLayout, QLabel)

from .results import SUMMARY_HEADER, write_summary_csv


class BatchResultsDialog(QDialog):
    def __init__(self, rows: list[list], skipped: list[str], folder: str | None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Batch results")
        self.resize(820, 420)
        self.rows = rows
        lay = QVBoxLayout(self)
        n_multi = sum(1 for r in rows if r[1] == "multiple minima")
        lay.addWidget(QLabel(f"{len(rows)} fits analysed; {n_multi} show multiple minima."
                             + (f" Figures saved in {folder}." if folder else "")))
        table = QTableWidget(len(rows), len(SUMMARY_HEADER))
        table.setHorizontalHeaderLabels(SUMMARY_HEADER)
        table.setEditTriggers(QTableWidget.NoEditTriggers)
        for r, row in enumerate(rows):
            for c, v in enumerate(row):
                table.setItem(r, c, QTableWidgetItem(str(v)))
        table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        lay.addWidget(table, 1)
        if skipped:
            note = QLabel("Skipped: " + "; ".join(skipped[:6]) + (" ..." if len(skipped) > 6 else ""))
            note.setObjectName("muted")
            note.setWordWrap(True)
            lay.addWidget(note)
        save = QPushButton("Save summary CSV...")
        save.clicked.connect(self._save)
        buttons = QDialogButtonBox(QDialogButtonBox.Close)
        buttons.rejected.connect(self.reject)
        buttons.addButton(save, QDialogButtonBox.ActionRole)
        lay.addWidget(buttons)

    def _save(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save summary", "batch-summary.csv", "CSV (*.csv)")
        if path:
            write_summary_csv(self.rows, path)
