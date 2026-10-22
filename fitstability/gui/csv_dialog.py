"""Import dialog for CSV/text spectra: preview, delimiter, columns, axis and units."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtWidgets import (QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout, QLabel,
                               QTableWidget, QTableWidgetItem, QVBoxLayout)

from ..io import read_spectrum, sniff_table

_DELIMS = [("Auto-detect", False), ("Comma", ","), ("Tab", "\t"), ("Semicolon", ";"), ("Whitespace", None)]


class CsvImportDialog(QDialog):
    """Choose which columns of a text file are energy and intensity.

    After ``exec()`` returns accepted, ``x``, ``y`` (binding energy ascending, intensity),
    ``count_scale`` and ``description`` are set.
    """

    def __init__(self, path, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Import spectrum - {Path(path).name}")
        self.resize(640, 520)
        self.path = path
        self.x = self.y = None
        self.count_scale = 1.0
        self.description = ""
        try:
            self.info = sniff_table(path)
        except ValueError as e:
            self.info = None
            self._error = str(e)
        else:
            self._error = ""

        lay = QVBoxLayout(self)
        self.table = QTableWidget()
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        lay.addWidget(QLabel("Preview (first rows):"))
        lay.addWidget(self.table, 1)

        form = QFormLayout()
        self.delim = QComboBox()
        for label, _ in _DELIMS:
            self.delim.addItem(label)
        self.xcol, self.ycol = QComboBox(), QComboBox()
        self.axis = QComboBox()
        self.axis.addItems(["Binding energy (eV)", "Kinetic energy (eV)"])
        self.hv = QDoubleSpinBox()
        self.hv.setRange(100.0, 20000.0)
        self.hv.setDecimals(2)
        self.hv.setValue(1486.68)
        self.hv.setSuffix(" eV  (photon energy, for kinetic energy)")
        self.units = QComboBox()
        self.units.addItems(["Counts", "Counts per second"])
        self.cps_factor = QDoubleSpinBox()
        self.cps_factor.setRange(0.0001, 1e7)
        self.cps_factor.setDecimals(4)
        self.cps_factor.setValue(1.0)
        self.cps_factor.setPrefix("dwell x scans = ")
        self.cps_factor.setSuffix(" s")
        form.addRow("Delimiter", self.delim)
        form.addRow("Energy column", self.xcol)
        form.addRow("Intensity column", self.ycol)
        form.addRow("Energy axis", self.axis)
        form.addRow("Photon energy", self.hv)
        form.addRow("Intensity units", self.units)
        form.addRow("Counts per second", self.cps_factor)
        lay.addLayout(form)
        self.msg = QLabel(self._error)
        self.msg.setObjectName("danger")
        self.msg.setWordWrap(True)
        lay.addWidget(self.msg)
        self.buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttons.accepted.connect(self._accept)
        self.buttons.rejected.connect(self.reject)
        lay.addWidget(self.buttons)

        self.delim.currentIndexChanged.connect(self._reload)
        self.axis.currentIndexChanged.connect(self._sync)
        self.units.currentIndexChanged.connect(self._sync)
        self._reload()

    def _sync(self):
        self.hv.setEnabled(self.axis.currentIndex() == 1)
        self.cps_factor.setEnabled(self.units.currentIndex() == 1)

    def _delimiter(self):
        return _DELIMS[self.delim.currentIndex()][1]

    def _reload(self):
        d = self._delimiter()
        info = self.info
        if info is None:
            self.buttons.button(QDialogButtonBox.Ok).setEnabled(False)
            return
        if d is not False:                                  # re-sniff the columns for a chosen delimiter
            from ..io import _split
            rows = [r for r in (_split(ln.strip(), d) for ln in Path(self.path).read_text(errors="replace").splitlines()[:60]
                                if ln.strip()) if len(r) >= 2]
            ncol = max((len(r) for r in rows), default=2)
            preview = rows[:12]
            header = None
        else:
            ncol, preview, header = info.n_columns, info.rows, info.header
        self.table.clear()
        self.table.setColumnCount(ncol)
        self.table.setRowCount(len(preview))
        self.table.setHorizontalHeaderLabels([(header[i] if header and header[i] else f"column {i + 1}") for i in range(ncol)])
        for r, row in enumerate(preview):
            for c in range(ncol):
                self.table.setItem(r, c, QTableWidgetItem(row[c] if c < len(row) else ""))
        names = [(header[i] if header and header[i] else f"column {i + 1}") for i in range(ncol)]
        for combo, default in ((self.xcol, 0), (self.ycol, 1 if ncol > 1 else 0)):
            combo.clear()
            combo.addItems(names)
            combo.setCurrentIndex(min(default, ncol - 1))
        # a header naming the axis hints at the right choice
        if header:
            for i, h in enumerate(header):
                low = (h or "").lower()
                if "kinetic" in low:
                    self.xcol.setCurrentIndex(i)
                    self.axis.setCurrentIndex(1)
        self._sync()

    def _accept(self):
        try:
            x, y = read_spectrum(self.path, self.xcol.currentIndex(), self.ycol.currentIndex(), self._delimiter())
        except ValueError as e:
            self.msg.setText(str(e))
            return
        if self.axis.currentIndex() == 1:
            x = self.hv.value() - x
            order = np.argsort(x)
            x, y = x[order], y[order]
        self.x, self.y = x, y
        self.count_scale = self.cps_factor.value() if self.units.currentIndex() == 1 else 1.0
        self.description = f"{len(x)} points, {x.min():.2f}-{x.max():.2f} eV"
        self.accept()
