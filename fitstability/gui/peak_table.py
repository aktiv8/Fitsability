"""Editable peak table: one row per peak, constraints as short text rules."""
from __future__ import annotations

import math
from dataclasses import replace

from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt, Signal
from PySide6.QtGui import QBrush, QColor

from ..model import PeakSpec
from . import protocol_model as pm

COLUMNS = ["Name", "Shape", "Position (eV)", "Position rule", "FWHM (eV)", "FWHM rule", "Area rule",
           "Pos min", "Pos max", "FWHM min", "FWHM max"]
_TIPS = {
    1: "Empty = pseudo-Voigt. Or a CasaXPS shape: GL(30), SGL(30), LA(1.1,1.9,7), LF(1,1,75,200)... "
       "or a loaded reference spectrum: %&<>LS-CE4",
    3: pm.RULE_HELP["pos"], 5: pm.RULE_HELP["fwhm"], 6: pm.RULE_HELP["area"],
    7: "Lower bound of the position (empty = nominal - 3 eV)", 8: "Upper bound of the position (empty = nominal + 3 eV)",
    9: "Lower bound of the FWHM", 10: "Upper bound of the FWHM",
}


class PeakTableModel(QAbstractTableModel):
    """Holds a list of :class:`PeakSpec`; edits replace the spec (they are immutable)."""

    peaksChanged = Signal()
    problem = Signal(str)          # last edit's error, '' when fine

    def __init__(self, parent=None):
        super().__init__(parent)
        self._peaks: list[PeakSpec] = []
        self.resolver = None        # callable(shape) -> ReferenceSpectrum | None, for CasaXPS library shapes (%&<>LS-...)

    # -- API used by the window -------------------------------------------------------------
    def peaks(self) -> list[PeakSpec]:
        return list(self._peaks)

    def set_peaks(self, peaks):
        self.beginResetModel()
        self._peaks = list(peaks)
        self.endResetModel()
        self.peaksChanged.emit()

    def add_peak(self, pos: float, fwhm: float = 1.2):
        self._insert(PeakSpec(pm.unique_name(self._peaks), float(pos), float(fwhm)))

    def duplicate(self, row: int):
        if 0 <= row < len(self._peaks):
            p = self._peaks[row]
            self._insert(replace(p, name=pm.unique_name(self._peaks), area_tie=None, fwhm_tie=None, pos_tie=None,
                                 pos=p.pos + 0.5))

    def remove(self, row: int):
        if 0 <= row < len(self._peaks):
            self.set_peaks(pm.remove_peak(self._peaks, row))

    def _insert(self, spec):
        n = len(self._peaks)
        self.beginInsertRows(QModelIndex(), n, n)
        self._peaks.append(spec)
        self.endInsertRows()
        self.peaksChanged.emit()

    # -- Qt model ---------------------------------------------------------------------------
    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self._peaks)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(COLUMNS)

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if orientation == Qt.Horizontal:
            if role == Qt.DisplayRole:
                return COLUMNS[section]
            if role == Qt.ToolTipRole:
                return _TIPS.get(section)
        elif role == Qt.DisplayRole:
            return section + 1
        return None

    def flags(self, index):
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable

    def _text(self, p: PeakSpec, c: int) -> str:
        if c == 0:
            return p.name
        if c == 1:
            return p.shape
        if c == 2:
            return f"{p.pos:.4f}"
        if c == 3:
            return pm.format_rule(p, "pos")
        if c == 4:
            return f"{p.fwhm:.4f}"
        if c == 5:
            return pm.format_rule(p, "fwhm")
        if c == 6:
            return pm.format_rule(p, "area")
        if c in (7, 8):
            b = p.pos_bounds
            return "" if b is None else f"{b[c - 7]:.4f}"
        return f"{p.fwhm_bounds[c - 9]:.4g}"

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid():
            return None
        p, c = self._peaks[index.row()], index.column()
        if role in (Qt.DisplayRole, Qt.EditRole):
            t = self._text(p, c)
            return ("pseudo-Voigt" if c == 1 and not t and role == Qt.DisplayRole else t)
        if role == Qt.ToolTipRole:
            return _TIPS.get(c)
        if role == Qt.ForegroundRole and c in (3, 5, 6):
            return QBrush(QColor("#888888")) if self._text(p, c) == "free" else None
        return None

    def setData(self, index, value, role=Qt.EditRole):
        if role != Qt.EditRole or not index.isValid():
            return False
        r, c = index.row(), index.column()
        p = self._peaks[r]
        text = str(value).strip()
        names = [q.name for q in self._peaks]
        try:
            if c == 0:
                if not text:
                    raise ValueError("a peak needs a name")
                if text != p.name and text in names:
                    raise ValueError(f"there is already a peak named {text!r}")
                self._peaks = pm.rename_peak(self._peaks, r, text)
                self.dataChanged.emit(self.index(0, 0), self.index(len(self._peaks) - 1, len(COLUMNS) - 1))
            else:
                self._peaks[r] = replace(p, **self._changes(p, c, text, names))
                self.dataChanged.emit(index, index)
                if c in (3, 5, 6):   # a rule changes how the other columns read
                    self.dataChanged.emit(self.index(r, 0), self.index(r, len(COLUMNS) - 1))
        except ValueError as e:
            self.problem.emit(str(e))
            return False
        self.problem.emit("")
        self.peaksChanged.emit()
        return True

    def _changes(self, p: PeakSpec, c: int, text: str, names) -> dict:
        def num():
            try:
                v = float(text)
            except ValueError:
                raise ValueError(f"{text!r} is not a number") from None
            if not math.isfinite(v):
                raise ValueError("value must be finite")
            return v
        if c == 1:
            if text.startswith("%&") or text.upper().startswith("LS-"):
                shape = text if text.startswith("%&") else "%&<>" + text
                ref = self.resolver(shape) if self.resolver else None
                if ref is None:
                    raise ValueError(f"the reference spectrum {shape.lstrip('%&<>')!r} is not loaded "
                                     "(File > Reference spectra...)")
                return dict(shape=shape, reference=ref, fix_fwhm=True, fwhm_tie=None, width_group="")
            return dict(shape="" if text.lower() in ("", "pv", "pseudo-voigt") else text, reference=None)
        if c == 2:
            return dict(pos=num())
        if c == 4:
            v = num()
            if v <= 0:
                raise ValueError("FWHM must be positive")
            return dict(fwhm=v)
        if c == 3:
            return pm.parse_rule("pos", text, names, p.name)
        if c == 5:
            return pm.parse_rule("fwhm", text, names, p.name)
        if c == 6:
            return pm.parse_rule("area", text, names, p.name)
        if c in (7, 8):
            if not text:
                return dict(pos_bounds=None)
            lo, hi = p.pos_bounds or (p.pos - p.pos_halfwindow, p.pos + p.pos_halfwindow)
            lo, hi = (num(), hi) if c == 7 else (lo, num())
            if lo >= hi:
                raise ValueError("position min must be below max")
            return dict(pos_bounds=(lo, hi))
        lo, hi = p.fwhm_bounds
        lo, hi = (num(), hi) if c == 9 else (lo, num())
        if not 0 < lo < hi:
            raise ValueError("FWHM bounds need 0 < min < max")
        return dict(fwhm_bounds=(lo, hi))
