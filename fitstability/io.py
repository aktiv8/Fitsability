"""Reading spectra and protocols, and building synthetic test spectra.

Spectra: delimited text/CSV with a chosen x (energy) and y (intensity) column. :func:`sniff_table`
previews a file and guesses its delimiter/header; :func:`read_spectrum` reads two columns (the
first two by default), skipping header and comment lines. VAMAS files are read by
:mod:`fitstability.vamas`.

Protocol files are JSON::

    {"background": "shirley",
     "count_scale": 1.0,
     "peaks": [{"name": "C1", "pos": 285.0, "fwhm": 1.25},
               {"name": "C2", "pos": 286.5, "fwhm": 1.25, "width_group": "main"},
               {"name": "Sat", "pos": 291.5, "fwhm": 3.5, "fix_pos": true},
               {"name": "D", "pos": 291.5, "fwhm": 1.0, "shape": "GL(30)",
                "area_tie": ["C1", 0.5], "fwhm_tie": ["C1", 1.2], "pos_tie": ["C1", 5.0]}]}

``background`` is ``"none"``, ``"linear"``, ``"shirley"`` or, for a protocol taken from a CasaXPS fit,
``{"casa": "Shirley", "avg": 1, "params": [...], "scale": 1.0}``. :func:`save_protocol` writes this.
"""
from __future__ import annotations

import csv
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .lineshapes import pseudo_voigt
from .model import PeakSpec, Protocol

_PEAK_FIELDS = {"name", "pos", "fwhm", "eta", "shape", "fix_pos", "fix_fwhm", "width_group", "area_tie",
                "fwhm_tie", "pos_tie", "pos_bounds", "fwhm_bounds", "area_bounds", "pos_halfwindow"}
_TUPLE_FIELDS = ("area_tie", "fwhm_tie", "pos_tie", "pos_bounds", "fwhm_bounds", "area_bounds")


# -- spectra ---------------------------------------------------------------------------------
@dataclass
class TableInfo:
    """A preview of a delimited text file."""
    delimiter: str | None            # None: whitespace
    n_columns: int
    header: list[str] | None         # column names if the file has a header row
    rows: list[list[str]] = field(default_factory=list)   # first rows after the header, for display
    n_numeric: int = 0               # rows with at least n_columns numeric cells


def _split(line: str, delimiter: str | None) -> list[str]:
    if delimiter is None:
        return line.split()
    return next(csv.reader([line], delimiter=delimiter))


def _is_num(t: str) -> bool:
    try:
        float(t)
        return True
    except ValueError:
        return False


def sniff_table(path, max_rows: int = 12) -> TableInfo:
    """Guess delimiter, header and column count of a text spectrum."""
    text = Path(path).read_text(errors="replace")
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith(("#", ";", "%"))]
    if not lines:
        raise ValueError(f"{path}: no data lines")
    sample = lines[:50]
    delimiter = None
    counts = {}
    for d in (",", "\t", ";"):
        counts[d] = sum(1 for ln in sample if d in ln)
    best = max(counts, key=counts.get)
    if counts[best] >= max(2, len(sample) // 2):
        delimiter = best
    rows = [_split(ln.strip(), delimiter) for ln in lines]
    rows = [[c.strip() for c in r] for r in rows]
    first_numeric = next((i for i, r in enumerate(rows) if len(r) >= 2 and _is_num(r[0]) and _is_num(r[1])), None)
    if first_numeric is None:
        raise ValueError(f"{path}: no numeric rows found")
    n_cols = max(len(r) for r in rows[first_numeric:first_numeric + 50])
    header = None
    if first_numeric > 0:
        cand = rows[first_numeric - 1]
        if len(cand) >= 2 and not all(_is_num(c) for c in cand if c):
            header = (cand + [""] * n_cols)[:n_cols]
    body = rows[first_numeric:]
    n_numeric = sum(1 for r in body if len(r) >= 2 and sum(_is_num(c) for c in r[:n_cols]) >= 2)
    return TableInfo(delimiter, n_cols, header, body[:max_rows], n_numeric)


def read_spectrum(path, xcol: int = 0, ycol: int = 1, delimiter: str | None | bool = False
                  ) -> tuple[np.ndarray, np.ndarray]:
    """Read columns ``xcol`` and ``ycol`` (0-based) from a text file; x is returned ascending.

    ``delimiter``: the default (False) auto-detects it as :func:`sniff_table` does; None means
    whitespace; or give the character. Rows that are not numeric in both columns are skipped.
    """
    if delimiter is False:
        delimiter = sniff_table(path).delimiter
    xs, ys = [], []
    for line in Path(path).read_text(errors="replace").splitlines():
        parts = _split(line.strip(), delimiter) if line.strip() else []
        if len(parts) <= max(xcol, ycol):
            continue
        try:
            xs.append(float(parts[xcol]))
            ys.append(float(parts[ycol]))
        except ValueError:
            continue  # header / comment
    if len(xs) < 10:
        raise ValueError(f"{path}: found only {len(xs)} numeric (x, y) rows in columns {xcol + 1} and {ycol + 1}")
    x, y = np.array(xs), np.array(ys)
    order = np.argsort(x)
    return x[order], y[order]


# -- protocols -------------------------------------------------------------------------------
def _protocol_from_cfg(cfg: dict, library=None) -> tuple[Protocol, float | None]:
    peaks = []
    for d in cfg["peaks"]:
        bad = set(d) - _PEAK_FIELDS
        if bad:
            raise ValueError(f"unknown peak field(s) {sorted(bad)} in {d.get('name')!r}")
        d = dict(d)
        for k in _TUPLE_FIELDS:
            if d.get(k) is not None:
                d[k] = tuple(d[k])
            elif k in d:
                del d[k]
        if "area_bounds" in d:
            lo, hi = d["area_bounds"]
            d["area_bounds"] = (0.0 if lo is None else lo, math.inf if hi is None else hi)
        if library is not None and str(d.get("shape", "")).startswith("%&"):
            d["reference"] = library.get(d["shape"])
        peaks.append(PeakSpec(**d))
    bg = cfg.get("background", "shirley")
    if isinstance(bg, dict):
        from .casa import CasaBackground
        bg = CasaBackground(bg["casa"], int(bg.get("avg", 1)), tuple(bg.get("params", ())), float(bg.get("scale", 1.0)))
    cs = cfg.get("count_scale")
    return Protocol(peaks, bg), (float(cs) if cs is not None else None)


def load_protocol(path, library=None) -> tuple[Protocol, float | None]:
    """Single-protocol JSON. Returns (protocol, count_scale or None if not given).

    ``library`` (a :class:`~fitstability.library.ReferenceLibrary`) supplies the reference spectra of any
    CasaXPS library-spectrum shapes (``%&<>LS-...``) the protocol uses."""
    return _protocol_from_cfg(json.loads(Path(path).read_text()), library)


def protocol_to_cfg(protocol: Protocol, count_scale: float | None = None) -> dict:
    """The JSON-able form of a protocol (defaults are omitted)."""
    defaults = PeakSpec("x", 0.0, 1.0)
    peaks = []
    for p in protocol.peaks:
        d = {"name": p.name, "pos": p.pos, "fwhm": p.fwhm}
        for k in sorted(_PEAK_FIELDS - {"name", "pos", "fwhm"}):
            v = getattr(p, k)
            if v == getattr(defaults, k):
                continue
            if k == "area_bounds":
                v = [v[0], None if math.isinf(v[1]) else v[1]]
            elif isinstance(v, tuple):
                v = list(v)
            d[k] = v
        peaks.append(d)
    bg = protocol.background
    if callable(bg):
        bg = {"casa": bg.kind, "avg": bg.avg, "params": list(bg.params), "scale": bg.scale}
    cfg = {"background": bg, "peaks": peaks}
    if count_scale is not None:
        cfg["count_scale"] = count_scale
    return cfg


def save_protocol(path, protocol: Protocol, count_scale: float | None = None) -> None:
    Path(path).write_text(json.dumps(protocol_to_cfg(protocol, count_scale), indent=2) + "\n", encoding="utf-8")


class ProtocolSet:
    """One protocol for every block, or a protocol per region.

    JSON is either a single protocol (applied to every block) or
    ``{"protocols": {"C 1s": {...}, "Ti 2p": {...}}}`` keyed by region name (case and spaces
    ignored, matched against species+transition then block id).
    """

    def __init__(self, default=None, by_region=None):
        self.default = default
        self.by_region = by_region or {}

    @staticmethod
    def _key(s: str) -> str:
        return re.sub(r"\s+", "", s).lower()

    def for_block(self, block) -> tuple[Protocol, float | None]:
        for name in (block.region, block.block_id):
            hit = self.by_region.get(self._key(name))
            if hit:
                return hit
        if self.default is None:
            raise KeyError(f"no protocol for {block.label}; protocols given for: {', '.join(self.by_region)}")
        return self.default


def load_protocols(path, library=None) -> ProtocolSet:
    cfg = json.loads(Path(path).read_text())
    if "protocols" in cfg:
        return ProtocolSet(by_region={ProtocolSet._key(k): _protocol_from_cfg(v, library) for k, v in cfg["protocols"].items()})
    return ProtocolSet(default=_protocol_from_cfg(cfg, library))


def synthetic_spectrum(peaks: list[tuple[float, float, float]], x=None, bg=(200.0, 800.0),
                       counts_scale: float = 1.0, seed: int | None = 0, eta: float = 0.3):
    """Synthetic C 1s-like spectrum: pseudo-Voigt peaks ``(pos, fwhm, area)`` on a
    Shirley-like step, with Poisson noise. Returns (x, y) with x ascending."""
    x = np.linspace(281.0, 296.0, 301) if x is None else np.asarray(x, float)
    y = np.zeros_like(x)
    for pos, fwhm, area in peaks:
        y += pseudo_voigt(x, pos, fwhm, area, eta)
    cum = np.cumsum(y)
    step = bg[0] + (bg[1] - bg[0]) * cum / cum[-1]
    clean = y + step
    if seed is None:
        return x, clean
    rng = np.random.default_rng(seed)
    return x, rng.poisson(clean * counts_scale) / counts_scale
