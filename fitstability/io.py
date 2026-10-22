"""Reading spectra and protocols, and building synthetic test spectra.

Spectra: two-column text/CSV (binding energy, intensity), any common delimiter, header and
comment lines skipped. VAMAS and vendor formats are out of scope for this standalone module;
export a two-column file (or call :func:`fitstability.FitProblem` directly) until it is
integrated with a reader.

Protocol files are JSON::

    {"background": "shirley",
     "count_scale": 1.0,
     "peaks": [{"name": "C1", "pos": 285.0, "fwhm": 1.25},
               {"name": "C2", "pos": 286.5, "fwhm": 1.25, "width_group": "main"},
               {"name": "Sat", "pos": 291.5, "fwhm": 3.5, "fix_pos": true}]}
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

from .lineshapes import pseudo_voigt
from .model import PeakSpec, Protocol

_PEAK_FIELDS = {"name", "pos", "fwhm", "eta", "fix_pos", "fix_fwhm", "width_group", "pos_bounds",
                "fwhm_bounds", "pos_halfwindow"}


def read_spectrum(path) -> tuple[np.ndarray, np.ndarray]:
    xs, ys = [], []
    for line in Path(path).read_text(errors="replace").splitlines():
        parts = [t for t in re.split(r"[,\s;\t]+", line.strip()) if t]
        if len(parts) < 2:
            continue
        try:
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
        except ValueError:
            continue  # header / comment
    if len(xs) < 10:
        raise ValueError(f"{path}: found only {len(xs)} numeric (x, y) rows")
    x, y = np.array(xs), np.array(ys)
    order = np.argsort(x)
    return x[order], y[order]


def _protocol_from_cfg(cfg: dict) -> tuple[Protocol, float | None]:
    peaks = []
    for d in cfg["peaks"]:
        bad = set(d) - _PEAK_FIELDS
        if bad:
            raise ValueError(f"unknown peak field(s) {sorted(bad)} in {d.get('name')!r}")
        d = dict(d)
        for k in ("pos_bounds", "fwhm_bounds"):
            if k in d:
                d[k] = tuple(d[k])
        peaks.append(PeakSpec(**d))
    cs = cfg.get("count_scale")
    return Protocol(peaks, cfg.get("background", "shirley")), (float(cs) if cs is not None else None)


def load_protocol(path) -> tuple[Protocol, float | None]:
    """Single-protocol JSON. Returns (protocol, count_scale or None if not given)."""
    return _protocol_from_cfg(json.loads(Path(path).read_text()))


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


def load_protocols(path) -> ProtocolSet:
    cfg = json.loads(Path(path).read_text())
    if "protocols" in cfg:
        return ProtocolSet(by_region={ProtocolSet._key(k): _protocol_from_cfg(v) for k, v in cfg["protocols"].items()})
    return ProtocolSet(default=_protocol_from_cfg(cfg))


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
