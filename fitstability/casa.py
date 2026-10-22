"""CasaXPS fits stored in VAMAS block comments, turned into fitstability protocols.

When CasaXPS saves a fitted spectrum as VAMAS it writes the fit into each block's comment lines::

    Calib M = 455.59 A = 458.6 BE ADD
    CASA region (*Cu 2p*) (*Shirley*) 546.68 560.13 25.39 5 0.32 0 0 -450 0 0 (*Cu 2p*) 63.549 0 25.39
    CASA comp (*Ti(0)*) (*LA(1.1,1.9,7)*) Area 27251 1e-20 4.4e6 -1 1 MFWHM 0.93 0.22 5.5 -1 1 Position 1032.4 1021.6 1036.4 -1 1 RSF ...
    CASA comp (*Ti(0)*) (*LA(1.1,1.9,7)*) Area 13625 1e-20 4.4e6  0 0.5 MFWHM 1.22 1e-05 3.2 0 1.318 Position 1026.3 0 0 0 -6.106 RSF ...

Each of Area / MFWHM / Position is ``value  lower  upper  tie-index  tie-factor``:

* tie index ``-1``: free, between the bounds (lower == upper: held fixed);
* tie index ``k >= 0``: tied to the k-th component of the block (0-based, in file order, across
  regions): Area = factor x Area[k]; MFWHM = factor x MFWHM[k]; Position = Position[k] + factor
  (kinetic energy, so a binding-energy offset of the opposite sign).

This was checked on the 253 VAMAS files found on the author's machine: every tied value in every
fit equals its reference combined that way, and no other constraint form occurs.

Conventions (shared with the author's eXPoSe SpectraDeck, from which the line shapes are taken):

* region limits and positions are kinetic energies in the frame set by the ``Calib`` line;
* areas are counts per second x eV, so for a counts ordinate they are multiplied by
  dwell x scans (``area_units="counts"``, the default; pass ``"cps"`` if the ordinate is already
  counts per second);
* the fit is restricted to the region window, with the region's own background (Shirley, Linear,
  Tougaard variants) recomputed from the data in that window.

Binding energies used here are the calibrated ones CasaXPS displays (raw axis + calibration).
Line shapes GL/SGL are exact; LA and LF are reconstructions (see ``casa_lineshapes``).
Components of an unsupported shape (CasaXPS library-spectrum components) make that region
unavailable, reported in :attr:`CasaFit.skipped`.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

from . import casa_lineshapes as cl
from .model import FitProblem, PeakSpec, Protocol, StartPoint

_TAG = re.compile(r"\(\*(.*?)\*\)")
_NUM = re.compile(r"^-?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?$")
_CALIB = re.compile(r"^Calib\s+M\s*=\s*([-\d.eE+]+)\s+A\s*=\s*([-\d.eE+]+)(.*)$")
_REGION = re.compile(r"^CASA region \(\*(.*?)\*\) \(\*(.*?)\*\)\s*(.*)$")
_COMP = re.compile(r"^CASA comp \(\*(.*?)\*\) \(\*(.*?)\*\)\s*(.*)$")
_INDEX = re.compile(r"\bINDEX\s+(-?\d+)\s*(?:\(\*(.*?)\*\))?")
_PARAM = re.compile(r"\b(Area|MFWHM|Position)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)")


class CasaFitError(ValueError):
    """The block has no usable CasaXPS fit."""


@dataclass
class CasaParam:
    value: float
    lo: float
    hi: float
    tie: int            # -1 = free, else 0-based index of the reference component in the block
    factor: float       # ratio (Area, MFWHM) or offset (Position, KE) to the reference

    @property
    def tied(self) -> bool:
        return self.tie >= 0

    @property
    def fixed(self) -> bool:
        return not self.tied and self.lo == self.hi


@dataclass
class CasaRegion:
    name: str
    background: str
    start_ke: float
    end_ke: float
    avg: int = 1
    params: tuple = ()          # the six numbers after the averaging width (offsets, Tougaard B/C/D)


@dataclass
class CasaComp:
    index: int                  # position in the block's component list
    name: str
    shape: str
    area: CasaParam
    fwhm: CasaParam
    pos: CasaParam
    region: str = ""
    state: int = -1             # CasaXPS INDEX: chemical-state group


@dataclass
class CasaFit:
    regions: list = field(default_factory=list)
    comps: list = field(default_factory=list)
    calib_shift: float = 0.0
    region_shift: float | None = None
    comp_shift: float | None = None
    skipped: list = field(default_factory=list)    # (region name, reason)

    def shift_of_regions(self) -> float:
        return self.calib_shift if self.region_shift is None else self.region_shift

    def shift_of_comps(self) -> float:
        return self.calib_shift if self.comp_shift is None else self.comp_shift

    def comps_of(self, region: CasaRegion) -> list:
        mine = [c for c in self.comps if c.region == region.name]
        if mine or len(self.regions) != 1:
            return mine
        return list(self.comps)


def _floats(text: str) -> list[float]:
    return [float(t) for t in text.split() if _NUM.match(t)]


def _calibration(lines):
    total = reg = comp = 0.0
    found = False
    for raw in lines or ():
        m = _CALIB.match(str(raw).strip())
        if not m:
            continue
        found = True
        d = float(m.group(2)) - float(m.group(1))
        flags = {w.upper() for w in m.group(3).split()}
        total += d
        reg += 0.0 if "REGIONS" in flags else d
        comp += 0.0 if "COMPS" in flags else d
    return (total, reg, comp) if found else None


def parse_fit(comment_lines) -> CasaFit | None:
    """The fit held in a block's comment lines, or None if there is none."""
    fit = CasaFit()
    cal = _calibration(comment_lines)
    if cal:
        fit.calib_shift, fit.region_shift, fit.comp_shift = cal
    current = ""
    for raw in comment_lines or ():
        line = str(raw).strip()
        m = _REGION.match(line)
        if m:
            nums = _floats(_TAG.sub(" ", m.group(3)))
            if len(nums) >= 2:
                avg = int(nums[3]) if len(nums) > 3 and 0 < nums[3] < 50 else 1
                fit.regions.append(CasaRegion(m.group(1), m.group(2), min(nums[:2]), max(nums[:2]), avg,
                                              tuple(nums[4:10])))
                current = m.group(1)
            continue
        m = _COMP.match(line)
        if m:
            rest = m.group(3)
            p = {k: CasaParam(float(v), float(lo), float(hi), int(float(ci)), float(cf))
                 for k, v, lo, hi, ci, cf in _PARAM.findall(rest)}
            if not {"Area", "MFWHM", "Position"} <= set(p):
                continue
            idx = _INDEX.search(rest)
            fit.comps.append(CasaComp(len(fit.comps), m.group(1), m.group(2), p["Area"], p["MFWHM"], p["Position"],
                                      current, int(idx.group(1)) if idx else -1))
    if len(fit.regions) > 1:
        _attribute_regions(fit)
    return fit if fit.regions or fit.comps else None


def _attribute_regions(fit: CasaFit) -> None:
    """When region lines come before their components, give each component to the one region whose
    window contains its position (only when exactly one does)."""
    by_name = {r.name: r for r in fit.regions}
    for c in fit.comps:
        cur = by_name.get(c.region)
        if cur is not None and cur.start_ke <= c.pos.value <= cur.end_ke:
            continue
        hits = [r for r in fit.regions if r.start_ke <= c.pos.value <= r.end_ke]
        if len(hits) == 1:
            c.region = hits[0].name


@dataclass
class CasaBackground:
    """A region's CasaXPS background as a picklable ``f(x, y)`` on ascending binding energy.

    ``scale`` converts y (counts) to the counts-per-second the background functions expect.
    """
    kind: str
    avg: int
    params: tuple
    scale: float = 1.0

    def __call__(self, x, y):
        ke = -np.asarray(x, float)[::-1]            # ascending kinetic energy (up to a constant)
        cps = np.asarray(y, float)[::-1] / self.scale
        bg = cl.background(self.kind, cps, self.avg, x=ke, params=self.params)
        if bg is None:
            raise CasaFitError(f"background type {self.kind!r} cannot be reconstructed")
        return (np.asarray(bg, float) * self.scale)[::-1]


@dataclass
class CasaModel:
    """One CasaXPS region ready to analyse: data window, protocol, and Casa's own fitted values."""
    region: str
    protocol: Protocol
    x: np.ndarray                # calibrated binding energy, ascending, inside the region window
    y: np.ndarray
    stored: StartPoint           # the parameters CasaXPS saved (a Casa fit, in this protocol's frame)
    window: tuple[float, float]
    notes: list = field(default_factory=list)
    count_scale: float = 1.0
    rel_rms: float = 0.0         # rms of (data - stored fit) / data range: how well Casa's saved fit reproduces
    rstd: float = 0.0            # the same, Poisson-weighted (CasaXPS "RSTD")

    def problem(self) -> FitProblem:
        return FitProblem(self.x, self.y, self.protocol, self.count_scale)


def _unique_names(comps) -> list[str]:
    seen: dict = {}
    for c in comps:
        seen[c.name] = seen.get(c.name, 0) + 1
    count: dict = {}
    out = []
    for c in comps:
        if seen[c.name] == 1:
            out.append(c.name)
        else:
            count[c.name] = count.get(c.name, 0) + 1
            out.append(f"{c.name} ({count[c.name]})")
    return out


RELAX_AREA_BOUNDS_ABOVE = 0.15      # rel. rms of the stored fit above which area units are not trusted


def casa_models(block, area_units: str = "counts", regions: list[str] | None = None) -> list[CasaModel]:
    """One :class:`CasaModel` per CasaXPS region of a VAMAS block that holds a fit.

    ``area_units``: ``"counts"`` (default) if the ordinate is counts, so stored areas (counts/s x eV)
    are multiplied by dwell x scans; ``"cps"`` if it is already counts per second.

    Each model carries how well CasaXPS's saved fit reproduces the data (``rel_rms``, ``rstd``). If it
    does not (relative rms above 0.15, or the envelope more than 3x the data), the stored areas are
    not trusted as being in these units: a note is added and the area bounds are opened to (0, inf),
    since multi-start draws its own random areas and only area *ratios* (ties) are needed.
    Raises :class:`CasaFitError` if the block has no fit. Regions that cannot be built (unsupported
    line shape or background, tie to a component outside the region) are listed in ``fit.skipped``
    of :func:`parse_fit`'s result and omitted here; see :func:`describe_skips`.
    """
    fit = parse_fit(block.comments)
    if fit is None or not fit.comps:
        raise CasaFitError(f"{block.label}: no CasaXPS fit found in the block comments")
    if area_units not in ("counts", "cps"):
        raise ValueError("area_units must be 'counts' or 'cps'")
    k = float(block.dwell_s * block.scans) if area_units == "counts" and block.dwell_s and block.scans else 1.0
    hv, calib = block.hv_ev, fit.calib_shift
    x_all = np.asarray(block.be, float) + calib
    y_all = np.asarray(block.y, float)
    rshift, cshift = fit.shift_of_regions(), fit.shift_of_comps()
    off_comp = hv + calib - cshift              # BE (calibrated) = off_comp - stored KE position
    off_reg = hv + calib - rshift
    models, fit.skipped = [], []
    for reg in fit.regions:
        comps = fit.comps_of(reg)
        if regions and reg.name not in regions:
            continue
        if not comps:
            continue
        try:
            models.append(_build(fit, reg, comps, x_all, y_all, off_comp, off_reg, k,
                                 1.0 if area_units == "counts" else (block.dwell_s * block.scans if block.dwell_s and block.scans else 1.0)))
        except CasaFitError as e:
            fit.skipped.append((reg.name, str(e)))
    block._casa_skipped = list(fit.skipped)      # surfaced by describe_skips
    if not models and not fit.skipped:
        raise CasaFitError(f"{block.label}: the CasaXPS fit has no region with components")
    return models


def describe_skips(block) -> list[str]:
    return [f"region {r!r} skipped: {why}" for r, why in getattr(block, "_casa_skipped", [])]


def _build(fit, reg, comps, x_all, y_all, off_comp, off_reg, k, count_scale=1.0) -> CasaModel:
    for c in comps:
        bad = cl.check_shape(c.shape)
        if not bad:
            try:
                cl.component_curve(np.linspace(-5.0, 5.0, 101), c.shape, 0.0, 1.0, 1.0)
            except Exception as e:  # noqa: BLE001 - any failure means the shape is not reproducible
                bad = [f"cannot be evaluated ({e})"]
        if bad:
            raise CasaFitError(f"{c.name} has line shape {c.shape!r}: " + "; ".join(bad))
    ids = {c.index: i for i, c in enumerate(comps)}
    names = _unique_names(comps)
    peaks, notes, approx = [], [], []
    for c in comps:
        pos_be = off_comp - c.pos.value
        kw = {}
        for kind, par in (("area", c.area), ("fwhm", c.fwhm), ("pos", c.pos)):
            if par.tied:
                if par.tie not in ids:
                    notes.append(f"{c.name}.{kind} is tied to a component outside region {reg.name!r}; held fixed")
                    kw[f"fix_{kind}"] = True
                else:
                    ref = names[ids[par.tie]]
                    kw[f"{kind}_tie"] = (ref, -par.factor if kind == "pos" else par.factor)
        if not c.pos.tied:
            lo, hi = sorted((off_comp - c.pos.lo, off_comp - c.pos.hi))
            kw["pos_bounds"] = (lo, hi)
            if c.pos.fixed:
                kw["fix_pos"] = True
        if c.fwhm.fixed:
            kw["fix_fwhm"] = True
        kw["fwhm_bounds"] = (max(c.fwhm.lo, 1e-3), max(c.fwhm.hi, c.fwhm.lo, 1e-3 * 2))
        kw["area_bounds"] = (max(c.area.lo, 0.0) * k, c.area.hi * k)
        peaks.append(PeakSpec(name=names[len(peaks)], pos=pos_be, fwhm=c.fwhm.value, shape=c.shape, **kw))
        if not cl.is_exact(c.shape) and c.shape.split("(")[0] not in approx:
            approx.append(c.shape.split("(")[0])
    if approx:
        notes.append(f"{'/'.join(approx)} line shapes are reconstructions of CasaXPS's, not its own code (GL/SGL are exact)")
    bg = CasaBackground(reg.background, reg.avg, reg.params, k)
    lo, hi = off_reg - reg.end_ke, off_reg - reg.start_ke
    sel = (x_all >= lo - 1e-6) & (x_all <= hi + 1e-6)
    if sel.sum() < 10:
        raise CasaFitError(f"only {int(sel.sum())} data points inside the region window")
    x, y = x_all[sel], y_all[sel]
    proto = Protocol(peaks, bg)
    try:
        bg(x, y)
    except CasaFitError:
        raise
    stored = proto.resolve(StartPoint(np.array([p.pos for p in peaks]),
                                      np.array([c.fwhm.value for c in comps]),
                                      np.array([c.area.value * k for c in comps])))
    model = CasaModel(reg.name, proto, x, y, stored, (float(lo), float(hi)), notes, count_scale)
    prob = model.problem()
    env = prob.model(stored)
    rng = max(float(y.max() - y.min()), 1e-30)
    model.rel_rms = float(np.sqrt(np.mean((y - env) ** 2)) / rng)
    model.rstd = float(np.sqrt(np.mean(((y - env) / prob.sigma) ** 2)))
    if model.rel_rms > RELAX_AREA_BOUNDS_ABOVE or float(env.max()) > 3.0 * float(y.max()):
        model.notes.append(f"CasaXPS's saved fit does not reproduce the data here (relative rms {model.rel_rms:.2f}); "
                           "its areas may be in other units, so area bounds are opened and only tie ratios are used")
        model.protocol = Protocol([PeakSpec(**{**p.__dict__, "area_bounds": (0.0, float("inf"))}) for p in peaks], bg)
    return model
