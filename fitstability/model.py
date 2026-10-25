"""Peak model definition and the least-squares fit problem.

A *protocol* is a list of :class:`PeakSpec` plus a background: everything an analyst
chooses before pressing "fit". The experiments in this package ask how sensitive the fit
outcome of one protocol is to its starting conditions, to noise, and to individual parameters.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from .lineshapes import casa_curve, make_background, peak_height_per_area, pseudo_voigt


@dataclass(frozen=True)
class PeakSpec:
    """One synthetic peak: nominal start values and the constraints on it.

    ``shape`` is empty for the built-in pseudo-Voigt (``eta`` = Lorentzian fraction) or a
    CasaXPS line-shape string such as ``"GL(30)"``, ``"LA(1.53,243)"`` or ``"LF(1,1,75,200)"``.
    A shape starting ``%&`` is a CasaXPS library spectrum: ``reference`` must hold its
    :class:`~fitstability.library.ReferenceSpectrum` (the FWHM is then unused).

    Constraints (a constrained parameter is not a free fit parameter):

    * ``fix_pos`` / ``fix_fwhm``: hold the value at its nominal start.
    * ``width_group``: every peak with the same non-empty label shares one FWHM parameter.
    * ``area_tie=(ref, f)``: area = f x area of peak ``ref``.
    * ``fwhm_tie=(ref, f)``: FWHM = f x FWHM of ``ref``.
    * ``pos_tie=(ref, d)``: position = position of ``ref`` + d (eV, in the axis direction).

    These are CasaXPS's constraints (a component tied to another by ratio or offset, e.g. a
    spin-orbit doublet). Bounds default to a window around the nominal values.
    """

    name: str
    pos: float
    fwhm: float
    eta: float = 0.3
    shape: str = ""
    fix_pos: bool = False
    fix_fwhm: bool = False
    width_group: str = ""
    area_tie: tuple[str, float] | None = None
    fwhm_tie: tuple[str, float] | None = None
    pos_tie: tuple[str, float] | None = None
    pos_bounds: tuple[float, float] | None = None
    fwhm_bounds: tuple[float, float] = (0.01, 6.0)
    area_bounds: tuple[float, float] = (0.0, float("inf"))
    pos_halfwindow: float = 3.0  # eV, used when pos_bounds is None
    # the reference spectrum of a CasaXPS library component (shape '%&<>LS-CE4'); see fitstability.library
    reference: object = field(default=None, compare=False, repr=False)


@dataclass
class StartPoint:
    """Starting values for every peak (full, un-tied): arrays aligned with the peak list."""

    pos: np.ndarray
    fwhm: np.ndarray
    area: np.ndarray


@dataclass
class FitResult:
    pos: np.ndarray
    fwhm: np.ndarray
    area: np.ndarray
    chi2: float          # sum of squared weighted residuals
    rstd: float          # sqrt(mean(weighted residual^2)), CasaXPS-style "RSTD"
    red_chi2: float
    success: bool
    n_eval: int
    theta: np.ndarray = field(repr=False, default_factory=lambda: np.empty(0))

    def as_start(self) -> StartPoint:
        return StartPoint(self.pos.copy(), self.fwhm.copy(), self.area.copy())


class Protocol:
    """Peaks + background; maps between full peak parameters and the free-parameter vector.

    ``background`` is ``"none"``, ``"linear"``, ``"shirley"`` or a callable ``f(x, y)`` returning
    the background on the (ascending) fit grid, so a fixed, picklable object can carry a
    CasaXPS background.
    """

    def __init__(self, peaks: list[PeakSpec], background="shirley"):
        if not peaks:
            raise ValueError("a protocol needs at least one peak")
        names = [p.name for p in peaks]
        if len(set(names)) != len(names):
            raise ValueError("peak names must be unique")
        self.peaks = list(peaks)
        self.background = background
        self._index = {n: i for i, n in enumerate(names)}
        for p in self.peaks:
            if p.shape.startswith("%&") and p.reference is None:
                raise ValueError(f"{p.name}: the reference spectrum {p.shape.lstrip('%&<>')!r} is not loaded "
                                 "(supply it with --library, or Reference spectra... in the GUI)")
            for kind in ("area_tie", "fwhm_tie", "pos_tie"):
                t = getattr(p, kind)
                if t is not None and t[0] not in self._index:
                    raise ValueError(f"{p.name}.{kind} refers to unknown peak {t[0]!r}")
                if t is not None and t[0] == p.name:
                    raise ValueError(f"{p.name}.{kind} refers to itself")
        for kind in ("area_tie", "fwhm_tie", "pos_tie"):  # no cycles
            for i in range(len(self.peaks)):
                seen, j = {i}, i
                while (t := getattr(self.peaks[j], kind)) is not None:
                    j = self._index[t[0]]
                    if j in seen:
                        raise ValueError(f"circular {kind} involving {self.peaks[i].name!r}")
                    seen.add(j)
        self._build_layout()

    # -- layout of the free-parameter vector -------------------------------------------
    def _build_layout(self):
        n = len(self.peaks)
        k = 0
        self.area_idx: list[int | None] = [None] * n
        self.pos_idx: list[int | None] = [None] * n
        self.fwhm_idx: list[int | None] = [None] * n
        for i, p in enumerate(self.peaks):
            if p.area_tie is None:
                self.area_idx[i] = k
                k += 1
        for i, p in enumerate(self.peaks):
            if not p.fix_pos and p.pos_tie is None:
                self.pos_idx[i] = k
                k += 1
        group_idx: dict[str, int] = {}
        for i, p in enumerate(self.peaks):
            if p.fix_fwhm or p.fwhm_tie is not None:
                continue
            if p.width_group:
                if p.width_group not in group_idx:
                    group_idx[p.width_group] = k
                    k += 1
                self.fwhm_idx[i] = group_idx[p.width_group]
            else:
                self.fwhm_idx[i] = k
                k += 1
        self.n_free = k
        self.width_groups = group_idx

    def is_free(self, peak: str, parameter: str) -> bool:
        """Is ``peak.parameter`` ("area", "pos" or "fwhm") an independent free fit parameter?
        A tied, fixed or (for fwhm) group-shared parameter is not."""
        i = self._index[peak]
        if parameter == "fwhm" and self.peaks[i].width_group:
            return False
        return {"area": self.area_idx, "pos": self.pos_idx, "fwhm": self.fwhm_idx}[parameter][i] is not None

    def free_names(self) -> list[str]:
        out = [""] * self.n_free
        for i, p in enumerate(self.peaks):
            if self.area_idx[i] is not None:
                out[self.area_idx[i]] = f"{p.name}.area"
            if self.pos_idx[i] is not None:
                out[self.pos_idx[i]] = f"{p.name}.pos"
        for g, j in self.width_groups.items():
            out[j] = f"width[{g}]"
        for i, p in enumerate(self.peaks):
            j = self.fwhm_idx[i]
            if j is not None and not out[j]:
                out[j] = f"{p.name}.fwhm"
        return out

    def pack(self, s: StartPoint) -> np.ndarray:
        th = np.zeros(self.n_free)
        for i in range(len(self.peaks)):
            if self.area_idx[i] is not None:
                th[self.area_idx[i]] = s.area[i]
            if self.pos_idx[i] is not None:
                th[self.pos_idx[i]] = s.pos[i]
            if self.fwhm_idx[i] is not None:
                th[self.fwhm_idx[i]] = s.fwhm[i]  # peaks sharing a width: the last start wins
        return th

    def unpack(self, th: np.ndarray) -> StartPoint:
        n = len(self.peaks)
        pos = [th[self.pos_idx[i]] if self.pos_idx[i] is not None else self.peaks[i].pos for i in range(n)]
        fwhm = [th[self.fwhm_idx[i]] if self.fwhm_idx[i] is not None else self.peaks[i].fwhm for i in range(n)]
        area = [th[self.area_idx[i]] if self.area_idx[i] is not None else 0.0 for i in range(n)]
        vals = {"area": area, "pos": pos, "fwhm": fwhm}
        done: set = set()

        def resolve(kind, i):
            t = getattr(self.peaks[i], kind + "_tie")
            if t is None or (kind, i) in done:
                return vals[kind][i]
            ref = resolve(kind, self._index[t[0]])
            vals[kind][i] = ref + t[1] if kind == "pos" else ref * t[1]
            done.add((kind, i))
            return vals[kind][i]

        for kind in vals:
            for i in range(n):
                resolve(kind, i)
        return StartPoint(np.array(pos, float), np.array(fwhm, float), np.array(area, float))

    def resolve(self, s: StartPoint) -> StartPoint:
        """``s`` with tied parameters recomputed from their references (and fixed ones reset)."""
        return self.unpack(self.pack(s))

    def bounds(self) -> tuple[np.ndarray, np.ndarray]:
        lo = np.full(self.n_free, -np.inf)
        hi = np.full(self.n_free, np.inf)
        for i, p in enumerate(self.peaks):
            if self.area_idx[i] is not None:
                lo[self.area_idx[i]], hi[self.area_idx[i]] = p.area_bounds
            if self.pos_idx[i] is not None:
                lo[self.pos_idx[i]], hi[self.pos_idx[i]] = p.pos_bounds or (p.pos - p.pos_halfwindow, p.pos + p.pos_halfwindow)
            j = self.fwhm_idx[i]
            if j is not None:  # shared widths: the tightest bounds of the group
                lo[j] = p.fwhm_bounds[0] if not np.isfinite(lo[j]) else max(lo[j], p.fwhm_bounds[0])
                hi[j] = p.fwhm_bounds[1] if not np.isfinite(hi[j]) else min(hi[j], p.fwhm_bounds[1])
        return lo, hi

    def nominal_start(self) -> StartPoint:
        n = len(self.peaks)
        return self.resolve(StartPoint(np.array([p.pos for p in self.peaks]),
                                       np.array([p.fwhm for p in self.peaks]), np.zeros(n)))

    def with_peak(self, name: str, **changes) -> "Protocol":
        return Protocol([replace(p, **changes) if p.name == name else p for p in self.peaks], self.background)


class FitProblem:
    """A protocol bound to a spectrum. ``count_scale`` converts the y units to counts
    (e.g. counts = CPS x dwell x scans) so the Poisson weighting is right; use 1.0 if y is
    already in counts."""

    def __init__(self, x, y, protocol: Protocol, count_scale: float = 1.0):
        self.x = np.asarray(x, float)
        self.y = np.asarray(y, float)
        if self.x.shape != self.y.shape or self.x.ndim != 1:
            raise ValueError("x and y must be 1-D arrays of equal length")
        self.protocol = protocol
        self.count_scale = float(count_scale)
        self.background = make_background(protocol.background, self.x, self.y)
        self.sigma = np.sqrt(np.maximum(self.y * self.count_scale, 1.0)) / self.count_scale

    def with_y(self, y) -> "FitProblem":
        """Same protocol and x on new intensities (used by Monte Carlo)."""
        return FitProblem(self.x, y, self.protocol, self.count_scale)

    def components(self, s: StartPoint) -> np.ndarray:
        return np.array([peak_curve(self.x, p, s.pos[i], s.fwhm[i], s.area[i])
                         for i, p in enumerate(self.protocol.peaks)])

    def model(self, s: StartPoint) -> np.ndarray:
        return self.background + self.components(s).sum(axis=0)

    def residuals(self, theta: np.ndarray) -> np.ndarray:
        s = self.protocol.unpack(theta)
        return (self.y - self.model(s)) / self.sigma

    def max_areas(self, pos: np.ndarray, fwhm: np.ndarray) -> np.ndarray:
        """Largest area each peak can have, alone, before it pokes above the data envelope.

        Following the paper: raise each peak at its nominal position until its maximum just
        touches the data envelope; the height is taken relative to the background beneath.
        """
        net = self.y - self.background
        out = np.empty(len(pos))
        for i, p in enumerate(self.protocol.peaks):
            env = max(float(np.interp(pos[i], self.x, net)), 0.0)
            out[i] = env / unit_height(p, pos[i], fwhm[i])
        return out


def peak_curve(x, spec: PeakSpec, pos, fwhm, area):
    """The peak on ``x`` (ascending): built-in pseudo-Voigt, a CasaXPS line shape, or a library spectrum."""
    if spec.reference is not None:
        return spec.reference.curve(x, pos, area)
    if spec.shape:
        return casa_curve(x, spec.shape, pos, fwhm, area)
    return pseudo_voigt(x, pos, fwhm, area, spec.eta)


def unit_height(spec: PeakSpec, pos: float, fwhm: float) -> float:
    """Maximum of the peak when its area is 1."""
    if spec.reference is not None:
        return float(spec.reference.unit.max())
    if not spec.shape:
        return float(peak_height_per_area(fwhm, spec.eta))
    grid = np.linspace(pos - 12 * fwhm, pos + 12 * fwhm, 1201)  # uniform, as the Casa shapes need
    return float(casa_curve(grid, spec.shape, pos, fwhm, 1.0).max())
