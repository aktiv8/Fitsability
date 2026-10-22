"""The three stability experiments.

* :func:`multistart` -- the paper's box-plot test: many random starting conditions, one
  protocol. Narrow final boxes => one minimum; wide boxes => multiple local minima.
* :func:`monte_carlo` -- refit noisy replicates of the best-fit curve: parameter uncertainty.
* :func:`uniqueness` -- fix one parameter on a grid, refit the rest: parameter correlation.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np

from .fit import fit, fit_nominal
from .model import FitProblem, FitResult, StartPoint


def _fit_one(args):
    problem, start = args
    return fit(problem, start)


def _run(problem, starts, n_jobs):
    tasks = [(problem, s) for s in starts]
    if n_jobs == 1 or len(tasks) < 2:
        return [_fit_one(t) for t in tasks]
    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
        return list(ex.map(_fit_one, tasks, chunksize=max(1, len(tasks) // (4 * n_jobs))))


# -- multi-start ---------------------------------------------------------------------------
@dataclass
class MultiStartResult:
    problem: FitProblem
    starts: list[StartPoint]
    results: list[FitResult]

    @property
    def names(self):
        return [p.name for p in self.problem.protocol.peaks]

    def start_areas(self) -> np.ndarray:
        return np.array([s.area for s in self.starts])

    def final(self, what: str = "area") -> np.ndarray:
        return np.array([getattr(r, what) for r in self.results])

    def chi2(self) -> np.ndarray:
        return np.array([r.chi2 for r in self.results])

    def best(self) -> FitResult:
        return self.results[int(np.argmin(self.chi2()))]


def random_start(problem: FitProblem, rng: np.random.Generator, randomize_pos: float = 0.0,
                 randomize_fwhm: float = 0.0, bounds_fraction: float = 0.0) -> StartPoint:
    """One random starting condition, built as in the paper.

    Areas: uniform(0, 1) x the largest area that peak could have at its nominal position
    and width before exceeding the data envelope. Positions/widths start at their nominal
    values unless ``randomize_pos`` (eV half-range) / ``randomize_fwhm`` (fractional
    half-range, e.g. 0.3 = +/-30 %) are set, or ``bounds_fraction`` f > 0: every free position and
    width is then drawn uniformly over the fraction f of its own allowed range on each side of
    its nominal value (f = 1: the whole range, e.g. a CasaXPS constraint window).
    """
    proto = problem.protocol
    nom = proto.nominal_start()
    pos, fwhm = nom.pos.copy(), nom.fwhm.copy()
    area = rng.uniform(0.0, 1.0, len(pos)) * problem.max_areas(pos, fwhm)
    if randomize_pos:
        pos = pos + rng.uniform(-randomize_pos, randomize_pos, len(pos))
    if randomize_fwhm:
        fwhm = fwhm * (1.0 + rng.uniform(-randomize_fwhm, randomize_fwhm, len(pos)))
    start = StartPoint(pos, fwhm, area)
    if bounds_fraction > 0:
        lo, hi = proto.bounds()
        th = proto.pack(start)
        nom = proto.pack(proto.nominal_start())
        area_cols = {j for j in proto.area_idx if j is not None}
        for j in range(proto.n_free):
            if j in area_cols:
                continue
            a = nom[j] - bounds_fraction * (nom[j] - lo[j])
            b = nom[j] + bounds_fraction * (hi[j] - nom[j])
            th[j] = rng.uniform(a, b) if np.isfinite(a) and np.isfinite(b) and b > a else nom[j]
        start = proto.unpack(th)
    return proto.resolve(start)   # tied/fixed parameters follow their references


def multistart(problem: FitProblem, n_starts: int = 32, seed: int | None = None,
               randomize_pos: float = 0.0, randomize_fwhm: float = 0.0, n_jobs: int = 1,
               bounds_fraction: float = 0.0) -> MultiStartResult:
    rng = np.random.default_rng(seed)
    starts = [random_start(problem, rng, randomize_pos, randomize_fwhm, bounds_fraction) for _ in range(n_starts)]
    return MultiStartResult(problem, starts, _run(problem, starts, n_jobs))


# -- Monte Carlo ---------------------------------------------------------------------------
@dataclass
class MonteCarloResult:
    problem: FitProblem
    best: FitResult
    results: list[FitResult]

    @property
    def names(self):
        return [p.name for p in self.problem.protocol.peaks]

    def final(self, what: str = "area") -> np.ndarray:
        return np.array([getattr(r, what) for r in self.results])


def monte_carlo(problem: FitProblem, n_trials: int = 100, seed: int | None = None,
                best: FitResult | None = None, n_jobs: int = 1) -> MonteCarloResult:
    """Poisson-noise replicates of the best-fit curve, each refitted from the best fit."""
    rng = np.random.default_rng(seed)
    best = best or fit_nominal(problem)
    curve = problem.model(best.as_start())
    cs = problem.count_scale
    tasks = []
    for _ in range(n_trials):
        y = rng.poisson(np.maximum(curve * cs, 0.0)) / cs
        tasks.append((problem.with_y(y), best.as_start()))
    if n_jobs == 1:
        results = [_fit_one(t) for t in tasks]
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as ex:
            results = list(ex.map(_fit_one, tasks, chunksize=max(1, n_trials // (4 * n_jobs))))
    return MonteCarloResult(problem, best, results)


# -- uniqueness ----------------------------------------------------------------------------
UNIQUENESS_METHOD_NOTE = (
    "Uniqueness plots (as in the paper): the data are first fitted with the stated protocol. One "
    "parameter of one peak is then fixed at a series of values and the fit repeated with all the other "
    "parameters free. Only parameters that are free in the protocol can be probed: a position that is "
    "fixed, or a FWHM that is fixed or tied to other peaks (the usual constraint), is not a fit parameter "
    "and is not scanned. U/V-shaped MSE curve: the parameter is significant. Flat curve: the other "
    "parameters compensate for it (correlation / over-parameterisation)."
)


def scannable_parameters(protocol, parameter: str | None = None) -> list[tuple[str, str]]:
    """(peak name, "pos"|"fwhm") pairs that are free in ``protocol`` and so can be scanned.

    A parameter is scannable only if it is an independent free fit parameter: not fixed, not
    tied to another peak (by offset or ratio, as for a spin-orbit doublet) and, for the FWHM,
    not shared through a ``width_group``.
    """
    out = []
    for p in protocol.peaks:
        for prm in ("pos", "fwhm"):
            if parameter in (None, prm) and unscannable_reason(protocol, p.name, prm) is None:
                out.append((p.name, prm))
    return out


def unscannable_reason(protocol, peak: str, parameter: str) -> str | None:
    """Why ``peak.parameter`` is not a free fit parameter, or None if it can be scanned."""
    p = next(q for q in protocol.peaks if q.name == peak)
    if parameter == "pos":
        if p.pos_tie:
            return f"{peak}.pos is tied to {p.pos_tie[0]} (offset {p.pos_tie[1]:+g} eV)"
        return f"{peak}.pos is fixed" if p.fix_pos else None
    if parameter != "fwhm":
        raise ValueError("parameter must be 'pos' or 'fwhm'")
    if p.fwhm_tie:
        return f"{peak}.fwhm is constrained (tied to {p.fwhm_tie[0]}, ratio {p.fwhm_tie[1]:g})"
    if p.fix_fwhm:
        return f"{peak}.fwhm is fixed"
    if p.width_group:
        return f"{peak}.fwhm is constrained (tied in width group {p.width_group!r})"
    return None


@dataclass
class UniquenessResult:
    peak: str
    parameter: str          # "pos" or "fwhm"
    values: np.ndarray
    mse: np.ndarray         # mean squared weighted residual per grid value
    best_value: float
    best_mse: float

    def classify(self, ratio: float = 1.5) -> str:
        """"unique" (U/V: MSE rises on both sides of an interior minimum), "correlated" (flat),
        or "grid too narrow" (a rising curve whose minimum is on the edge of the grid)."""
        lo = self.mse.min()
        rise = self.mse.max() / max(lo, 1e-30)
        k = int(np.argmin(self.mse))
        if rise < ratio:
            return "correlated"
        if k in (0, len(self.mse) - 1):
            return "grid too narrow"
        return "unique" if min(self.mse[0], self.mse[-1]) / max(lo, 1e-30) >= ratio else "grid too narrow"


def uniqueness(problem: FitProblem, peak: str, parameter: str = "pos", half_range: float = 0.8,
               n_points: int = 9, best: FitResult | None = None) -> UniquenessResult:
    """The paper's uniqueness plot for one free parameter (see :data:`UNIQUENESS_METHOD_NOTE`).

    Raises ``ValueError`` for a parameter that is not free in the protocol (fixed position,
    fixed or tied FWHM): such a parameter does not enter the fit and cannot be probed.
    """
    if parameter not in ("pos", "fwhm"):
        raise ValueError("parameter must be 'pos' or 'fwhm'")
    proto = problem.protocol
    names = [p.name for p in proto.peaks]
    if peak not in names:
        raise ValueError(f"no peak named {peak!r}; peaks: {', '.join(names)}")
    why = unscannable_reason(proto, peak, parameter)
    if why:
        raise ValueError(f"{why}, so it is not a free fit parameter and has no uniqueness plot")
    i = names.index(peak)
    best = best or fit_nominal(problem)
    centre = float(getattr(best, parameter)[i])
    lo = centre - half_range
    if parameter == "fwhm":
        lo = max(lo, 0.05)
    values = np.linspace(lo, centre + half_range, n_points)
    key = "fix_pos" if parameter == "pos" else "fix_fwhm"
    mse = []
    for v in values:
        prob2 = FitProblem(problem.x, problem.y, proto.with_peak(peak, **{key: True, parameter: float(v)}),
                           problem.count_scale)
        s = best.as_start()
        getattr(s, parameter)[i] = v
        mse.append(fit(prob2, s).chi2 / len(problem.x))
    mse = np.array(mse)
    k = int(np.argmin(mse))
    return UniquenessResult(peak, parameter, values, mse, float(values[k]), float(mse[k]))
