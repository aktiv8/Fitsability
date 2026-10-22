"""Numeric summaries and an advisory verdict for the box-plot test.

The paper judges box width by eye. Here the spread is made scale-free so the same
thresholds work for any spectrum: each peak's 5-95 percentile range of final area is divided
by the total fitted area of the best run. Thresholds are this package's own choices, not
the paper's, and are configurable. A verdict is a prompt to look, never a pass mark.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .experiments import MultiStartResult


@dataclass
class Verdict:
    label: str                   # "single minimum", "mostly single minimum", "multiple minima"
    max_spread: float            # worst per-peak (p95-p5)/total area
    frac_at_best: float          # share of runs whose chi2 is within `chi2_tol` of the best
    n_distinct_minima: int       # chi2 clusters at the given tolerance
    per_peak_spread: dict[str, float]
    best_rstd: float
    notes: list[str]

    def __str__(self):
        lines = [f"Verdict: {self.label}",
                 f"  worst peak spread (5-95 %, fraction of total area): {self.max_spread:.3f}",
                 f"  runs within tolerance of best chi2: {100 * self.frac_at_best:.0f} %",
                 f"  distinct chi2 minima found: {self.n_distinct_minima}",
                 f"  best-fit RSTD: {self.best_rstd:.3f}"]
        lines += [f"  note: {n}" for n in self.notes]
        return "\n".join(lines)


def assess(ms: MultiStartResult, narrow: float = 0.02, wide: float = 0.10, chi2_tol: float = 0.01) -> Verdict:
    area = ms.final("area")
    chi2 = ms.chi2()
    best = ms.best()
    total = max(float(best.area.sum()), 1e-30)
    spread = (np.percentile(area, 95, axis=0) - np.percentile(area, 5, axis=0)) / total
    per_peak = dict(zip(ms.names, map(float, spread)))
    worst = float(spread.max())
    at_best = chi2 <= chi2.min() * (1 + chi2_tol)
    frac = float(np.mean(at_best))
    # cluster chi2 values: sort, start a new cluster on a relative gap larger than tol
    order = np.sort(chi2)
    clusters = 1 + int(np.sum(np.diff(order) > chi2_tol * np.maximum(order[:-1], 1e-30)))
    if worst <= narrow and frac >= 0.9:
        label = "single minimum"
    elif worst >= wide or frac < 0.5:
        label = "multiple minima"
    else:
        label = "mostly single minimum"
    notes = []
    if not all(r.success for r in ms.results):
        notes.append(f"{sum(not r.success for r in ms.results)} of {len(ms.results)} fits did not report convergence")
    pegged = _pegged(ms)
    if pegged:
        notes.append("parameters at a bound in the best fit: " + ", ".join(pegged))
    notes.append("narrow boxes show one minimum, not a good fit: check the residuals (RSTD) as well")
    return Verdict(label, worst, frac, clusters, per_peak, best.rstd, notes)


def _pegged(ms: MultiStartResult, eps: float = 1e-3) -> list[str]:
    proto = ms.problem.protocol
    lo, hi = proto.bounds()
    th = ms.best().theta
    names = proto.free_names()
    out = []
    for j, v in enumerate(th):
        span = hi[j] - lo[j] if np.isfinite(hi[j] - lo[j]) else 1.0
        if (np.isfinite(lo[j]) and v - lo[j] <= eps * span) or (np.isfinite(hi[j]) and hi[j] - v <= eps * span):
            if not names[j].endswith(".area"):
                out.append(names[j])
    return out


def summary_table(ms: MultiStartResult) -> str:
    """Plain-text per-peak statistics of the final fitted areas."""
    area = ms.final("area")
    start = ms.start_areas()
    rows = [f"{'peak':<8}{'start med':>11}{'start IQR':>11}{'final med':>11}{'final IQR':>11}{'final min':>11}{'final max':>11}"]
    for k, n in enumerate(ms.names):
        q = lambda a, p: np.percentile(a[:, k], p)
        rows.append(f"{n:<8}{q(start, 50):>11.3g}{q(start, 75) - q(start, 25):>11.3g}"
                    f"{q(area, 50):>11.3g}{q(area, 75) - q(area, 25):>11.3g}"
                    f"{area[:, k].min():>11.3g}{area[:, k].max():>11.3g}")
    return "\n".join(rows)
