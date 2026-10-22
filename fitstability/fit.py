"""Weighted least-squares fit of a :class:`FitProblem` from one starting point."""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .model import FitProblem, FitResult, StartPoint


def fit(problem: FitProblem, start: StartPoint, max_refits: int = 5, rtol: float = 1e-8) -> FitResult:
    """Minimise chi-squared (Poisson-weighted) from ``start``.

    The optimiser is restarted from its own result until chi-squared stops improving, as an
    analyst would by pressing "fit" repeatedly; a single call often stops early on a flat
    valley floor.
    """
    proto = problem.protocol
    lo, hi = proto.bounds()
    theta = np.clip(proto.pack(start), lo, hi)
    n_eval = 0
    chi2_prev = np.inf
    res = None
    for _ in range(max_refits):
        res = least_squares(problem.residuals, theta, bounds=(lo, hi), method="trf", x_scale="jac")
        n_eval += res.nfev
        theta = res.x
        chi2 = float(np.sum(res.fun ** 2))
        if chi2_prev - chi2 <= rtol * max(chi2, 1e-30):
            break
        chi2_prev = chi2
    s = proto.unpack(theta)
    r = res.fun
    dof = max(len(r) - proto.n_free, 1)
    chi2 = float(np.sum(r ** 2))
    return FitResult(
        pos=s.pos, fwhm=s.fwhm, area=s.area, chi2=chi2,
        rstd=float(np.sqrt(np.mean(r ** 2))), red_chi2=chi2 / dof,
        success=bool(res.success), n_eval=n_eval, theta=theta,
    )


def fit_nominal(problem: FitProblem, area_fraction: float = 0.5, **kw) -> FitResult:
    """Fit from the nominal positions/widths with areas at a fraction of their maximum."""
    s = problem.protocol.nominal_start()
    s.area = area_fraction * problem.max_areas(s.pos, s.fwhm)
    return fit(problem, s, **kw)
