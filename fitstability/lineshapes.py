"""Peak line shapes and backgrounds.

All peak shapes are area-normalised: ``peak(x, pos, fwhm, area, eta)`` integrates to ``area``.
``eta`` is the Lorentzian fraction of a pseudo-Voigt (0 = Gaussian, 1 = Lorentzian).

The background is computed once from the data (linear or Shirley) and held fixed during the
fits, so the fitted peak areas are the only quantities that depend linearly on the model.
This keeps the multi-start experiment about the peak parameters, which is what the box-plot
test probes.
"""
from __future__ import annotations

import numpy as np

_SQRT_LN2 = np.sqrt(np.log(2.0))
_GAUSS_NORM = 2.0 * _SQRT_LN2 / np.sqrt(np.pi)  # peak height of unit-area Gaussian, per unit 1/FWHM


def pseudo_voigt(x, pos, fwhm, area, eta=0.3):
    """Area-normalised pseudo-Voigt: ``area * (eta*L + (1-eta)*G)``."""
    x = np.asarray(x, dtype=float)
    d = x - pos
    g = _GAUSS_NORM / fwhm * np.exp(-4.0 * np.log(2.0) * (d / fwhm) ** 2)
    lor = (2.0 / (np.pi * fwhm)) / (1.0 + 4.0 * (d / fwhm) ** 2)
    return area * ((1.0 - eta) * g + eta * lor)


def peak_height_per_area(fwhm, eta=0.3):
    """Peak height of a unit-area pseudo-Voigt of the given FWHM."""
    return (1.0 - eta) * _GAUSS_NORM / fwhm + eta * 2.0 / (np.pi * fwhm)


def linear_background(x, y, n_avg=3):
    """Straight line between the mean of the first and last ``n_avg`` points."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    x0, x1 = x[:n_avg].mean(), x[-n_avg:].mean()
    y0, y1 = y[:n_avg].mean(), y[-n_avg:].mean()
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def shirley_background(x, y, n_avg=3, tol=1e-6, max_iter=100):
    """Iterative Shirley background (Proctor-Sherwood form), independent of the fit.

    Works for either energy ordering. The background steps from its low-BE level to its
    high-BE level in proportion to the net peak area integrated from the low-BE side.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    y_lo = ys[:n_avg].mean()
    y_hi = ys[-n_avg:].mean()
    bg = np.full_like(ys, y_lo)
    for _ in range(max_iter):
        net = ys - bg
        # cumulative net-signal integral from the low-BE end up to each point
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (net[1:] + net[:-1]) * np.diff(xs))])
        total = cum[-1]
        if total == 0:
            break
        # B(E) = y_lo + (y_hi - y_lo) * integral_{E_low}^{E} net / total
        new_bg = y_lo + (y_hi - y_lo) * cum / total
        if np.max(np.abs(new_bg - bg)) < tol * max(1.0, np.abs(ys).max()):
            bg = new_bg
            break
        bg = new_bg
    out = np.empty_like(bg)
    out[order] = bg
    return out


BACKGROUNDS = {
    "none": lambda x, y: np.zeros_like(np.asarray(y, float)),
    "linear": linear_background,
    "shirley": shirley_background,
}


def casa_curve(x, shape, pos, fwhm, area):
    """A CasaXPS line shape (``"GL(30)"``, ``"LA(1.53,243)"``, ...) on ascending binding energy ``x``.

    CasaXPS evaluates shapes on a kinetic-energy axis, where the asymmetry of LA/LF points the
    other way, so the grid is mirrored (KE = -BE up to a constant, which shapes do not see),
    evaluated, and mirrored back. ``x`` must be uniformly spaced.
    """
    from .casa_lineshapes import component_curve
    ke = -np.asarray(x, float)[::-1]
    return component_curve(ke, shape, -pos, fwhm, area)[::-1]


def make_background(kind, x, y):
    if callable(kind):
        return np.asarray(kind(np.asarray(x, float), np.asarray(y, float)), float)
    try:
        return BACKGROUNDS[kind](x, y)
    except KeyError:
        raise ValueError(f"unknown background {kind!r}; choose from {sorted(BACKGROUNDS)}") from None
