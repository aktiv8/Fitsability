"""Figures: paper-style box plots, fit overlays, Monte Carlo and uniqueness plots.

Every ``plot_*`` takes an optional ``fig``. Given a matplotlib ``Figure`` (e.g. a Qt canvas's) it
clears it and draws there without touching pyplot; without one it makes a pyplot figure, as the
command line does. Line colours follow ``matplotlib.rcParams`` so a dark theme just works.
"""
from __future__ import annotations

import matplotlib
import numpy as np

from .experiments import MonteCarloResult, MultiStartResult, UniquenessResult
from .model import FitProblem, StartPoint


def _figure(fig, figsize):
    if fig is None:
        import matplotlib.pyplot as plt
        return plt.figure(figsize=figsize, constrained_layout=True)
    fig.clear()
    try:
        fig.set_layout_engine("constrained")
    except Exception:  # noqa: BLE001 - older matplotlib
        pass
    return fig


def _fg() -> str:
    return matplotlib.rcParams["text.color"]


def box_stats(data: np.ndarray, label: str) -> dict:
    """Matplotlib ``bxp`` stats in the paper's convention: whiskers at 5/95 %, min/max flier."""
    return dict(label=label, med=np.median(data), q1=np.percentile(data, 25), q3=np.percentile(data, 75),
                whislo=np.percentile(data, 5), whishi=np.percentile(data, 95),
                fliers=[data.min(), data.max()], mean=data.mean())


def draw_boxes(ax, matrix: np.ndarray, names: list[str], ylabel: str = "Peak area"):
    colors = matplotlib.colormaps["tab10"](np.arange(len(names)) % 10)
    stats = [box_stats(matrix[:, k], n) for k, n in enumerate(names)]
    bp = ax.bxp(stats, showmeans=True, showfliers=True, patch_artist=False,
                meanprops=dict(marker="s", markerfacecolor="none", markersize=5),
                flierprops=dict(marker="*", markersize=7))
    for k, c in enumerate(colors):
        for key in ("boxes", "medians", "whiskers", "caps", "fliers", "means"):
            items = bp[key]
            per = len(items) // len(names)
            for art in items[k * per:(k + 1) * per]:
                art.set_color(c)
                if key in ("fliers", "means"):
                    art.set_markeredgecolor(c)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, 1.05 * float(np.max(matrix)))
    ax.margins(x=0.05)
    if max(len(n) for n in names) > 6:                  # long (e.g. CasaXPS) component names
        ax.tick_params(axis="x", labelrotation=40)
        for lab in ax.get_xticklabels():
            lab.set_horizontalalignment("right")


def plot_multistart(ms: MultiStartResult, title: str = "", path: str | None = None, verdict=None, fig=None):
    """Start vs final areas (paper's Fig. 3/4 layout) plus the best fit with residuals."""
    fig = _figure(fig, (11, 7))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 2])
    ax_s = fig.add_subplot(gs[0, 0])
    ax_f = fig.add_subplot(gs[0, 1])
    start, final = ms.start_areas(), ms.final("area")
    draw_boxes(ax_s, start, ms.names)
    draw_boxes(ax_f, final, ms.names)
    top = 1.05 * float(max(start.max(), final.max()))   # same scale on both panels
    ax_s.set_ylim(0, top)
    ax_f.set_ylim(0, top)
    ax_s.set_title("Random starting areas", fontsize=9)
    lab = f"Final areas ({len(ms.results)} fits)" + (f" - {verdict.label}" if verdict else "")
    ax_f.set_title(lab, fontsize=9)
    ax_r = fig.add_subplot(gs[1, :])
    plot_fit(ms.problem, ms.best().as_start(), ax=ax_r)
    if title:
        fig.suptitle(title, fontsize=9)
    if path:
        fig.savefig(path, dpi=150)
    return fig


def plot_fit(problem: FitProblem, s: StartPoint, ax=None, residual_offset: float | None = None, fig=None):
    """Data, components, background, envelope, and the (offset) weighted residual.

    Draws into ``ax``; else into a new axes of ``fig``; else into a new pyplot figure.
    """
    if ax is None:
        ax = _figure(fig, (7, 4)).add_subplot(111)
    x = problem.x
    comps = problem.components(s)
    ax.plot(x, problem.y, "-", color=_fg(), lw=1.2, label="data")
    ax.plot(x, problem.model(s), "r-", lw=1.0, label="fit")
    ax.plot(x, problem.background, "-", color="grey", lw=0.8, label="background")
    for k, p in enumerate(problem.protocol.peaks):
        ax.plot(x, problem.background + comps[k], lw=0.9, label=p.name)
    resid = (problem.y - problem.model(s))
    off = residual_offset if residual_offset is not None else 1.08 * problem.y.max()
    ax.plot(x, off + resid, color="olive", lw=0.8, label="residual")
    ax.set_xlabel("Binding energy (eV)")
    ax.set_ylabel("Intensity")
    if x[0] < x[-1]:
        ax.invert_xaxis()
    ax.legend(ncol=2, fontsize=7, frameon=False)
    return ax


def plot_data(x, y, title: str = "", fig=None, background=None, markers=None):
    """The spectrum (and optional background) with dotted markers at ``markers`` = [(energy, label)],
    for a data set with no fitted model yet."""
    ax = _figure(fig, (7, 4)).add_subplot(111)
    ax.plot(x, y, "-", color=_fg(), lw=1.2, label="data")
    if background is not None:
        ax.plot(x, background, "-", color="grey", lw=0.8, label="background")
        ax.legend(frameon=False, fontsize=8)
    top = float(np.max(y))
    for e, lab in markers or []:
        ax.axvline(e, color="tab:blue", ls=":", lw=0.9)
        ax.text(e, top, lab, rotation=90, va="top", ha="right", fontsize=7, color="tab:blue")
    ax.set_xlabel("Binding energy (eV)")
    ax.set_ylabel("Intensity")
    if x[0] < x[-1]:
        ax.invert_xaxis()
    if title:
        ax.set_title(title, fontsize=10)
    return ax.figure


def plot_monte_carlo(mc: MonteCarloResult, title: str = "", path: str | None = None, fig=None):
    fig = _figure(fig, (11, 4))
    a, b = fig.subplots(1, 2)
    draw_boxes(a, mc.final("area"), mc.names)
    a.set_title(f"Monte Carlo areas ({len(mc.results)} noisy replicates)")
    pos = mc.final("pos")
    draw_boxes(b, pos, mc.names, ylabel="Peak position (eV)")
    b.set_ylim(pos.min() - 0.05, pos.max() + 0.05)
    b.set_title("Monte Carlo positions")
    if title:
        fig.suptitle(title)
    if path:
        fig.savefig(path, dpi=150)
    return fig


def plot_uniqueness(results: list[UniquenessResult], title: str = "", path: str | None = None, fig=None):
    """Paper-style uniqueness plots: MSE against the fixed value of one free parameter, one panel per
    probed parameter (positions and widths have different ranges, so they are never overlaid).
    U/V = significant, flat = correlated with other parameters."""
    n = max(len(results), 1)
    cols = min(3, n)
    rows = -(-n // cols)
    fig = _figure(fig, (4.2 * cols, 3.3 * rows + 0.4))
    axes = np.atleast_1d(fig.subplots(rows, cols, squeeze=False)).ravel()
    for ax, u in zip(axes, results):
        ax.plot(u.values, u.mse, "o-", color="tab:blue")
        ax.axvline(u.best_value, color="grey", ls=":", lw=0.8)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(f"fixed {'position' if u.parameter == 'pos' else 'FWHM'} (eV)", fontsize=8)
        ax.set_ylabel("MSE (weighted)", fontsize=8)
        ax.set_title(f"{u.peak}: {u.classify()}", fontsize=8)
        ax.tick_params(labelsize=7)
    for ax in axes[len(results):]:
        ax.set_visible(False)
    fig.suptitle(title or "Uniqueness (free parameters only)", fontsize=9)
    if path:
        fig.savefig(path, dpi=150)
    return fig
