"""Text and CSV forms of the experiment results (no Qt)."""
from __future__ import annotations

import csv

import numpy as np

from ..experiments import UNIQUENESS_METHOD_NOTE, MonteCarloResult, MultiStartResult, UniquenessResult
from ..verdict import Verdict, summary_table


def multistart_csv(ms: MultiStartResult, path) -> None:
    """One row per random start: its starting areas and everything it converged to."""
    names = ms.names
    start, area, pos, fwhm = ms.start_areas(), ms.final("area"), ms.final("pos"), ms.final("fwhm")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        head = ["run", "chi2", "rstd", "converged"]
        for n in names:
            head += [f"{n} start area", f"{n} area", f"{n} position", f"{n} fwhm"]
        w.writerow(head)
        for i, r in enumerate(ms.results):
            row = [i + 1, f"{r.chi2:.6g}", f"{r.rstd:.6g}", int(r.success)]
            for k in range(len(names)):
                row += [f"{start[i, k]:.6g}", f"{area[i, k]:.6g}", f"{pos[i, k]:.6g}", f"{fwhm[i, k]:.6g}"]
            w.writerow(row)


def multistart_report(title: str, ms: MultiStartResult, verdict: Verdict, casa_rstd: float | None = None,
                      notes: list[str] | None = None) -> str:
    lines = [f"Multi-start box-plot test: {title}", ""]
    for n in notes or []:
        lines.append(f"note: {n}")
    if casa_rstd is not None:
        lines.append(f"CasaXPS saved fit: RSTD {casa_rstd:.3f}; best of {len(ms.results)} starts: RSTD {ms.best().rstd:.3f}")
    lines += [summary_table(ms), "", str(verdict)]
    return "\n".join(lines)


def monte_carlo_report(title: str, mc: MonteCarloResult) -> str:
    a = mc.final("area")
    lines = [f"Monte Carlo ({len(mc.results)} Poisson replicates): {title}", "",
             f"{'peak':<24}{'area':>12}{'sd':>12}{'sd %':>8}"]
    for k, n in enumerate(mc.names):
        m, sd = a[:, k].mean(), a[:, k].std(ddof=1) if len(a) > 1 else float("nan")
        lines.append(f"{n[:23]:<24}{m:>12.4g}{sd:>12.3g}{100 * sd / max(abs(m), 1e-30):>8.1f}")
    return "\n".join(lines)


def uniqueness_report(title: str, results: list[UniquenessResult], skipped: list[str]) -> str:
    lines = [f"Uniqueness: {title}", "", UNIQUENESS_METHOD_NOTE, ""]
    lines += [f"not scanned: {s}" for s in skipped]
    for u in results:
        ratio = u.mse.max() / max(u.mse.min(), 1e-30)
        lines.append(f"{u.peak}.{u.parameter}: MSE max/min = {ratio:.2f} -> {u.classify()}")
    if not results:
        lines.append("no free parameters of the requested kind in this protocol, nothing to scan")
    return "\n".join(lines)


def summary_row(title: str, ms: MultiStartResult, verdict: Verdict, casa_rstd: float | None) -> list:
    return [title, verdict.label, f"{verdict.max_spread:.3f}", f"{100 * verdict.frac_at_best:.0f}",
            verdict.n_distinct_minima, f"{ms.best().rstd:.3f}", "" if casa_rstd is None else f"{casa_rstd:.3f}"]


SUMMARY_HEADER = ["block", "verdict", "worst spread", "% runs at best", "distinct minima", "best RSTD", "Casa RSTD"]


def write_summary_csv(rows: list[list], path) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(SUMMARY_HEADER)
        w.writerows(rows)
