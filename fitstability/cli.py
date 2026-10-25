"""Command line: ``fitstability list|boxplot|montecarlo|uniqueness|demo``.

Spectra come from VAMAS (.vms/.vamas; one or many blocks) or two-column text/CSV.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib
import numpy as np

from . import (FitProblem, UNIQUENESS_METHOD_NOTE, assess, fit_nominal, load_protocols, monte_carlo, multistart,
               read_spectrum, read_vamas, scannable_parameters, summary_table, synthetic_spectrum, uniqueness,
               unscannable_reason)
from .model import PeakSpec, Protocol
from .casa import CasaFitError, casa_models, describe_skips
from .library import ReferenceLibrary, calibrate_origins
from .vamas import VamasError

VAMAS_SUFFIXES = {".vms", ".vamas", ".vam"}


def _items(args):
    """Yield (tag, title, FitProblem, info) for every selected block (or the one text spectrum).

    With a protocol file the same protocol(s) are applied; without one the CasaXPS fit stored in
    each VAMAS block is used (every region of it)."""
    src = Path(args.source)
    lib = _library(args)
    protos = load_protocols(args.protocol, lib) if args.protocol else None
    if src.suffix.lower() in VAMAS_SUFFIXES:
        vf = read_vamas(src)
        blocks = vf.select(args.block)
        if len(blocks) > 1 and args.block is None:
            print(f"{src.name}: {len(blocks)} blocks; running all (use --block to choose, `fitstability list` to see them)")
        for b in blocks:
            tag0 = re.sub(r"[^A-Za-z0-9]+", "_", f"{b.index}_{b.region}").strip("_")
            if protos is None:
                try:
                    models = casa_models(b, args.count_units, library=lib)
                except CasaFitError as e:
                    print(f"skipping {b.label}: {e}")
                    continue
                for msg in describe_skips(b):
                    print(f"{b.label}: {msg}")
                for m in models:
                    tag = tag0 + (("_" + re.sub(r"[^A-Za-z0-9]+", "_", m.region).strip("_")) if len(models) > 1 else "")
                    yield tag, f"{src.name} {b.label}" + (f" / {m.region}" if len(models) > 1 else ""), m.problem(), m
                continue
            try:
                proto, cs = protos.for_block(b)
            except KeyError as e:
                print(f"skipping {b.label}: {e.args[0]}")
                continue
            scale = args.count_scale or cs or b.count_scale(args.count_units)
            yield tag0, f"{src.name} {b.label}", FitProblem(b.be, b.y, proto, scale), None
    else:
        if protos is None:
            sys.exit("a text/CSV spectrum needs a protocol file (the Casa fit is read from VAMAS files)")
        x, y = read_spectrum(src)
        proto, cs = protos.for_block(_Anon())
        yield "", src.name, FitProblem(x, y, proto, args.count_scale or cs or 1.0), None


def _library(args):
    paths = getattr(args, "library", None)
    if not paths:
        return None
    lib = ReferenceLibrary.load(paths, args.library_axis, args.library_hv)
    if not lib:
        raise ValueError("no reference spectra found in " + ", ".join(paths))
    return lib


class _Anon:
    region = block_id = ""
    label = "spectrum"


def _out_path(args, tag: str, kind: str) -> str:
    stem = Path(args.source).with_suffix("")
    if args.out and Path(args.out).suffix:        # explicit file
        p = Path(args.out)
        if tag:
            p = p.with_name(f"{p.stem}.{tag}{p.suffix}")
    else:
        base = Path(args.out) if args.out else stem.parent
        name = f"{stem.name}{'.' + tag if tag else ''}.{kind}.png"
        p = base / name
    p.parent.mkdir(parents=True, exist_ok=True)
    return str(p)


def cmd_library_list(args):
    lib = _library(args)
    if lib is None:
        sys.exit("give --library PATH (a file or folder of reference spectra)")
    print(f"{len(lib)} reference spectra:")
    for n in lib.names():
        r = lib.get(n)
        o = f"origin {r.origin:.3f} eV" if r.calibrated else f"origin = maximum ({r.max_energy:.2f} eV)"
        print(f"  {n:<24}{r.be.min():.1f}-{r.be.max():.1f} eV  {o}  [{r.source}]")


def cmd_library_calibrate(args):
    import json
    lib = _library(args)
    if lib is None:
        sys.exit("give --library PATH (a file or folder of reference spectra)")
    vf = read_vamas(args.source)
    blocks = vf.select(args.block)
    fits = calibrate_origins(lib, blocks, args.count_units)
    if not fits:
        sys.exit(f"no region of {Path(args.source).name} uses the loaded reference spectra: {', '.join(lib.names())}")
    print(f"Origins fitted to CasaXPS's saved fits in {Path(args.source).name} ({len(blocks)} blocks searched)\n")
    print(f"{'reference':<22}{'regions':>8}{'origin (eV)':>14}{'spread':>9}{'rel rms before':>16}{'after':>9}")
    for n, f in sorted(fits.items()):
        flag = "   <- spread is large: not a fixed offset" if f.spread > 0.15 and len(f.values) > 2 else ""
        print(f"{n:<22}{len(f.values):>8}{f.median:>14.3f}{f.spread:>9.3f}{f.rel_before:>16.3f}{f.rel_after:>9.3f}{flag}")
    if args.write:
        dest = Path(args.write)
        dest = dest / "origins.json" if dest.is_dir() else dest
        existing = json.loads(dest.read_text(encoding="utf-8")) if dest.exists() else {}
        existing.update({n: round(f.median, 4) for n, f in fits.items()})
        dest.write_text(json.dumps(existing, indent=2) + "\n", encoding="utf-8")
        print(f"\nWrote {dest}")
    else:
        print("\nRe-run with --write FOLDER to save these origins beside the reference spectra.")


def cmd_list(args):
    src = Path(args.source)
    if src.suffix.lower() not in VAMAS_SUFFIXES:
        sys.exit("list is for VAMAS files")
    vf = read_vamas(src)
    print(f"{src.name}: {vf.mode}, {len(vf.blocks)} blocks")
    for b in vf.blocks:
        print("  " + b.describe())


def cmd_boxplot(args):
    from .plots import plot_multistart
    for tag, title, prob, info in _items(args):
        print(f"\n=== {title} ===")
        for note in (info.notes if info else []):
            print(f"note: {note}")
        ms = multistart(prob, args.n, args.seed, args.randomize_pos, args.randomize_fwhm, args.jobs,
                        args.bounds_fraction)
        v = assess(ms)
        if info:
            print(f"CasaXPS saved fit: RSTD {info.rstd:.3f}; best of {args.n} starts: RSTD {ms.best().rstd:.3f}")
        print(summary_table(ms))
        print()
        print(v)
        path = _out_path(args, tag, "boxplot")
        plot_multistart(ms, title=title, path=path, verdict=v)
        print(f"Figure: {path}")


def cmd_montecarlo(args):
    from .plots import plot_monte_carlo
    for tag, title, prob, info in _items(args):
        print(f"\n=== {title} ===")
        mc = monte_carlo(prob, args.n, args.seed, n_jobs=args.jobs)
        a = mc.final("area")
        print(f"{'peak':<8}{'area':>12}{'sd':>12}{'sd %':>8}")
        for k, n in enumerate(mc.names):
            m, sd = a[:, k].mean(), a[:, k].std(ddof=1)
            print(f"{n:<8}{m:>12.4g}{sd:>12.3g}{100 * sd / max(abs(m), 1e-30):>8.1f}")
        path = _out_path(args, tag, "montecarlo")
        plot_monte_carlo(mc, title=title, path=path)
        print(f"Figure: {path}")


def cmd_uniqueness(args):
    from .plots import plot_uniqueness
    print(UNIQUENESS_METHOD_NOTE)
    wanted = {"pos": "pos", "fwhm": "fwhm", "both": None}[args.parameter]
    for tag, title, prob, info in _items(args):
        print(f"\n=== {title} ===")
        proto = prob.protocol
        names = args.peak or [p.name for p in proto.peaks]
        params = ("pos", "fwhm") if wanted is None else (wanted,)
        targets = []
        for n in names:
            for prm in params:
                why = unscannable_reason(proto, n, prm)
                if why:
                    print(f"  not scanned: {why}")
                else:
                    targets.append((n, prm))
        if not targets:
            print("  no free parameters of the requested kind in this protocol, nothing to scan")
            continue
        best = fit_nominal(prob)
        res = [uniqueness(prob, n, prm, args.half_range, args.points, best) for n, prm in targets]
        for u in res:
            print(f"  {u.peak}.{u.parameter}: MSE max/min = {u.mse.max() / max(u.mse.min(), 1e-30):.2f} -> {u.classify()}")
        path = _out_path(args, tag, "uniqueness")
        plot_uniqueness(res, title=title + " (free parameters only)", path=path)
        print(f"Figure: {path}")


def cmd_demo(args):
    """Reproduce the paper's contrast on a synthetic C 1s-like spectrum."""
    from .plots import plot_multistart
    x, y = synthetic_spectrum([(285.0, 1.2, 12000), (286.5, 1.3, 9000), (288.5, 1.2, 3000), (291.5, 3.0, 600)],
                              counts_scale=0.2, seed=1)
    outdir = Path(args.out or "fitstability_demo")
    outdir.mkdir(parents=True, exist_ok=True)
    four = [PeakSpec("C1", 285.0, 1.25), PeakSpec("C2", 286.5, 1.25), PeakSpec("C3", 288.5, 1.25), PeakSpec("Sat", 291.5, 3.5)]
    ten_pos = [284.5, 285.2, 285.9, 286.6, 287.3, 288.0, 288.7, 289.4, 290.1]
    ten = [PeakSpec(f"C{i + 1}", p, 1.25) for i, p in enumerate(ten_pos)] + [PeakSpec("Sat", 291.5, 3.5)]
    for label, peaks in (("four_components", four), ("ten_components", ten)):
        prob = FitProblem(x, y, Protocol(peaks, "shirley"), count_scale=0.2)
        ms = multistart(prob, args.n, args.seed, n_jobs=args.jobs)
        v = assess(ms)
        print(f"== {label} ==\n{v}\n")
        plot_multistart(ms, title=label, path=str(outdir / f"{label}.png"), verdict=v)
    print(f"Figures in {outdir}")


def main(argv=None):
    matplotlib.use("Agg")
    ap = argparse.ArgumentParser(prog="fitstability", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def add_library(p):
        p.add_argument("--library", action="append", metavar="PATH",
                       help="reference spectra for CasaXPS library components (%%&<>LS-...): a file or folder of VAMAS or "
                            "two-column text files (repeatable); an origins.json beside them gives calibrated origins")
        p.add_argument("--library-axis", choices=["binding", "kinetic"], default="binding",
                       help="energy axis of text-file references (VAMAS files state their own)")
        p.add_argument("--library-hv", type=float, default=None, help="photon energy, for kinetic-energy text references")

    def common(p):
        p.add_argument("source", help="VAMAS file (.vms) or two-column text/CSV: binding energy, intensity")
        p.add_argument("protocol", nargs="?", default=None,
                       help="protocol JSON: one protocol, or {\"protocols\": {region: protocol}}. Omit for a VAMAS "
                            "file saved from CasaXPS with a fit: the fit (shapes, bounds, ties, background) is the protocol")
        p.add_argument("--block", help="VAMAS blocks: 'all' (default), numbers like 1,3,5-7, or a region name such as 'C 1s'")
        p.add_argument("--count-units", choices=["counts", "cps"], default="counts",
                       help="VAMAS ordinate units for the Poisson weights: counts (default) or counts/s (x dwell x scans)")
        p.add_argument("--count-scale", type=float, default=0.0, help="counts per y unit; overrides --count-units")
        add_library(p)
        p.add_argument("--seed", type=int, default=None)
        p.add_argument("--jobs", type=int, default=1)
        p.add_argument("--out", help="figure file (one block) or folder; default: next to the source")

    p = sub.add_parser("list", help="list the blocks in a VAMAS file")
    p.add_argument("source")
    p.set_defaults(func=cmd_list)

    p = sub.add_parser("library-list", help="list the reference spectra found by --library")
    add_library(p)
    p.set_defaults(func=cmd_library_list)

    p = sub.add_parser("library-calibrate", help="fit the origin of each reference spectrum to CasaXPS's saved fits")
    p.add_argument("source", help="a VAMAS file fitted in CasaXPS with those library components")
    add_library(p)
    p.add_argument("--block", help="VAMAS blocks to use (default: all)")
    p.add_argument("--count-units", choices=["counts", "cps"], default="counts")
    p.add_argument("--write", metavar="FOLDER_OR_FILE", help="write/update origins.json here")
    p.set_defaults(func=cmd_library_calibrate)

    p = sub.add_parser("boxplot", help="multi-start box-plot test (the paper's method)")
    common(p)
    p.add_argument("-n", type=int, default=32, help="number of random starts (paper: 32)")
    p.add_argument("--randomize-pos", type=float, default=0.0, help="also randomise positions, +/- eV")
    p.add_argument("--randomize-fwhm", type=float, default=0.0, help="also randomise widths, +/- fraction")
    p.add_argument("--bounds-fraction", type=float, default=0.0,
                   help="also draw every free position and width over this fraction (0-1) of its allowed range "
                        "(e.g. a CasaXPS constraint window) on each side of its nominal value")
    p.set_defaults(func=cmd_boxplot)

    p = sub.add_parser("montecarlo", help="Poisson-noise Monte Carlo of the best fit")
    common(p)
    p.add_argument("-n", type=int, default=100)
    p.set_defaults(func=cmd_montecarlo)

    p = sub.add_parser("uniqueness", help="uniqueness plots of FREE parameters (fixed/tied ones are skipped)",
                       description=UNIQUENESS_METHOD_NOTE)
    common(p)
    p.add_argument("--peak", action="append", help="peak name (repeatable; default all)")
    p.add_argument("--parameter", choices=["pos", "fwhm", "both"], default="pos",
                       help="pos (paper's example, default); fwhm only where it is free and untied; both")
    p.add_argument("--half-range", type=float, default=0.8)
    p.add_argument("--points", type=int, default=9)
    p.set_defaults(func=cmd_uniqueness)

    p = sub.add_parser("demo", help="synthetic four- vs ten-component comparison")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--out", help="output folder")
    p.add_argument("-n", type=int, default=32)
    p.set_defaults(func=cmd_demo)

    args = ap.parse_args(argv)
    try:
        args.func(args)
    except (VamasError, ValueError, OSError) as e:
        sys.exit(f"error: {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
