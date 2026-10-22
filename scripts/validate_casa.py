"""Validate the Casa-fit reader over many VAMAS files.

For every block with a CasaXPS fit, every region is rebuilt as a fitstability model and CasaXPS's
stored parameters are evaluated on the data. If frames, units, ties, line shapes and
backgrounds are right the envelope follows the data to within noise (small ``rel rms``, the rms
residual as a fraction of the data range; ``RSTD`` near 1 for a good fit).

    python scripts/validate_casa.py "dir/**/*.vms" @list.txt --csv out.csv --worst 20

``@file`` reads paths one per line. Files without a CasaXPS fit are skipped and identical files (same bytes) are counted once.
Paths on a cloud-synced folder that are not downloaded (OneDrive placeholders) are read, and so
downloaded, if you list them: pass only local files. The report
groups regions by line-shape family and background type, lists every reason a region or file
could not be reproduced, and the worst-reproduced regions.
"""
import argparse
import collections
import csv
import glob
import hashlib
import re
import sys

import numpy as np

from fitstability import read_vamas
from fitstability.casa import CasaFitError, background_check, casa_models, describe_skips, parse_fit

TOL_OK = 0.03          # rel rms below which a region counts as reproduced ...
TOL_RSTD = 2.0         # ... or Poisson-weighted RSTD below this (rel rms is noise-dominated for weak spectra)


def reproduced(r) -> bool:
    return r["rel"] < TOL_OK or r["rstd"] < TOL_RSTD


def shape_family(shape: str) -> str:
    """Base name of a Casa shape string: 'LA(1.1,1.9,7)' -> 'LA', 'SB(5)DS(1,2)' -> 'SB+DS'."""
    names = re.findall(r"([A-Za-z%&<>\-]+[A-Za-z0-9\-]*)\(", shape)
    if names:
        return "+".join(names)
    return "library" if shape.startswith("%") else shape[:12]


def collect(paths):
    seen, files = set(), []
    for p in paths:
        try:
            data = open(p, "rb").read()
        except OSError:
            continue
        if b"CASA comp" not in data:                    # no CasaXPS fit in it: nothing to validate
            continue
        h = hashlib.sha1(data).hexdigest()
        if h not in seen:
            seen.add(h)
            files.append(p)
    return files


def run(files):
    rows, skips, unread = [], collections.Counter(), collections.Counter()
    run.bgchecks = []
    for f in files:
        short = f.replace("\\", "/").split("/")[-1]
        try:
            vf = read_vamas(f)
        except Exception as e:  # noqa: BLE001
            unread[str(e)[:80]] += 1
            continue
        for b in vf.blocks:
            fit = parse_fit(b.comments)
            if fit is None or not fit.comps:
                continue
            fams = sorted({shape_family(c.shape) for c in fit.comps})
            run.bgchecks += background_check(b)
            bgs = sorted({r.background for r in fit.regions}) or ["?"]
            try:
                models = casa_models(b)
            except CasaFitError as e:
                skips[re.sub(r"[-\d.]+", "#", str(e).split(": ", 1)[-1])[:90]] += 1
                continue
            for s in describe_skips(b):
                skips[re.sub(r"'[^']*'", "'x'", s)[:110]] += 1
            for m in models:
                prob = m.problem()
                rows.append(dict(file=short, block=b.label, region=m.region, n=len(m.protocol.peaks),
                                 shapes="/".join(fams), strings="|".join(sorted({q.shape for q in m.protocol.peaks})),
                                 background=m.protocol.background.kind,
                                 rel=m.rel_rms, rstd=m.rstd, notes=len(m.notes), points=len(prob.x)))
    return rows, skips, unread


def report(rows, skips, unread, worst):
    print(f"{len(rows)} regions reproduced to some degree")
    if unread:
        print("\nFiles that could not be read:")
        for k, v in unread.most_common():
            print(f"  {v:>3}  {k}")
    if skips:
        print("\nRegions/blocks skipped (not reproducible):")
        for k, v in skips.most_common():
            print(f"  {v:>4}  {k}")
    checks = getattr(run, "bgchecks", [])
    if checks:
        print("\nBackground sanity (every region, no line shapes needed): reconstructed background meets the data at "
              "both ends and stays below it")
        g = collections.defaultdict(list)
        for c in checks:
            g[c.kind].append(c.ok)
        for k, v in sorted(g.items(), key=lambda kv: -len(kv[1])):
            print(f"  {k[:30]:<32}{len(v):>6} regions  {100 * np.mean(v):>5.0f} % pass")
    if not rows:
        return

    def table(title, key):
        groups = collections.defaultdict(list)
        for r in rows:
            groups[r[key]].append(r)
        print(f"\nBy {title}: regions, median rel rms, median RSTD, share reproduced (rel < {TOL_OK:.0%} or RSTD < {TOL_RSTD:g})")
        for k, v in sorted(groups.items(), key=lambda kv: -len(kv[1])):
            print(f"  {k[:34]:<36}{len(v):>5}  {np.median([r['rel'] for r in v]):>7.3f}  {np.median([r['rstd'] for r in v]):>7.2f}"
                  f"  {100 * np.mean([reproduced(r) for r in v]):>5.0f} %")
    table("line shape family", "shapes")
    table("background", "background")
    one = [r for r in rows if "|" not in r["strings"]]
    groups = collections.defaultdict(list)
    for r in one:
        groups[r["strings"]].append(r["rel"])
    bad = sorted(((np.median(v), k, len(v)) for k, v in groups.items() if len(v) >= 2 and np.median(v) > TOL_OK),
                 reverse=True)[:12]
    if bad:
        print(f"\nExact shape strings (regions using only that string, 2+ regions) with median rel rms above {TOL_OK:.0%}:")
        for med, k, n in bad:
            print(f"  {k[:30]:<32}{n:>4} regions  median {med:.3f}")
    rel = np.array([r["rel"] for r in rows])
    print(f"\nOverall: median rel rms {np.median(rel):.3f}, 90th percentile {np.percentile(rel, 90):.3f}; "
          f"{100 * np.mean([reproduced(r) for r in rows]):.1f} % reproduced (rel < {TOL_OK:.0%} or RSTD < {TOL_RSTD:g}); "
          f"{int((rel > 0.15).sum())} regions above 15 % rel rms")
    print(f"\nWorst {worst} regions:")
    print(f"  {'rel rms':>8} {'RSTD':>8}  {'n':>3}  {'background':<16}{'shapes':<14}file / block")
    for r in sorted(rows, key=lambda r: -r["rel"])[:worst]:
        print(f"  {r['rel']:>8.3f} {r['rstd']:>8.1f}  {r['n']:>3}  {r['background'][:15]:<16}{r['shapes'][:13]:<14}{r['file'][:40]} {r['block'][:30]}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sources", nargs="+", help="globs (recursive) or @listfile")
    ap.add_argument("--csv", help="write every region to this CSV")
    ap.add_argument("--worst", type=int, default=15)
    args = ap.parse_args(argv)
    paths = []
    for s in args.sources:
        if s.startswith("@"):
            paths += [ln.strip() for ln in open(s[1:], encoding="utf-8") if ln.strip()]
        else:
            paths += glob.glob(s, recursive=True)
    files = collect(sorted(set(paths)))
    print(f"{len(paths)} paths, {len(files)} distinct files")
    rows, skips, unread = run(files)
    report(rows, skips, unread, args.worst)
    if args.csv and rows:
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    sys.exit(main())
