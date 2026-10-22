"""Validate the Casa-fit reader: rebuild CasaXPS's stored fit on its own data and compare.

For every block with a CasaXPS fit in the given VAMAS files, every region is turned into a
fitstability model and CasaXPS's stored parameters are evaluated on the data. If the reader's
frames, units, ties, line shapes and backgrounds are right the envelope follows the data to
within noise: the Poisson-weighted RSTD is near 1 and the largest relative deviation is small.

    python scripts/validate_casa.py "path/to/*.vms" ...
"""
import glob
import sys

import numpy as np

from fitstability import read_vamas
from fitstability.casa import CasaFitError, casa_models, describe_skips, parse_fit


def main(patterns):
    rows = []
    for pat in patterns:
        for f in sorted(glob.glob(pat, recursive=True)):
            try:
                vf = read_vamas(f)
            except Exception as e:  # noqa: BLE001
                print(f"{f}: not read ({e})")
                continue
            for b in vf.blocks:
                if parse_fit(b.comments) is None:
                    continue
                try:
                    models = casa_models(b)
                except CasaFitError as e:
                    print(f"{f.split('/')[-1].split(chr(92))[-1]} {b.label}: {e}")
                    continue
                for s in describe_skips(b):
                    print(f"{f.split('/')[-1].split(chr(92))[-1]} {b.label}: {s}")
                for m in models:
                    prob = m.problem()
                    s = m.stored
                    env = prob.model(s)
                    r = (prob.y - env) / prob.sigma
                    rstd = float(np.sqrt(np.mean(r ** 2)))
                    rel = float(np.sqrt(np.mean((prob.y - env) ** 2)) / max(prob.y.max() - prob.y.min(), 1e-30))
                    rows.append((f.replace("\\", "/").split("/")[-1], b.label, m.region, len(m.protocol.peaks), rstd, rel))
    print(f"{'file':<42}{'block':<26}{'region':<14}{'n':>3}{'RSTD':>8}{'rel rms':>9}")
    for r in rows:
        print(f"{r[0][:41]:<42}{r[1][:25]:<26}{r[2][:13]:<14}{r[3]:>3}{r[4]:>8.2f}{r[5]:>9.3f}")
    if rows:
        rel = np.array([r[5] for r in rows])
        print(f"\n{len(rows)} regions; relative rms median {np.median(rel):.3f}, 90th pct {np.percentile(rel, 90):.3f}, max {rel.max():.3f}")


if __name__ == "__main__":
    main(sys.argv[1:] or ["*.vms"])
