"""End-to-end demo of library-spectrum support on synthetic data (no real references needed).

Builds a Casa-fitted VAMAS file whose components are library spectra, writes the reference spectra
to a folder, then runs the CLI: list, calibrate (recovering the hidden origin), and a box plot.
"""
import pathlib
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tests"))
from test_library import make_block, make_reference  # noqa: E402

from fitstability import write_vamas  # noqa: E402


def main():
    d = pathlib.Path(tempfile.mkdtemp(prefix="fitstability_lib_"))
    b1, _, _ = make_block()
    b2, _, _ = make_block(pos4=897.9, pos3=886.3, areas=(500.0, 700.0))
    b2.index = 2
    write_vamas(d / "fit.vms", [b1, b2])
    (d / "refs").mkdir()
    for n, sh in (("LS-CE4", 0.0), ("LS-CE3", -3.0)):
        r = make_reference(n, sh)
        np.savetxt(d / "refs" / f"{n}.csv", np.c_[r.be, r.y], delimiter=",", header="be,counts", comments="")
    exe = pathlib.Path(sys.executable).with_name("fitstability.exe")
    for args in (["library-list", "--library", str(d / "refs")],
                 ["library-calibrate", str(d / "fit.vms"), "--library", str(d / "refs"), "--write", str(d / "refs")],
                 ["library-list", "--library", str(d / "refs")],
                 ["boxplot", str(d / "fit.vms"), "--library", str(d / "refs"), "-n", "6", "--seed", "1", "--out", str(d)]):
        r = subprocess.run([str(exe)] + args, capture_output=True, text=True)
        print(f"$ fitstability {args[0]}\n{(r.stdout + r.stderr)[-1400:]}\n")
    print("files in", d)


if __name__ == "__main__":
    main()
