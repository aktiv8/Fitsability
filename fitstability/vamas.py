"""VAMAS (ISO 14976) reader for XPS spectra, plus a minimal writer for tests and demos.

Supported: experiment modes NORM and MAP (single analysis position, as CasaXPS exports of
saved spectra use), scan mode REGULAR, technique XPS. Several blocks (regions, repeats, depth
levels, tilt angles) per file are returned together so a protocol can be run over each.

Conventions and quirks, checked against Thermo K-Alpha and CasaXPS-written files:

* ``1e37`` means "not specified" and is returned as ``None``.
* With several corresponding variables (Counts and Transmission on the K-Alpha) the ordinate
  values are interleaved and their count is the total, so each block has n / n_vars points.
* The abscissa is read as given: Kinetic Energy is converted to binding energy as
  ``hv - KE`` (no work function), Binding Energy is used as is. An unlabelled axis is taken as kinetic
  (Casa writes text-converted files so, with hv = 0 and negative energies: BE = -axis). Axes are returned ascending.
* MAP experiment mode adds three header lines and four per-block lines (x, y, field of view);
  they are skipped.
* Whether the ordinate is total counts or counts per second is not reliably stated in the
  file; see :meth:`Block.count_scale`.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

NOT_SPECIFIED = 1e37
FORMAT_ID = "VAMAS Surface Chemical Analysis Standard Data Transfer Format 1988 May 4"
SUPPORTED_MODES = ("NORM", "MAP")


class VamasError(ValueError):
    """The file is not a VAMAS file this reader can handle."""


def _opt(text: str):
    try:
        v = float(text)
    except ValueError:
        raise VamasError(f"expected a number, found {text!r}") from None
    return None if abs(v) >= NOT_SPECIFIED * 0.99 else v


@dataclass
class Block:
    index: int                    # 1-based position in the file
    block_id: str
    sample_id: str = ""
    species: str = ""
    transition: str = ""
    hv_ev: float = 0.0
    be: np.ndarray = field(default_factory=lambda: np.zeros(0))   # ascending binding energy (eV)
    y: np.ndarray = field(default_factory=lambda: np.zeros(0))
    y_name: str = "counts"
    dwell_s: float | None = None
    scans: int | None = None
    pass_energy_ev: float | None = None
    exp_values: list = field(default_factory=list)
    comments: list = field(default_factory=list)
    timestamp: str = ""
    # reference spectra stored as blocks named "%&<>LS-..." in the same file (CasaXPS library components)
    embedded_library: object = field(default=None, repr=False, compare=False)

    @property
    def region(self) -> str:
        """Short region name, e.g. "C 1s" (species + transition), falling back to the block id."""
        return f"{self.species} {self.transition}".strip() or self.block_id

    @property
    def label(self) -> str:
        return f"#{self.index} {self.region}" + (f" [{self.block_id}]" if self.block_id and self.block_id != self.region else "")

    def count_scale(self, units: str = "counts") -> float:
        """Factor turning ``y`` into counts, for the Poisson weights.

        ``"counts"``: ``y`` is already counts (scale 1). ``"cps"``: ``y`` is counts per
        second, so counts = y x dwell x scans; falls back to 1.0 if dwell or scans are not
        in the file. VAMAS does not say reliably which the ordinate is, so the caller picks.
        """
        if units == "counts":
            return 1.0
        if units == "cps":
            return float(self.dwell_s * self.scans) if self.dwell_s and self.scans else 1.0
        raise ValueError("units must be 'counts' or 'cps'")

    def describe(self) -> str:
        pts = len(self.be)
        rng = f"{self.be.min():.2f}-{self.be.max():.2f} eV" if pts else "empty"
        extra = []
        if self.dwell_s:
            extra.append(f"dwell {self.dwell_s:g} s")
        if self.scans:
            extra.append(f"{self.scans} scans")
        if self.pass_energy_ev:
            extra.append(f"PE {self.pass_energy_ev:g}")
        return f"{self.label}: {pts} points, {rng}" + (f" ({', '.join(extra)})" if extra else "")


@dataclass
class VamasFile:
    path: str
    mode: str
    blocks: list

    def regions(self) -> list[str]:
        seen: dict = {}
        for b in self.blocks:
            seen.setdefault(b.region, None)
        return list(seen)

    def select(self, spec: str | None = None) -> list[Block]:
        """Pick blocks by ``spec``: None/"all" = every block; a number (1-based) or comma/range
        list like "1,3,5-7"; or a region/block-id substring ("C 1s", "c1s"; spaces ignored)."""
        if spec is None or spec.strip().lower() == "all":
            return list(self.blocks)
        s = spec.strip()
        if re.fullmatch(r"[\d,\s\-]+", s):
            wanted: set = set()
            for part in s.split(","):
                part = part.strip()
                if "-" in part:
                    a, b = part.split("-")
                    wanted.update(range(int(a), int(b) + 1))
                elif part:
                    wanted.add(int(part))
            picked = [b for b in self.blocks if b.index in wanted]
            missing = wanted - {b.index for b in picked}
            if missing:
                raise VamasError(f"no block number(s) {sorted(missing)}; the file has 1-{len(self.blocks)}")
            return picked
        key = _norm(s)
        picked = [b for b in self.blocks if key in _norm(b.region) or key in _norm(b.block_id)]
        if not picked:
            raise VamasError(f"no block matches {spec!r}; regions in file: {', '.join(self.regions())}")
        return picked


def _norm(text: str) -> str:
    return re.sub(r"\s+", "", text).lower()


class _Lines:
    def __init__(self, lines):
        self.lines, self.pos = lines, 0

    def next(self) -> str:
        if self.pos >= len(self.lines):
            raise VamasError("file ends early (truncated VAMAS file?)")
        v = self.lines[self.pos]
        self.pos += 1
        return v.strip()

    def int(self) -> int:
        t = self.next()
        try:
            return int(float(t))
        except ValueError:
            raise VamasError(f"expected an integer, found {t!r}") from None

    def skip(self, n: int):
        for _ in range(n):
            self.next()


def read_vamas(path) -> VamasFile:
    try:
        text = Path(path).read_text(encoding="latin-1")
    except OSError as exc:
        raise VamasError(f"cannot open {path}: {exc}") from None
    return parse(text, str(path))


def parse(text: str, path: str = "") -> VamasFile:
    lines = text.splitlines()
    if not lines or not lines[0].strip().upper().startswith("VAMAS"):
        raise VamasError("not a VAMAS file (first line should start with 'VAMAS Surface Chemical Analysis')")
    ln = _Lines(lines)
    ln.skip(5)                       # format id, institution, instrument, operator, experiment id
    ln.skip(ln.int())                # file comments
    mode, scan = ln.next().upper(), ln.next().upper()
    if mode not in SUPPORTED_MODES:
        raise VamasError(f"experiment mode {mode} is not supported (NORM and MAP only)")
    if scan != "REGULAR":
        raise VamasError(f"scan mode {scan} is not supported (REGULAR only)")
    if mode == "MAP":
        ln.skip(3)                   # analysis positions, x and y coordinate counts
    ln.int()                         # number of spectral regions
    nev = ln.int()
    ln.skip(2 * nev)                 # experimental variable labels and units
    ln.skip(4)                       # inclusion list, manually entered items, two upgrade counts
    nblocks = ln.int()
    blocks = [_read_block(ln, mode, nev, i + 1) for i in range(nblocks)]
    _attach_embedded_library(blocks)
    return VamasFile(path, mode, blocks)


def _attach_embedded_library(blocks) -> None:
    """CasaXPS saves the library spectra a fit uses as extra blocks named ``%&<>LS-...`` in the same file;
    make them available to every block (see :func:`fitstability.casa.casa_models`)."""
    marked = [b for b in blocks if b.block_id.strip().startswith("%&")]
    if not marked:
        return
    from .library import ReferenceLibrary, ReferenceSpectrum
    lib = ReferenceLibrary()
    for b in marked:
        try:
            lib.add(ReferenceSpectrum(b.block_id.strip(), b.be, b.y, None, f"block {b.index} of this file"))
        except ValueError:
            continue
    for b in blocks:
        b.embedded_library = lib


def _read_block(ln: _Lines, mode: str, nev: int, index: int) -> Block:
    b = Block(index=index, block_id=ln.next(), sample_id=ln.next())
    y, mo, d, h, mi, s = (ln.next() for _ in range(6))
    ln.next()                        # hours from GMT
    b.timestamp = f"{y}-{mo.zfill(2)}-{d.zfill(2)} {h.zfill(2)}:{mi.zfill(2)}:{s.zfill(2)}"
    b.comments = [ln.next() for _ in range(ln.int())]
    technique = ln.next()
    if technique.upper() != "XPS":
        raise VamasError(f"block {index} ({b.block_id!r}) uses technique {technique}; only XPS is supported")
    if mode == "MAP":
        ln.skip(2)                   # x, y coordinates
    b.exp_values = [_opt(ln.next()) for _ in range(nev)]
    ln.next()                        # analysis source label
    b.hv_ev = _opt(ln.next()) or 0.0
    ln.skip(5)                       # strength, beam widths, incidence polar and azimuth
    if mode == "MAP":
        ln.skip(2)                   # field of view x, y
    ln.next()                        # analyser mode
    b.pass_energy_ev = _opt(ln.next())
    ln.skip(7)                       # magnification, work function, bias, widths, take-off angles
    b.species, b.transition = ln.next(), ln.next()
    ln.next()                        # charge of detected particle
    abscissa = ln.next().lower()
    ln.next()                        # abscissa unit
    start, step = float(ln.next()), float(ln.next())
    ncv = ln.int()
    names = []
    for _ in range(ncv):
        names.append(ln.next())
        ln.next()
    ln.next()                        # signal mode
    b.dwell_s = _opt(ln.next())
    sc = _opt(ln.next())
    b.scans = int(sc) if sc is not None else None
    ln.next()                        # signal time correction
    ln.skip(3)                       # sample-normal polar, azimuth, rotation
    for _ in range(ln.int()):        # additional numerical parameters (label, unit, value)
        ln.skip(3)
    n = ln.int()
    ln.skip(2 * ncv)                 # min and max of each corresponding variable
    if ncv < 1 or n % ncv:
        raise VamasError(f"block {index} has {n} values for {ncv} variables")
    data = np.array([float(ln.next()) for _ in range(n)]).reshape(-1, ncv)
    axis = start + step * np.arange(data.shape[0])
    if "binding" in abscissa:
        be = axis
    elif "kinetic" in abscissa or not abscissa.strip():
        be = b.hv_ev - axis     # an unlabelled axis (Casa-converted text files) is a kinetic-energy axis, hv possibly 0
    else:
        raise VamasError(f"block {index}: unknown energy axis {abscissa!r}")
    order = np.argsort(be)
    col = next((i for i, nm in enumerate(names) if re.search(r"count|intens", nm, re.I)), 0)
    b.be = be[order]
    b.y = data[order, col]
    b.y_name = names[col]
    return b


# -- writer (tests, demos) -----------------------------------------------------------------
def _num(v) -> str:
    return "1e+37" if v is None else f"{v:.10g}"


def write_vamas(path, blocks: list[Block], mode: str = "NORM") -> None:
    """Write blocks as a minimal NORM (or MAP) XPS VAMAS file with a single Counts column."""
    map_mode = mode == "MAP"
    out = [FORMAT_ID, "Not Specified", "fitstability", "Not Specified", "Not Specified", "0", mode, "REGULAR"]
    if map_mode:
        out += ["1", "0", "0"]
    out += [str(len({b.region for b in blocks})), "1", "Exp Variable", "d", "0", "0", "0", "0", str(len(blocks))]
    for b in blocks:
        be = np.asarray(b.be, float)
        step = float(be[1] - be[0])
        ke = b.hv_ev - be[::-1]
        out += [b.block_id, b.sample_id, "2025", "1", "1", "0", "0", "0", "0", str(len(b.comments))] + list(b.comments)
        out += ["XPS"]
        if map_mode:
            out += ["16384", "16384"]
        out += [_num(b.exp_values[0] if b.exp_values else None), "KA1486", _num(b.hv_ev)] + ["1e+37"] * 5
        if map_mode:
            out += ["1e+37"] * 2
        out += ["FAT", _num(b.pass_energy_ev)] + ["1e+37"] * 7
        out += [b.species, b.transition, "-1", "Kinetic Energy", "eV", _num(float(ke[0])), _num(step), "1",
                b.y_name, "d", "pulse counting", _num(b.dwell_s), _num(b.scans), "1e+37", "1e+37", "1e+37", "1e+37",
                "0", str(len(be)), "0", "1"]
        out += [f"{v:.10g}" for v in np.asarray(b.y, float)[::-1]]
    out.append("end of experiment")
    Path(path).write_bytes(("\r\n".join(out) + "\r\n").encode("latin-1"))  # bytes: text mode would double the CR on Windows
