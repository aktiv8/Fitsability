"""Reference spectra for CasaXPS library-spectrum components.

A CasaXPS component whose line shape is a *library spectrum* (``%&<>LS-CE4``) is not a formula: it is
a measured or computed reference spectrum, shifted to ``Position`` and scaled to ``Area``. The
reference itself lives in CasaXPS's library, not in the VAMAS file, so it has to be supplied here.

Supply them as files or folders (:meth:`ReferenceLibrary.load`):

* VAMAS files: every block is a reference, named by its block id (and sample id, if different). A CasaXPS
  data file that carries its own library spectra as blocks named ``%&<>LS-...`` supplies just those;
* two-column text/CSV files: named by the file name (``LS-CE4.csv`` -> ``LS-CE4``).

Names match the way CasaXPS writes them: ``LS-CE4``, ``%&<>LS-CE4`` and ``CE4`` are the same
reference, case and punctuation ignored. References should already be background-subtracted.

**Origin.** CasaXPS pins a library spectrum at ``Position`` by some reference energy of the spectrum that
the VAMAS file does not state. Each reference carries an ``origin``: the energy on its own axis that sits at
the component's Position. ``None`` means the maximum of the reference, which is what CasaXPS does: checked
on the Ce 3d references LS-CE4 and LS-CE3 against 101 saved CasaXPS fits, the fitted origins were the
maxima to within 0.02 eV with no spread. For a reference where that is not so, :func:`calibrate_origins`
finds the origin by fitting it to CasaXPS's own saved fits (which are in the file); the result is written to
``origins.json`` beside the references and read back on the next load.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

REFERENCE_SUFFIXES = {".vms", ".vamas", ".vam", ".csv", ".txt", ".dat", ".asc", ".xy"}
ORIGINS_FILE = "origins.json"


def reference_key(name: str) -> str:
    """The comparison key of a library-spectrum name: ``'%&<>LS-CE4'``, ``'LS-CE4'`` and ``'ce4'`` agree."""
    s = re.sub(r"^%&<>|^%&", "", str(name).strip()).lstrip("$").upper()
    s = re.sub(r"^LS\d*[-_ ]", "", s)
    return re.sub(r"[^A-Z0-9]", "", s)


@dataclass(frozen=True, eq=False)
class ReferenceSpectrum:
    """A background-subtracted reference on an ascending binding-energy axis.

    ``curve(x, pos, area)`` is the reference moved so that its ``origin`` sits at ``pos`` and scaled to
    unit-area x ``area``; it is zero outside the reference's own energy range.
    """

    name: str
    be: np.ndarray
    y: np.ndarray
    origin: float | None = None
    source: str = ""
    unit: np.ndarray = field(init=False, repr=False, compare=False, default=None)

    def __post_init__(self):
        be = np.asarray(self.be, float)
        y = np.asarray(self.y, float)
        if be.ndim != 1 or be.shape != y.shape or len(be) < 5:
            raise ValueError(f"reference {self.name!r}: need matching 1-D energy and intensity arrays (>= 5 points)")
        order = np.argsort(be)
        be, y = be[order], y[order]
        area = float(np.trapezoid(y, be))
        if not np.isfinite(area) or area <= 0:
            raise ValueError(f"reference {self.name!r} has no positive area")
        object.__setattr__(self, "be", be)
        object.__setattr__(self, "y", y)
        object.__setattr__(self, "unit", y / area)

    @property
    def max_energy(self) -> float:
        return float(self.be[int(np.argmax(self.y))])

    @property
    def effective_origin(self) -> float:
        return self.max_energy if self.origin is None else float(self.origin)

    @property
    def calibrated(self) -> bool:
        return self.origin is not None

    def with_origin(self, origin: float | None) -> "ReferenceSpectrum":
        return replace(self, origin=origin)

    def curve(self, x, pos: float, area: float) -> np.ndarray:
        u = np.asarray(x, float) - pos + self.effective_origin
        return area * np.interp(u, self.be, self.unit, left=0.0, right=0.0)


class ReferenceLibrary:
    """A set of named :class:`ReferenceSpectrum` objects."""

    def __init__(self, items=()):
        self._items: dict[str, ReferenceSpectrum] = {}
        for r in items:
            self.add(r)

    def __len__(self):
        return len(self._items)

    def __bool__(self):
        return bool(self._items)

    def names(self) -> list[str]:
        return sorted(r.name for r in self._items.values())

    def add(self, ref: ReferenceSpectrum) -> None:
        self._items[reference_key(ref.name)] = ref

    def get(self, name: str) -> ReferenceSpectrum | None:
        return self._items.get(reference_key(name))

    def over(self, other: "ReferenceLibrary | None") -> "ReferenceLibrary":
        """This library with ``other``'s references as fall-backs: a name found here wins."""
        if other is None or not other:
            return self
        out = ReferenceLibrary()
        out._items = {**other._items, **self._items}
        return out

    def missing(self, shapes) -> list[str]:
        """The distinct library-spectrum names among ``shapes`` that are not loaded."""
        out = []
        for s in shapes:
            if str(s).startswith("%&") and self.get(s) is None and s.lstrip("%&<>") not in out:
                out.append(s.lstrip("%&<>"))
        return out

    def set_origins(self, origins: dict[str, float]) -> None:
        for name, o in origins.items():
            ref = self.get(name)
            if ref is not None:
                self.add(ref.with_origin(float(o)))

    def origins(self) -> dict[str, float]:
        return {r.name: r.origin for r in self._items.values() if r.origin is not None}

    # -- loading ----------------------------------------------------------------------------
    @classmethod
    def load(cls, paths, axis: str = "binding", hv: float | None = None) -> "ReferenceLibrary":
        """Load references from files and folders (folders are searched recursively).

        ``axis``/``hv`` apply to text files: ``"kinetic"`` converts KE to BE as hv - KE. VAMAS axes are
        read as the file states them. An ``origins.json`` in a folder, or beside a single file, supplies
        calibrated origins.
        """
        lib = cls()
        origin_files = []
        if isinstance(paths, (str, Path)):
            paths = [paths]
        for p in paths:
            p = Path(p)
            if p.is_dir():
                files = sorted(f for f in p.rglob("*") if f.suffix.lower() in REFERENCE_SUFFIXES)
                origin_files.append(p / ORIGINS_FILE)
            elif p.is_file():
                files = [p]
                origin_files.append(p.parent / ORIGINS_FILE)
            else:
                raise FileNotFoundError(f"reference spectra not found: {p}")
            for f in files:
                lib._load_file(f, axis, hv)
        for of in dict.fromkeys(origin_files):
            if of.is_file():
                try:
                    lib.set_origins(json.loads(of.read_text(encoding="utf-8")))
                except (OSError, ValueError) as e:
                    raise ValueError(f"cannot read {of}: {e}") from None
        return lib

    def _load_file(self, f: Path, axis: str, hv: float | None) -> None:
        if f.suffix.lower() in (".vms", ".vamas", ".vam"):
            from .vamas import VamasError, read_vamas
            try:
                vf = read_vamas(f)
            except VamasError:
                return                                   # not a spectrum file: ignore silently
            blocks = vf.blocks
            marked = [b for b in blocks if b.block_id.strip().startswith("%&")]
            if marked:                       # a CasaXPS data file that carries its own library spectra (as %&<>LS-... blocks)
                blocks = marked
            for b in blocks:
                names = [n for n in dict.fromkeys([b.block_id.strip(), b.sample_id.strip()]) if n and n != "Not Specified"]
                if marked:
                    names = names[:1]
                for n in names:
                    try:
                        self.add(ReferenceSpectrum(n, b.be, b.y, None, f"{f.name} block {b.index}"))
                    except ValueError:
                        break
            return
        from .io import read_spectrum
        try:
            x, y = read_spectrum(f)
        except ValueError:
            return
        if axis == "kinetic":
            if not hv:
                raise ValueError("a kinetic-energy reference needs the photon energy (hv)")
            x = hv - x
        try:
            self.add(ReferenceSpectrum(f.stem, x, y, None, f.name))
        except ValueError:
            pass

    def save_origins(self, path) -> None:
        Path(path).write_text(json.dumps(self.origins(), indent=2) + "\n", encoding="utf-8")


# -- calibration -------------------------------------------------------------------------------
@dataclass
class OriginFit:
    name: str
    values: list[float]              # best origin found in each region that uses the reference
    rel_before: float
    rel_after: float

    @property
    def median(self) -> float:
        return float(np.median(self.values))

    @property
    def spread(self) -> float:
        """Interquartile range of the per-region best origins (eV): small means one fixed convention."""
        return float(np.subtract(*np.percentile(self.values, [75, 25]))) if len(self.values) > 1 else 0.0


def calibrate_origins(library: ReferenceLibrary, blocks, area_units: str = "counts", search: float | None = None,
                      step: float = 0.02, coarse: float = 0.1, passes: int = 2, progress=None,
                      cancel=None) -> dict[str, OriginFit]:
    """Find each reference's origin from CasaXPS's own saved fits.

    For every region of ``blocks`` that uses library components, the stored positions and areas are held
    and the origin of each reference involved is scanned, keeping the value with the smallest rms residual
    against the data: first over the whole energy range of the reference in ``coarse`` eV steps (the
    convention is not known to be near the maximum), then +/-``coarse`` around the best in ``step`` steps.
    ``search`` (eV) limits the first scan to +/- that distance about the reference's maximum instead.
    The per-region results are collected per reference; their median is the origin, and a large spread
    means the offset is not a constant of the reference (check the files, or the reference).
    """
    from .casa import CasaFitError, casa_models
    from .model import FitProblem, Protocol

    per_ref: dict[str, list[float]] = {}
    before: dict[str, list[float]] = {}
    after: dict[str, list[float]] = {}
    blocks = list(blocks)
    for bi, b in enumerate(blocks):
        if cancel and cancel():
            from .experiments import Cancelled
            raise Cancelled()
        if progress:
            progress(bi, len(blocks))
        try:
            models = casa_models(b, area_units, library=library)
        except CasaFitError:
            continue
        for m in models:
            peaks = list(m.protocol.peaks)
            idx = [i for i, p in enumerate(peaks) if p.reference is not None]
            if not idx:
                continue
            prob0 = m.problem()
            rng = max(float(prob0.y.max() - prob0.y.min()), 1e-30)

            def rel(pk):
                prob = FitProblem(prob0.x, prob0.y, Protocol(pk, m.protocol.background), prob0.count_scale)
                return float(np.sqrt(np.mean((prob.y - prob.model(m.stored)) ** 2)) / rng)

            def scan(n, ii, grid):
                best = (np.inf, None)
                for o in grid:
                    trial = list(peaks)
                    for j in ii:
                        trial[j] = replace(peaks[j], reference=peaks[j].reference.with_origin(float(o)))
                    r = rel(trial)
                    if r < best[0]:
                        best = (r, float(o))
                return best[1]

            names = {}
            for i in idx:
                names.setdefault(peaks[i].reference.name, []).append(i)
            r0 = rel(peaks)
            cur = {}
            for n, ii in names.items():
                ref = peaks[ii[0]].reference
                if search is None:
                    lo, hi = float(ref.be.min()), float(ref.be.max())
                else:
                    lo, hi = ref.max_energy - search, ref.max_energy + search
                cur[n] = scan(n, ii, np.arange(lo, hi + coarse / 2, coarse))
                for j in ii:
                    peaks[j] = replace(peaks[j], reference=peaks[j].reference.with_origin(cur[n]))
            for _ in range(passes):
                for n, ii in names.items():
                    cur[n] = scan(n, ii, cur[n] + np.arange(-coarse, coarse + step / 2, step))
                    for j in ii:
                        peaks[j] = replace(peaks[j], reference=peaks[j].reference.with_origin(cur[n]))
            r1 = rel(peaks)
            for n in names:
                per_ref.setdefault(n, []).append(cur[n])
                before.setdefault(n, []).append(r0)
                after.setdefault(n, []).append(r1)
    return {n: OriginFit(n, v, float(np.median(before[n])), float(np.median(after[n]))) for n, v in per_ref.items()}
