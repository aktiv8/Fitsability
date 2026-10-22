# fitstability

Peak-fit stability diagnostics for XPS, after Moeini, Haack, Fairley, Fernandez, Gengenbach, Easton and
Linford, *Box plots: a simple graphical tool for visualizing overfitting in peak fitting*,
J. Electron Spectrosc. Relat. Phenom. 250 (2021) 147094. GUI-free core, CLI and matplotlib plots.

## GUI

```bash
.venv\Scripts\python -m pip install -e ".[gui]"     # adds PySide6
.venv\Scripts\fitstability-gui                       # or: python -m fitstability.gui [file]
```

A desktop window (light/dark) with the same engine as the command line:

- **Open** a VAMAS (`.vms`) or CSV/text file with the file dialog (Ctrl+O), drag it onto the window, or pick
  from *Open recent*. A VAMAS file lists its blocks, marking those that hold a CasaXPS fit; a CSV opens an
  import dialog (preview, delimiter, energy and intensity columns, binding or kinetic energy, counts or
  counts/s).
- **Peak model:** for a Casa-fitted block the fit is loaded into an editable table (shape, position, FWHM,
  bounds, and constraints as short rules: `free`, `fixed`, `C1 + 5.2`, `C1 * 1.2`, `group main`). For CSV
  data, or any block, add peaks by hand or load/save the protocol as JSON. *Reset to Casa fit* undoes edits.
- **Run:** *Fit once*, *Box plot test*, *Monte Carlo*, *Uniqueness* (F5-F9), in the background with a progress
  bar and Cancel. Results appear in their own tabs with the verdict and notes; editing the model marks them stale.
  Uniqueness lists the tied or fixed parameters it did not scan.
- **Export:** the current figure (PNG/PDF/SVG), the multi-start CSV, a text report; *Run > Box plot on every
  block with a Casa fit* analyses a whole file and saves a summary CSV and one figure per block.

## Methods

| Command | Method | Reads as |
|---|---|---|
| `boxplot` | **Multi-start box plot (the paper's method).** N random starts (default 32): areas = U(0,1) × the largest area the peak can take before touching the data envelope; positions and widths start nominal (optionally randomised). Poisson-weighted χ² fit, restarted until converged. | Narrow final boxes: one minimum. Wide boxes: multiple local minima. |
| `montecarlo` | Poisson-noise replicates of the best fit, refitted. | Parameter uncertainty. |
| `uniqueness` | Fix one position or FWHM on a grid, refit the rest, plot MSE. | U/V: significant. Flat: correlated with other parameters. |

### Uniqueness plots: free parameters only (a point of this method)

As in the paper, the data are fitted with the stated protocol, one parameter of one peak is then **fixed** at a
series of values, and the fit is repeated with all the other parameters free. The curve of MSE against the fixed
value is U/V-shaped when the parameter is significant, and flat when the others compensate for it.

It follows that **only parameters that are free in the protocol can be probed**. A fixed position, a fixed FWHM,
or a FWHM tied to other peaks (the common constraint) is not a fit parameter, so it is not scanned and
`uniqueness()` raises an explanatory error for it. With the usual tied-width protocol, uniqueness is
therefore a test of peak *positions* (the paper's own example); FWHM enters only where it is free and
independent. The CLI prints which parameters it skipped and why.

Narrow boxes show a single minimum, not a good fit, so read them with the residual panel and RSTD.
The `single minimum` / `multiple minima` verdict uses this package's own thresholds (`verdict.assess`),
since the paper judges width by eye.

## Input: VAMAS and text files

`.vms`/`.vamas` files (NORM and MAP modes; Kinetic or Binding Energy axis; checked against K-Alpha and
CasaXPS files) load every block of the file, so several regions, repeats, depth levels or angles run together:

```bash
fitstability list data.vms                                   # blocks in the file
fitstability boxplot data.vms protocol.json --block "C 1s"   # a region name, numbers (1,3,5-7) or all
```

A protocol file is one protocol for every selected block, or `{"protocols": {"C 1s": {...}, "O 1s": {...}}}` to
give each region its own. One figure per block is written next to the source. VAMAS does not say reliably
whether the ordinate is counts or counts per second, so the Poisson weighting assumes counts; use
`--count-units cps` (counts = y x dwell x scans) or `--count-scale` if yours is not. Python:
`read_vamas(path).select("C 1s")` gives blocks with `.be` and `.y` ready for `FitProblem`.

## Protocol from a CasaXPS fit (VAMAS saved from Casa)

Omit the protocol file and the fit CasaXPS stored in each block's comments *is* the protocol:

```bash
fitstability boxplot "Titanium Metal Depth Profile.vms" --block 11 --bounds-fraction 0.5
fitstability uniqueness "Titanium Metal Depth Profile.vms" --block 11 --parameter both
```

What is read (`fitstability.casa`): every `CASA region` (window, background type, end-point averaging,
offsets, Tougaard parameters) and `CASA comp` line (line shape, and for Area / FWHM / Position the start
value, bounds and tie). Ties are Casa's: a component tied to another by ratio (Area, FWHM) or offset
(Position), as for spin-orbit doublets. Checked against 253 VAMAS files: every tied value equals its
reference combined that way, and no other constraint form occurs. Lower = upper holds a parameter fixed.

- **Frames and units:** positions and windows are kinetic energies; they are converted to the calibrated
  binding energies CasaXPS shows (the `Calib` line). Areas are counts/s x eV and are multiplied by
  dwell x scans for a counts ordinate (`--count-units cps` if the ordinate is already counts/s).
- **Line shapes and backgrounds:** GL/SGL are exact; LA and LF (and DS) are reconstructions of CasaXPS's
  from the author's eXPoSe SpectraDeck (`casa_lineshapes.py`, vendored). Shirley, Linear and Tougaard
  (U2/U3/W) backgrounds are recomputed from the data in the region window. Library-spectrum
  components, spline backgrounds and E-Tougaard cannot be reproduced: that region is skipped and the reason printed.
- **Check:** each region reports how well CasaXPS's saved fit reproduces the data (`rel_rms`, `RSTD`).
  `scripts/validate_casa.py` does this over many files: on 313 regions in the author's files the median
  relative rms is 0.6 % of the data range (90th percentile 2.2 %). A few files do not reproduce (areas in
  other units, an edited file): the area bounds are then opened and only the tie ratios are used.
- **Starts:** by default the Casa fit's own positions and widths are the nominal start, so only the areas
  are random (the paper's protocol). Casa's saved values are usually already an optimum, so for a real test
  also spread positions and widths: `--bounds-fraction f` draws each free position and width uniformly
  over fraction f of its Casa constraint window on each side of its saved value.
- **Uniqueness** probes only what is free in Casa: tied positions, tied or fixed FWHM are listed as not
  scanned, with the tie that removes them.
- **Not handled:** overlapping regions that CasaXPS merges into one segmented background (each region is
  analysed on its own), and VAMAS experiment modes other than NORM and MAP.

## Use

```bash
.venv\Scripts\python -m pip install -e .[dev]
.venv\Scripts\fitstability demo --out examples/demo          # 4- vs 10-component synthetic comparison
.venv\Scripts\fitstability boxplot spectrum.csv examples/protocol_c1s.json -n 32 --jobs 4
.venv\Scripts\fitstability montecarlo spectrum.csv examples/protocol_c1s.json
.venv\Scripts\fitstability uniqueness spectrum.csv examples/protocol_c1s.json --peak C2
.venv\Scripts\python -m pytest
```

Python API: `FitProblem(x, y, Protocol([PeakSpec(...)], "shirley"), count_scale)`, then `multistart`,
`assess`, `monte_carlo`, `uniqueness`; figures in `fitstability.plots`.

## Scope and differences from the paper

- Hand-written protocols use an area-normalised pseudo-Voigt (`eta`, fixed per peak) and a linear or Shirley
  background; a protocol read from a Casa fit uses Casa's own line shapes and backgrounds (above). The
  background is computed from the data and held fixed during the fits, as in CasaXPS.
- Input is VAMAS or two-column text/CSV. `count_scale` converts y to counts for the Poisson weights (counts = y x scale).
- Tied widths use `width_group`; `fix_pos` / `fix_fwhm` hold values at their nominal start.
- Uniqueness supports position and FWHM (free parameters only), not area.
