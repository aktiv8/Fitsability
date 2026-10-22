"""The fitstability main window.

Left: the data (open a VAMAS or CSV file, pick a block, the fit window). Centre: plots
(spectrum and model, box plot, Monte Carlo, uniqueness). Right: the peak model (from the
CasaXPS fit in the file, or edited by hand) and the run controls.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
from PySide6.QtCore import QSettings, Qt, QTimer
from PySide6.QtGui import QAction, QFont, QKeySequence
from PySide6.QtWidgets import (QAbstractItemView, QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog,
                               QFormLayout, QGridLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel, QMainWindow, QMessageBox,
                               QPlainTextEdit, QProgressBar, QPushButton, QScrollArea, QSpinBox, QSplitter, QTableView,
                               QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from .. import __version__
from .. import plots
from ..casa import CasaFitError, casa_models, describe_skips, parse_fit
from ..experiments import (UNIQUENESS_METHOD_NOTE, monte_carlo, multistart, scannable_parameters, uniqueness,
                           unscannable_reason)
from ..fit import fit, fit_nominal
from ..io import load_protocol, save_protocol
from ..model import FitProblem, PeakSpec, Protocol
from ..vamas import VamasError, VamasFile, read_vamas
from ..verdict import assess
from . import protocol_model as pm
from . import results, theme
from .batch_dialog import BatchResultsDialog
from .csv_dialog import CsvImportDialog
from .peak_table import PeakTableModel
from .plot_pane import PlotPane
from .workers import Worker

VAMAS_SUFFIXES = {".vms", ".vamas", ".vam"}
OPEN_FILTER = "Spectra (*.vms *.vamas *.csv *.txt *.dat *.asc);;VAMAS (*.vms *.vamas);;CSV / text (*.csv *.txt *.dat *.asc);;All files (*)"
TAB_SPECTRUM, TAB_BOX, TAB_MC, TAB_UNIQ = range(4)


def _spin(lo, hi, value, step=1, decimals=None, suffix=""):
    w = QSpinBox() if decimals is None else QDoubleSpinBox()
    w.setRange(lo, hi)
    if decimals is not None:
        w.setDecimals(decimals)
    w.setSingleStep(step)
    w.setValue(value)
    if suffix:
        w.setSuffix(suffix)
    return w


class MainWindow(QMainWindow):
    def __init__(self, theme_name: str | None = None, settings: QSettings | None = None):
        super().__init__()
        self.settings = settings or QSettings("fitstability", "fitstability-gui")
        self.theme_name = theme_name or str(self.settings.value("theme", "light"))
        self.setWindowTitle("fitstability")
        self.resize(1480, 880)
        self.setAcceptDrops(True)

        # -- data state ---------------------------------------------------------------------
        self.path: str | None = None
        self.vamas: VamasFile | None = None
        self.block = None
        self.base_x = self.base_y = None          # the whole spectrum (calibrated binding energy ascending)
        self.window = None                        # (lo, hi) fit window
        self.csv_count_scale = 1.0
        self.casa = []                            # CasaModel list for the current block
        self.casa_model = None                    # the one in use (None: hand-made model)
        self.casa_bg = None                       # callable background from the Casa region
        self.display_start = None                 # parameters to draw in the spectrum tab
        self.notes: list[str] = []
        self._loading = False
        # -- results ------------------------------------------------------------------------
        self.ms = self.verdict = self.mc = None
        self.uniq: list = []
        self.uniq_skipped: list[str] = []
        self.worker: Worker | None = None
        self._queue: list = []
        self._title = ""

        self._build()
        self.apply_theme(self.theme_name)
        self._restore_settings()
        self._refresh_timer = QTimer(self, singleShot=True, interval=120)
        self._refresh_timer.timeout.connect(self.refresh_spectrum)
        self._update_actions()
        self.refresh_spectrum()

    # ======================================================================================
    # construction
    # ======================================================================================
    def _build(self):
        split = QSplitter(Qt.Horizontal)
        split.addWidget(self._build_left())
        split.addWidget(self._build_centre())
        split.addWidget(self._build_right())
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)
        split.setStretchFactor(2, 0)
        split.setSizes([330, 700, 520])
        self.splitter = split
        self.setCentralWidget(split)
        self.progress = QProgressBar()
        self.progress.setMaximumWidth(220)
        self.progress.setVisible(False)
        self.status = self.statusBar()
        self.status.addPermanentWidget(self.progress)
        self.status.showMessage("Open a VAMAS or CSV spectrum (Ctrl+O), or drop a file on the window.")
        self._build_menus()

    def _build_left(self):
        box = QGroupBox("Data")
        lay = QVBoxLayout(box)
        row = QHBoxLayout()
        self.btn_open = QPushButton("Open spectrum...")
        self.btn_open.setObjectName("primary")
        self.btn_open.clicked.connect(self.open_dialog)
        row.addWidget(self.btn_open)
        self.btn_reimport = QPushButton("Import options...")
        self.btn_reimport.setToolTip("Change the columns, axis or units of the CSV file")
        self.btn_reimport.clicked.connect(self.reimport_csv)
        row.addWidget(self.btn_reimport)
        lay.addLayout(row)
        self.lbl_file = QLabel("No file loaded")
        self.lbl_file.setWordWrap(True)
        self.lbl_file.setObjectName("muted")
        lay.addWidget(self.lbl_file)

        lay.addWidget(QLabel("Blocks in the file:"))
        self.blocks = QTableWidget(0, 4)
        self.blocks.setHorizontalHeaderLabels(["#", "Region", "Points", "Casa fit"])
        self.blocks.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.blocks.setSelectionMode(QAbstractItemView.SingleSelection)
        self.blocks.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.blocks.verticalHeader().setVisible(False)
        hdr = self.blocks.horizontalHeader()
        hdr.setStretchLastSection(False)
        hdr.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        hdr.setSectionResizeMode(1, QHeaderView.Stretch)
        hdr.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        hdr.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        self.blocks.itemSelectionChanged.connect(self._block_selected)
        lay.addWidget(self.blocks, 1)

        form = QFormLayout()
        self.units = QComboBox()
        self.units.addItems(["Counts", "Counts per second"])
        self.units.setToolTip("What the file's intensity is. Sets the Poisson weighting and the Casa area scale.")
        self.units.currentIndexChanged.connect(self._units_changed)
        form.addRow("Intensity is", self.units)
        self.win_lo = _spin(-1e5, 1e5, 0.0, 0.5, 3, " eV")
        self.win_hi = _spin(-1e5, 1e5, 0.0, 0.5, 3, " eV")
        for w in (self.win_lo, self.win_hi):
            w.editingFinished.connect(self._window_changed)
        form.addRow("Fit from", self.win_lo)
        form.addRow("Fit to", self.win_hi)
        self.bg = QComboBox()
        self.bg.addItems(["Shirley", "Linear", "None"])
        self.bg.currentIndexChanged.connect(self._model_changed)
        form.addRow("Background", self.bg)
        lay.addLayout(form)
        return box

    def _build_centre(self):
        self.tabs = QTabWidget()
        self.panes = [PlotPane(self.theme_name) for _ in range(4)]
        for pane, name in zip(self.panes, ["Spectrum && model", "Box plot", "Monte Carlo", "Uniqueness"]):
            self.tabs.addTab(pane, name)
        self.panes[TAB_BOX].hint.setText("The paper's test: many random starting conditions, one protocol. Narrow final "
                                         "boxes = one minimum; wide boxes = multiple local minima.")
        self.panes[TAB_MC].hint.setText("Poisson-noise replicates of the best fit, each refitted: parameter uncertainty.")
        self.panes[TAB_UNIQ].hint.setText(UNIQUENESS_METHOD_NOTE)
        for p, msg in ((TAB_BOX, "Run the box plot test"), (TAB_MC, "Run Monte Carlo"), (TAB_UNIQ, "Run uniqueness")):
            self.panes[p].show_message(msg)
        return self.tabs

    def _build_right(self):
        panel = QWidget()
        lay = QVBoxLayout(panel)
        lay.setContentsMargins(0, 0, 0, 0)

        model_box = QGroupBox("Peak model")
        ml = QVBoxLayout(model_box)
        top = QHBoxLayout()
        top.addWidget(QLabel("Casa region:"))
        self.region = QComboBox()
        self.region.currentIndexChanged.connect(self._region_selected)
        top.addWidget(self.region, 1)
        self.btn_reset = QPushButton("Reset to Casa fit")
        self.btn_reset.clicked.connect(lambda: self._use_casa(self.region.currentIndex()))
        top.addWidget(self.btn_reset)
        ml.addLayout(top)
        self.table = PeakTableModel(self)
        self.view = QTableView()
        self.view.setModel(self.table)
        self.view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.view.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.view.horizontalHeader().setDefaultSectionSize(84)
        self.view.setColumnWidth(0, 130)
        self.view.setMinimumHeight(190)
        self.table.peaksChanged.connect(self._model_changed)
        self.table.problem.connect(self._table_problem)
        ml.addWidget(self.view, 1)
        for group in ((("Add peak", self.add_peak), ("Duplicate", self.duplicate_peak), ("Remove", self.remove_peak)),
                      (("Load protocol (JSON)...", self.load_protocol_dialog),
                       ("Save protocol (JSON)...", self.save_protocol_dialog))):
            row = QHBoxLayout()
            for text, fn in group:
                b = QPushButton(text)
                b.clicked.connect(fn)
                row.addWidget(b)
            ml.addLayout(row)
        self.lbl_valid = QLabel("")
        self.lbl_valid.setWordWrap(True)
        self.lbl_valid.setObjectName("danger")
        ml.addWidget(self.lbl_valid)
        self.lbl_notes = QLabel("")
        self.lbl_notes.setWordWrap(True)
        self.lbl_notes.setObjectName("muted")
        ml.addWidget(self.lbl_notes)
        lay.addWidget(model_box, 3)

        run_box = QGroupBox("Run")
        rl = QGridLayout(run_box)
        self.n_starts = _spin(2, 2000, 32)
        self.n_starts.setToolTip("Random starting conditions (the paper used 32)")
        self.bounds_fraction = _spin(0.0, 1.0, 0.0, 0.05, 2)
        self.bounds_fraction.setToolTip("Also draw every free position and width over this fraction of its allowed range "
                                        "(for a Casa fit: of its constraint window) on each side of its saved value. "
                                        "0 = only the areas are random, as in the paper.")
        self.rand_pos = _spin(0.0, 10.0, 0.0, 0.1, 2, " eV")
        self.rand_fwhm = _spin(0.0, 1.0, 0.0, 0.05, 2)
        self.rand_fwhm.setToolTip("Fractional half-range, 0.3 = +/-30 %")
        seed_row = QHBoxLayout()
        self.use_seed = QCheckBox("fixed")
        self.use_seed.setChecked(True)
        self.seed = _spin(0, 2 ** 31 - 1, 1)
        seed_row.addWidget(self.use_seed)
        seed_row.addWidget(self.seed, 1)
        self.mc_trials = _spin(5, 5000, 100)
        self.u_param = QComboBox()
        self.u_param.addItems(["Position", "FWHM", "Both"])
        self.u_param.setToolTip("Only free parameters can be scanned; tied or fixed ones are listed as not scanned")
        self.u_range = _spin(0.05, 20.0, 0.8, 0.1, 2, " eV")
        self.u_points = _spin(3, 41, 9)
        self.jobs = _spin(1, max(1, os.cpu_count() or 1), 1)
        self.jobs.setToolTip("Parallel processes (1 = run in the background of this window)")
        seed_w = QWidget()
        seed_w.setLayout(seed_row)
        seed_row.setContentsMargins(0, 0, 0, 0)
        cells = [("Random starts", self.n_starts), ("Spread in bounds", self.bounds_fraction),
                 ("Randomise pos +/-", self.rand_pos), ("Randomise FWHM +/-", self.rand_fwhm),
                 ("Seed", seed_w), ("Processes", self.jobs),
                 ("Monte Carlo trials", self.mc_trials), ("Uniqueness scan", self.u_param),
                 ("Scan half-range", self.u_range), ("Scan points", self.u_points)]
        for i, (label, w) in enumerate(cells):
            r, c = divmod(i, 2)
            rl.addWidget(QLabel(label), r, 2 * c)
            rl.addWidget(w, r, 2 * c + 1)
        lay.addWidget(run_box)

        btns = QHBoxLayout()
        self.btn_fit = QPushButton("Fit once")
        self.btn_box = QPushButton("Box plot test")
        self.btn_box.setObjectName("primary")
        self.btn_mc = QPushButton("Monte Carlo")
        self.btn_uniq = QPushButton("Uniqueness")
        self.btn_cancel = QPushButton("Cancel")
        self.btn_cancel.setEnabled(False)
        self.btn_fit.clicked.connect(self.run_fit)
        self.btn_box.clicked.connect(self.run_box)
        self.btn_mc.clicked.connect(self.run_mc)
        self.btn_uniq.clicked.connect(self.run_uniq)
        self.btn_cancel.clicked.connect(self.cancel)
        for b in (self.btn_fit, self.btn_box, self.btn_mc, self.btn_uniq, self.btn_cancel):
            btns.addWidget(b)
        lay.addLayout(btns)

        self.out = QPlainTextEdit()
        self.out.setReadOnly(True)
        self.out.setFont(QFont("Consolas", 9))
        self.out.setPlaceholderText("Results appear here.")
        self.out.setMinimumHeight(170)
        lay.addWidget(self.out, 2)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(panel)
        scroll.setMinimumWidth(480)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        return scroll

    def _build_menus(self):
        mb = self.menuBar()
        f = mb.addMenu("&File")
        self.act_open = self._act(f, "&Open spectrum...", self.open_dialog, QKeySequence.Open)
        self.recent_menu = f.addMenu("Open &recent")
        f.addSeparator()
        self._act(f, "Load protocol (JSON)...", self.load_protocol_dialog)
        self._act(f, "Save protocol (JSON)...", self.save_protocol_dialog)
        f.addSeparator()
        self.act_fig = self._act(f, "Export current &figure...", self.export_figure, QKeySequence("Ctrl+E"))
        self.act_csv = self._act(f, "Export multi-start &CSV...", self.export_csv)
        self.act_report = self._act(f, "Save text &report...", self.export_report)
        f.addSeparator()
        self._act(f, "&Quit", self.close, QKeySequence.Quit)
        r = mb.addMenu("&Run")
        self._act(r, "Fit once", self.run_fit, QKeySequence("F5"))
        self._act(r, "Box plot test", self.run_box, QKeySequence("F6"))
        self._act(r, "Monte Carlo", self.run_mc, QKeySequence("F7"))
        self._act(r, "Uniqueness", self.run_uniq, QKeySequence("F8"))
        self._act(r, "All three", self.run_all, QKeySequence("F9"))
        r.addSeparator()
        self.act_batch = self._act(r, "Box plot on every block with a Casa fit...", self.run_batch)
        self._act(r, "Cancel", self.cancel, QKeySequence("Esc"))
        v = mb.addMenu("&View")
        self._act(v, "Light theme", lambda: self.apply_theme("light"))
        self._act(v, "Dark theme", lambda: self.apply_theme("dark"))
        h = mb.addMenu("&Help")
        self._act(h, "About the method", self.about_method)
        self._act(h, "About fitstability", self.about)
        self._rebuild_recent()

    def _act(self, menu, text, fn, shortcut=None):
        a = QAction(text, self)
        a.triggered.connect(lambda checked=False: fn())
        if shortcut is not None:
            a.setShortcut(shortcut)
        menu.addAction(a)
        return a

    # ======================================================================================
    # theme, settings, recent files
    # ======================================================================================
    def apply_theme(self, name: str):
        self.theme_name = name
        theme.apply(QApplication.instance(), name)
        for p in self.panes:
            p.set_theme(name)
        self.settings.setValue("theme", name)
        if self.has_data():
            self.refresh_spectrum()
            if self.ms is not None:
                self._draw_box()
            if self.mc is not None:
                self._draw_mc()
            if self.uniq:
                self._draw_uniq()

    def _restore_settings(self):
        geo = self.settings.value("geometry")
        if geo is not None:
            self.restoreGeometry(geo)
        for key, w in (("n_starts", self.n_starts), ("mc_trials", self.mc_trials), ("u_points", self.u_points)):
            v = self.settings.value(key)
            if v is not None:
                w.setValue(int(v))

    def closeEvent(self, ev):
        if self.worker is not None and self.worker.isRunning():
            self.worker.request_cancel()
            self.worker.wait(3000)
        self.settings.setValue("geometry", self.saveGeometry())
        for key, w in (("n_starts", self.n_starts), ("mc_trials", self.mc_trials), ("u_points", self.u_points)):
            self.settings.setValue(key, w.value())
        super().closeEvent(ev)

    def _recent(self) -> list[str]:
        v = self.settings.value("recent", [])
        return [v] if isinstance(v, str) else list(v or [])

    def _remember(self, path: str):
        r = [p for p in self._recent() if p != path]
        self.settings.setValue("recent", [path] + r[:7])
        self.settings.setValue("last_dir", str(Path(path).parent))
        self._rebuild_recent()

    def _rebuild_recent(self):
        self.recent_menu.clear()
        recent = [p for p in self._recent() if Path(p).exists()]
        self.recent_menu.setEnabled(bool(recent))
        for p in recent:
            self._act(self.recent_menu, p, lambda p=p: self.open_file(p))

    # ======================================================================================
    # opening data
    # ======================================================================================
    def dragEnterEvent(self, ev):
        if ev.mimeData().hasUrls():
            ev.acceptProposedAction()

    def dropEvent(self, ev):
        urls = ev.mimeData().urls()
        if urls:
            self.open_file(urls[0].toLocalFile())

    def open_dialog(self):
        start = str(self.settings.value("last_dir", ""))
        path, _ = QFileDialog.getOpenFileName(self, "Open spectrum", start, OPEN_FILTER)
        if path:
            self.open_file(path)

    def open_file(self, path: str, csv_options: dict | None = None):
        """Load a VAMAS or text spectrum. ``csv_options`` (x, y, count_scale) skips the import dialog (tests)."""
        path = str(path)
        if not Path(path).exists():
            self._error(f"{path} does not exist")
            return
        try:
            if Path(path).suffix.lower() in VAMAS_SUFFIXES:
                self._open_vamas(path)
            else:
                self._open_table(path, csv_options)
        except (VamasError, ValueError, OSError) as e:
            self._error(f"Cannot open {Path(path).name}:\n{e}")

    def _open_vamas(self, path):
        vf = read_vamas(path)
        self._reset_results()
        self.path, self.vamas = path, vf
        self.lbl_file.setText(f"{Path(path).name}\n{vf.mode}, {len(vf.blocks)} blocks")
        self.blocks.blockSignals(True)
        self.blocks.setRowCount(len(vf.blocks))
        first_fit = None
        for r, b in enumerate(vf.blocks):
            fit_ = parse_fit(b.comments)
            has = bool(fit_ and fit_.comps)
            if has and first_fit is None:
                first_fit = r
            for c, text in enumerate([str(b.index), b.region, str(len(b.be)), f"{len(fit_.comps)} peaks" if has else "-"]):
                self.blocks.setItem(r, c, QTableWidgetItem(text))
        self.blocks.blockSignals(False)
        self.blocks.setEnabled(True)
        self.units.setEnabled(True)
        self.btn_reimport.setEnabled(False)
        self._remember(path)
        self.blocks.selectRow(first_fit if first_fit is not None else 0)
        self._block_selected()
        self._update_actions()

    def _open_table(self, path, options=None):
        if options is None:
            dlg = CsvImportDialog(path, self)
            if dlg.exec() != dlg.Accepted:
                return
            x, y, cs, desc = dlg.x, dlg.y, dlg.count_scale, dlg.description
        else:
            x, y, cs = options["x"], options["y"], options.get("count_scale", 1.0)
            desc = f"{len(x)} points"
        self._reset_results()
        self.path, self.vamas, self.block = path, None, None
        self.casa, self.casa_model, self.casa_bg = [], None, None
        self.csv_count_scale = cs
        self.base_x, self.base_y = np.asarray(x, float), np.asarray(y, float)
        self.lbl_file.setText(f"{Path(path).name}\n{desc}")
        self.blocks.setRowCount(0)
        self.blocks.setEnabled(False)
        self.units.setEnabled(False)
        self.btn_reimport.setEnabled(True)
        self._set_window(float(self.base_x.min()), float(self.base_x.max()))
        self._set_casa_ui([])
        self.display_start = None
        self.notes = []
        self.bg.setEnabled(True)
        self._remember(path)
        self._after_data_change("Model: add peaks in the table (or load a protocol JSON), then run.")

    def reimport_csv(self):
        if self.path and self.vamas is None:
            peaks = self.table.peaks()
            self._open_table(self.path)
            if peaks and not self.table.peaks():
                self.table.set_peaks(peaks)

    def _block_selected(self):
        if self.vamas is None:
            return
        rows = self.blocks.selectionModel().selectedRows()
        if not rows:
            return
        self._reset_results()
        self._load_block(self.vamas.blocks[rows[0].row()])

    def _units_name(self) -> str:
        return "counts" if self.units.currentIndex() == 0 else "cps"

    def _load_block(self, b):
        self.block = b
        cf = parse_fit(b.comments)
        calib = cf.calib_shift if cf else 0.0
        self.base_x = np.asarray(b.be, float) + calib
        self.base_y = np.asarray(b.y, float)
        self.casa, self.casa_model, self.casa_bg = [], None, None
        if cf and cf.comps:
            try:
                self.casa = casa_models(b, self._units_name())
            except (CasaFitError, ValueError) as e:
                self.casa = []
                self.notes = [str(e)]
            skips = describe_skips(b)
        else:
            skips = []
        self._set_casa_ui(self.casa)
        if self.casa:
            self._use_casa(0)
            self.notes += skips
        else:
            self._set_window(float(self.base_x.min()), float(self.base_x.max()))
            self.display_start = None
            self.notes = skips or (["This block has no CasaXPS fit: add peaks in the table or load a protocol."]
                                   if not cf else self.notes)
            self.bg.setEnabled(True)
            self._after_data_change()
        self._show_notes()

    def _set_casa_ui(self, models):
        self.region.blockSignals(True)
        self.region.clear()
        for m in models:
            self.region.addItem(m.region)
        self.region.blockSignals(False)
        self.region.setEnabled(bool(models))
        self.btn_reset.setEnabled(bool(models))

    def _region_selected(self):
        if self.casa and self.region.currentIndex() >= 0:
            self._use_casa(self.region.currentIndex())

    def _use_casa(self, i: int):
        if not (0 <= i < len(self.casa)):
            return
        m = self.casa[i]
        self._loading = True
        self.casa_model = m
        self.casa_bg = m.protocol.background
        self.table.set_peaks(m.protocol.peaks)
        self._set_window(*m.window)
        self.bg.setEnabled(False)
        self._loading = False
        self.display_start = m.stored
        self.notes = list(m.notes)
        self.notes.insert(0, f"CasaXPS saved fit: RSTD {m.rstd:.2f} ({m.protocol.background.kind} background, {len(m.protocol.peaks)} peaks)")
        self._show_notes()
        self._after_data_change()

    def _units_changed(self):
        if self.block is not None:
            self._reset_results()
            self._load_block(self.block)

    # ======================================================================================
    # model / window plumbing
    # ======================================================================================
    def has_data(self) -> bool:
        return self.base_x is not None

    def _set_window(self, lo, hi):
        for w, v in ((self.win_lo, lo), (self.win_hi, hi)):
            w.blockSignals(True)
            w.setValue(v)
            w.blockSignals(False)
        self.window = (lo, hi)

    def _window_changed(self):
        self.window = (self.win_lo.value(), self.win_hi.value())
        self._model_changed()

    def _table_problem(self, msg: str):
        if msg:
            self.status.showMessage(msg, 6000)

    def _model_changed(self, *args):
        """Peaks, background or window edited: the drawn model and old results no longer match."""
        if self._loading:
            return
        self.display_start = None
        msg = pm.validate(self.table.peaks(), "shirley")
        self.lbl_valid.setText(msg)
        for p in self.panes[1:]:
            p.set_stale(True)
        self._update_actions()
        if hasattr(self, "_refresh_timer"):
            self._refresh_timer.start()

    def _after_data_change(self, message: str = ""):
        msg = pm.validate(self.table.peaks(), "shirley")
        self.lbl_valid.setText(msg if self.table.peaks() else "")
        self._update_actions()
        self.refresh_spectrum()
        if message:
            self.status.showMessage(message)

    def _reset_results(self):
        self.ms = self.verdict = self.mc = None
        self.uniq, self.uniq_skipped = [], []
        self.out.clear()
        for p, msg in ((TAB_BOX, "Run the box plot test"), (TAB_MC, "Run Monte Carlo"), (TAB_UNIQ, "Run uniqueness")):
            self.panes[p].show_message(msg)
        self._update_actions()

    def _show_notes(self):
        self.lbl_notes.setText("\n".join(self.notes))

    def current_background(self):
        return self.casa_bg if self.casa_bg is not None else self.bg.currentText().lower()

    def count_scale(self) -> float:
        if self.casa_model is not None:
            return self.casa_model.count_scale
        if self.block is not None:
            return self.block.count_scale(self._units_name())
        return self.csv_count_scale

    def build_problem(self) -> FitProblem:
        """The fit problem for the current data, window and peak table. Raises ValueError with a readable reason."""
        if not self.has_data():
            raise ValueError("open a spectrum first")
        peaks = self.table.peaks()
        bg = self.current_background()
        msg = pm.validate(peaks, bg)
        if msg:
            raise ValueError(msg)
        lo, hi = self.window
        sel = (self.base_x >= lo) & (self.base_x <= hi)
        if sel.sum() < 10:
            raise ValueError("the fit window holds fewer than 10 points")
        x, y = self.base_x[sel], self.base_y[sel]
        if any(p.shape for p in peaks):
            dx = np.diff(x)
            if dx.size and dx.max() > 1.02 * dx.min():
                raise ValueError("Casa line shapes need an evenly spaced energy axis; use the pseudo-Voigt (empty shape)")
        return FitProblem(x, y, Protocol(peaks, bg), self.count_scale())

    def _start_for(self, problem):
        """Best single fit to start Monte Carlo/uniqueness from: Casa's saved values if untouched, else nominal."""
        if self.display_start is not None and len(self.display_start.pos) == len(problem.protocol.peaks):
            return fit(problem, problem.protocol.resolve(self.display_start))
        return fit_nominal(problem)

    # -- peak table buttons -----------------------------------------------------------------
    def add_peak(self):
        if not self.has_data():
            return
        lo, hi = self.window
        sel = (self.base_x >= lo) & (self.base_x <= hi)
        x, y = self.base_x[sel], self.base_y[sel]
        pos = float(x[int(np.argmax(y))]) if len(x) else (lo + hi) / 2
        self.table.add_peak(pos)

    def duplicate_peak(self):
        rows = self.view.selectionModel().selectedRows()
        if rows:
            self.table.duplicate(rows[0].row())

    def remove_peak(self):
        rows = self.view.selectionModel().selectedRows()
        if rows:
            self.table.remove(rows[0].row())

    def load_protocol_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, "Load protocol", str(self.settings.value("last_dir", "")),
                                              "Protocol JSON (*.json)")
        if path:
            self.load_protocol(path)

    def load_protocol(self, path: str):
        try:
            proto, cs = load_protocol(path)
        except (ValueError, OSError, KeyError) as e:
            self._error(f"Cannot load {Path(path).name}:\n{e}")
            return
        self.casa_model, self.casa_bg = None, None
        self.csv_count_scale = cs if cs is not None else self.csv_count_scale
        self._loading = True
        self.table.set_peaks(proto.peaks)
        if isinstance(proto.background, str):
            self.bg.setCurrentText(proto.background.capitalize())
        else:
            self.casa_bg = proto.background
        self.bg.setEnabled(self.casa_bg is None)
        self._loading = False
        self.display_start = None
        self.notes = [f"Protocol loaded from {Path(path).name}"]
        self._show_notes()
        self._after_data_change()

    def save_protocol_dialog(self):
        peaks = self.table.peaks()
        msg = pm.validate(peaks, self.current_background())
        if msg:
            self._error(msg)
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save protocol", str(self.settings.value("last_dir", "")),
                                              "Protocol JSON (*.json)")
        if path:
            save_protocol(path, Protocol(peaks, self.current_background()), self.count_scale())
            self.status.showMessage(f"Saved {path}", 5000)

    # ======================================================================================
    # drawing
    # ======================================================================================
    def refresh_spectrum(self):
        pane = self.panes[TAB_SPECTRUM]
        if not self.has_data():
            pane.show_message("Open a spectrum to begin")
            return
        lo, hi = self.window if self.window else (self.base_x.min(), self.base_x.max())
        try:
            prob = self.build_problem()
        except (ValueError, KeyError) as e:
            sel = (self.base_x >= lo) & (self.base_x <= hi)
            if sel.sum() >= 2:
                pane.draw(plots.plot_data, self.base_x[sel], self.base_y[sel],
                          markers=[(p.pos, p.name) for p in self.table.peaks()])
            else:
                pane.show_message(str(e))
            return
        ds = self.display_start
        if ds is not None and len(ds.pos) == len(prob.protocol.peaks):
            pane.draw(plots.plot_fit, prob, prob.protocol.resolve(ds))
        else:
            pane.draw(plots.plot_data, prob.x, prob.y, background=prob.background,
                      markers=[(p.pos, p.name) for p in self.table.peaks()])

    def _draw_box(self):
        self.panes[TAB_BOX].draw(plots.plot_multistart, self.ms, title=self._title, verdict=self.verdict)

    def _draw_mc(self):
        self.panes[TAB_MC].draw(plots.plot_monte_carlo, self.mc, title=self._title)

    def _draw_uniq(self):
        self.panes[TAB_UNIQ].draw(plots.plot_uniqueness, self.uniq, title=self._title + " (free parameters only)")

    # ======================================================================================
    # running experiments
    # ======================================================================================
    def _busy(self) -> bool:
        return self.worker is not None and self.worker.isRunning()

    def _update_actions(self):
        ready = self.has_data() and not self._busy() and not pm.validate(self.table.peaks(), self.current_background())
        for b in (self.btn_fit, self.btn_box, self.btn_mc, self.btn_uniq):
            b.setEnabled(ready)
        self.btn_cancel.setEnabled(self._busy())
        self.btn_open.setEnabled(not self._busy())
        self.act_fig.setEnabled(self.has_data())
        self.act_csv.setEnabled(self.ms is not None)
        self.act_report.setEnabled(self.ms is not None or self.mc is not None or bool(self.uniq))
        self.act_batch.setEnabled(self.vamas is not None and not self._busy())

    def _seed(self):
        return self.seed.value() if self.use_seed.isChecked() else None

    def _title_for(self) -> str:
        name = Path(self.path).name if self.path else ""
        if self.block is not None:
            name += f" {self.block.label}"
        if self.casa_model is not None and len(self.casa) > 1:
            name += f" / {self.casa_model.region}"
        return name.strip()

    def _problem_or_warn(self):
        try:
            return self.build_problem()
        except (ValueError, KeyError) as e:
            self._error(str(e))
            return None

    def _start(self, label: str, job, on_done):
        if self._busy():
            return
        self.worker = Worker(job, self)
        self.progress.setRange(0, 0)
        self.progress.setVisible(True)
        self.status.showMessage(f"{label}...")
        self.worker.progressed.connect(self._progressed)
        self.worker.succeeded.connect(lambda res: self._finished(on_done, res))
        self.worker.cancelled.connect(self._cancelled)
        self.worker.failed.connect(self._failed)
        self.worker.finished.connect(self._worker_done)
        self._update_actions()
        self.worker.start()

    def _progressed(self, done, total):
        if total > 0:
            self.progress.setRange(0, total)
            self.progress.setValue(done)

    def _finished(self, on_done, res):
        self.status.showMessage("Done", 5000)
        try:
            on_done(res)
        except Exception as e:  # noqa: BLE001 - a bug in result handling must not kill the window
            self._error(f"{type(e).__name__}: {e}")
            self._queue.clear()

    def _cancelled(self):
        self._queue.clear()
        self.status.showMessage("Cancelled", 4000)

    def _failed(self, msg):
        self._queue.clear()
        self._error(msg)

    def _worker_done(self):
        self.progress.setVisible(False)
        self._update_actions()
        if self._queue:
            self._queue.pop(0)()

    def cancel(self):
        if self._busy():
            self.worker.request_cancel()
            self.status.showMessage("Cancelling...")

    def run_fit(self):
        prob = self._problem_or_warn()
        if prob is None:
            return
        start_ds = self.display_start

        def job(progress, cancel):
            if start_ds is not None and len(start_ds.pos) == len(prob.protocol.peaks):
                return fit(prob, prob.protocol.resolve(start_ds))
            return fit_nominal(prob)

        def done(r):
            self.display_start = r.as_start()
            self.refresh_spectrum()
            names = [p.name for p in prob.protocol.peaks]
            lines = [f"Single fit: RSTD {r.rstd:.3f}, reduced chi2 {r.red_chi2:.3f}, converged: {r.success}",
                     f"{'peak':<24}{'position':>10}{'FWHM':>8}{'area':>12}"]
            lines += [f"{n[:23]:<24}{r.pos[k]:>10.3f}{r.fwhm[k]:>8.3f}{r.area[k]:>12.4g}" for k, n in enumerate(names)]
            self.out.setPlainText("\n".join(lines))
            self.tabs.setCurrentIndex(TAB_SPECTRUM)
        self._start("Fitting", job, done)

    def run_box(self):
        prob = self._problem_or_warn()
        if prob is None:
            return
        n, seed, bf = self.n_starts.value(), self._seed(), self.bounds_fraction.value()
        rp, rf, jobs = self.rand_pos.value(), self.rand_fwhm.value(), self.jobs.value()
        self._title = self._title_for()
        casa_rstd = self.casa_model.rstd if self.casa_model is not None and self.display_start is not None else None

        def job(progress, cancel):
            return multistart(prob, n, seed, rp, rf, jobs, bf, progress, cancel)

        def done(ms):
            self.ms = ms
            self.verdict = assess(ms)
            self._draw_box()
            self.panes[TAB_BOX].set_stale(False)
            self.out.setPlainText(results.multistart_report(self._title, ms, self.verdict, casa_rstd, self.notes))
            self.tabs.setCurrentIndex(TAB_BOX)
            self.status.showMessage(f"Box plot test: {self.verdict.label}", 8000)
            self._update_actions()
        self._start("Fitting from random starts", job, done)

    def run_mc(self):
        prob = self._problem_or_warn()
        if prob is None:
            return
        n, seed, jobs = self.mc_trials.value(), self._seed(), self.jobs.value()
        self._title = self._title_for()

        def job(progress, cancel):
            return monte_carlo(prob, n, seed, self._start_for(prob), jobs, progress, cancel)

        def done(mc):
            self.mc = mc
            self._draw_mc()
            self.out.setPlainText(results.monte_carlo_report(self._title, mc))
            self.tabs.setCurrentIndex(TAB_MC)
            self._update_actions()
        self._start("Monte Carlo", job, done)

    def run_uniq(self):
        prob = self._problem_or_warn()
        if prob is None:
            return
        which = {0: "pos", 1: "fwhm", 2: None}[self.u_param.currentIndex()]
        hr, npts = self.u_range.value(), self.u_points.value()
        proto = prob.protocol
        targets = scannable_parameters(proto, which)
        skipped = [r for r in (unscannable_reason(proto, p.name, k) for p in proto.peaks
                               for k in (("pos", "fwhm") if which is None else (which,))) if r]
        if not targets:
            self.out.setPlainText(results.uniqueness_report(self._title_for(), [], skipped))
            self.status.showMessage("No free parameters of that kind to scan", 6000)
            return
        self._title = self._title_for()

        def job(progress, cancel):
            best = self._start_for(prob)
            total, out = len(targets) * npts, []
            for i, (name, prm) in enumerate(targets):
                base = i * npts
                out.append(uniqueness(prob, name, prm, hr, npts, best,
                                      progress=lambda d, t, base=base: progress(base + d, total), cancel=cancel))
            return out

        def done(res):
            self.uniq, self.uniq_skipped = res, skipped
            self._draw_uniq()
            self.panes[TAB_UNIQ].set_stale(False)
            self.out.setPlainText(results.uniqueness_report(self._title, res, skipped))
            self.tabs.setCurrentIndex(TAB_UNIQ)
            self._update_actions()
        self._start("Uniqueness scans", job, done)

    def run_all(self):
        if self._busy():
            return
        self._queue = [self.run_mc, self.run_uniq]
        self.run_box()

    # -- batch over a VAMAS file ------------------------------------------------------------
    def run_batch(self):
        if self.vamas is None or self._busy():
            return
        folder = QFileDialog.getExistingDirectory(self, "Folder for the figures (Cancel = no figures)",
                                                  str(Path(self.path).parent))
        n, seed, bf = self.n_starts.value(), self._seed(), self.bounds_fraction.value()
        rp, rf = self.rand_pos.value(), self.rand_fwhm.value()
        units, vf, stem = self._units_name(), self.vamas, Path(self.path).stem
        blocks = [b for b in vf.blocks if (f := parse_fit(b.comments)) and f.comps]
        if not blocks:
            self._error("No block in this file has a CasaXPS fit")
            return

        def job(progress, cancel):
            from matplotlib.figure import Figure
            rows, skipped = [], []
            for i, b in enumerate(blocks):
                if cancel():
                    from ..experiments import Cancelled
                    raise Cancelled()
                try:
                    models = casa_models(b, units)
                except (CasaFitError, ValueError) as e:
                    skipped.append(f"{b.label}: {e}")
                    continue
                skipped += describe_skips(b)
                for m in models:
                    ms = multistart(m.problem(), n, seed, rp, rf, 1, bf, None, cancel)
                    v = assess(ms)
                    label = f"{b.label} / {m.region}" if len(models) > 1 else b.label
                    rows.append(results.summary_row(label, ms, v, m.rstd))
                    if folder:
                        fig = Figure(figsize=(11, 7))
                        plots.plot_multistart(ms, title=f"{stem} {label}", verdict=v, fig=fig)
                        safe = "".join(ch if ch.isalnum() else "_" for ch in f"{stem}_{b.index}_{m.region}")
                        fig.savefig(str(Path(folder) / f"{safe}.boxplot.png"), dpi=130)
                progress(i + 1, len(blocks))
            return rows, skipped, folder

        def done(res):
            rows, skipped, fld = res
            BatchResultsDialog(rows, skipped, fld, self).exec()
        self._start(f"Batch over {len(blocks)} blocks", job, done)

    # ======================================================================================
    # exports and dialogs
    # ======================================================================================
    def export_figure(self):
        pane = self.panes[self.tabs.currentIndex()]
        path, _ = QFileDialog.getSaveFileName(self, "Export figure", str(self.settings.value("last_dir", "")),
                                              "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
        if path:
            pane.fig.savefig(path, dpi=200)
            self.status.showMessage(f"Saved {path}", 5000)

    def export_csv(self):
        if self.ms is None:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export multi-start results", "multistart.csv", "CSV (*.csv)")
        if path:
            results.multistart_csv(self.ms, path)
            self.status.showMessage(f"Saved {path}", 5000)

    def export_report(self):
        text = self.out.toPlainText()
        parts = []
        if self.ms is not None:
            parts.append(results.multistart_report(self._title, self.ms, self.verdict, None, self.notes))
        if self.mc is not None:
            parts.append(results.monte_carlo_report(self._title, self.mc))
        if self.uniq:
            parts.append(results.uniqueness_report(self._title, self.uniq, self.uniq_skipped))
        path, _ = QFileDialog.getSaveFileName(self, "Save report", "fitstability-report.txt", "Text (*.txt)")
        if path:
            Path(path).write_text("\n\n".join(parts) or text, encoding="utf-8")
            self.status.showMessage(f"Saved {path}", 5000)

    def about_method(self):
        QMessageBox.information(
            self, "About the method",
            "Box plot test (Moeini et al., J. Electron Spectrosc. Relat. Phenom. 250 (2021) 147094): fit the same "
            "spectrum from many random starting conditions. Narrow boxes of the final peak areas mean one minimum; "
            "wide boxes mean multiple local minima, i.e. an over-parameterised fit.\n\n"
            "A narrow box shows only that the minimum is unique, not that the fit is good: read it with the residual.\n\n"
            + UNIQUENESS_METHOD_NOTE)

    def about(self):
        QMessageBox.about(self, "fitstability", f"fitstability {__version__}\nPeak-fit stability diagnostics for XPS.")

    def _error(self, msg: str):
        self.status.showMessage(msg.splitlines()[0], 8000)
        if not os.environ.get("FITSTABILITY_NO_DIALOGS"):
            QMessageBox.warning(self, "fitstability", msg)
