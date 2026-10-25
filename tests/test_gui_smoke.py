"""Offscreen smoke tests of the main window: open files, run the experiments, edit the table."""
import os
import time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("FITSTABILITY_NO_DIALOGS", "1")

import numpy as np
import pytest

pytest.importorskip("PySide6")
from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import QApplication

from fitstability import PeakSpec, Protocol, save_protocol, synthetic_spectrum, write_vamas
from fitstability.gui.main_window import TAB_BOX, TAB_MC, TAB_UNIQ, MainWindow

from test_casa import _make_block   # a synthetic Casa-fitted VAMAS block with a tied doublet


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture()
def win(app, tmp_path):
    w = MainWindow("light", QSettings(str(tmp_path / "s.ini"), QSettings.IniFormat))
    w.n_starts.setValue(6)
    w.mc_trials.setValue(8)
    w.u_points.setValue(5)
    w.show()
    yield w
    w.close()


def wait(app, win, timeout=120):
    t0 = time.time()
    while win._busy() or win._queue:
        app.processEvents()
        time.sleep(0.01)
        assert time.time() - t0 < timeout, "worker did not finish"
    app.processEvents()


def casa_file(tmp_path):
    f = tmp_path / "casa.vms"
    write_vamas(f, [_make_block(0.0, "GL(30)")])
    return str(f)


def test_open_vamas_loads_casa_fit_into_table(app, win, tmp_path):
    win.open_file(casa_file(tmp_path))
    names = [p.name for p in win.table.peaks()]
    assert names == ["main", "partner", "free"]
    assert win.blocks.rowCount() == 1 and win.blocks.item(0, 3).text() == "3 peaks"
    assert win.region.count() == 1 and win.casa_model is not None
    assert win.btn_box.isEnabled()
    assert win.panes[0].fig.axes                       # spectrum drawn with the Casa model
    assert "tied" not in win.lbl_valid.text()
    # the table shows Casa's constraints as rules
    assert win.table.index(1, 3).data() == "main + 5"
    assert win.table.index(1, 5).data() == "main * 1.2"
    assert win.table.index(1, 6).data() == "main * 0.5"


def test_box_plot_run_populates_tab_and_verdict(app, win, tmp_path):
    win.open_file(casa_file(tmp_path))
    win.run_box()
    wait(app, win)
    assert win.ms is not None and win.verdict.label == "single minimum"
    assert win.panes[TAB_BOX].has_result and win.panes[TAB_BOX].fig.axes
    assert "single minimum" in win.out.toPlainText()
    assert win.tabs.currentIndex() == TAB_BOX
    assert win.act_csv.isEnabled()
    out = tmp_path / "ms.csv"
    from fitstability.gui import results
    results.multistart_csv(win.ms, out)
    assert len(out.read_text().splitlines()) == 7


def test_all_three_run_in_sequence_and_uniqueness_skips_ties(app, win, tmp_path):
    win.open_file(casa_file(tmp_path))
    win.u_param.setCurrentIndex(2)                      # both position and FWHM
    win.run_all()
    wait(app, win)
    assert win.ms is not None and win.mc is not None and win.uniq
    scanned = {(u.peak, u.parameter) for u in win.uniq}
    assert ("main", "pos") in scanned and ("partner", "pos") not in scanned and ("free", "fwhm") not in scanned
    assert any("partner.pos is tied" in s for s in win.uniq_skipped)
    assert "not scanned" in win.out.toPlainText()
    assert win.panes[TAB_MC].has_result and win.panes[TAB_UNIQ].has_result


def test_editing_table_marks_results_stale_and_validates(app, win, tmp_path):
    win.open_file(casa_file(tmp_path))
    win.run_box()
    wait(app, win)
    assert not win.panes[TAB_BOX].banner.isVisible()
    idx = win.table.index(0, 2)
    assert win.table.setData(idx, "531.3")
    assert win.panes[TAB_BOX].banner.isVisible()       # stale banner
    assert win.display_start is None
    # a bad rule is rejected and the table is unchanged
    assert not win.table.setData(win.table.index(1, 3), "Nope + 1")
    assert win.table.index(1, 3).data() == "main + 5"
    # a cycle blocks running
    assert win.table.setData(win.table.index(0, 6), "partner * 2")
    assert "circular" in win.lbl_valid.text() and not win.btn_box.isEnabled()


def test_csv_open_with_manual_peaks_and_protocol_json(app, win, tmp_path):
    x, y = synthetic_spectrum([(285.0, 1.2, 12000), (286.5, 1.3, 9000)], counts_scale=0.2, seed=1)
    csvf = tmp_path / "c.csv"
    csvf.write_text("BE,counts\n" + "\n".join(f"{a},{b}" for a, b in zip(x, y)))
    win.open_file(str(csvf), csv_options=dict(x=x, y=y, count_scale=0.2))
    assert win.table.rowCount() == 0 and not win.btn_box.isEnabled()
    win.add_peak()
    win.add_peak()
    assert win.table.rowCount() == 2
    win.table.setData(win.table.index(0, 2), "285.0")
    win.table.setData(win.table.index(1, 2), "286.5")
    assert win.btn_box.isEnabled()
    win.run_fit()
    wait(app, win)
    assert win.display_start is not None and "Single fit" in win.out.toPlainText()
    # save and reload the protocol
    pj = tmp_path / "p.json"
    save_protocol(pj, Protocol(win.table.peaks(), "shirley"), 0.2)
    win.table.set_peaks([])
    win.load_protocol(str(pj))
    assert win.table.rowCount() == 2 and win.count_scale() == pytest.approx(0.2)
    win.run_box()
    wait(app, win)
    assert win.ms is not None and win.verdict.label in ("single minimum", "mostly single minimum", "multiple minima")


def test_cancel_stops_a_long_run(app, win, tmp_path):
    win.open_file(casa_file(tmp_path))
    win.n_starts.setValue(2000)
    win.run_box()
    app.processEvents()
    assert win._busy()
    win.cancel()
    wait(app, win)
    assert win.ms is None and win.btn_box.isEnabled()


def test_real_casa_file_if_present(app, win):
    f = ("C:/Users/david/Documents/Brazil Training Course September 2025/September 2nd - Advanced Chemical State Analysis/"
         "Data/Titanium Metal Depth Profile - INSTRUCTORS.vms")
    if not os.path.exists(f):
        pytest.skip("real example not present")
    win.open_file(f)
    assert win.blocks.rowCount() == 24
    assert win.blocks.selectionModel().selectedRows()[0].row() >= 0
    assert win.table.rowCount() >= 2 and win.casa_model is not None
    win.n_starts.setValue(8)
    win.bounds_fraction.setValue(0.5)
    win.run_box()
    wait(app, win)
    assert win.ms is not None and win.panes[TAB_BOX].has_result


def test_batch_and_theme_switch(app, win, tmp_path, monkeypatch):
    win.open_file(casa_file(tmp_path))
    win.apply_theme("dark")
    assert win.theme_name == "dark" and win.panes[0].theme_name == "dark"
    win.apply_theme("light")
    from PySide6.QtWidgets import QFileDialog
    monkeypatch.setattr(QFileDialog, "getExistingDirectory", staticmethod(lambda *a, **k: str(tmp_path)))
    shown = {}
    import fitstability.gui.main_window as mw
    monkeypatch.setattr(mw.BatchResultsDialog, "exec", lambda self: shown.setdefault("rows", self.rows))
    win.run_batch()
    wait(app, win)
    assert len(shown["rows"]) == 1
    assert list(tmp_path.glob("*.boxplot.png"))


# -- CasaXPS library-spectrum components ---------------------------------------------------------------
def library_setup(tmp_path):
    from test_library import make_block, make_reference
    blk, _, _ = make_block()
    f = tmp_path / "ce.vms"
    write_vamas(f, [blk])
    refs = tmp_path / "refs"
    refs.mkdir()
    for n, sh in (("LS-CE4", 0.0), ("LS-CE3", -3.0)):
        r = make_reference(n, sh)
        np.savetxt(refs / f"{n}.csv", np.c_[r.be, r.y], delimiter=",", header="be,counts", comments="")
    return str(f), str(refs)


def test_library_spectra_missing_then_supplied_then_calibrated(app, win, tmp_path):
    f, refs = library_setup(tmp_path)
    win.open_file(f)
    assert win.table.rowCount() == 0 and "not loaded" in win.lbl_notes.text() and "LS-CE4" in win.lbl_notes.text()
    win._add_reference_paths([refs])                      # the dialog buttons call this
    assert len(win.library) == 2 and win.table.rowCount() == 2
    assert any("not calibrated" in n for n in win.notes)
    assert win.act_calibrate.isEnabled()
    win.calibrate_references()
    wait(app, win)
    assert win.library.get("LS-CE4").origin == pytest.approx(899.1, abs=0.03)
    assert "origin (eV)" in win.out.toPlainText()
    assert win.casa_model is not None and win.casa_model.reproduced and win.casa_model.rel_rms < 0.01
    win.run_box()
    wait(app, win)
    assert win.verdict is not None and win.panes[TAB_BOX].has_result


def test_library_peak_in_table_resolves_by_name_and_survives_protocol_json(app, win, tmp_path):
    f, refs = library_setup(tmp_path)
    win.open_file(f)
    win._add_reference_paths([refs])
    assert not win.table.setData(win.table.index(0, 1), "%&<>LS-NOPE")                 # unknown reference refused
    assert win.table.setData(win.table.index(0, 1), "LS-CE3")                          # short spelling resolves
    assert win.table.peaks()[0].shape == "%&<>LS-CE3" and win.table.peaks()[0].reference is not None
    pj = tmp_path / "p.json"
    from fitstability import Protocol
    save_protocol(pj, Protocol(win.table.peaks(), win.current_background()), None)
    win.table.set_peaks([])
    win.load_protocol(str(pj))
    assert win.table.rowCount() == 2 and win.table.peaks()[0].reference is not None
