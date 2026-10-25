"""CasaXPS library-spectrum components: reference spectra supplied by the user."""
import json
import pickle

import numpy as np
import pytest

from fitstability import Block, PeakSpec, Protocol, multistart, write_vamas
from fitstability.casa import casa_models, describe_skips
from fitstability.library import ReferenceLibrary, ReferenceSpectrum, calibrate_origins, reference_key

from test_casa import _comp_line

HV = 1486.68
K = 0.1 * 5


def make_reference(name, shift=0.0):
    """A multiplet-like background-free reference: a main line, a shoulder and a satellite."""
    be = np.arange(870.0, 925.0, 0.05)
    y = (1.0 * np.exp(-0.5 * ((be - (882.5 + shift)) / 0.7) ** 2) + 0.45 * np.exp(-0.5 * ((be - (888.8 + shift)) / 1.2) ** 2)
         + 0.7 * np.exp(-0.5 * ((be - (898.3 + shift)) / 0.9) ** 2))
    return ReferenceSpectrum(name, be, y)


def make_block(origin_true=899.1, pos4=897.6, pos3=886.0, areas=(900.0, 300.0), names=("LS-CE4", "LS-CE3")):
    """A Casa fit with two library components (stored in kinetic energy), data generated from the references."""
    ref4, ref3 = make_reference(names[0]), make_reference(names[1], shift=-3.0)
    ref4, ref3 = ref4.with_origin(origin_true), ref3.with_origin(ref3.max_energy + 0.4)
    x_raw = np.linspace(870.0, 925.0, 551)                     # raw binding energy
    cps = 30.0 + 0.6 * (x_raw - 870.0) + ref4.curve(x_raw, pos4, areas[0]) + ref3.curve(x_raw, pos3, areas[1])
    lo_be, hi_be = 872.0, 922.0
    comments = [f"CASA region (*Ce 3d*) (*Linear*) {HV - hi_be:.4f} {HV - lo_be:.4f} 1 1 0 0 0 0 0 0 (*Ce 3d*) 140 0 1",
                _comp_line("Ce 4+", f"%&<>{names[0]}", areas[0], 1, HV - pos4, fwhm_b=(1, 1), pos_b=(HV - pos4 - 3, HV - pos4 + 3)),
                _comp_line("Ce 3+", f"%&<>{names[1]}", areas[1], 1, HV - pos3, fwhm_b=(1, 1), pos_b=(HV - pos3 - 3, HV - pos3 + 3))]
    blk = Block(index=1, block_id="Ce 3d/2", sample_id="s", species="Ce", transition="3d", hv_ev=HV, be=x_raw,
                y=cps * K, dwell_s=0.1, scans=5, pass_energy_ev=40.0, comments=comments)
    return blk, ref4, ref3


def test_reference_key_matches_casa_spellings():
    keys = {reference_key(n) for n in ("%&<>LS-CE4", "LS-CE4", "ce4", "LS_CE4", "$LS21-CE4")}
    assert keys == {"CE4"}
    assert reference_key("%&<>LS-CE(III)") == "CEIII" and reference_key("LS-CE3") != reference_key("LS-CE4")


def test_reference_curve_is_unit_area_and_positioned_at_origin():
    r = make_reference("X").with_origin(898.3)
    x = np.linspace(850, 950, 20001)
    c = r.curve(x, 900.0, 1234.0)
    assert np.trapezoid(c, x) == pytest.approx(1234.0, rel=1e-3)
    window = np.abs(x - 900.0) < 0.3
    assert x[window][np.argmax(c[window])] == pytest.approx(900.0, abs=0.06)       # the satellite maximum sits at pos
    assert np.all(r.curve(np.array([100.0, 2000.0]), 900.0, 5.0) == 0)             # zero outside its range
    with pytest.raises(ValueError):
        ReferenceSpectrum("bad", [1, 2, 3], [1, 1, 1])


def test_library_loads_vamas_blocks_and_text_files(tmp_path):
    r4, r3 = make_reference("LS-CE4"), make_reference("LS-CE3", -3.0)
    write_vamas(tmp_path / "refs.vms", [Block(index=1, block_id="LS-CE4", species="Ce", transition="3d", hv_ev=HV, be=r4.be, y=r4.y),
                                        Block(index=2, block_id="LS-CE3", species="Ce", transition="3d", hv_ev=HV, be=r3.be, y=r3.y)])
    np.savetxt(tmp_path / "LS-SB3D.csv", np.c_[r4.be, r4.y], delimiter=",", header="be,counts", comments="")
    (tmp_path / "origins.json").write_text(json.dumps({"LS-CE4": 899.1}))
    lib = ReferenceLibrary.load(tmp_path)
    assert lib.names() == ["LS-CE3", "LS-CE4", "LS-SB3D"]
    assert lib.get("%&<>LS-CE4").origin == 899.1 and lib.get("LS-CE3").origin is None
    assert lib.missing(["%&<>LS-CE4", "%&<>LS-ZZ", "GL(30)"]) == ["LS-ZZ"]
    assert np.allclose(lib.get("ce4").be, r4.be, atol=1e-6)
    with pytest.raises(FileNotFoundError):
        ReferenceLibrary.load(tmp_path / "nope")


def test_kinetic_energy_text_reference_needs_hv(tmp_path):
    r = make_reference("LS-CE4")
    np.savetxt(tmp_path / "LS-CE4.txt", np.c_[HV - r.be, r.y])
    with pytest.raises(ValueError, match="photon energy"):
        ReferenceLibrary.load(tmp_path / "LS-CE4.txt", axis="kinetic")
    lib = ReferenceLibrary.load(tmp_path / "LS-CE4.txt", axis="kinetic", hv=HV)
    assert np.allclose(lib.get("LS-CE4").be, r.be, atol=1e-6)


def test_library_component_is_skipped_with_names_until_references_are_supplied():
    blk, r4, r3 = make_block()
    assert casa_models(blk) == []
    msg = " ".join(describe_skips(blk))
    assert "LS-CE4" in msg and "LS-CE3" in msg and "not loaded" in msg
    assert casa_models(blk, library=ReferenceLibrary([r4])) == []
    msg = " ".join(describe_skips(blk))
    assert "LS-CE3" in msg and "LS-CE4" not in msg


def test_model_with_library_reproduces_saved_fit_and_runs_multistart():
    blk, r4, r3 = make_block()
    (m,) = casa_models(blk, library=ReferenceLibrary([r4, r3]))
    assert m.reproduced and m.rel_rms < 0.01, m.notes
    p = {q.name: q for q in m.protocol.peaks}
    assert p["Ce 4+"].reference is r4 and p["Ce 4+"].fix_fwhm and p["Ce 4+"].shape == "%&<>LS-CE4"
    assert p["Ce 4+"].pos == pytest.approx(897.6, abs=1e-6)
    ms = multistart(m.problem(), n_starts=4, seed=0)
    assert ms.best().rstd < 3
    assert pickle.loads(pickle.dumps(m.protocol)).peaks[0].reference.name == "LS-CE4"   # survives parallel runs


def test_protocol_without_reference_object_is_rejected_with_a_reason():
    with pytest.raises(ValueError, match="not loaded"):
        Protocol([PeakSpec("a", 900.0, 1.0, shape="%&<>LS-CE4")])


def test_uncalibrated_origin_is_flagged_and_calibration_recovers_the_truth(tmp_path):
    origin_true = 899.1
    blk, r4, r3 = make_block(origin_true=origin_true)
    naive = ReferenceLibrary([make_reference("LS-CE4"), make_reference("LS-CE3", -3.0)])      # origin unknown: assumed at the maximum
    (m0,) = casa_models(blk, library=naive)
    assert any("not calibrated" in n for n in m0.notes) and m0.rel_rms > 0.02
    fits = calibrate_origins(naive, [blk])
    assert set(fits) == {"LS-CE4", "LS-CE3"}
    assert fits["LS-CE4"].median == pytest.approx(origin_true, abs=0.03)
    assert fits["LS-CE3"].median == pytest.approx(make_reference("LS-CE3", -3.0).max_energy + 0.4, abs=0.03)
    assert fits["LS-CE4"].rel_after < 0.5 * fits["LS-CE4"].rel_before and fits["LS-CE4"].rel_after < 0.01
    naive.set_origins({n: f.median for n, f in fits.items()})
    (m1,) = casa_models(blk, library=naive)
    assert m1.rel_rms < 0.01 and not any("not calibrated" in n for n in m1.notes)
    naive.save_origins(tmp_path / "origins.json")
    assert set(json.loads((tmp_path / "origins.json").read_text())) == {"LS-CE4", "LS-CE3"}
