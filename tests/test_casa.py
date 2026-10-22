"""CasaXPS fit reader: frames, units, ties, shapes, backgrounds."""
import numpy as np
import pytest

from fitstability import (Block, assess, fit_nominal, multistart, read_vamas, scannable_parameters,
                          unscannable_reason, write_vamas)
from fitstability import casa_lineshapes as cl
from fitstability.casa import CasaFitError, casa_models, describe_skips, parse_fit

HV = 1486.68
DWELL, SCANS = 0.1, 5
K = DWELL * SCANS


def _comp_line(name, shape, area, fwhm, pos_ke, area_c=(-1, 1), fwhm_c=(-1, 1), pos_c=(-1, 1), state=-1,
               area_b=(1e-20, 1e7), fwhm_b=(0.3, 4.0), pos_b=None):
    pos_b = pos_b or (pos_ke - 3, pos_ke + 3)
    return (f"CASA comp (*{name}*) (*{shape}*) Area {area} {area_b[0]} {area_b[1]} {area_c[0]} {area_c[1]} "
            f"MFWHM {fwhm} {fwhm_b[0]} {fwhm_b[1]} {fwhm_c[0]} {fwhm_c[1]} "
            f"Position {pos_ke} {pos_b[0]} {pos_b[1]} {pos_c[0]} {pos_c[1]} RSF 1 MASS 1 INDEX {state} (*{name}*) CONST (**)")


def _make_block(calib=0.0, shape="LA(1.1,1.9,7)", background="Linear"):
    """A block holding a synthetic CasaXPS doublet fit, generated in the kinetic-energy frame the way
    CasaXPS stores it (so the test does not share the reader's own frame conversions).

    Stored positions are calibrated-frame KEs: BE_calibrated = hv - stored_KE - comp_shift + calib, and
    with no Regions/Comps flags comp_shift = calib, so BE_calibrated = hv - stored_KE."""
    be_axis = np.linspace(525.0, 545.0, 401)                      # raw binding energy, ascending
    ke_asc = (HV - be_axis)[::-1]
    pos0 = HV - 531.0                                             # main peak at calibrated BE 531.0
    c0 = dict(area=900.0, fwhm=1.3, pos_ke=pos0)                  # area in cps*eV
    c1 = dict(area=0.5 * c0["area"], fwhm=1.2 * c0["fwhm"], pos_ke=pos0 - 5.0)   # tied partner, 5 eV higher in BE
    free = dict(area=300.0, fwhm=1.8, pos_ke=HV - 536.5)
    # on the raw axis the shift is undone: raw KE = stored KE + calib
    comps_cps = sum(cl.component_curve(ke_asc, shape, c["pos_ke"] + calib, c["fwhm"], c["area"]) for c in (c0, c1, free))
    bg = np.linspace(40.0, 90.0, len(ke_asc))                     # cps, higher at the high-BE (low-KE) end
    cps_asc = comps_cps + bg[::-1]
    counts = (cps_asc * K)[::-1]                                  # back to ascending-BE order
    lo_be, hi_be = 526.0, 544.0                                   # window in calibrated BE, ends well clear of the peaks
    comments = []
    if calib:
        comments.append(f"Calib M = 100 A = {100 + calib} BE ADD")
    comments.append(f"CASA region (*O 1s*) (*{background}*) {HV - hi_be:.4f} {HV - lo_be:.4f} 0.5 1 0 0 0 0 0 0 (*O 1s*) 16 0 0.5")
    comments.append(_comp_line("main", shape, c0["area"], c0["fwhm"], c0["pos_ke"], state=0))
    comments.append(_comp_line("partner", shape, c1["area"], c1["fwhm"], c1["pos_ke"], area_c=(0, 0.5), fwhm_c=(0, 1.2),
                               pos_c=(0, -5.0), state=0, pos_b=(0, 0)))
    comments.append(_comp_line("free", shape, free["area"], free["fwhm"], free["pos_ke"], fwhm_b=(1.8, 1.8)))
    return Block(index=1, block_id="O1s Scan", species="O", transition="1s", hv_ev=HV, be=be_axis, y=counts,
                 dwell_s=DWELL, scans=SCANS, pass_energy_ev=50.0, comments=comments)


def _model(tmp_path, calib=0.0, shape="GL(30)"):
    f = tmp_path / "c.vms"
    write_vamas(f, [_make_block(calib, shape)])
    (m,) = casa_models(read_vamas(f).blocks[0])
    return m


@pytest.mark.parametrize("calib", [0.0, 1.25])
@pytest.mark.parametrize("shape", ["GL(30)", "LA(1.1,1.9,7)"])
def test_stored_fit_reproduces_data(tmp_path, calib, shape):
    m = _model(tmp_path, calib, shape)
    assert m.rel_rms < 0.01, m.notes
    assert not [n for n in m.notes if "does not reproduce" in n]
    assert m.window == pytest.approx((526.0, 544.0), abs=0.1)        # the calibrated BE window CasaXPS shows


def test_ties_bounds_and_frames(tmp_path):
    m = _model(tmp_path)
    p = {q.name: q for q in m.protocol.peaks}
    assert p["main"].pos == pytest.approx(531.0, abs=1e-6)                       # BE = hv - KE
    assert p["main"].pos_bounds == pytest.approx((528.0, 534.0), abs=1e-6)       # KE window flipped
    assert p["partner"].area_tie == ("main", 0.5) and p["partner"].fwhm_tie == ("main", 1.2)
    assert p["partner"].pos_tie == ("main", 5.0)                                 # -(KE offset): higher in BE
    assert p["free"].fix_fwhm and not p["free"].fix_pos                          # lower == upper -> fixed
    assert m.stored.area[0] == pytest.approx(900.0 * K)                          # cps x dwell x scans
    s = m.protocol.nominal_start()
    assert s.pos[1] == pytest.approx(536.0, abs=1e-6) and s.fwhm[1] == pytest.approx(1.56)


def test_uniqueness_skips_casa_ties(tmp_path):
    m = _model(tmp_path)
    got = scannable_parameters(m.protocol)
    assert ("main", "pos") in got and ("main", "fwhm") in got and ("free", "pos") in got
    assert ("partner", "pos") not in got and ("partner", "fwhm") not in got and ("free", "fwhm") not in got
    assert "tied to main" in unscannable_reason(m.protocol, "partner", "pos")
    assert "tied to main" in unscannable_reason(m.protocol, "partner", "fwhm")
    assert "fixed" in unscannable_reason(m.protocol, "free", "fwhm")


def test_multistart_on_casa_fit_keeps_ties(tmp_path):
    m = _model(tmp_path)
    ms = multistart(m.problem(), n_starts=6, seed=0)
    a = ms.final("area")
    assert np.allclose(a[:, 1], 0.5 * a[:, 0])                                    # tie holds in every fit
    assert np.allclose(ms.final("pos")[:, 1] - ms.final("pos")[:, 0], 5.0)
    assert np.allclose(ms.start_areas()[:, 1], 0.5 * ms.start_areas()[:, 0])      # and in the displayed starts
    assert assess(ms).label == "single minimum"


def test_bounds_fraction_draws_inside_bounds(tmp_path):
    m = _model(tmp_path)
    ms = multistart(m.problem(), n_starts=8, seed=1, bounds_fraction=0.5)
    lo, hi = m.protocol.bounds()
    for s in ms.starts:
        th = m.protocol.pack(s)
        assert np.all(th >= lo - 1e-9) and np.all(th <= hi + 1e-9)
    assert np.ptp([s.pos[0] for s in ms.starts]) > 0.1                            # positions really were randomised


def test_fit_from_stored_start_is_no_worse(tmp_path):
    m = _model(tmp_path, shape="LA(1.1,1.9,7)")
    r = fit_nominal(m.problem())
    assert r.rstd <= m.rstd * 1.05 + 1e-9


def test_la_orientation_matches_casa_frame():
    """An LA tail built by CasaXPS on the kinetic-energy axis must come out on the high-BE side."""
    from fitstability.lineshapes import casa_curve
    x = np.linspace(525, 545, 2001)                               # ascending BE
    v = casa_curve(x, "LA(1.1,3.0,7)", 535.0, 1.2, 1.0)
    ke = (HV - x)[::-1]
    ref = cl.component_curve(ke, "LA(1.1,3.0,7)", HV - 535.0, 1.2, 1.0)[::-1]
    assert np.allclose(v, ref, atol=1e-9)
    assert abs(float(np.sum(v[x < 535.0])) - float(np.sum(v[x > 535.0]))) > 0.05 * float(v.sum())   # really asymmetric


def test_no_fit_and_unsupported_shape():
    b = _make_block()
    b.comments = []
    with pytest.raises(CasaFitError):
        casa_models(b)
    assert parse_fit(["hello", "Calib M = 1 A = 2 BE ADD"]) is None
    b2 = _make_block()
    b2.comments = [c.replace("LA(1.1,1.9,7)", "%&LS2CE2P") for c in b2.comments]
    assert casa_models(b2) == []                                   # region omitted, not silently wrong ...
    assert any("skipped" in s for s in describe_skips(b2))          # ... and the reason is reported
