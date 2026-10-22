import numpy as np
import pytest

from fitstability import (FitProblem, PeakSpec, Protocol, assess, fit_nominal, monte_carlo, multistart,
                          synthetic_spectrum, uniqueness)
from fitstability.lineshapes import pseudo_voigt, shirley_background
from fitstability.io import load_protocol, read_spectrum

TRUE = [(285.0, 1.2, 12000.0), (286.5, 1.3, 9000.0), (288.5, 1.2, 3000.0), (291.5, 3.0, 600.0)]


@pytest.fixture(scope="module")
def spectrum():
    return synthetic_spectrum(TRUE, counts_scale=0.2, seed=1)


def four():
    return Protocol([PeakSpec("C1", 285.0, 1.25), PeakSpec("C2", 286.5, 1.25),
                     PeakSpec("C3", 288.5, 1.25), PeakSpec("Sat", 291.5, 3.5)])


def ten():
    pos = [284.5, 285.2, 285.9, 286.6, 287.3, 288.0, 288.7, 289.4, 290.1]
    return Protocol([PeakSpec(f"C{i}", p, 1.25) for i, p in enumerate(pos)] + [PeakSpec("Sat", 291.5, 3.5)])


def test_pseudo_voigt_area():
    x = np.linspace(-2000, 2000, 400001)  # Lorentzian tails are long
    for eta in (0.0, 0.3, 1.0):
        assert np.trapezoid(pseudo_voigt(x, 0, 1.5, 7.0, eta), x) == pytest.approx(7.0, rel=2e-3)


def test_shirley_bounded_and_monotone():
    x = np.linspace(280, 296, 200)
    y = 100 + 50 * np.exp(-((x - 288) / 1.5) ** 2)
    bg = shirley_background(x, y)
    assert np.all(np.diff(bg) >= -1e-9)
    assert bg[0] == pytest.approx(y[:3].mean(), rel=1e-6)


def test_layout_ties_and_fixed():
    p = Protocol([PeakSpec("a", 1, 1, width_group="g"), PeakSpec("b", 2, 1, width_group="g"),
                  PeakSpec("c", 3, 1, fix_pos=True, fix_fwhm=True)])
    assert p.n_free == 3 + 2 + 1  # 3 areas, pos of a,b, one shared width
    assert p.fwhm_idx[0] == p.fwhm_idx[1] and p.fwhm_idx[2] is None
    lo, hi = p.bounds()
    assert np.all(lo[:3] == 0)
    s = p.unpack(p.pack(p.nominal_start()))
    assert s.pos[2] == 3 and s.fwhm[2] == 1


def test_fit_recovers_truth(spectrum):
    x, y = spectrum
    r = fit_nominal(FitProblem(x, y, four(), count_scale=0.2))
    assert r.rstd < 1.3
    assert np.allclose(r.area[:3], [t[2] for t in TRUE[:3]], rtol=0.08)
    # the weak broad satellite is partly degenerate with the background: only require it present
    assert r.area[3] > 0


def test_wellposed_fit_is_single_minimum(spectrum):
    x, y = spectrum
    ms = multistart(FitProblem(x, y, four(), 0.2), n_starts=16, seed=3)
    v = assess(ms)
    assert v.label == "single minimum", str(v)


def test_overparameterised_fit_shows_multiple_minima(spectrum):
    x, y = spectrum
    ms = multistart(FitProblem(x, y, ten(), 0.2), n_starts=16, seed=3)
    v = assess(ms)
    assert v.label == "multiple minima", str(v)


def test_multistart_reproducible_with_seed(spectrum):
    x, y = spectrum
    prob = FitProblem(x, y, four(), 0.2)
    a = multistart(prob, 4, seed=7).final("area")
    b = multistart(prob, 4, seed=7).final("area")
    assert np.array_equal(a, b)


def test_monte_carlo_spread_small_for_wellposed(spectrum):
    x, y = spectrum
    mc = monte_carlo(FitProblem(x, y, four(), 0.2), n_trials=20, seed=2)
    a = mc.final("area")
    assert np.all(a[:, :3].std(axis=0, ddof=1) / a[:, :3].mean(axis=0) < 0.05)


def test_uniqueness_distinguishes_wellposed_from_correlated(spectrum):
    x, y = spectrum
    u = uniqueness(FitProblem(x, y, four(), 0.2), "C2", "pos", 0.6, 7)
    assert u.mse.max() / u.mse.min() > 1.5
    # ten unconstrained components: a mid-envelope peak is compensated for by its neighbours
    u10 = uniqueness(FitProblem(x, y, ten(), 0.2), "C4", "pos", 0.6, 7)
    assert u10.mse.max() / u10.mse.min() < u.mse.max() / u.mse.min()


def test_io_roundtrip(tmp_path, spectrum):
    x, y = spectrum
    f = tmp_path / "s.csv"
    f.write_text("BE,counts\n" + "\n".join(f"{a},{b}" for a, b in zip(x, y)))
    xr, yr = read_spectrum(f)
    assert np.allclose(xr, x) and np.allclose(yr, y)
    pj = tmp_path / "p.json"
    pj.write_text('{"background":"linear","count_scale":0.2,"peaks":[{"name":"A","pos":285,"fwhm":1.2,"pos_bounds":[284,286]}]}')
    proto, cs = load_protocol(pj)
    assert proto.background == "linear" and cs == 0.2 and proto.peaks[0].pos_bounds == (284, 286)
    pj.write_text('{"peaks":[{"name":"A","pos":285,"fwhm":1.2,"bogus":1}]}')
    with pytest.raises(ValueError):
        load_protocol(pj)


# -- uniqueness follows the paper: free parameters only -------------------------------------
from fitstability import scannable_parameters, unscannable_reason  # noqa: E402


def test_scannable_excludes_constrained_fwhm():
    p = Protocol([PeakSpec("a", 1, 1), PeakSpec("b", 2, 1, width_group="g"), PeakSpec("c", 3, 1, width_group="g"),
                  PeakSpec("d", 4, 1, fix_pos=True, fix_fwhm=True)])
    assert ("a", "fwhm") in scannable_parameters(p)
    assert ("b", "fwhm") not in scannable_parameters(p) and ("c", "fwhm") not in scannable_parameters(p)
    assert ("d", "pos") not in scannable_parameters(p) and ("d", "fwhm") not in scannable_parameters(p)
    assert "tied" in unscannable_reason(p, "b", "fwhm")
    assert "fixed" in unscannable_reason(p, "d", "pos")
    assert unscannable_reason(p, "b", "pos") is None


def test_uniqueness_refuses_tied_fwhm(spectrum):
    x, y = spectrum
    proto = Protocol([PeakSpec("C1", 285.0, 1.25, width_group="m"), PeakSpec("C2", 286.5, 1.25, width_group="m"),
                      PeakSpec("C3", 288.5, 1.25, width_group="m"), PeakSpec("Sat", 291.5, 3.5, fix_fwhm=True)])
    prob = FitProblem(x, y, proto, 0.2)
    with pytest.raises(ValueError, match="not a free fit parameter"):
        uniqueness(prob, "C2", "fwhm")
    with pytest.raises(ValueError, match="not a free fit parameter"):
        uniqueness(prob, "Sat", "fwhm")
    assert uniqueness(prob, "C2", "pos", 0.5, 5).classify() == "unique"  # positions are still free


def test_uniqueness_free_fwhm_scan(spectrum):
    x, y = spectrum
    u = uniqueness(FitProblem(x, y, four(), 0.2), "C1", "fwhm", 0.4, 5)
    assert u.classify() == "unique"


# -- VAMAS ----------------------------------------------------------------------------------
from fitstability import Block, read_vamas, write_vamas, load_protocols  # noqa: E402


def _blocks(spectrum):
    x, y = spectrum
    out = []
    for i, (reg, tr) in enumerate([("C", "1s"), ("C", "1s"), ("O", "1s")]):
        out.append(Block(index=i + 1, block_id=f"{reg}{tr} Scan", sample_id="s", species=reg, transition=tr,
                         hv_ev=1486.68, be=x, y=y * (1 + i), dwell_s=0.05, scans=4, pass_energy_ev=50.0))
    return out


@pytest.mark.parametrize("mode", ["NORM", "MAP"])
def test_vamas_roundtrip_multiblock(tmp_path, spectrum, mode):
    x, y = spectrum
    f = tmp_path / "t.vms"
    write_vamas(f, _blocks(spectrum), mode)
    vf = read_vamas(f)
    assert vf.mode == mode and len(vf.blocks) == 3
    assert vf.regions() == ["C 1s", "O 1s"]
    b = vf.blocks[1]
    assert np.allclose(b.be, x, atol=1e-6) and np.allclose(b.y, y * 2, atol=1e-6)
    assert b.dwell_s == 0.05 and b.scans == 4 and b.pass_energy_ev == 50.0
    assert b.count_scale("cps") == pytest.approx(0.2) and b.count_scale() == 1.0


def test_vamas_select(tmp_path, spectrum):
    f = tmp_path / "t.vms"
    write_vamas(f, _blocks(spectrum))
    vf = read_vamas(f)
    assert [b.index for b in vf.select(None)] == [1, 2, 3]
    assert [b.index for b in vf.select("c1s")] == [1, 2]
    assert [b.index for b in vf.select("1,3")] == [1, 3]
    assert [b.index for b in vf.select("2-3")] == [2, 3]
    with pytest.raises(Exception):
        vf.select("N 1s")
    with pytest.raises(Exception):
        vf.select("9")


def test_vamas_rejects_non_vamas(tmp_path):
    f = tmp_path / "x.vms"
    f.write_text("hello\n")
    with pytest.raises(Exception, match="not a VAMAS"):
        read_vamas(f)


def test_protocolset_per_region(tmp_path, spectrum):
    f = tmp_path / "t.vms"
    write_vamas(f, _blocks(spectrum))
    vf = read_vamas(f)
    pj = tmp_path / "p.json"
    pj.write_text('{"protocols": {"C 1s": {"peaks":[{"name":"A","pos":285,"fwhm":1.2}]},'
                  ' "o1s": {"peaks":[{"name":"B","pos":531,"fwhm":1.5}], "background":"linear"}}}')
    ps = load_protocols(pj)
    assert ps.for_block(vf.blocks[0])[0].peaks[0].name == "A"
    assert ps.for_block(vf.blocks[2])[0].background == "linear"


REAL = "C:/Users/david/Documents/Brazil Training Course September 2025/September 2nd - Advanced Chemical State Analysis/Data/"


@pytest.mark.skipif(not __import__("os").path.isdir(REAL), reason="real VAMAS examples not present")
def test_real_vamas_files_parse():
    import glob
    files = glob.glob(REAL + "*.vms")
    assert files
    for f in files:
        vf = read_vamas(f)
        assert vf.blocks and all(len(b.be) > 10 and np.all(np.diff(b.be) > 0) for b in vf.blocks)
