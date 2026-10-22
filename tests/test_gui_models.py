"""GUI-independent pieces: rule parsing/formatting, protocol JSON round-trip, CSV sniffing, callbacks."""
import json
from dataclasses import replace

import numpy as np
import pytest

from fitstability import PeakSpec, Protocol, FitProblem, multistart, monte_carlo, uniqueness, synthetic_spectrum
from fitstability.experiments import Cancelled
from fitstability.gui import protocol_model as pm
from fitstability.io import (load_protocol, protocol_to_cfg, read_spectrum, save_protocol, sniff_table)

NAMES = ["C1", "C2", "Ti 2p Ti(IV) (1)"]


@pytest.mark.parametrize("kind,text,expect", [
    ("pos", "free", dict(fix_pos=False, pos_tie=None)),
    ("pos", "FIXED", dict(fix_pos=True, pos_tie=None)),
    ("pos", "C1 + 5.2", dict(fix_pos=False, pos_tie=("C1", 5.2))),
    ("pos", "C1 - 6.1", dict(fix_pos=False, pos_tie=("C1", -6.1))),
    ("pos", "Ti 2p Ti(IV) (1) - 5.7753", dict(fix_pos=False, pos_tie=("Ti 2p Ti(IV) (1)", -5.7753))),
    ("fwhm", "C1 * 1.2", dict(fix_fwhm=False, fwhm_tie=("C1", 1.2), width_group="")),
    ("fwhm", "group main", dict(fix_fwhm=False, fwhm_tie=None, width_group="main")),
    ("fwhm", "fixed", dict(fix_fwhm=True, fwhm_tie=None, width_group="")),
    ("area", "C1 * 0.5", dict(area_tie=("C1", 0.5))),
    ("area", "", dict(area_tie=None)),
])
def test_parse_rule(kind, text, expect):
    assert pm.parse_rule(kind, text, NAMES, "C2") == expect


@pytest.mark.parametrize("kind,text", [
    ("pos", "C1 * 2"), ("pos", "Nope + 1"), ("pos", "C2 + 1"), ("fwhm", "C1 + 1"), ("fwhm", "C1 * -1"),
    ("area", "C1 * -0.5"), ("area", "gibberish"),
])
def test_parse_rule_rejects(kind, text):
    with pytest.raises(ValueError):
        pm.parse_rule(kind, text, NAMES, "C2")


def test_format_parse_round_trip():
    p = PeakSpec("D", 290.0, 1.0, pos_tie=("C1", -5.7753), fwhm_tie=("C1", 1.807), area_tie=("C1", 0.5))
    q = replace(p, **pm.parse_rule("pos", pm.format_rule(p, "pos"), ["C1", "D"], "D"),
                **pm.parse_rule("fwhm", pm.format_rule(p, "fwhm"), ["C1", "D"], "D"),
                **pm.parse_rule("area", pm.format_rule(p, "area"), ["C1", "D"], "D"))
    assert (q.pos_tie, q.fwhm_tie, q.area_tie) == (p.pos_tie, p.fwhm_tie, p.area_tie)


def test_rename_and_remove_update_ties():
    peaks = [PeakSpec("A", 1, 1), PeakSpec("B", 2, 1, area_tie=("A", 0.5), pos_tie=("A", 1.0))]
    r = pm.rename_peak(peaks, 0, "Main")
    assert r[0].name == "Main" and r[1].area_tie == ("Main", 0.5) and r[1].pos_tie == ("Main", 1.0)
    assert pm.validate(r) == ""
    gone = pm.remove_peak(r, 0)
    assert len(gone) == 1 and gone[0].area_tie is None and gone[0].pos_tie is None
    assert pm.validate(gone) == ""


def test_validate_reports_cycles_and_empty():
    assert pm.validate([]) != ""
    a = PeakSpec("A", 1, 1, area_tie=("B", 1.0))
    b = PeakSpec("B", 2, 1, area_tie=("A", 1.0))
    assert "circular" in pm.validate([a, b])


def test_protocol_json_round_trip_with_ties_and_shapes(tmp_path):
    peaks = [PeakSpec("C1", 285.0, 1.2, shape="GL(30)", pos_bounds=(283.0, 287.0), area_bounds=(0.0, 5000.0)),
             PeakSpec("C2", 290.0, 1.2, shape="GL(30)", area_tie=("C1", 0.5), fwhm_tie=("C1", 1.2), pos_tie=("C1", 5.0)),
             PeakSpec("S", 292.0, 3.0, fix_pos=True, width_group="w", eta=0.1)]
    f = tmp_path / "p.json"
    save_protocol(f, Protocol(peaks, "linear"), 0.5)
    json.loads(f.read_text())                                  # valid JSON (no Infinity)
    proto, cs = load_protocol(f)
    assert cs == 0.5 and proto.background == "linear"
    for a, b in zip(peaks, proto.peaks):
        assert a == b
    cfg = protocol_to_cfg(Protocol(peaks[1:2] + [PeakSpec("C1", 285.0, 1.2)], "none"))
    assert "eta" not in cfg["peaks"][0]                        # defaults omitted


def test_casa_background_round_trips(tmp_path):
    from fitstability.casa import CasaBackground
    bg = CasaBackground("Shirley", 3, (1.0, 2.0, 0.0, 0.0, 0.0, 0.0), 2.5)
    f = tmp_path / "p.json"
    save_protocol(f, Protocol([PeakSpec("A", 1, 1)], bg))
    proto, _ = load_protocol(f)
    assert isinstance(proto.background, CasaBackground) and proto.background == bg


def test_sniff_and_read_table(tmp_path):
    x = np.linspace(280, 296, 40)
    for name, text in {
        "comma.csv": "BE,Counts,Other\n" + "\n".join(f"{a},{a * 2},{a * 3}" for a in x),
        "semi.txt": "# comment\nBE;Counts\n" + "\n".join(f"{a};{a * 2}" for a in x),
        "tab.dat": "\n".join(f"{a}\t{a * 2}" for a in x),
        "ws.asc": "\n".join(f"{a}   {a * 2}" for a in x),
    }.items():
        f = tmp_path / name
        f.write_text(text)
        info = sniff_table(f)
        assert info.n_numeric == 40 and info.n_columns >= 2
        xs, ys = read_spectrum(f)
        assert np.allclose(xs, x) and np.allclose(ys, x * 2)
    info = sniff_table(tmp_path / "comma.csv")
    assert info.header == ["BE", "Counts", "Other"] and info.delimiter == ","
    xs, ys = read_spectrum(tmp_path / "comma.csv", 0, 2)
    assert np.allclose(ys, x * 3)
    with pytest.raises(ValueError):
        read_spectrum(tmp_path / "comma.csv", 0, 7)


def _problem():
    x, y = synthetic_spectrum([(285.0, 1.2, 12000), (286.5, 1.3, 9000)], counts_scale=0.2, seed=1)
    return FitProblem(x, y, Protocol([PeakSpec("C1", 285.0, 1.25), PeakSpec("C2", 286.5, 1.25)]), 0.2)


def test_progress_and_cancel_callbacks():
    prob = _problem()
    seen = []
    multistart(prob, 5, seed=0, progress=lambda d, t: seen.append((d, t)))
    assert seen == [(i, 5) for i in range(1, 6)]
    seen.clear()
    monte_carlo(prob, 4, seed=0, progress=lambda d, t: seen.append(d))
    assert seen == [1, 2, 3, 4]
    seen.clear()
    uniqueness(prob, "C1", "pos", 0.4, 5, progress=lambda d, t: seen.append((d, t)))
    assert seen[-1] == (5, 5)
    n = []
    with pytest.raises(Cancelled):
        multistart(prob, 10, seed=0, progress=lambda d, t: n.append(d), cancel=lambda: len(n) >= 2)
    assert len(n) == 2
    with pytest.raises(Cancelled):
        uniqueness(prob, "C1", "pos", 0.4, 5, cancel=lambda: True)


def test_plots_draw_into_given_figure_without_pyplot():
    import sys
    from matplotlib.figure import Figure
    from fitstability import assess
    from fitstability.plots import plot_data, plot_fit, plot_monte_carlo, plot_multistart, plot_uniqueness
    prob = _problem()
    ms = multistart(prob, 4, seed=0)
    for draw in (lambda f: plot_multistart(ms, verdict=assess(ms), fig=f), lambda f: plot_monte_carlo(monte_carlo(prob, 4, seed=0), fig=f),
                 lambda f: plot_uniqueness([uniqueness(prob, "C1", "pos", 0.4, 5)], fig=f),
                 lambda f: plot_fit(prob, ms.best().as_start(), fig=f), lambda f: plot_data(prob.x, prob.y, fig=f, markers=[(285.0, "C1")])):
        fig = Figure()
        draw(fig)
        assert fig.axes
