"""Editing a protocol as text rules: parsing and formatting, with no Qt dependency.

A peak's constraints are shown and typed as short rules, one per parameter:

========== ==========================================================================
position   ``free`` | ``fixed`` | ``<peak> + 5.2`` / ``<peak> - 5.2`` (offset to another peak)
FWHM       ``free`` | ``fixed`` | ``<peak> * 1.2`` (ratio) | ``group <label>`` (shared width)
area       ``free`` | ``<peak> * 0.5`` (ratio)
========== ==========================================================================

These are the same constraints CasaXPS stores, so a Casa fit round-trips through the table.
"""
from __future__ import annotations

import re
from dataclasses import replace

from ..model import PeakSpec, Protocol

_NUM = r"-?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?|-?\.\d+"
_TIE = re.compile(rf"^(?P<ref>.*\S)\s+(?P<op>[+\-*])\s*(?P<val>{_NUM})$")

RULE_HELP = {
    "pos": "free | fixed | <peak> + 5.2 (offset in eV; use - for below)",
    "fwhm": "free | fixed | <peak> * 1.2 (ratio) | group <label> (shared width)",
    "area": "free | <peak> * 0.5 (ratio)",
}


def format_rule(p: PeakSpec, kind: str) -> str:
    if kind == "pos":
        if p.pos_tie:
            ref, d = p.pos_tie
            return f"{ref} {'+' if d >= 0 else '-'} {abs(d):g}"
        return "fixed" if p.fix_pos else "free"
    if kind == "fwhm":
        if p.fwhm_tie:
            return f"{p.fwhm_tie[0]} * {p.fwhm_tie[1]:g}"
        if p.width_group:
            return f"group {p.width_group}"
        return "fixed" if p.fix_fwhm else "free"
    if kind == "area":
        return f"{p.area_tie[0]} * {p.area_tie[1]:g}" if p.area_tie else "free"
    raise ValueError(kind)


def parse_rule(kind: str, text: str, names: list[str], own: str = "") -> dict:
    """The :class:`PeakSpec` field changes a rule implies. Raises ``ValueError`` with a readable message."""
    t = (text or "").strip()
    low = t.lower()
    if kind == "pos":
        if low in ("", "free"):
            return dict(fix_pos=False, pos_tie=None)
        if low == "fixed":
            return dict(fix_pos=True, pos_tie=None)
        ref, val = _tie(t, names, own, ("+", "-"), "position, e.g. 'C1 + 5.2'")
        return dict(fix_pos=False, pos_tie=(ref, val))
    if kind == "fwhm":
        if low in ("", "free"):
            return dict(fix_fwhm=False, fwhm_tie=None, width_group="")
        if low == "fixed":
            return dict(fix_fwhm=True, fwhm_tie=None, width_group="")
        m = re.match(r"^group[\s:]+(\S.*)$", t, re.I)
        if m:
            return dict(fix_fwhm=False, fwhm_tie=None, width_group=m.group(1).strip())
        ref, val = _tie(t, names, own, ("*",), "FWHM, e.g. 'C1 * 1.2'")
        if val <= 0:
            raise ValueError("a FWHM ratio must be positive")
        return dict(fix_fwhm=False, fwhm_tie=(ref, val), width_group="")
    if kind == "area":
        if low in ("", "free"):
            return dict(area_tie=None)
        ref, val = _tie(t, names, own, ("*",), "area, e.g. 'C1 * 0.5'")
        if val < 0:
            raise ValueError("an area ratio cannot be negative")
        return dict(area_tie=(ref, val))
    raise ValueError(kind)


def _tie(text, names, own, ops, example):
    m = _TIE.match(text)
    if not m or m.group("op") not in ops:
        raise ValueError(f"cannot read {text!r}: type free, fixed, or a tie of the form {example}")
    ref = m.group("ref").strip()
    if ref not in names:
        raise ValueError(f"no peak named {ref!r}")
    if ref == own:
        raise ValueError("a peak cannot be tied to itself")
    val = float(m.group("val"))
    return ref, (-val if m.group("op") == "-" and ops == ("+", "-") else val)


def rename_peak(peaks: list[PeakSpec], index: int, new: str) -> list[PeakSpec]:
    """Rename a peak and update every tie that refers to it."""
    old = peaks[index].name
    out = []
    for i, p in enumerate(peaks):
        ch = {}
        if i == index:
            ch["name"] = new
        for k in ("area_tie", "fwhm_tie", "pos_tie"):
            t = getattr(p, k)
            if t is not None and t[0] == old:
                ch[k] = (new, t[1])
        out.append(replace(p, **ch) if ch else p)
    return out


def remove_peak(peaks: list[PeakSpec], index: int) -> list[PeakSpec]:
    """Remove a peak; a peak tied to it keeps its current value and becomes free/fixed (ties dropped)."""
    gone = peaks[index].name
    out = []
    for i, p in enumerate(peaks):
        if i == index:
            continue
        ch = {}
        for k in ("area_tie", "fwhm_tie", "pos_tie"):
            t = getattr(p, k)
            if t is not None and t[0] == gone:
                ch[k] = None
        out.append(replace(p, **ch) if ch else p)
    return out


def unique_name(peaks: list[PeakSpec], base: str = "P") -> str:
    names = {p.name for p in peaks}
    n = len(peaks) + 1
    while f"{base}{n}" in names:
        n += 1
    return f"{base}{n}"


def validate(peaks: list[PeakSpec], background="shirley") -> str:
    """'' if the peaks form a valid protocol, else the reason."""
    if not peaks:
        return "add at least one peak"
    try:
        Protocol(peaks, background)
    except ValueError as e:
        return str(e)
    return ""
