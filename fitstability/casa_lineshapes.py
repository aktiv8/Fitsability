# Vendored unchanged from the author's eXPoSe SpectraDeck (Mutli-Data-Browser/lineshapes.py, 2026-10-22).
# CasaXPS line shapes and backgrounds reconstructed from stored fits; see the module docstring for validation notes.
r"""Line shapes and backgrounds for reconstructing CasaXPS fits.

All curves are computed on an **ascending kinetic-energy grid** (the frame in
which CasaXPS stores positions) and scaled to the stored *area* (intensity
integrated over energy).

* ``GL(m)``  product Gaussian/Lorentzian, ``m`` = % Lorentzian
* ``SGL(m)`` sum Gaussian/Lorentzian
* ``LA(a,b,m)`` (also ``LA(a,m)`` = ``LA(a,a,m)``, or ``LA(m)``) asymmetric
  Lorentzian, Gaussian-broadened
* ``LF(a,b,w,m)`` as LA with a finite tail of width ``w``

GL and SGL are exact. **LA and LF are reconstructions**: their Gaussian
broadening scale is still an empirical calibration (below), but the base
kernel itself is no longer a third-party guess -- CasaXPS's own *Cookbook*
(Casa Software Ltd., 2026, pp.65-66) gives ``LA(x: alpha,beta,w) = N *
INTEGRAL(lg(tau; alpha,beta) * g(x-tau; w) dtau)`` with ``lg(x: alpha,beta)
= l(x)**alpha`` for ``x <= 0`` and ``l(x)**beta`` for ``x > 0`` -- matching
this module's formula exactly, and matching the independent third-party
write-up this was originally reconstructed from (Major et al., *Surf.
Interface Anal.* 53 (2021) 689, Eq. 6; Major, Shah, Fernandez, Fairley (Casa
Software Ltd.) and Linford, *Vacuum Technology & Coating*, March 2020,
Eq. 3). The Gaussian broadening's own scale was calibrated on real CasaXPS
fits (``tests/`` and the validation notes in the README give the residuals).
The stored areas, positions and widths are CasaXPS's own numbers.

**The shared Lorentzian width.** Raising a Lorentzian to a power ``p != 1``
shrinks its own half-max distance by ``sqrt(2**(1/p) - 1)`` relative to
``p = 1``, so using the file's ``fwhm`` unscaled as each side's own width (as
this module did before) makes a strongly asymmetric component (e.g.
``LA(1.2,5,8)``, CasaXPS's sharp metallic-tail cutoff) reconstruct far too
tall for its stored area -- confirmed on two independent real files
(titanium and vanadium metallic-tail examples in ``tests/``), the
reconstructed peak standing up to ~10-18 % of the peak height above the raw
data at the component's own maximum, hidden by ``residual_rms`` because that
metric averages over the whole curve. The fix (from KherveFitting's own
``LA``/``LAxG`` implementation, ``libraries/Peak_Functions.py``,
github.com/KherveFitting/KherveFitting, retrieved 2026-09-25) is **one
Lorentzian width shared by both sides**, chosen so the two sides' own
half-max distances add up to the stored ``fwhm`` (not each side independently
reproducing it, which was tried first and made CasaXPS's own worked LA
examples with mild asymmetry measurably worse while helping the severe
cases): ``F = 2*fwhm / (sqrt(2**(1/a)-1) + sqrt(2**(1/b)-1))``, which is
exactly ``fwhm`` for an unmodified Lorentzian (``a == b == 1``) -- the GL/SGL
shapes and the 2-argument ``LA(m)``/``LF(...)`` shorthand (which default
``a`` and ``b`` to 1) are unaffected; a symmetric but non-unity pair (e.g.
``LA(2,2,m)``) is rescaled too, since raising *either* side to a power other
than 1 narrows it the same way.
Checked against every real ``LA``/``LF``-fitted file available (titanium,
copper, both vanadium exports, MXene): titanium's worst-case overshoot fell
from 18.3 % to 8.7 %, this file's vanadium example from 9.85 % to 4.70 %, and
6 of MXene's 8 affected regions improved too (the other 2 regressed by under
one percentage point). ``TestLAAsymmetryAccuracy`` in ``tests/test_casafit.py``
is a canary against this getting worse.

**The 1-argument ``LA(m)`` shorthand is not ``LA(1,1,m)``.** Every finding
above this paragraph, and two rejected leads that used to be recorded here
(a "corrected direction" retune, and a claim that no ``GAUSS_K`` can round out
a symmetric ``LA(m)`` component's own apex), were all tested by feeding the
shorthand's ``m`` (0-100) straight into the same formula the explicit
``LA(a,b,m)`` form uses -- which is wrong. CasaXPS's own definition (Fairley
et al., "Practical guide to understanding goodness-of-fit metrics ... using
nylon as an example," *J. Vac. Sci. Technol. A* 41(1) 2023, supplementary
information, Eq. 4-6; also *Cookbook* 2026 p.19, Figure 15's caption) is
``LA(x: alpha,beta,n) = N * INTEGRAL(lg(tau; alpha,beta) * g(x-tau; f_G(n))
dtau)`` (a generalised Lorentzian convolved with a Gaussian whose width is
controlled by ``n``) and, separately, ``LA(x,m) = LA(x: 1,1, 1401 -
(m/100)*1401)``: the shorthand's ``m`` and the explicit form's ``n`` run in
**opposite directions** (``LA(0)`` -> ``n = 1401``, ``LA(100)`` -> ``n =
0``), not the same variable. ``parse_shape`` used to set ``m = ps[0]`` for
the shorthand directly; it now applies this conversion, so every real
``LA(m)``-shorthand file (PET's ``LA(50)`` among them) is now evaluated at
its true ``n = 1401 - 14.01*m``, and the explicit 3-argument form is
unaffected (it was already using ``n`` correctly).

**The 2-argument ``LA(a,n)`` form is a separate case, and was not
distinguished from the 1-argument shorthand above at all** -- a second,
independent bug, found while reading the *Cookbook* itself (2026 p.66):
*"Abbreviation for LA lineshape: LA(x: alpha, w) = LA(x: alpha, alpha, w)"*
-- ``b`` is implicitly ``a`` (a symmetric-but-non-unity Lorentzian power),
and the trailing number is the explicit form's own ``n``, used directly --
*not* run through the 1-argument shorthand's ``m``-to-``n`` conversion
above. ``parse_shape`` used to fall into the same ``else`` branch as the
1-argument case for *any* fewer-than-3-argument string, taking only
``ps[0]`` through the wrong (1-argument) formula regardless of whether a
second number was present at all -- so a real ``LA(1.53,243)`` file was
silently evaluated as ``a=b=1, n=1401-14.01*1.53 approx 1380`` instead of
the correct ``a=b=1.53, n=243``, both the exponent and the Gaussian width
wrong.
Confirmed on the one real file with this form found across the full local
corpus (~700 files, ``D:\Temp`` and ``~\Downloads``):
``~\Downloads\assigned.vms``, a 5-sample Mo 3d/S 2s fit (MoS2/MoO2/MoO3/WS2
mixed system, ~10 ``LA(1.53,243)`` components per region). Comparing
``casafit.curves()``'s reconstruction against this file's own raw data,
before/after the fix, ``residual_rms`` / ``chi2_red`` per region:
``RW_WS2_MoS2_thicker`` 8.86% / 139.3 -> 5.32% / 54.9,
``RW_Tha_MoS2`` 7.23% / 930.2 -> 3.98% / 329.8,
``RW_sonic_MoS2`` 7.73% / 2339.2 -> 4.39% / 917.5,
``RW_WS2_MoS2`` 8.59% / 129.8 -> 5.38% / 52.6,
``RW_Nb_MoS2`` 10.81% / 455.3 -> 6.84% / 209.7 -- both metrics roughly halve
across every region.

This was found the same way the ``DS`` calibration below was: two new
synthetic files built specifically to sweep ``m`` in a controlled way against
one shared raw peak, with the components' FWHM deliberately linked equal so
it is the only free parameter left to compensate for a "wrong" ``m`` --
``D:\Temp\for claude files\synthetic PMMA.vms`` (``LA(m)`` shorthand, ``m =
0..100``, 11 CasaXPS refits) and ``D:\Temp\for claude files\synthetic PMMA.
3 param LAvms.vms`` (explicit ``LA(1,1,m)``, ``m = 1000..0``), both refitting
the same raw data. Comparing the reconstructed envelope to the raw data
directly does not work on this pair (their components' stored ``Area`` is
~100-450x too small for the data's actual peak height -- an unrelated
data-quality artifact of the file, not a shape question); instead, since FWHM
is the fit's only free width parameter, how the reported (linked) FWHM must
drift as ``m`` sweeps -- via the standard pseudo-Voigt width approximation,
``FWHM_V = 0.5346*FWHM_L + sqrt(0.2166*FWHM_L**2 + FWHM_G**2)`` -- is a clean
signature of the true broadening-vs-``m`` relationship: the correct ``(K, P)``
should make the implied total width stay constant across the whole sweep,
since it is the same underlying peak throughout. Feeding both files' raw
``m`` into one formula gave opposite-signed answers (this is what surfaced
the shorthand/explicit scale bug above); on the corrected, unified ``n``
scale, a joint fit -- minimising that width-consistency spread (today's
formula: 8.9 % coefficient of variation across the 21-point sweep) while
never letting any real ``LA``/``LF``-fitted region in the existing corpus get
more than ~0.25 percentage points worse than it already was -- converges on
``GAUSS_K["LA"] = 0.41``, ``GAUSS_P["LA"] = 0.60`` (replacing ``0.20`` and an
implicit exponent of 1; used the same way ``GAUSS_P["DS"]`` already is,
below), bringing the sweep's spread down to ~7.9 %. **The real corpus this
was checked against is bigger than the "six real fits" this constant was
first calibrated on**, and bigger than what the earlier "not adopted"
attempts above were checked against: 25 ``LA``/``LF``-fitted regions across 7
files (titanium `fitting example` and the `Metal Depth Profile` instructors'
file -- the latter has both a clean 2-component metal region and an
8-component mixed metal+oxide region in one file, the busiest region in the
corpus and the one that ruled out a more aggressive retune tried first --
copper, vanadium (both the plain and `For Instructors` exports), two MXene
files, PET). Measured effect: titanium `fitting example`'s worst-case
overshoot (the number the old ``TestLAAsymmetryAccuracy`` ceiling comment was
built around) falls from 8.70 % to 6.27 %; every other region in the corpus
moves by at most ~0.24 percentage points either way; PET's shorthand
``LA(50)`` component (now correctly ``n = 700.5``, not ``50``) moves its
C 1s/O 1s overshoot by ~0.02 pp. The worst case across the whole corpus stays
PET's O 1s, 9.46 % to 9.43 % (very slightly better). Do not retune
``GAUSS_K``/``GAUSS_P`` again without new evidence beyond what produced this
result -- but do re-run this comparison (this file's own module docstring
here, and ``TestLAAsymmetryAccuracy``/``TestPETFile``/``TestLASweepShape``)
before touching either constant, since the previous two "not adopted"
conclusions in this file were reached on a mis-scaled ``m`` and turned out to
be wrong, and a first attempt at this same retune (``K=0.90, P=0.80``, fit
only against a 5-file subset of the corpus) passed every test available at
the time but pushed the `Metal Depth Profile` file's 8-component region's
overshoot from 8.8 % to 16.3 % -- always check the *busiest*, most
overlapping real region available, not just the cleanest ones.

**The 1-argument ``LA(m)`` shorthand's residual overshoot is not a
``GAUSS_K``/``GAUSS_P`` calibration gap -- it is intrinsic to ``a=b=1`` and
cannot be tuned away.** Checked directly against real data, not a synthetic
sweep, thanks to a file built specifically for the comparison:
``D:\Temp\for claude files\PtCl2_new\PtCl2_refitted.vms`` refits the *same*
Cl 2p / Pt 4d / Pt 4f data twice -- once with the shorthand (``LA(70)``,
``LA(90)``), once with the explicit form (``LA(1.52,43)``, ``LA(1.2,43)``) --
so the two reconstructions can be measured against the same raw counts.
Sweeping ``GAUSS_K["LA"]`` from today's ``0.41`` down to ``0`` (no Gaussian
broadening at all -- a bare, unconvolved ``a=b=1`` Lorentzian) moves the
shorthand's overshoot by well under half a percentage point (Cl 2p 4.51% ->
4.50%, Pt 4d 4.82% -> 4.79%, Pt 4f 4.42% -> 4.32%): the convolution width is
not the cause. Printing the residual point by point shows two distinct
effects. First, a slow, low-level excess that persists tens of eV past the
peak (e.g. Pt 4f's high-BE tail, ~700-1800 counts above a ~7000-count flat
continuum from BE 78-88 eV) -- but this is present, at similar magnitude, in
the *explicit*-form reconstruction of the same data too (which plots well),
so it is a generic property of summing large-area Lorentzian tails, not
something specific to the shorthand. Second, and the larger effect: a
**valley-bridging excess between the doublet's two main components**
(Pt 4f's biggest single residual, ~2200 counts on a ~50000-count peak, sits
in the dip between the two spin-orbit lines, not in either peak's own tail).
Raising a Lorentzian's power above 1 narrows it and suppresses exactly this
kind of inter-peak bridging (this is the whole reason the explicit
``LA(1.2,43)``/``LA(1.52,43)`` refit of the same components measures lower
overshoot, 3.49%/2.70%, than the ``a=b=1`` shorthand's 4.42%/4.51%) -- but
CasaXPS's own shorthand definition (Fairley et al. 2023, *Cookbook* 2026, both
cited above) is ``a=b=1`` by construction, so there is no calibration
constant that can narrow it without changing what ``LA(m)`` means. **Do not
retune ``GAUSS_K``/``GAUSS_P`` to chase this again without evidence beyond
what produced this result** (the sweep above already covers the full range
either constant could plausibly take).

**The LF finite tail was reconstructed with the wrong functional form
entirely -- found and fixed from a primary source, not a retune.** The two
"tried and rejected" investigations that used to be recorded in this
paragraph (a polynomial ``(1 - r^2)^2`` hard cutoff at ``r = |x-pos|/fwhm/w``,
and a same-shaped logistic-sigmoid alternative from a third-party
reconstruction) were both variations on the same wrong idea: a taper
*multiplying* the plain ``LA(a,b)`` curve down to (or towards) zero at some
radius. Both failed for the same reason -- ``component_curve`` rescales a
component so its analytic integral equals CasaXPS's stored area
(``test_a_narrow_window_does_not_inflate_the_peak`` protects this), so any
tail removed by a multiplicative taper has to come back as extra height at
the peak. The actual formula, confirmed directly from CasaXPS's own primary
description with runnable reference code -- Major, Shah, Avval, Fernandez,
Fairley, Linford, "Advanced Line Shapes in X-Ray Photoelectron Spectroscopy
II. The Finite Lorentzian (LF) Line Shape" (with MATLAB code), *Vacuum
Technology & Coating*, April 2020, cross-checked against Part I (the same
authors, March 2020, which independently confirms the plain ``LA(alpha,beta)``
equation and its exponent-side convention -- low-KE side -> alpha, high-KE
side -> beta -- exactly matching this module's ``t < 0 -> a`` split) -- is
not a taper on the *curve* at all: it is a smoothly *varying exponent*.
``ex(x) = 3.0 - (3.0 - a) / (1 + 4*((x-pos)/w)**2)`` on the low-KE side (and
the analogous expression with ``b`` on the high-KE side), with the base
Lorentzian raised to ``ex(x)`` instead of a fixed ``a``/``b``: at ``x = pos``,
``ex = a`` (unchanged); as ``|x-pos|`` grows past ``w`` (used directly, in
eV -- **not** scaled by ``fwhm``, unlike the removed taper), ``ex`` rises
towards a fixed ceiling of ``3.0`` (a value the source states outright, not a
fitted constant), so the tail decays by an ever-steeper power law rather than
being clipped to exactly zero. The source's own reference code also uses the
base Lorentzian's width as a plain, fixed ``F`` (``F = 1; %peak width``, used
directly, with no rescaling formula given anywhere in either article) --
implemented here as plain ``fwhm``, not the ``_shared_width`` rescale plain
``LA`` uses (checked separately: removing ``_shared_width`` from plain,
explicit ``LA(a,b,n)`` -- no ``w`` -- regresses badly on the wider corpus,
e.g. ``assigned.vms``'s W 4f overshoot jumping from 4.6% to 54%, so that fix
stays; it is specifically the ``w>0`` LF branch that had the wrong mechanism).
Checked against every real ``LF``-fitted region available (5, across
``Titanium Metal Depth Profile - INSTRUCTORS.vms``'s Ti 2p and
``D:\Temp\for claude files\PtCl2\PtCl2_quantified.vms``'s Cl 2p / Pt 4d /
Pt 4f x2): ``residual_rms`` improves on every single region, substantially --
Pt 4d 0.1499 -> 0.0686 (more than halved), Ti 2p 0.0410 -> 0.0307, Cl 2p
0.0629 -> 0.0523, Pt 4f 0.0592 -> 0.0491, Pt 4f (area2) 0.0718 -> 0.0607.
Overshoot (a single worst-point metric, unlike ``residual_rms``'s whole-curve
average) moves in mixed directions region to region -- expected, and not a
sign the fix is wrong, since a better whole-curve fit can still have a
slightly different single peak excursion. The tail's visible extent in
``plots.draw_fit`` (``_COMPONENT_VISIBLE_FLOOR``) is unaffected -- still a
display choice, not a change to the shape or its stored area.

One further real form was found but not adopted: ``PtCl2_quantified.vms``'s
Pt 4d region uses a genuine 5-argument string, ``LF(0.5,0.6,45,180,1)`` /
``LF(0.5,0.6,45,180,3)`` -- a 5th number ``parse_shape`` still drops. The
obvious guess (the 5th number overrides the ``3.0`` ceiling per component)
was tried and made the reconstruction measurably worse (``residual_rms``
0.069 -> 0.091), so it is wrong; the 5th argument is left unreconstructed and
honestly unflagged rather than guessed again, the same precedent as this
module's other "not enough evidence" cases (``QF``, ``H``/``F``, the Tougaard
region line's unused leading numbers).

**A(a,b,n)GL(m) / A(a,b,n)SGL(m): the Gelius asymmetric shape.** Sourced
from CasaXPS's own "Peak Fitting in XPS" (Casa Software Ltd, 2006, p.20 --
read directly, not via a third party): a symmetric ``GL(m)``/``SGL(m)`` base
plus a one-sided tail on the low-KE (high binding energy) side,
``A(x) = base(x) + w(a,b)*[AW(x) - G0(x)]`` for ``x`` below the position,
plain ``base(x)`` above it, where ``w(a,b) = b*(0.7 + 0.3/(a+0.01))``,
``AW(x) = exp(-(u/(fwhm + a*u))**2)`` with ``u = 2*sqrt(ln2)*|x-pos|``, and
``G0`` the plain Gaussian of the same ``fwhm`` (``AW`` at ``a = 0``). This
independent reading matches KherveFitting's own ``A_GL``/``A_SGL``
(``Peak_Functions.py``, retrieved 2026-09-27), which cites the same primary
source and, once its own BE-ascending axis is un-mirrored back to this
module's ascending-KE grid, is algebraically identical -- good agreement
between an independent derivation and a third party's, from the same source.

**The tail term does not always decay to zero.** ``AW(x) - G0(x)`` approaches
a nonzero plateau ``w*exp(-1/a**2)`` as ``x`` moves far from the peak
(confirmed both analytically and numerically), rather than the power-law
decay every other reconstructed shape here has -- CasaXPS's own whitepaper
warns of exactly this ("intensity from the nominal ... region ... can cause
... to move beyond the region defining the background ... not suitable for
quantification without intensity calibration"). For small ``a`` (the one
real value seen anywhere, in the whitepaper's own worked poly(propylene)
example, ``A(0.15,0.7,20)SGL(12)``, is ``0.15``) the plateau is negligible
(``exp(-1/0.15**2)`` ~ 1e-19) and ``component_curve``'s wide-grid area
normalisation is unaffected; a fit using a much larger ``a`` would make that
normalisation window-dependent and unreliable. Not guarded against here --
no real file with a large ``a`` has been seen to know whether CasaXPS's own
usage ever goes there.

**The Gaussian convolution width (the shape string's ``n``, 0-499) is
uncalibrated.** The primary source's own formula has no ``n`` in it at all
(a separate line says only that the shape is "convoluted with a Gaussian
with width characterized by an integer 0 <= n <= 499", the same phrasing
used for ``DS``); lacking a real ``A``/``SGL``-Gelius-fitted file anywhere in
the available corpora (checked: the Brazil Training Course, both Avantage/
Kratos reference sets, the CasaXPS VAMAS corpus and the MXene/PET files all
have zero components using this shape), ``GAUSS_K["A"]`` and
``GAUSS_P["A"]`` are left at the same neutral ``(25/n)**1`` reference point
DS/LA already use (``K = 1.0, P = 1.0`` -- i.e. ``gw = fwhm`` at ``n = 25``),
not a calibrated value. Retune when a real fitted file becomes available,
the same way DS's own calibration was redone once real ``(a, n)`` coverage
appeared.

**TLA: formula and parameter order now confirmed from a primary source; the
Gaussian-convolution calibration is not.** Originally implemented (flagged,
at the user's explicit direction) from KherveFitting's own reconstruction
alone (``Peak_Functions.py``, ``TLA``/``_tla_unit_profile``/
``_tla_width_scale``), which itself cited an unlocatable "CasaXPS Cookbook
(2026), Line shapes section" -- no way at the time to cross-check it against
a primary CasaXPS source, the same "third party might be wrong" risk the
``LA(m)`` shorthand bug above was a cautionary tale about. That Cookbook has
since been located (Casa Software Ltd., 2026, p.66) and checked term by
term: ``TLA(x: alpha,mu,w) = N * INTEGRAL(T(tau; alpha,mu) * g(x-tau; w)
dtau)``, ``T(x: alpha,mu) = [1/(1+4x^2)]**alpha * [pi/2 - atan(2x) + pi/mu]``
for ``mu > 0, alpha > 0`` -- an exact match to this module's own
``_tla_values`` (Lorentzian raised to ``alpha`` = ``sp["a"]``, modulated by
the same arctangent tail factor with ``mu`` = ``sp["b"]``, larger ``mu`` =
weaker asymmetry), **and** confirms the shape string's parameter order is
``TLA(alpha, mu, n)`` as this module already assumed (KherveFitting's own
function signature groups ``(mu, wg, alpha)`` for its *fitting grid
columns* only, which said nothing about the on-disk string's own order --
this was a guess that turned out right, not something the Cookbook itself
was needed to get the code working, but it removes the doubt). What the
Cookbook does *not* give is a numeric Gaussian-convolution width: the width
``w``/``n`` has no closed form connecting the shape string's ``0-499`` code
to an actual eV width for a given ``fwhm``, so ``GAUSS_K["TLA"] =
GAUSS_P["TLA"] = 1.0`` remains the same neutral, uncalibrated reference
point DS/``A`` above use (``gw = fwhm`` at ``n = 25``) -- retune only when a
real ``TLA``-fitted file becomes available, the same way DS's and the
2-argument ``LA`` bug above were. ``_tla_width_scale`` measuring the base
shape's FWHM numerically (no closed form exists for that either) is
unaffected by any of this. Until a real file calibrates the convolution
width, treat a ``TLA``-reconstructed component's overall breadth as
illustrative even though its underlying asymmetric-Lorentzian shape and
parameter order are now on solid ground.

**Tail suffix -- formula now known, still unvalidated.** CasaXPS also lets a
``GL``/``SGL`` shape string carry a trailing ``T(k)`` tail modifier
(``GL(30)T(1.5)``), used for asymmetric metallic peaks. ``parse_shape``
recognises and strips this suffix so the base shape's own parameters parse
correctly, but the tail itself is not reconstructed. Unlike when this was
first written, the formula is no longer unpublished: the *Cookbook* (2026,
p.67) gives it in full -- ``E(x: j,k) = j*exp(k*x)`` for ``x < 0`` else
``0``, and ``SGL(x,m,n)T(k) = (SGL(x,m) + (1-SGL(x,m))*E(x: 1,k)) * G(n)``
(``T(k)`` is the ``j=1`` case of a more general ``PHI(j,k)``; ``GL`` has the
same two variants). The Cookbook also documents two further mechanisms not
implemented here: a general ``ST(mu,gamma)`` asymmetry prefix applicable to
*any* line shape (p.66: ``S(mu)F(x) = N*F(x)*[INTEGRAL(F(phi)dphi, x, inf) +
A/mu]``, ``ST(mu,gamma)F(x) = S(mu)F convolved with a Gaussian of width
gamma`` -- distinct from both ``T(k)``/``PHI(j,k)`` above and the Gelius
``A(a,b,n)`` shape below), and a plain 3-argument extension of ``GL``/
``SGL`` themselves, ``GL(x,m,n) = GL(x,m) * G(n)`` (p.67) -- a second
Gaussian convolution on top of the base pseudo-Voigt that this module's own
``GL``/``SGL``/``VOIGT`` branch of ``parse_shape`` would today silently
truncate to plain ``GL(m)``, reading only ``ps[0]``, if a string like that
ever appeared. None of ``PHI(j,k)``, ``ST(mu,gamma)``, or a 2-argument
``GL(m,n)``/``SGL(m,n)`` turned up anywhere in the ~700-file local corpus
scan that found the real 2-argument ``LA`` file above, nor in any of the
individually-named real files elsewhere in this module -- so none of those
three is implemented; an affected component still draws as its plain base
shape and is flagged ``approximate`` by ``is_exact`` (honestly
unreconstructed, the same treatment ``QF``/``H``/``F`` get below), rather
than silently wrong.

**``GL(m)T(k)`` / ``SGL(m)T(k)`` is now drawn** (``_raw_values``): the
Cookbook formula above for a one-argument shape, no Gaussian convolution,
with ``x`` in units of the component's FWHM and ``x < 0`` the low-KE
(high-binding-energy) side, ``F = base + (1 - base) * exp(k x)`` there and
the plain base above the position, so it is the base peak with an
exponential tail added on the high-BE side (continuous at the position,
larger ``k`` fades faster). **Evidence, honestly:** the formula is
Casa's own (Cookbook p.67), and HarwellXPS's Lineshape Tester (Armoury 0.9.5,
read as text only, never run) independently implements ``GL(m)T(k)`` as "GL
with an exponential tail on the high binding energy side; larger k fades
faster" and says it measured it against Casa's own curves on ``GL(m)`` only
(1 Oct 2026), not on other shapes. **Now validated**: CasaXPS's own export
of an Au 4d doublet fitted with ``SGL(70)T(1)``, ``T(1.5)``, ``T(2)`` and
``GL(70)T(1)``, ``T(1.5)``, ``T(2)`` (``D:\\Temp\\GL SGL and T tail.csv``) is
reproduced by the table of ``_gl_table`` to 0.03-0.18 % of the peak height
(the formula above drawn on the old analytic shape was 1.2 % out for ``SGL`` and
2.2 % for ``GL``), with the stated Area used as it is; ``is_exact`` is
therefore True for a one-number ``GL``/``SGL`` with a positive ``k``. A
``T()`` on any other shape is still drawn as the plain base.

**LF's ``w`` is in eV, not in line widths (checked, rejected).** The tester's
help text says LF's tails "steepen beyond about w line widths", while this
module takes ``w`` in eV straight from the MATLAB listing (whose own ``F`` is
1, so the two readings only coincide at FWHM = 1). Tested on every real LF
region (5: Ti 2p of the titanium depth profile; Cl 2p, Pt 4d and Pt 4f x2 of
PtCl2_quantified), evaluating ``w * fwhm`` instead of ``w``: ``residual_rms``
got **worse on all five** (Ti 2p 0.0303 -> 0.0312, Cl 2p 0.0523 -> 0.0539,
Pt 4d 0.0686 -> 0.1401, Pt 4f 0.0491 -> 0.0530, Pt 4f area2 0.0607 ->
0.0652), so ``w`` stays in eV and the tester's wording is loose.

**Invalid parameters** (``check_shape``): the tester refuses shapes with
non-finite numbers, LA exponents <= 0 or a negative Gaussian index, a negative
LF ``m``, a DS asymmetry outside [0, 1) or a negative ``n``, and a ``T(k)``
with ``k <= 0``; ``check_shape`` reports the same (plus a few from the
Cookbook's own ranges) and ``casafit.curves`` leaves such a component out of
the envelope and says so, rather than drawing whatever the arithmetic gives.

**Unrecognised shape names.** A shape name CasaXPS writes that this module
does not implement (``QF``, or CasaXPS's own undocumented ``H``/``F``
families -- seen, unreconstructed, on a real file with three refits of one
C 1s region: ``H(0.09,250)SGL(90)`` and ``F(0.09,32,150)SGL(90)`` alongside
the ``DS`` fit below) used to be silently coerced into an exact ``GL(mix)``
using its first numeric token as a 0-100 % mix -- for ``DS(0.09,500)`` (see
below) this drew an almost-pure Gaussian and reported it as *exact*, 12 % of
peak height off and 68 % over the real peak at its own maximum, with no
warning. ``parse_shape`` now keeps the real name (so ``is_exact`` correctly
reports it unreconstructed) and ``_raw_values`` draws it as a plain
unmodified Lorentzian of the component's own ``fwhm`` (the ``a = b = 1``,
``w = 0`` default of the LA/LF branch below, no Gaussian broadening) purely
so something renders -- not a guess at the real shape, just a placeholder
that is honestly flagged. A compound name CasaXPS writes as a base shape
followed by a second one in its own parentheses (``H(0.09,250)SGL(90)``, or
the documented ``DS(a,n)GL(m)``/``DS(a,n)SGL(m)`` blend) parses as the
leading shape with the remainder kept verbatim in ``parse_shape``'s
``suffix`` key rather than failing to parse at all; the blend itself is not
reconstructed (no real file with it has been available), so a ``DS`` with a
suffix draws as plain ``DS``, same precedent as the ``T(k)`` tail.

**DS(a, n): Doniach-Sunjic.** CasaXPS's asymmetric-tail shape for metallic /
graphitic peaks (e.g. HOPG's C 1s). Reconstructed from the published kernel
(Doniach & Sunjic 1970; the form here matches an independent reconstruction
retrieved 2026-09-26 from public papers, itself unvalidated) --
``t = 2(x-pos)/fwhm``, ``DS_raw(t; a) = cos(pi*a/2 + (1-a)*atan(t)) /
(1+4t^2)**((1-a)/2)`` -- convolved with a Gaussian through the same
``gauss_conv`` + area-normalisation pipeline LA/LF already use. This raw
kernel's own asymmetry direction is exact (a fast, power-law-t^-(2-a) decay
on the low-binding-energy side, KE > pos; the well-known long t^-(1-a) tail
on the high-BE side, KE < pos) -- only the Gaussian convolution width is a
calibration.

The convolution width uses the same ``GAUSS_K`` dial as LA/LF but, unlike
them, its own exponent on ``(25/n)`` (``GAUSS_P``, 1 for every other shape):
``gw = fwhm * GAUSS_K["DS"] * (25/n) ** GAUSS_P["DS"]``. First calibrated on
one ``DS(0.09,500)`` component (a plain ``K``, implicit exponent 1) it
overshot the real data specifically on the low-BE side once more real
``(a, n)`` variations of the *same* underlying spectrum became available
(``D:\Temp\for claude files\DS Variations.vms``: 7 CasaXPS refits of one C1s
scan, ``a`` in 0.05-0.10, ``n`` in 200/400/500) -- a symmetric Gaussian
widens a naturally sharp edge (the fast-decaying low-BE side) far more
visibly than it perturbs an already-broad tail, so too large a ``gw`` there
reads exactly as "the low-BE side extends past the data". A per-region
residual-minimising search over ``gw`` on those 7 regions showed the true
optimum scales more weakly with ``n`` than ``25/n`` (exponent 1): fitting
``gw* = fwhm * K' * (25/n)**p`` by ordinary least squares on
``log(gw*/fwhm)`` vs ``log(25/n)`` gives **``GAUSS_K["DS"] = 1.8445``,
``GAUSS_P["DS"] = 0.5737``** (R^2 = 0.886 over the 7 points -- real but
imperfect; ``a`` has a secondary effect on the true optimum, ~5.5 % relative
spread within one ``n`` group, not modelled -- a third parameter is not
worth fitting on 7 points). Against today's numbers (``K=5.6``, implicit
``p=1``) this drops the mean ``residual_rms`` across the 7 regions from
1.49 % to 1.30 % of peak height and the mean low-BE-side overshoot from
5.27 % to 3.83 % (worst case 6.68 % to 6.69 %; one region regresses
~0.03 pp, the rest improve, e.g. ``DS(0.05,200)``'s overshoot 5.68 % to
3.41 %). A ``gw`` with no ``n`` dependence at all was tried and rejected: it
removes the low-BE overshoot almost entirely but then under-broadens the
peak apex at ``n=200``, visibly worsening the whole-curve fit there --
``n`` has to stay in the formula, just with a weaker exponent. **This is
still one calibration spectrum** (now with rich ``(a, n)`` coverage rather
than one point, not an independent second measurement): retune again only
with a ``DS`` fit on genuinely different real data, not on general
principle.

numpy is needed for evaluation; everything else here is plain Python.
"""

from __future__ import annotations

import math
import re

GAUSS_K = {"LF": 0.60, "DS": 1.8445, "A": 1.0, "TLA": 1.0}
GAUSS_P = {"DS": 0.5737, "A": 1.0, "TLA": 1.0}  # exponent on
                                    # (25/m); every other shape is 1 (below)
_LN2_4 = 2.772588722239781               # 4 ln 2


def gaussian_index(kind, fwhm, gauss_fwhm):
    """The CasaXPS Gaussian index ``n`` whose convolution width, for a
    ``kind`` (``"LF"``, ``"DS"``) component of width ``fwhm``,
    equals ``gauss_fwhm`` (same unit): the inverse of the width rule in
    :func:`component_curve`, ``gw = fwhm * GAUSS_K * (25/n) ** GAUSS_P``.
    ``0.0`` (no convolution) when either width is missing or not positive.
    ``LA`` is not a width rule any more (it is the program's own table, see
    ``_la_table``): a source that states an LA's Gaussian width in eV is
    ``LAG(a,b,width/fwhm)`` instead."""
    k = GAUSS_K.get(kind, 0.0)
    if not (k > 0 and fwhm and fwhm > 0 and gauss_fwhm and gauss_fwhm > 0):
        return 0.0
    return 25.0 * (fwhm * k / gauss_fwhm) ** (1.0 / GAUSS_P.get(kind, 1.0))


_TAIL_SUFFIX_RE = re.compile(r"^(.*\))\s*T\(\s*([^()]*?)\s*\)\s*$", re.IGNORECASE)
_SHAPE_RE = re.compile(r"^\s*([A-Za-z]+)\s*\(([^()]*)\)\s*(.*)$")
_A_BASE_RE = re.compile(r"(?i)^(SGL|GL)\(([^()]*)\)")


def parse_shape(text) -> dict:
    """``{"kind", "a", "b", "w", "m", "mix", "tail", "suffix", "params"}``
    from a shape string such as ``GL(30)``, ``LA(1.1,1.9,7)``,
    ``LF(1.1,1.2,75,200)``, ``DS(0.09,500)``, a ``GL``/``SGL`` shape with a
    CasaXPS tail suffix (``GL(30)T(1.5)``, used for asymmetric metallic
    peaks) or a compound shape naming a second one in its own parentheses
    (``H(0.09,250)SGL(90)``, or the documented ``DS(a,n)GL(m)``/
    ``DS(a,n)SGL(m)`` blend). ``tail`` is True when the ``T(k)`` suffix was
    present and ``tk`` is its ``k`` (None when unreadable); the tail is drawn
    for ``GL``/``SGL`` only (see the module docstring).
    ``suffix`` holds a compound name's second shape verbatim, unparsed. A
    name this module does not implement (``QF``, CasaXPS's own ``H``/``F``
    families, ...) keeps its real ``kind`` and stores its raw numeric tokens
    in ``params`` rather than being coerced into a fabricated ``GL(mix)`` --
    see the module docstring. An unreadable string (no ``NAME(...)`` at all)
    is a symmetric Lorentzian-Gaussian mix ("GL(30)")."""
    text = str(text or "")
    tail, tk = False, None
    tm = _TAIL_SUFFIX_RE.match(text)
    if tm:
        text = tm.group(1)
        tail = True
        try:
            tk = float(tm.group(2))       # T(k): the decay (None if unreadable)
        except ValueError:
            tk = None
    m = _SHAPE_RE.match(text)
    if not m:
        return {"kind": "GL", "mix": 30.0, "a": 1.0, "b": 1.0, "w": 0.0,
                "m": 0.0, "tail": tail, "tk": tk, "suffix": "", "params": []}
    kind = m.group(1).upper()
    ps = []
    for tok in m.group(2).split(","):
        try:
            ps.append(float(tok))
        except ValueError:
            pass
    suffix = m.group(3).strip()
    out = {"kind": kind, "a": 1.0, "b": 1.0, "w": 0.0, "m": 0.0, "mix": 0.0,
           "tail": tail, "tk": tk, "suffix": suffix, "params": ps}
    if kind in ("LF", "LSF", "QF", "QSF"):
        out.update(a=ps[0] if ps else 1.0, b=ps[1] if len(ps) > 1 else 1.0,
                   w=ps[2] if len(ps) > 2 else 0.0,
                   m=ps[3] if len(ps) > 3 else 0.0)
    elif kind == "LA":
        if len(ps) >= 3:
            out.update(a=ps[0], b=ps[1], m=ps[2])
        elif len(ps) == 2:
            # CasaXPS Cookbook 2026, p.66: "Abbreviation for LA lineshape:
            # LA(x: alpha, w) = LA(x: alpha, alpha, w)" -- the 2-argument
            # form's own w/n is used directly, unlike the 1-argument
            # shorthand below (see module docstring: this used to be
            # conflated with that shorthand, a real bug on a real file).
            out.update(a=ps[0], b=ps[0], m=ps[1])
        else:
            # The 1-argument shorthand's m (0-100) is NOT the explicit form's
            # own m/n -- CasaXPS defines LA(m) = LA(1,1, 1401 - (m/100)*1401)
            # (Cookbook 2026 p.19/66; Fairley et al., JVST A 41(1) 2023
            # supplementary, Eq. 6): the two run in opposite directions
            # (LA(0) -> n=1401, LA(100) -> n=0).
            raw = min(max(ps[0] if ps else 0.0, 0.0), 100.0)
            out.update(m=float(_LA_MAX_N - math.floor(raw / 100.0 * _LA_MAX_N)))
    elif kind in ("LS", "QA", "QS"):
        # LS(a, b, n) and QA / QS alike; the two-number form is (a, a, n)
        if len(ps) >= 3:
            out.update(a=ps[0], b=ps[1], m=ps[2])
        elif len(ps) == 2:
            out.update(a=ps[0], b=ps[0], m=ps[1])
    elif kind == "LAG":
        # LAG(a, b, r): this module's own name (not a CasaXPS shape) for
        # KherveFitting's LA x G, the shared-width asymmetric Lorentzian
        # convolved with a Gaussian of FWHM r * fwhm (see
        # readers/khervefitting_kfit.py); a and b as for LA.
        out.update(a=ps[0] if ps else 1.0, b=ps[1] if len(ps) > 1 else 1.0,
                   m=ps[2] if len(ps) > 2 else 0.0)
    elif kind == "DS":
        out.update(a=ps[0] if ps else 0.0, m=ps[1] if len(ps) > 1 else 0.0)
    elif kind in ("GL", "SGL", "VOIGT"):
        out["mix"] = ps[0] if ps else 30.0
    elif kind in ("AV", "AVP"):
        # AV(mix, CT, TM, ET, step) / AVP(...): an Avantage fitted peak
        # (readers/thermo_vgp.py). mix = the Lorentzian percentage of the
        # Gaussian/Lorentzian *sum* (AV) or *product* (AVP) function; CT, TM,
        # ET the constant tail height, tail mixing ratio and exponential tail
        # of Avantage's tail function (Data Processing Reference Manual v5,
        # p.73); step = the width of one channel (eV), in which ET is
        # measured. See _avantage_values.
        out["mix"] = ps[0] if ps else 30.0
        out["tail_terms"] = tuple(ps[1:4])
        out["step"] = ps[4] if len(ps) > 4 else 0.0
    elif kind == "A":
        # A(a,b,n)GL(m) / A(a,b,n)SGL(m): the Gelius asymmetric shape (see
        # module docstring). a, b control the tail; n (kept in "m", like DS's
        # own Gaussian-conv code) is the Gaussian convolution width; the
        # trailing GL(m)/SGL(m) is the symmetric base and its own mix.
        out.update(a=ps[0] if ps else 0.2, b=ps[1] if len(ps) > 1 else 0.4,
                   m=ps[2] if len(ps) > 2 else 0.0)
        bm = _A_BASE_RE.match(suffix)
        out["base"] = bm.group(1).upper() if bm else "GL"
        out["mix"] = float(bm.group(2)) if bm and bm.group(2) else 30.0
    elif kind == "TLA":
        # TLA(alpha, mu, n): see the module docstring for how uncertain this
        # parameter order is. a = alpha (Lorentzian power), b = mu (tail
        # strength), m = n (Gaussian-conv code, same convention as DS/A).
        if len(ps) == 2:
            # TLA(mu, n): the program reads it as TLA(1, mu, n)
            out.update(a=1.0, b=ps[0], m=ps[1])
        else:
            out.update(a=ps[0] if ps else 1.0,
                       b=ps[1] if len(ps) > 1 else 20.0,
                       m=ps[2] if len(ps) > 2 else 0.0)
    return out


def is_exact(shape) -> bool:
    """True for shapes that are exact rather than reconstructed: a plain
    ``GL(m)``/``SGL(m)`` and an ``LA`` with equal exponents. A second argument (``GL(m,n)``, the extra Gaussian
    convolution of Cookbook p.67) or anything after the parentheses is not
    drawn, so such a component is not exact whatever its first number."""
    ps = parse_shape(shape)
    if ps["kind"] == "LA":
        # the program's own table (see _la_table), checked against CasaXPS's
        # exported component curves to the printed precision for a == b (the
        # LA(m) and LA(a,m) forms); an unequal pair uses the same code but has
        # no exported curves to compare with, so it stays flagged
        return (not ps["tail"] and not ps["suffix"]
                and abs(ps["a"] - ps["b"]) < 1e-9)
    if ps["kind"] == "DS" and len(ps["params"]) == 3 and not ps["suffix"]:
        # DS(a, n, m): checked against CasaXPS's exported curves (a 0-0.7,
        # n 2-800, m 1-20; a cut-off of 0.5 or less has no positive area and
        # CasaXPS draws nothing, as the table does at an Area of 0)
        return _table_for(ps) is not None
    if ps["kind"] in ("LSF", "QF", "QSF"):
        # the program's table with the boxcar (LSF, QSF) and / or the 16 x^4
        # base (QF, QSF), checked against CasaXPS's exported curves (QF / QSF:
        # a 0.5-1, b 1-12, w 20-50, n 0-100, k 1-12; LSF: a 0.8-1, b 1-1.2,
        # w 10-50, n 0-500)
        return (len(ps["params"]) in (4, 5) and not ps["tail"]
                and not ps["suffix"] and ps["a"] > 0 and ps["b"] > 0
                and ps["w"] > 0)
    if ps["kind"] in ("LS", "QA", "QS"):
        # the program's table with the "Moving Average" kernel (LS, QS) and / or
        # the 16 x^4 base (QA, QS), checked against CasaXPS's exported curves
        # (LS: a, b 0.5-1.5, n 10-500; QA / QS: a 0.5-1, b 1-5, n 0-300)
        return (len(ps["params"]) in (2, 3) and not ps["tail"]
                and not ps["suffix"] and ps["a"] > 0 and ps["b"] > 0)
    if ps["kind"] == "E":
        # checked against CasaXPS's exported curves (a = 1, 0.5 and 0.1)
        bm = _A_BASE_RE.match(ps.get("suffix") or "")
        return len(ps["params"]) == 3 and bool(bm) and not ps["tail"]
    if ps["kind"] in ("GL", "SGL") and ps["suffix"] and not ps["tail"]:
        # GL(m)PHI(j,w) / GL(m,n)PHI(j,w) (and SL): checked against CasaXPS's
        # exported curves for n of 30 and 80
        return (len(ps["params"]) in (1, 2)
                and bool(_PHI_RE.match(ps["suffix"]))
                and _table_for(ps) is not None)
    if ps["kind"] in ("GL", "SGL") and ps["tail"]:
        # GL(m)T(k) / SGL(m)T(k): the program's table with its low-KE tail
        # (see _gl_table), checked against CasaXPS's exported curves
        tk = ps.get("tk")
        return (tk is not None and math.isfinite(tk) and tk > 0
                and len(ps["params"]) == 1 and not ps["suffix"])
    return (ps["kind"] in ("GL", "SGL") and not ps["tail"]
            and len(ps["params"]) <= 1 and not ps["suffix"])


def shape_notes(shape) -> list:
    """Caveats about a shape that is drawn, but not as its source defines it
    (the component stays in the envelope): [] for the shapes that are drawn
    as they are. An Avantage peak (``AV`` / ``AVP``) with tail terms but no
    channel width cannot draw its tail and says so."""
    ps = parse_shape(shape)
    if ps["kind"] in ("AV", "AVP"):
        t = ps.get("tail_terms") or ()
        has_tail = len(t) == 3 and not (abs(t[0]) < 1e-9
                                        and abs(t[1] - 1.0) < 1e-9)
        if has_tail and not _avantage_tail(ps):
            return ["Avantage tail terms are not drawn (the data step is "
                    "unknown)"]
    return []


def area_for_height(shape, fwhm, height):
    """The ``area`` to give :func:`component_curve` so that a ``GL`` / ``SGL``
    / ``AV`` / ``AVP`` component peaks at ``height``: those shapes are 1 at
    their centre and are not broadened further, so it is ``height`` times the
    area :func:`component_curve` measures for a unit-height peak (over its
    own extent, tail included)."""
    np = _np()
    sp = parse_shape(shape)
    fwhm = fwhm if fwhm and fwhm > 0 else 1.0
    span = _NORM_HALF_WIDTH * fwhm
    wide = np.linspace(-span, span, _NORM_POINTS)
    s = float(_raw_values(wide, sp, 0.0, fwhm).sum()) * (wide[1] - wide[0])
    return height * s


def window_integral(shape, pos, fwhm, lo, hi, step=0.004):
    """The integral over ``[lo, hi]`` (eV) of the unit-height shape centred on
    ``pos`` (for the shapes :func:`area_for_height` covers)."""
    np = _np()
    sp = parse_shape(shape)
    fwhm = fwhm if fwhm and fwhm > 0 else 1.0
    x = np.arange(lo, hi, step)
    if len(x) == 0:
        return 0.0
    return float(_raw_values(x, sp, pos, fwhm).sum()) * step


def check_shape(shape) -> list:
    """What is wrong with a shape string's parameters, as short sentences
    ([] when nothing is). Pure, never raises: file content is not rejected at
    parse time; ``casafit.curves`` uses this to leave a component out of the
    envelope instead of drawing arithmetic nonsense. The rules are the ones
    HarwellXPS's Lineshape Tester refuses (non-finite numbers, LA exponents
    <= 0 or a negative Gaussian index, a negative LF ``m``, a DS asymmetry
    outside [0, 1) or a negative ``n``, a ``T(k)`` with ``k <= 0``) plus the
    Cookbook's own ranges for ``GL``/``SGL`` (0-100) and ``TLA`` (alpha, mu >
    0). A name this module does not implement has nothing to check."""
    ps = parse_shape(shape)
    kind, vals, bad = ps["kind"], ps["params"], []
    if any(not math.isfinite(v) for v in vals):
        return ["a parameter is not a finite number"]
    n = len(vals)
    if kind == "LA":
        if n not in (1, 2, 3):
            bad.append("LA needs one, two or three parameters")
        elif n >= 2 and (ps["a"] <= 0 or ps["b"] <= 0):
            bad.append("LA's exponents must be positive")
        if n >= 2 and ps["m"] < 0:
            bad.append("LA's Gaussian index cannot be negative")
    elif kind in ("LS", "QA", "QS"):
        if n not in (2, 3):
            bad.append(f"{kind} needs two or three parameters")
        elif ps["a"] <= 0 or ps["b"] <= 0:
            bad.append(f"{kind}'s exponents must be positive")
        elif ps["m"] < 0:
            bad.append(f"{kind}'s kernel width cannot be negative")
    elif kind in ("LF", "LSF", "QF", "QSF"):
        if n < 4:
            bad.append(f"{kind} needs four (or five) parameters")
        else:
            if ps["a"] <= 0 or ps["b"] <= 0:
                bad.append(f"{kind}'s exponents must be positive")
            if ps["w"] < 0:
                bad.append(f"{kind}'s tail width cannot be negative")
            if ps["m"] < 0:
                bad.append(f"{kind}'s Gaussian parameter m cannot be negative")
            if n >= 5 and not 0 <= vals[4] <= _LF_MAX_K:
                bad.append(f"{kind}'s repeat count must be 0 to {_LF_MAX_K}")
    elif kind == "DS":
        if n < 2:
            bad.append("DS needs two parameters")
        else:
            if not 0 <= ps["a"] < 1:
                bad.append("DS's asymmetry must be at least 0 and below 1")
            if ps["m"] < 0:
                bad.append("DS's Gaussian index cannot be negative")
            if n >= 3 and vals[2] < 0:
                bad.append("DS's cut-off width cannot be negative")
    elif kind in ("GL", "SGL"):
        if n < 1:
            bad.append(f"{kind} needs one parameter")
        elif not 0 <= vals[0] <= 100:
            bad.append(f"{kind}'s mix must be between 0 and 100")
    elif kind == "A":
        if n < 3:
            bad.append("A() needs three parameters")
        elif ps["a"] < 0 or ps["b"] < 0 or ps["m"] < 0:
            bad.append("A()'s parameters cannot be negative")
    elif kind == "TLA":
        if n < 3:
            bad.append("TLA needs three parameters")
        elif ps["a"] <= 0 or ps["b"] <= 0 or ps["m"] < 0:
            bad.append("TLA's alpha and mu must be positive and n not negative")
    elif kind in ("AV", "AVP"):
        if n < 1:
            bad.append(f"{kind} needs a mix")
        elif not 0 <= vals[0] <= 100:
            bad.append(f"{kind}'s mix must be between 0 and 100")
        if n >= 4 and not (vals[1] >= 0 and 0 <= vals[2] <= 1
                           and vals[3] >= 0):
            bad.append(f"{kind}'s tail terms must have CT >= 0, TM between 0 "
                       "and 1 and ET >= 0")
        if n >= 5 and vals[4] < 0:
            bad.append(f"{kind}'s channel width cannot be negative")
    if ps["tail"]:
        tk = ps["tk"]
        if tk is None or not math.isfinite(tk) or tk <= 0:
            bad.append("T(k) needs a positive decay k")
        elif kind not in ("GL", "SGL"):
            bad.append("T(k) is only reconstructed for GL and SGL")
    return bad


def _np():
    import numpy
    return numpy


def gauss_conv(x, y, gfwhm):
    """``y`` convolved with a Gaussian of FWHM ``gfwhm`` (same units as
    ``x``, a uniform grid); edges are extended flat."""
    np = _np()
    if not gfwhm > 0 or len(x) < 3:
        return y
    dx = abs(x[1] - x[0])
    sig = gfwhm / 2.3548200450309493
    n = int(min(200, max(1, round(4 * sig / dx))))
    k = np.exp(-0.5 * (np.arange(-n, n + 1) * dx / sig) ** 2)
    k /= k.sum()
    pad = np.concatenate([np.full(n, y[0]), y, np.full(n, y[-1])])
    return np.convolve(pad, k, mode="valid")


# --- the CasaXPS shape tables (read from the program, see below) -------------
#
# CasaXPS does not evaluate ``LA`` / ``LF`` / ``DS`` as formulas: it tabulates
# each shape once and looks the table up (the routines that fill it, and the
# one that measures it, are in the Ghidra analysis of CasaXPS.exe: ``LA``
# 0x14040dad0, ``LF`` 0x14040ead0, ``DS`` 0x14040c5b0, measure 0x140416980).
# Every constant below is read from there:
#
# * A table of ``npts`` samples (601 for LA / LF, 1203 for DS) on a grid of
#   ``_TABLE_STEP`` (the class default, 0x140a6b0a8 = 0.025); a second, coarse
#   array keeps every ``_TABLE_RATIO``-th point, so the table reaches ``npts``
#   * ``_TABLE_RATIO`` / 2 steps from the centre.
# * ``LA(a, b, n)``: ``y = 1 / (1 + 4 x^2) ** a`` for the low-KE (high-BE) side
#   and ``** b`` for the other side. When ``n * step`` exceeds 2 the Gaussian
#   kernel is held at ``floor(2/step) + 1`` samples and the *Lorentzian's* x
#   axis is stretched by ``n * step / 2`` instead, which is how a large ``n``
#   gives a wide Gaussian.
# * ``GL(m)`` / ``SGL(m)``: the plain analytic shape (``GL`` the product
#   ``exp(-(1-m) 4 ln2 x^2) / (1 + 4 m x^2)``, ``SGL`` the sum), no kernel, area
#   from the whole table. **Its width is measured like every other table's**:
#   the product form's half maximum is not at x = 0.5, so the old assumption
#   of a FWHM of exactly 1 was out by up to 1.3 % of the peak height.
# * ``GL(m)T(k)`` / ``SGL(m)T(k)``: that shape with, on the low-KE (high-BE)
#   side, ``y + (1 - y) exp(k x)`` (x in the table's units, k > 0); the same
#   measured width and area.
# * ``GL(m)PHI(j, w)`` / ``SL(j, w)`` (the two names make the same call, the
#   negative-k slot of the builder): on the low-KE side ``y + j exp(13.8 x / w)
#   (1 - y)`` (13.8 = ln 1e6, so ``w`` is the distance over which the added
#   tail falls to a millionth; ``j`` may exceed 1); ``GL(m, n)`` adds a Gaussian
#   kernel of ``n`` samples (0-1999) after the modifier, and the area is that
#   of the whole table.
# * ``A(a, b, n)GL(m)`` / ``SGL(m)`` (the Gelius asymmetric shape): the same
#   shape on the high-KE side; on the low-KE side
#   ``GL(m) + b (0.7 + 0.3 / (a + 0.01)) (AW(x) - G0(x))`` with
#   ``G0 = exp(-4 ln2 x^2)`` and ``AW = exp(-4 ln2 x^2 / (1 - 2 a sqrt(ln2) x)^2)``
#   (zero where ``|1 - a x|`` is below 1e-100); ``n`` (0 or 2-1199) a Gaussian
#   kernel; area from the right half doubled.
# * ``TLA(a, mu, n)`` (type 7 of the same routine, so the same kernel clamp and
#   stretch): ``y = (pi/2 - atan(2x) + pi/mu) / (1 + 4 x^2) ** a`` on both
#   sides (``TLA(mu, n)`` is ``a`` = 1; ``mu`` 0 drops the tail term).
# * ``LS(a, b, n)`` / ``LS(a, n)`` (b = a): type 1 of the same builder, which
#   takes the kernel by name, ``Moving Average``: a boxcar of ``2 * floor(n/2) +
#   1`` equal samples (the "square Lorentzian" of the literature) in place of
#   the Gaussian; the stretch / clamp at ``n * step > 2`` is the same.
# * ``LSF(a, b, w, n, k)`` / ``QF`` / ``QSF`` (builder 0x14040ead0 types 1, 4,
#   5; ``QF`` and ``QSF`` checked against 13 exported curves, ``LSF`` against 8):
#   ``LSF`` is ``LF`` with the "Moving Average" boxcar in place of the Gaussian;
#   ``QF`` has ``16 x^4`` where ``LF`` has ``4 x^2`` in the base,
#   ``1 / (1 + 16 x^4) ** p``, with the same exponent p; ``QSF`` is ``QF`` with
#   the boxcar. ``k`` = 0 applies the kernel no times.
# * ``QA(a, b, n)`` / ``QA(a, n)`` / ``QS(...)`` (types 4 and 5 of the same
#   builder; checked against 11 exported curves, 0.04-0.63 % of the peak):
#   ``LA`` with ``16 x^4`` where it has ``4 x^2``, ``1 / (1 + 16 x^4)
#   ** a`` (``** b`` on the other side), the same stretch / clamp; ``QA`` has the
#   Gaussian, ``QS`` the boxcar.
# * ``LF(a, b, w, n, k)``: the same with the exponent
#   ``3 - (3 - a) / (1 + 4 x^2 / w^2)`` (``b`` on the other side), no stretch,
#   and the Gaussian applied ``k`` times (``k`` = 1 for the four-number form).
# * ``DS(a, n)``: the Doniach-Sunjic form on the low-KE side
#   ``cos(pi a / 2 + (1 - a) atan(2x)) / (1 + 4x^2) ** ((1 - a) / 2)``, with x
#   shifted by ``-1 / tan(pi / (2 - a))``; with a ``GL(m)`` / ``SGL(m)``
#   suffix the other side is that Gaussian-Lorentzian shape, otherwise the
#   same form on both sides. One Gaussian kernel of ``n`` samples.
# * ``DS(a, n, m)`` (builder type 3; checked against 20 exported curves,
#   a 0-0.7, n 2-800, m 1-20, to 0.1-1.3 % of the peak): the Doniach-Sunjic form on both sides, zeroed
#   below ``m`` (table units) to the low-KE side of the centre, then
#   ``y - (y_cut / I_cut) * I`` where ``I(x)`` is the integral of that table from
#   ``x`` up to the table's high-KE end (Savitzky-Golay window 5, order 2
#   integration kernel, then a running sum), so the shape reaches 0 at the cut.
# * ``H(a, n)GL(m)`` / ``H(a, n)SGL(m)``: the ``DS(a, n)GL(m)`` shape, but its
#   area is taken from the **right half only, doubled** (the symmetric
#   Gaussian-Lorentzian side), so the stated Area is that of the symmetric
#   peak and the table's own area is larger (1.4 times for a = 0.09).
# * ``F(a, b, n)GL(m)`` / ``F(a, b, n)SGL(m)``: on the high-KE side the
#   Gaussian-Lorentzian shape ``GL(m)`` / ``SGL(m)``; on the low-KE side that
#   shape mixed with a power law, ``b/100 / (1 + 4x^2) ** ((1 - a) / 2) +
#   (1 - b/100) * GL(m)``; same right-half area rule.
# * ``E(a, b, n)GL(m)`` / ``SGL(m)`` (builder 0x14040f950; checked against
#   CasaXPS's exported curves, 12 fits including a = 0.5 and 0.1): as
#   ``F``, but the low-KE term is the Gaussian-Lorentzian shape itself raised
#   to a power, ``b/100 * GL(m) ** ((1 - a) / 2) + (1 - b/100) * GL(m)``.
# * A Gaussian kernel of ``2 * floor(n/2) + 1`` samples (``n`` + 1 when even),
#   ``exp(-(3 j / half)^2)`` for j = -half..half, normalised, applied by an
#   edge-extended correlation.
# * The table is divided by its trapezoid area, its height and half-maximum
#   points are found (parabola through the top three samples, linear
#   interpolation at the crossings) and its width is taken as **twice the
#   smaller of the two half-widths** (this is what CasaXPS scales a component's
#   FWHM by: a sum or a larger-side width fitted the real files worse).
#
# ``LA(m)`` is ``LA(1, 1, 1401 - floor(m / 100 * 1401))`` (the program's own
# clamp: m is limited to 0..100 and n to 0..1401).
_TABLE_STEP = 0.025
_TABLE_RATIO = 40
_KERNEL_SPAN = 3.0                       # kernel runs to +-3 "half" widths
_LA_MAX_N = 1401
_LA_POINTS = 601
_DS_POINTS = 1203
_LF_MAX_W = 1000.0
_LF_MAX_N = 1401
_LF_MAX_K = 20


def _table_size(npts):
    return (npts - 1) * _TABLE_RATIO + 1


def _la_n(m):
    """The whole-number Gaussian index CasaXPS uses for ``m`` (its parser
    truncates ``n`` and caps it at 1401); never negative."""
    try:
        v = int(math.floor(float(m)))
    except (TypeError, ValueError, OverflowError):
        return 0
    return max(0, min(_LA_MAX_N, v))


def _measure_table(y, npts):
    """``(width, offset)`` of a table in grid units: twice the smaller
    half-width at half maximum, and where the peak sits relative to the
    centre. ``None`` when no half-maximum point is found."""
    np = _np()
    c = (len(y) - 1) // 2
    half_n = (npts - 1) // 2
    fine = y[c - half_n:c + half_n + 1]
    coarse = y[c - half_n * _TABLE_RATIO:c + half_n * _TABLE_RATIO + 1:
               _TABLE_RATIO]
    mid = half_n
    p = int(np.argmax(fine))
    if fine[p] <= fine[mid]:
        p = mid
    pk, top = float(p), float(fine[p])
    if 0 < p < len(fine) - 1:
        d2 = (fine[p - 1] + fine[p + 1] - 2.0 * fine[p]) * 0.5
        d1 = (fine[p + 1] - fine[p - 1]) * 0.5
        if abs(d2) > 1e-200:
            pk = p - d1 / (2.0 * d2)
            top = float(fine[p]) - d1 * d1 / (4.0 * d2)
    half = top * 0.5

    def side(sign):
        for arr, step in ((fine, _TABLE_STEP),
                          (coarse, _TABLE_STEP * _TABLE_RATIO)):
            seq = arr[p:] if sign > 0 else arr[:p + 1][::-1]
            below = np.nonzero(seq < half)[0]
            if len(below) == 0 or below[0] == 0:
                continue
            k = int(below[0])
            y1, y0 = float(seq[k]), float(seq[k - 1])     # y0 above, y1 below
            frac = (y0 - half) / (y0 - y1) if y0 != y1 else 1.0
            # (the coarse array is read from the same index, as the program
            # does, so it only matters for a table whose half maximum lies
            # beyond the fine range)
            return abs(p + sign * (k - 1 + frac) - pk) * step
        return None

    right, left = side(+1), side(-1)
    if right is None or left is None or min(left, right) <= 0:
        return None
    return 2.0 * min(left, right), (p - mid) * _TABLE_STEP


def _gauss_smooth(y, n, repeats=1):
    """``y`` correlated with the Gaussian kernel of ``n`` samples (made odd),
    ``repeats`` times, the ends extended flat. ``n`` < 2 changes nothing."""
    np = _np()
    if n <= 1:
        return y
    if (n // 2) * 2 == n:
        n += 1
    half = n // 2
    j = np.arange(n) - half
    kern = np.exp(-(_KERNEL_SPAN * j / half) ** 2)
    kern /= kern.sum()
    for _ in range(max(1, int(repeats))):
        pad = np.concatenate([np.full(half, y[0]), y, np.full(half, y[-1])])
        y = np.convolve(pad, kern, mode="valid")
    return y


def _boxcar_smooth(y, n, repeats=1):
    """``y`` correlated with the "Moving Average" kernel of ``n`` samples (made
    odd): all samples equal, normalised; the ends extended flat."""
    np = _np()
    if n <= 1:
        return y
    if (n // 2) * 2 == n:
        n += 1
    half = n // 2
    kern = np.full(n, 1.0 / n)
    for _ in range(max(1, int(repeats))):
        pad = np.concatenate([np.full(half, y[0]), y, np.full(half, y[-1])])
        y = np.convolve(pad, kern, mode="valid")
    return y


def _grid(npts):
    np = _np()
    size = _table_size(npts)
    return np.arange(size), (size - 1) // 2


def _finish_table(y, npts, right_half=False, signed=False):
    """``(y, width, offset, npts)``: the table divided by its trapezoid area,
    with its measured width and peak offset (width 1 when it cannot be
    measured). ``right_half``: the area is that of the right half from the
    centre, doubled (``H`` and ``F``). ``signed``: divide by the area even when
    it is negative (``DS(a, n, m)`` with a very short cut-off: the program does
    the same, its table is then negative and a fit can only give the component
    an Area of 0)."""
    if right_half:
        c = (len(y) - 1) // 2
        area = (float(y[c:].sum()) - 0.5 * float(y[c] + y[-1])) \
            * _TABLE_STEP * 2.0
    else:
        area = (float(y.sum()) - 0.5 * float(y[0] + y[-1])) * _TABLE_STEP
    y = y / area if (area > 0 or (signed and area < 0)) else y
    m = _measure_table(y, npts)
    return (y, m[0] if m else 1.0, m[1] if m else 0.0, npts)


_TABLE_CACHE = {}


def _cached(key, build):
    hit = _TABLE_CACHE.get(key)
    if hit is None:
        if len(_TABLE_CACHE) > 96:
            _TABLE_CACHE.clear()
        hit = _TABLE_CACHE[key] = build()
    return hit


def _la_table(a, b, n, boxcar=False, quartic=False):
    """The table of ``LA(a, b, n)``, ``LS`` (``boxcar``), ``QA`` (``quartic``,
    the ``16 x^4`` base) or ``QS`` (both); see the block comment above."""
    key = (("QS" if boxcar else "QA") if quartic
           else ("LS" if boxcar else "LA"),
           round(float(a), 9), round(float(b), 9), int(n))

    def build():
        np = _np()
        kind_, a_, b_, n_ = key
        scale, ne = 1.0, n_
        if 2.0 < n_ * _TABLE_STEP:
            scale = 2.0 / (n_ * _TABLE_STEP)
            ne = int(math.floor(2.0 / _TABLE_STEP)) + 1
        idx, c = _grid(_LA_POINTS)
        x = _TABLE_STEP * (idx - c) / scale
        lor = 1.0 / (1.0 + (16.0 * x ** 4 if kind_ in ("QA", "QS")
                            else 4.0 * x * x))
        y = np.where(idx > c, lor ** b_, lor ** a_)
        smooth = _boxcar_smooth if kind_ in ("LS", "QS") else _gauss_smooth
        return _finish_table(smooth(y, ne), _LA_POINTS)
    return _cached(key, build)


def _gl_curve(x, mode, mix):
    """The ``GL(m)`` (``mode`` 0, product) or ``SGL(m)`` (1, sum) shape."""
    np = _np()
    if mode == 0:
        return (np.exp(-(1.0 - mix) * _LN2_4 * x * x)
                / (1.0 + mix * x * x * 4.0))
    return (1.0 - mix) * np.exp(-_LN2_4 * x * x) + mix / (1.0 + 4.0 * x * x)


_PHI_DECAY = 13.8                        # 0x140a81650, ln(1e6)
_PHI_RE = re.compile(r"(?i)^(PHI|SL)\(\s*([^(),]+?)\s*,\s*([^(),]+?)\s*\)$")
_GL_MAX_N = 1999


def _gl_table(mode, mix, tail=0.0, phi=(0.0, 0.0), n=0):
    """The table of ``GL(m)`` / ``SGL(m)``; ``mix`` = m / 100. ``tail`` > 0
    is the ``T(k)`` modifier and ``phi`` = (j, w) the ``PHI(j, w)`` / ``SL(j,
    w)`` modifier, both on the low-KE (high-BE) side only: ``y + (1 - y)
    exp(k x)`` and ``y + j exp(13.8 x / w) (1 - y)``. ``n`` > 1 is the
    Gaussian kernel (``GL(m, n)``), applied after the modifier."""
    key = ("GL", int(mode), round(float(mix), 9), round(float(tail), 9),
           round(float(phi[0]), 9), round(float(phi[1]), 9), int(n))

    def build():
        np = _np()
        _, mode_, mix_, tail_, j_, w_, n_ = key
        idx, c = _grid(_LA_POINTS)
        x = _TABLE_STEP * (idx - c)
        y = _gl_curve(x, mode_, mix_)
        if tail_ > 0.0:
            y = np.where(idx < c, y + (1.0 - y) * np.exp(x * tail_), y)
        elif j_ > 0.0 and w_ > 1e-20:
            fall = np.exp(np.minimum(x, 0.0) * _PHI_DECAY / w_)
            y = np.where(idx < c, y + j_ * fall * (1.0 - y), y)
        return _finish_table(_gauss_smooth(y, n_), _LA_POINTS)
    return _cached(key, build)


def _e_table(a, b, n, mode, mix):
    """The table of ``E(a, b, n)GL(m)`` (``mode`` 0) / ``SGL(m)`` (1)."""
    key = ("E", round(float(a), 9), round(float(b), 9), int(n), int(mode),
           round(float(mix), 9))

    def build():
        np = _np()
        _, a_, b_, n_, mode_, mix_ = key
        idx, c = _grid(_LA_POINTS)
        x = _TABLE_STEP * (idx - c)
        gl = _gl_curve(x, mode_, mix_)
        w = b_ / 100.0
        left = w * gl ** ((1.0 - a_) * 0.5) + (1.0 - w) * gl
        y = np.where(idx > c, gl, left)
        return _finish_table(_gauss_smooth(y, n_), _LA_POINTS, True)
    return _cached(key, build)


def _a_table(a, b, n, mode, mix):
    """The table of ``A(a, b, n)GL(m)`` (``mode`` 0) / ``SGL(m)`` (1)."""
    key = ("A", round(float(a), 9), round(float(b), 9), int(n), int(mode),
           round(float(mix), 9))

    def build():
        np = _np()
        _, a_, b_, n_, mode_, mix_ = key
        idx, c = _grid(_LA_POINTS)
        x = _TABLE_STEP * (idx - c)
        gl = _gl_curve(x, mode_, mix_)
        g0 = np.exp(-_LN2_4 * x * x)
        den = 1.0 - a_ * 2.0 * math.sqrt(math.log(2.0)) * x
        ok = np.abs(1.0 - a_ * x) > 1e-100
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            aw = np.where(ok, np.exp(-_LN2_4 * x * x / (den * den)), 0.0)
        aw = np.nan_to_num(aw, nan=0.0, posinf=0.0)
        w = b_ * (0.7 + 0.3 / (a_ + 0.01))
        y = np.where(idx > c, gl, gl + w * (aw - g0))
        return _finish_table(_gauss_smooth(y, n_), _LA_POINTS, True)
    return _cached(key, build)


def _tla_table(a, mu, n):
    """The table of ``TLA(a, mu, n)``; see the block comment above."""
    key = ("TLA", round(float(a), 9), round(float(mu), 9), int(n))

    def build():
        np = _np()
        _, a_, mu_, n_ = key
        scale, ne = 1.0, n_
        if 2.0 < n_ * _TABLE_STEP:
            scale = 2.0 / (n_ * _TABLE_STEP)
            ne = int(math.floor(2.0 / _TABLE_STEP)) + 1
        idx, c = _grid(_LA_POINTS)
        x = _TABLE_STEP * (idx - c) / scale
        tail = math.pi / 2.0 - np.arctan(2.0 * x)
        tail = tail + math.pi / mu_ if mu_ > 1e-200 else 1.0 + 0.0 * x
        y = tail / (1.0 + 4.0 * x * x) ** a_
        return _finish_table(_gauss_smooth(y, ne), _LA_POINTS)
    return _cached(key, build)


def _lf_table(a, b, w, n, k, quartic=False, boxcar=False):
    """The table of ``LF(a, b, w, n, k)`` (``w`` in the table's own units);
    ``quartic`` gives ``QF`` (``16 x^4`` in the base), ``boxcar`` the ``S``
    variants (``LSF``, ``QSF``)."""
    key = ("LF", round(float(a), 9), round(float(b), 9), round(float(w), 9),
           int(n), int(k), bool(quartic), bool(boxcar))

    def build():
        np = _np()
        _, a_, b_, w_, n_, k_, quartic_, boxcar_ = key
        idx, c = _grid(_LA_POINTS)
        x = _TABLE_STEP * (idx - c)
        suppress = 1.0 / (1.0 + 4.0 * x * x / (w_ * w_))
        side = np.where(idx > c, b_, a_)
        base = 1.0 + (16.0 * x ** 4 if quartic_ else 4.0 * x * x)
        y = 1.0 / base ** (3.0 - (3.0 - side) * suppress)
        if k_ > 0:                       # the program applies the kernel k times
            smooth = _boxcar_smooth if boxcar_ else _gauss_smooth
            y = smooth(y, n_, k_)
        return _finish_table(y, _LA_POINTS)
    return _cached(key, build)


def _sg_integration_kernel():
    """The Savitzky-Golay (window 5, order 2) kernel of the program's
    "integration" filter: ``c0 - c1/2 + c2/3`` of the least-squares
    polynomial coefficients, i.e. the area under the local quadratic over the
    sample interval before the window centre."""
    np = _np()
    t = np.arange(-2, 3, dtype=float)
    design = np.vander(t, 3, increasing=True)          # 1, t, t^2
    w = np.linalg.pinv(design)                          # rows: c0, c1, c2 kernels
    return w[0] - w[1] / 2.0 + w[2] / 3.0


def _ds_cut(y, m):
    """The ``DS(a, n, m)`` block: ``y`` zeroed below ``m`` table units to the
    low-KE side of the centre, then ``y - (y_cut / I_cut) * I``."""
    np = _np()
    size = len(y)
    c = (size - 1) // 2
    cut = c - int(math.floor(m / _TABLE_STEP))
    cut = 0 if cut < 0 else (size - 1 if cut >= size else cut)
    a = np.where(np.arange(size) < cut, 0.0, y)
    # integral from each point up to the high-KE end, by integrating the
    # reversed table: correlate with the kernel, then a running sum x step
    rev = a[::-1]
    kern = _sg_integration_kernel()
    pad = np.concatenate([np.full(2, rev[0]), rev, np.full(2, rev[-1])])
    filt = np.correlate(pad, kern, mode="valid")
    running = np.concatenate([[0.0], np.cumsum(filt[1:]) * _TABLE_STEP])
    integ = running[::-1]
    if integ[cut] == 0.0:
        return a
    out = a - (a[cut] / integ[cut]) * integ
    out[:cut] = 0.0
    return out


def _ds_table(a, n, mode, mix, right_half=False, cutoff=0.0):
    """The table of ``DS(a, n)`` (``mode`` 2), ``DS(a, n)GL(m)`` (0) or
    ``DS(a, n)SGL(m)`` (1), ``mix`` = m / 100; ``right_half`` makes it the
    ``H`` shape."""
    key = ("DS", round(float(a), 9), int(n), int(mode), round(float(mix), 9),
           bool(right_half), round(float(cutoff), 9))

    def build():
        np = _np()
        _, a_, n_, mode_, mix_, half_, cut_ = key
        idx, c = _grid(_DS_POINTS)
        shift = 0.0
        if 1e-8 < a_ < 0.9999999:
            shift = -1.0 / math.tan(math.pi / (2.0 - a_))
        x = _TABLE_STEP * (idx - c) - shift
        ds = (np.cos(math.pi * a_ * 0.5 + (1.0 - a_) * np.arctan(2.0 * x))
              / (1.0 + 4.0 * x * x) ** ((1.0 - a_) * 0.5))
        if mode_ in (2, 3):
            y = _ds_cut(ds, cut_) if mode_ == 3 else ds
        elif mode_ == 0:
            right = (np.exp(-(1.0 - mix_) * _LN2_4 * x * x)
                     / (1.0 + mix_ * x * x * 4.0))
            y = np.where(idx > c, right, ds)
        else:
            right = ((1.0 - mix_) * np.exp(-_LN2_4 * x * x)
                     + mix_ / (1.0 + 4.0 * x * x))
            y = np.where(idx > c, right, ds)
        return _finish_table(_gauss_smooth(y, n_), _DS_POINTS, half_,
                             signed=(mode_ == 3))
    return _cached(key, build)


def _f_table(a, b, n, mode, mix):
    """The table of ``F(a, b, n)GL(m)`` (``mode`` 0) / ``SGL(m)`` (1); ``b``
    is the per-cent weight of the power-law term, ``mix`` = m / 100."""
    key = ("F", round(float(a), 9), round(float(b), 9), int(n), int(mode),
           round(float(mix), 9))

    def build():
        np = _np()
        _, a_, b_, n_, mode_, mix_ = key
        idx, c = _grid(_LA_POINTS)
        x = _TABLE_STEP * (idx - c)
        if mode_ == 0:
            gl = (np.exp(-(1.0 - mix_) * _LN2_4 * x * x)
                  / (1.0 + mix_ * x * x * 4.0))
        else:
            gl = ((1.0 - mix_) * np.exp(-_LN2_4 * x * x)
                  + mix_ / (1.0 + 4.0 * x * x))
        w = b_ / 100.0
        left = w / (1.0 + 4.0 * x * x) ** ((1.0 - a_) * 0.5) + (1.0 - w) * gl
        y = np.where(idx > c, gl, left)
        return _finish_table(_gauss_smooth(y, n_), _LA_POINTS, True)
    return _cached(key, build)


def _table_values(table, x, pos, fwhm, area):
    """A component of the given ``area`` on ``x``: the table looked up at
    ``(x - pos) / fwhm`` in units of its own measured width, so its area is
    the table's own (1) times ``area``. The **position is the table's centre,
    not its maximum**: CasaXPS's exported LF and DS curves peak 0.05-0.10 eV
    away from the stated Position, and re-centring on the maximum (the
    measured ``offset``, kept in the table but not used) put them 9 % of the
    peak out against 0.5-1.3 % without it."""
    np = _np()
    y, width, _off, npts = table
    u = np.asarray((np.asarray(x, dtype=float) - pos) / fwhm * width)
    c = (len(y) - 1) // 2
    half = (npts - 1) // 2
    # The program keeps each table twice: a fine array (``npts`` samples
    # ``_TABLE_STEP`` apart) and a coarse one (``_TABLE_RATIO`` times wider
    # spacing, so reaching ``_TABLE_RATIO`` times further), and reads the fine
    # one only inside its own range. Beyond it the coarse array is read with
    # linear interpolation, which rounds off a steep shoulder that falls there
    # (the edge of a wide boxcar): CasaXPS's own LSF(1,1,50,500) shows exactly
    # that, 1.8 % of the peak at +-7.5 table units, against 0.1 % elsewhere.
    fine = y[c - half:c + half + 1]
    coarse = y[c - half * _TABLE_RATIO:c + half * _TABLE_RATIO + 1:_TABLE_RATIO]
    fine_u = _TABLE_STEP * np.arange(-half, half + 1)
    coarse_u = _TABLE_STEP * _TABLE_RATIO * np.arange(-half, half + 1)
    inside = np.abs(u) <= _TABLE_STEP * half
    # a position beyond the coarse array (hundreds of FWHM away) counts as
    # zero (the program's own lookup past its table was not read)
    out = np.where(inside, np.interp(u, fine_u, fine),
                   np.interp(u, coarse_u, coarse, left=0.0, right=0.0))
    return out * (area * width / fwhm)


def _table_for(sp):
    """The table of a parsed shape, or ``None`` for a shape that is not (or
    not yet) read from the program."""
    kind = sp["kind"]
    if kind == "LA":
        return _la_table(sp["a"], sp["b"], _la_n(sp["m"]))
    if kind in ("LS", "QA", "QS") and sp["a"] > 0 and sp["b"] > 0 \
            and len(sp["params"]) in (2, 3):
        return _la_table(sp["a"], sp["b"], _la_n(sp["m"]), kind in ("LS", "QS"),
                         kind in ("QA", "QS"))
    if kind in ("GL", "SGL") and len(sp["params"]) in (1, 2) \
            and 0.0 <= sp["mix"] <= 100.0:
        # GL(m), GL(m, n), each with an optional T(k) or PHI(j, w) / SL(j, w)
        n = 0
        if len(sp["params"]) == 2:
            n = int(math.floor(sp["params"][1]))
            if not 0 <= n <= _GL_MAX_N:
                return None
        mode, mix = (0 if kind == "GL" else 1), sp["mix"] / 100.0
        tk, phi = 0.0, (0.0, 0.0)
        if sp["tail"]:
            tk = sp.get("tk")
            if tk is None or not math.isfinite(tk) or tk < 0 or sp["suffix"]:
                return None
        elif sp["suffix"]:
            pm = _PHI_RE.match(sp["suffix"])
            if not pm:
                return None
            try:
                phi = (float(pm.group(2)), float(pm.group(3)))
            except ValueError:
                return None
            if not (math.isfinite(phi[0]) and math.isfinite(phi[1])
                    and phi[0] >= 0 and phi[1] > 0):
                return None
        return _gl_table(mode, mix, tk, phi, n)
    if kind == "A" and 0.0 <= sp["a"] <= 1.0 and 0.0 <= sp["b"] <= 1.0 \
            and 0.0 <= sp["mix"] <= 100.0 and len(sp["params"]) == 3:
        n = _la_n(sp["m"])
        n = n if 2 <= n <= 1199 else 0
        return _a_table(sp["a"], sp["b"], n, 0 if sp["base"] == "GL" else 1,
                        sp["mix"] / 100.0)
    if kind == "TLA" and sp["a"] > 0 and sp["b"] >= 0:
        return _tla_table(sp["a"], sp["b"], _la_n(sp["m"]))
    if kind in ("LF", "LSF", "QF", "QSF") and sp["w"] > 0 \
            and sp["a"] > 0 and sp["b"] > 0 and len(sp["params"]) in (4, 5):
        ps = sp["params"]
        k = int(ps[4]) if len(ps) > 4 else 1
        return _lf_table(sp["a"], sp["b"], min(sp["w"], _LF_MAX_W),
                         min(_la_n(sp["m"]), _LF_MAX_N),
                         min(max(k, 0), _LF_MAX_K),
                         kind in ("QF", "QSF"), kind in ("LSF", "QSF"))
    if kind in ("H", "F", "E"):
        # H(a, n)GL(m) / H(a, n)SGL(m); F(a, b, n)GL(m) / F(a, b, n)SGL(m);
        # E(a, b, n)GL(m) / E(a, b, n)SGL(m)
        ps = sp["params"]
        if len(ps) != (2 if kind == "H" else 3):
            return None
        n = _la_n(ps[-1])
        n = n if 2 <= n <= 1400 else 0
        bm = _A_BASE_RE.match(sp.get("suffix") or "")
        if not (bm and bm.group(2)):
            return None
        try:
            mix = float(bm.group(2)) / 100.0
        except ValueError:
            return None
        mode = 0 if bm.group(1).upper() == "GL" else 1
        if not (0.0 <= mix <= 1.0 and 0.0 <= ps[0] <= 1.0):
            return None
        if kind == "H":
            return _ds_table(ps[0], n, mode, mix, True)
        if not 0.0 <= ps[1] <= 100.0:
            return None
        if kind == "E":
            return _e_table(ps[0], ps[1], n, mode, mix)
        return _f_table(ps[0], ps[1], n, mode, mix)
    if kind == "DS" and len(sp["params"]) == 3 and not sp["suffix"]:
        # DS(a, n, m): the finite-extent form; a in [0, 1], n 0 or 2-1400, m >= 0
        a, n, m = sp["params"]
        n = int(math.floor(n))
        # (the program also accepts m = 0, which cuts at the centre and leaves
        # a table of negative area; that case is left on the older path)
        if 0.0 <= a <= 1.0 and m > 0.0 and math.isfinite(m) \
                and (n == 0 or 2 <= n <= 1400):
            return _ds_table(a, n, 3, 0.0, False, m)
        return None
    if kind == "DS" and len(sp["params"]) in (1, 2):
        # DS(a, n)GL(m) / DS(a, n)SGL(m) / DS(a, n) / DS(a); n of 1 or past
        # 1400 is refused by the program, read here as none
        n = _la_n(sp["m"])
        n = n if 2 <= n <= 1400 else 0
        suffix = sp.get("suffix") or ""
        bm = _A_BASE_RE.match(suffix)
        if bm and bm.group(2):
            try:
                mix = float(bm.group(2)) / 100.0
            except ValueError:
                return None
            if not 0.0 <= mix <= 1.0:
                return None
            return _ds_table(sp["a"], n, 0 if bm.group(1).upper() == "GL" else 1,
                             mix)
        if not suffix:
            return _ds_table(sp["a"], n, 2, 0.0)
    return None


def _half_width_term(p):
    """``sqrt(2**(1/p) - 1)``: how far (in units of the shared Lorentzian
    width ``F``) ``1/(1+4t^2)**p`` must go to reach half its maximum. ``1``
    at ``p = 1``, so a symmetric shape's ``F`` below is exactly ``fwhm``."""
    if p <= 0:
        return 1.0
    return math.sqrt(max(1e-12, 2.0 ** (1.0 / p) - 1.0))


def _shared_width(fwhm, a, b):
    """The one Lorentzian width shared by both sides of an ``LA``/``LF``
    shape, chosen so the two sides' own half-max distances add up to the
    file's stated ``fwhm`` (KherveFitting's ``LA`` formula; see the module
    docstring). Reduces to ``fwhm`` only for an unmodified Lorentzian
    (``a == b == 1``); a symmetric but non-unity pair is rescaled too."""
    denom = _half_width_term(a) + _half_width_term(b)
    return 2.0 * fwhm / denom if denom > 0 else fwhm


def _gl_sgl_values(x, kind, mix, pos, fwhm):
    """The plain symmetric GL(mix) product or SGL(mix) sum, no broadening.
    Shared by the standalone GL/SGL shapes and, as the base term, by the
    Gelius asymmetric shape (``kind == "A"``) below."""
    np = _np()
    t = (x - pos) / fwhm
    lor = 1.0 / (1.0 + 4.0 * t * t)
    m = min(max(mix, 0.0), 100.0) / 100.0
    gau = np.exp(-_LN2_4 * t * t)
    if kind == "GL":
        return (lor ** m) * (gau ** (1.0 - m))
    return m * lor + (1.0 - m) * gau


def _gelius_tail(x, pos, fwhm, a, b):
    """The Gelius asymmetric shape's one-sided tail term,
    ``w(a,b)*[AW(x) - G0(x)]`` on the low-KE (high binding energy) side,
    zero elsewhere (see the module docstring: CasaXPS "Peak Fitting in
    XPS", 2006, p.20)."""
    np = _np()
    t = x - pos
    w = b * (0.7 + 0.3 / (a + 0.01))
    u = 2.0 * math.sqrt(math.log(2.0)) * np.where(t < 0, -t, 0.0)
    aw = np.exp(-(u / (fwhm + a * u)) ** 2)
    g0 = np.exp(-(u / fwhm) ** 2)
    return np.where(t < 0, w * (aw - g0), 0.0)


def _tla_width_scale(alpha, mu):
    """FWHM (in Lorentzian-width units) of the unconvolved CasaXPS ``TLA``
    base shape, found numerically since it has no closed form (mirrors
    KherveFitting's own ``_tla_width_scale``; see the module docstring for
    how uncertain this reconstruction is)."""
    np = _np()
    t = np.linspace(-40.0, 40.0, 8001)
    lor = (1.0 / (1.0 + 4.0 * t * t)) ** alpha
    tail = math.pi / 2.0 - np.arctan(2.0 * t) + math.pi / mu
    y = lor * tail
    above = np.nonzero(y >= float(y.max()) / 2.0)[0]
    if len(above) >= 2:
        return float(t[above[-1]] - t[above[0]])
    return 1.0


def _tla_values(x, pos, fwhm, alpha, mu):
    """The CasaXPS ``TLA`` raw shape (before Gaussian convolution, applied
    afterwards by ``component_curve`` like every other reconstructed shape
    here): a Lorentzian raised to ``alpha``, modulated by an arctangent tail
    factor. See the module docstring -- this shape is substantially less
    certain than the others in this module."""
    np = _np()
    alpha = max(float(alpha), 0.05)
    mu = max(float(mu), 1e-3)
    f_lor = fwhm / _tla_width_scale(alpha, mu)
    u = (x - pos) / f_lor
    lor = (1.0 / (1.0 + 4.0 * u * u)) ** alpha
    tail = math.pi / 2.0 - np.arctan(2.0 * u) + math.pi / mu
    return lor * tail


def _voigt_fwhm_split(fwhm, fraction):
    """``(Gaussian FWHM, Lorentzian FWHM)`` for a total Voigt FWHM and an
    L/G fraction (per cent Lorentzian of the total width), via the analytic
    Olivero-Longbothum inverse -- an exact match to KherveFitting's own
    ``voigt_fwhm_split`` (``Peak_Functions.py``, retrieved 2026-09-27): the
    forward relation ``F = 0.5346*f_l + sqrt(0.2166*f_l**2 + f_g**2)`` is
    linear in the overall scale, so splitting a stated total ``fwhm`` by a
    fraction ``r`` is analytic (``f_g = (1-r)*scale``, ``f_l = r*scale``,
    ``scale = fwhm/k``)."""
    r = min(max(fraction / 100.0, 0.0), 0.999)
    k = 0.5346 * r + math.sqrt(0.2166 * r * r + (1.0 - r) ** 2)
    scale = fwhm / k if k > 0 else fwhm
    return (1.0 - r) * scale, r * scale


def _voigt_values(x, pos, fwhm, fraction):
    """A true Voigt (Gaussian convolved with Lorentzian), split from its
    total FWHM and L/G fraction the same way KherveFitting's own
    ``voigt_simple`` does. Not a CasaXPS shape string -- ``"VOIGT(fraction)"``
    is this module's own encoding, used only by the KherveFitting ``.kfit``
    reader (CasaXPS itself has no true-Voigt shape string; its own
    Voigt-like shapes are the ``GL``/``SGL`` pseudo-Voigt already above).
    Computed as a plain Lorentzian of the split Lorentzian FWHM, Gaussian-
    convolved by the split Gaussian FWHM -- the same convolution
    KherveFitting evaluates analytically via the Faddeeva function
    (``lmfit``'s ``voigt``), not an approximation of a different kind."""
    np = _np()
    f_g, f_l = _voigt_fwhm_split(fwhm, fraction)
    f_l = f_l if f_l > 0 else 1e-6
    t = (x - pos) / f_l
    lor = 1.0 / (1.0 + 4.0 * t * t)
    return gauss_conv(x, lor, f_g)


def _avantage_tail(sp) -> bool:
    """True when an ``AV`` / ``AVP`` shape has a tail that can be drawn: the
    terms are not the "no tail" values ``CT = 0, TM = 1`` and the channel width
    is known."""
    t = sp.get("tail_terms") or ()
    if len(t) != 3 or not sp.get("step", 0.0) > 0:
        return False
    return not (abs(t[0]) < 1e-9 and abs(t[1] - 1.0) < 1e-9)


def _avantage_values(x, sp, pos, fwhm):
    """An Avantage peak, unit height at ``pos``, from the Data Processing
    Reference Manual (v5, pp.72-73): ``GL`` is the Gaussian/Lorentzian
    *product* ``1 / ([1 + M t^2] exp{(1 - M) ln2 t^2})`` with ``t = 2 (x - pos)
    / FWHM`` (``AVP``) or the *sum* ``M L + (1 - M) G`` (``AV``), ``M`` the
    L/G mix / 100. A tail is added on the low-kinetic-energy (high binding
    energy) side only: ``GL + (1 - GL) T`` with ``T = TM CT + (1 - TM)
    exp(-ET |dx|)``, ``dx`` the distance from the centre **in channels**
    (data points, ``step`` eV each); the other side is ``GL``."""
    np = _np()
    t = (x - pos) / fwhm
    m = min(max(sp["mix"], 0.0), 100.0) / 100.0
    if sp["kind"] == "AVP":
        g = np.exp(-_LN2_4 * (1.0 - m) * t * t) / (1.0 + 4.0 * m * t * t)
    else:
        g = m / (1.0 + 4.0 * t * t) + (1.0 - m) * np.exp(-_LN2_4 * t * t)
    if not _avantage_tail(sp):
        return g
    ct, tm, et = sp["tail_terms"]
    d = np.where(x < pos, pos - x, 0.0) / sp["step"]
    tail = tm * ct + (1.0 - tm) * np.exp(-et * d)
    return np.where(x < pos, g + (1.0 - g) * tail, g)


def _raw_values(x, sp, pos, fwhm):
    """The lineshape before any Gaussian broadening (GL/SGL/VOIGT have their
    own convolution already folded in): the GL product, the SGL sum, the
    true Voigt, the Gelius asymmetric shape, the TLA shape, the raw
    asymmetric Lorentzian, or (when ``w > 0``) its LF finite-tail variant, on
    grid ``x``. Used both for the values a caller asked for and, on a
    separate wide grid, for area normalisation."""
    np = _np()
    if sp["kind"] in ("GL", "SGL"):
        base = _gl_sgl_values(x, sp["kind"], sp["mix"], pos, fwhm)
        tk = sp.get("tk")
        if sp.get("tail") and tk is not None and math.isfinite(tk) and tk > 0:
            # Cookbook p.67: base + (1 - base) * E(1, k), E = exp(k x) for
            # x < 0 (x in FWHM, low KE = high binding energy), else 0
            t = np.minimum((x - pos) / fwhm, 0.0)
            return base + (1.0 - base) * np.where(t < 0, np.exp(tk * t), 0.0)
        return base
    if sp["kind"] == "VOIGT":
        return _voigt_values(x, pos, fwhm, sp["mix"])
    if sp["kind"] in ("AV", "AVP"):
        return _avantage_values(x, sp, pos, fwhm)
    if sp["kind"] == "DS":
        t = 2.0 * (x - pos) / fwhm
        a = sp["a"]
        return np.cos(np.pi * a / 2.0 + (1.0 - a) * np.arctan(t)) \
            / (1.0 + 4.0 * t * t) ** ((1.0 - a) / 2.0)
    if sp["kind"] == "A":
        base = _gl_sgl_values(x, sp["base"], sp["mix"], pos, fwhm)
        return base + _gelius_tail(x, pos, fwhm, sp["a"], sp["b"])
    if sp["kind"] == "TLA":
        return _tla_values(x, pos, fwhm, sp["a"], sp["b"])
    if sp["kind"] == "LF" and sp["w"] > 0:
        # The finite tail (see the module docstring): CasaXPS's own formula
        # (Major, Shah, Avval, Fernandez, Fairley, Linford, "Advanced Line
        # Shapes in XPS II: The Finite Lorentzian (LF) Line Shape", VT&C,
        # April 2020, MATLAB listing) is not a hard cutoff -- the base
        # Lorentzian's own exponent smoothly rises from a/b at the centre
        # towards a fixed ceiling of 3.0 as |x - pos| grows relative to w
        # (used directly, in eV -- not scaled by fwhm), so the tail decays
        # by a power law rather than being clipped to exactly zero. The base
        # Lorentzian itself uses plain fwhm here, not the shared-width
        # rescale plain LA/LF(w=0) use below (the source's own MATLAB uses a
        # single fixed "F = 1; %peak width" with no rescaling formula).
        t = (x - pos) / fwhm
        lor = 1.0 / (1.0 + 4.0 * t * t)
        u = (x - pos) / sp["w"]
        suppress = 1.0 / (1.0 + 4.0 * u * u)
        exL = 3.0 - (3.0 - sp["a"]) * suppress
        exR = 3.0 - (3.0 - sp["b"]) * suppress
        return np.where(t < 0, lor ** exL, lor ** exR)
    F = _shared_width(fwhm, sp["a"], sp["b"])
    t = (x - pos) / F
    lor = 1.0 / (1.0 + 4.0 * t * t)
    return np.where(t < 0, lor ** sp["a"], lor ** sp["b"])


# How far a component's own area is measured, in FWHM either side of its
# position: independent of the (possibly narrower) CasaXPS fit-region window
# it is drawn on, so a tail cut off by a tight region box is not mistaken for
# missing area and used to inflate the visible peak.
_NORM_HALF_WIDTH = 100.0
_NORM_POINTS = 4001


def component_curve(ke, shape, pos, fwhm, area):
    """The component on the ascending KE grid ``ke`` (numpy array), scaled so
    that its integral over its own full extent — not just over ``ke``, which
    may be a CasaXPS region window narrower than the shape's tail — is
    ``area``. A component whose tail runs past the displayed window is
    therefore not inflated to make up area that simply isn't in the window."""
    np = _np()
    ke = np.asarray(ke, dtype=float)
    sp = parse_shape(shape)
    fwhm = fwhm if fwhm and fwhm > 0 else 1.0
    table = _table_for(sp)
    if table is not None:
        # LA / LF / DS: a table lookup, not a formula (see _la_table)
        return _table_values(table, ke, pos, fwhm, area)
    v = _raw_values(ke, sp, pos, fwhm)
    if sp["kind"] in ("GL", "SGL", "VOIGT", "AV", "AVP"):
        conv = v
    else:
        if sp["kind"] == "LAG":
            gw = fwhm * max(sp["m"], 0.0)          # the width, stated directly
        else:
            gw = fwhm * GAUSS_K.get(sp["kind"], 0.0) * (
                (25.0 / sp["m"]) ** GAUSS_P.get(sp["kind"], 1.0)
                if sp["m"] else 0.0)
        conv = gauss_conv(ke, v, gw)
    span = _NORM_HALF_WIDTH * fwhm
    n_points = _NORM_POINTS
    if sp["kind"] == "LF" and sp["w"] > 0 and 20.0 * sp["w"] > span:
        # A real LF's w (an eV-scale width, not a multiple of fwhm -- see
        # _raw_values) can exceed _NORM_HALF_WIDTH*fwhm for a narrow peak
        # with a wide tail (w up to ~280 eV seen on real files against a
        # fwhm of a few eV); widen the normalisation grid so it still spans
        # well past where the exponent has settled near its 3.0 ceiling --
        # and grow the point count in proportion, or the (much sparser)
        # wide grid under-resolves the still-narrow peak at its centre and
        # silently mis-normalises the whole component (checked: at a fixed
        # 4001 points, widening the span alone from 130 to 1500 eV for a
        # fwhm=1.3 component leaves only ~2 points across the peak's own
        # FWHM, giving a component ~1.3% too tall).
        newspan = 20.0 * sp["w"]
        n_points = min(200001, int(_NORM_POINTS * newspan / span))
        span = newspan
    tk = sp.get("tk")
    if (sp.get("tail") and sp["kind"] in ("GL", "SGL") and tk and tk > 0
            and 15.0 * fwhm / tk > span):
        # a slow T(k) tail (small k) reaches past the default grid: widen it,
        # with the points in proportion, so its area is not cut off
        newspan = 15.0 * fwhm / tk
        n_points = min(200001, int(_NORM_POINTS * newspan / span))
        span = newspan
    wide = np.linspace(pos - span, pos + span, n_points)
    wstep = wide[1] - wide[0]
    s = float(_raw_values(wide, sp, pos, fwhm).sum()) * wstep
    return conv * (area / s if s > 0 else 0.0)


def shirley(y, avg=1, iters=100, st_offset=0.0, end_offset=0.0):
    """Iterative Shirley background of ``y`` on an ascending-KE slice (the
    low-KE end is the high-binding-energy side, which sits higher). ``avg``
    points are averaged at each end for the end levels, then ``st_offset``/
    ``end_offset`` (percent, CasaXPS's own "St Offset"/"End Offset" region
    parameters -- the region line's ``params[0]``/``params[1]``) reduce those
    two levels by that percentage before the iteration: ``I_used = I_natural *
    (1 - offset / 100)``, confirmed against CasaXPS's own "Peak Fitting in
    XPS" whitepaper (Casa Software Ltd, 2006, p.4: "these intensities can be
    further adjusted using the St Offset and End Offset parameters to reduce
    the intensities required to calculate the background beneath the peaks...
    A value of '0' means no offset while '100' means the background at the
    end point is zero") and, independently, against a real file
    (``D:\\Temp\\for claude files\\PtCl2_new\\PtCl2_refitted.vms``): its
    ``Pt 4f`` region's ``St Offset`` of 11.427384% lowers ``residual_rms``
    against the real raw counts (0.03663 -> 0.03552) when applied, rather
    than left at the raw data mean as before. ``st_offset`` is the low-KE/
    high-BE end (``y[:avg]``, CasaXPS's "E1"), ``end_offset`` the high-KE/
    low-BE end (``y[-avg:]``, "E2")."""
    np = _np()
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3:
        return y.copy()
    k = max(1, min(int(avg), n // 2 or 1))
    lo = float(y[:k].mean()) * (1.0 - st_offset / 100.0)
    hi = float(y[-k:].mean()) * (1.0 - end_offset / 100.0)
    b = np.full(n, lo)
    for _ in range(iters):
        cum = np.cumsum(np.maximum(y - b, 0.0))
        tot = cum[-1] if cum[-1] > 0 else 1.0
        nb = lo + (hi - lo) * cum / tot
        done = float(np.abs(nb - b).max()) < 1e-9 * max(1.0, abs(hi - lo))
        b = nb
        if done:
            break
    return b


def linear_bg(y, avg=1, st_offset=0.0, end_offset=0.0):
    """Linear background between the two averaged end levels, with the same
    CasaXPS ``St Offset``/``End Offset`` percentage reduction as
    :func:`shirley` (see its docstring)."""
    np = _np()
    y = np.asarray(y, dtype=float)
    n = len(y)
    k = max(1, min(int(avg), n // 2 or 1))
    lo = float(y[:k].mean()) * (1.0 - st_offset / 100.0)
    hi = float(y[-k:].mean()) * (1.0 - end_offset / 100.0)
    return lo + (hi - lo) * np.linspace(0.0, 1.0, n)


def tougaard_u2(x, y, b, c, avg=1, end_offset=0.0):
    """The two-parameter universal Tougaard background of ``y`` on the
    ascending-KE grid ``x``: the level at the high-KE end (mean of ``avg``
    points) plus the inelastic tail of everything above each point,
    ``sum K(T) (y - end) dE`` over ``T = x' - x``, with the cross-section
    ``K(T) = B T / (C + T^2)^2``. CasaXPS stores C with a minus sign (its
    general form is ``B T / ((C - T^2)^2 + D T^2)`` and the two-parameter
    version has D = 0). ``b`` is B (eV^2), ``c`` the positive C (eV^2).
    Checked on real CasaXPS fits (see ``tests/test_casafit.py``). ``end_offset``
    is CasaXPS's "End Offset" percentage (see :func:`shirley`), applied to this
    function's only anchor, ``base``; there is no "St Offset" parameter here
    because this construction has no separate start-anchor for one to act on."""
    np = _np()
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3:
        return y.copy()
    k = max(1, min(int(avg), n // 2 or 1))
    base = float(y[-k:].mean()) * (1.0 - end_offset / 100.0)
    dx = float(np.abs(np.diff(x)).mean())
    yb = y - base
    bg = np.empty(n)
    for i in range(n):
        t = x[i + 1:] - x[i]
        bg[i] = float(np.sum(b * t / ((c + t * t) ** 2) * yb[i + 1:]))
    return base + bg * dx


def tougaard_3param(x, y, b, c, d, avg=1, end_offset=0.0):
    """The three-parameter universal Tougaard background of ``y`` on the
    ascending-KE grid ``x``: same construction as :func:`tougaard_u2` (the
    level at the high-KE end plus the inelastic tail of everything above each
    point) but with the three-parameter cross section
    ``K(T) = B T / ((C - T^2)^2 + D T^2)`` (Tougaard, *Surf. Interface Anal.*
    25, 137 (1997); CasaXPS's own "Peak Fitting in XPS" names this the
    ``U 4 Tougaard`` cross section, with ``U Poly``/``U Si``/``U SiO2``/
    ``U Ge``/``U Al`` as its built-in material presets). Unlike
    :func:`tougaard_u2`, CasaXPS stores ``B``, ``C`` and ``D`` here with their
    natural sign (no minus on ``C``). Checked against a real CasaXPS
    ``U Poly Tougaard`` fit (``D:\\Temp\\for claude files\\PET``): the file's
    own ``params[2:5]`` for both its C 1s and O 1s regions are ``396, 551,
    436`` -- an exact match to Tougaard's published Polymers row -- and
    reconstructing background + the file's own components against its own raw
    data with these values gives residuals in the same few-percent range as
    this module's pre-existing ``LA``-shape reconstruction noise, not a
    background-shape mismatch. The two further numbers CasaXPS's region line
    carries ahead of ``B`` (``params[0]``, ``params[1]``) are now confirmed,
    not guessed: CasaXPS's own "Peak Fitting in XPS" whitepaper (Casa
    Software Ltd, 2006, p.4) names them the generic per-region "St Offset"/
    "End Offset" percentages every background type carries (see
    :func:`shirley`), not a per-cross-section ``T0`` energy-loss cutoff --
    that earlier guess was already ruled out by ``params[1]`` as a literal
    ``T0`` collapsing this same file's C 1s background to flat (19.6 eV
    exceeds the whole 15 eV fit window). ``end_offset`` (``params[1]``) is
    applied to ``base``, this function's only anchor; there is no "St Offset"
    parameter here for the same reason :func:`tougaard_u2` has none."""
    np = _np()
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3:
        return y.copy()
    k = max(1, min(int(avg), n // 2 or 1))
    base = float(y[-k:].mean()) * (1.0 - end_offset / 100.0)
    dx = float(np.abs(np.diff(x)).mean())
    yb = y - base
    bg = np.empty(n)
    for i in range(n):
        t = x[i + 1:] - x[i]
        bg[i] = float(np.sum(b * t / ((c - t * t) ** 2 + d * t * t) * yb[i + 1:]))
    return base + bg * dx


def tougaard_w(x, y, b, c):
    """The W Tougaard background of ``y`` on the ascending-KE grid ``x``:
    reconstructed from KherveFitting's own
    ``calculate_w_tougaard_background`` (github.com/KherveFitting/
    KherveFitting, ``libraries/Peak_Functions.py``, read directly 2026-09-27)
    -- a genuinely different cross section from both :func:`tougaard_u2` and
    :func:`tougaard_3param`, not a renamed duplicate: ``K(T) = B' T /
    (C + T^2)`` (no ``D`` term, and ``C`` is not squared), where
    ``B' = B * y[0] / y[-1]`` -- an endpoint-intensity-ratio adjustment to
    ``B`` KherveFitting's own comment calls "Adjust B based on endpoint
    intensities" that no other Tougaard variant in this module has. Ported
    faithfully, including that KherveFitting's own version does **not**
    subtract and restore a flat baseline the way :func:`tougaard_u2`/
    :func:`tougaard_3param` do (so no ``avg`` parameter here: there is no
    baseline level to average).

    **Two things are left unconfirmed, same tier as the Gelius/TLA shapes**:
    no real CasaXPS-fitted file using literally ``"W Tougaard"`` as its
    region type is available anywhere in the corpora this module was
    checked against, so this is a reconstruction from a credible secondary
    source only, self-consistency tested, not validated against a real
    residual; and KherveFitting's own ``y[0]``/``y[-1]`` endpoint indices are
    taken as-is with no confirmation of which physical end (high or low
    binding energy) they are meant to be on *this* module's own
    ascending-kinetic-energy convention.

    Neither CasaXPS "St Offset"/"End Offset" percentage (see :func:`shirley`)
    is applied here: unlike :func:`tougaard_u2`/:func:`tougaard_3param`, this
    construction has no flat baseline term at all for either to act on."""
    np = _np()
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3:
        return y.copy()
    b_adj = b * (y[0] / y[-1]) if y[-1] else b
    dx = float(np.abs(np.diff(x)).mean())
    bg = np.empty(n)
    for i in range(n):
        t = x[i + 1:] - x[i]
        bg[i] = float(np.sum(b_adj * t / (c + t * t) * y[i + 1:]))
    return bg * dx


_TOUGAARD_U2 = re.compile(r"^u\s*2\s*tougaard")
_TOUGAARD_3P = re.compile(r"^(?:u\s*(?:poly|sio2|si|ge|al|4)\s*)?tougaard\b")
_TOUGAARD_W = re.compile(r"^w\s*tougaard\b")


def background(kind, y, avg=1, x=None, params=()):
    """Background under ``y`` for a CasaXPS type name ('Shirley', 'Linear',
    'None', 'Tougaard', 'U 2 Tougaard', 'U Poly Tougaard', 'U Si Tougaard',
    'U SiO2 Tougaard', 'U Ge Tougaard', 'U Al Tougaard', 'U 4 Tougaard',
    'W Tougaard' ...); None for a type this module cannot reproduce ('E
    Tougaard', a Spline background -- checked against KherveFitting's own
    readable background source too: it has no ``E`` Tougaard function at all
    either, and its own Spline background is anchor-point/interactive, not a
    closed form of the handful of numbers a CasaXPS region line's ``params``
    can carry, so there is nothing to reconstruct it from). ``x`` (ascending
    KE) is needed for every Tougaard variant; ``params`` (the region line's
    six numbers after the averaging width) is needed for every type except
    'None': ``params[0]``/``params[1]`` are CasaXPS's "St Offset"/
    "End Offset" percentages (see :func:`shirley`), read here and passed to
    every background function that has an anchor for them to act on, and
    ``params[2:5]`` are the Tougaard B/C/D as each variant's own docstring
    describes.

    Plain ``'Tougaard'`` (no ``U ...`` prefix) shares :func:`tougaard_3param`
    with the ``U ...`` family: KherveFitting's own ``calculate_
    tougaard_background`` uses the identical ``B T / ((C - T^2)^2 + D T^2)``
    cross section (confirmed by reading both side by side) with generic
    default constants instead of a named material preset -- a second,
    independent confirmation of that function (Tougaard's own paper was the
    first), though still not a validation against a real CasaXPS file using
    literally ``"Tougaard"`` as its type (none found), so it carries the
    same caveat :func:`tougaard_3param`'s own docstring already states.
    ``'W Tougaard'`` is a distinct cross section (:func:`tougaard_w`, its own
    caveats there); its ``C`` sign convention (which family it follows,
    :func:`tougaard_u2`'s negated storage or :func:`tougaard_3param`'s
    natural one) is unconfirmed, so it is read with the natural sign, the
    more common convention among CasaXPS's Tougaard types."""
    t = str(kind or "").strip().lower()
    st = params[0] if len(params) > 0 else 0.0
    en = params[1] if len(params) > 1 else 0.0
    if t.startswith("shirley"):
        return shirley(y, avg, st_offset=st, end_offset=en)
    if _TOUGAARD_U2.match(t):
        if x is None or len(params) < 4 or not params[2]:
            return None
        return tougaard_u2(x, y, params[2], abs(params[3]) or 1643.0, avg,
                           end_offset=en)
    if _TOUGAARD_W.match(t):
        if x is None or len(params) < 4 or not params[2] or not params[3]:
            return None
        return tougaard_w(x, y, params[2], params[3])
    if _TOUGAARD_3P.match(t):
        if x is None or len(params) < 5 or not params[2] or not params[3]:
            return None
        return tougaard_3param(x, y, params[2], params[3], params[4], avg,
                               end_offset=en)
    if t.startswith("linear"):
        return linear_bg(y, avg, st_offset=st, end_offset=en)
    if t in ("none", "", "offset"):
        np = _np()
        return np.zeros(len(y))
    return None
