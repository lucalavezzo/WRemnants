"""Plot the SCETlib NP form factors: CS γ_ν^NP(b_T) and TMD F_eff(b_T, y).

A pure plotting library: takes physical λ values (the two parameter dicts the
model uses) and draws the two NP functions. Where the λ come from (a new
continuous-λ fit, an old template-based fit, or hand-picked values) is the
caller's job. The companion reader :mod:`fitresult_lambdas` turns a fitresults
HDF5 into the λ sets / toy ensembles this module consumes; ``main()`` glues the
two together, but the plot functions stay reader-agnostic.

The curves call the same form factors the fit integrates
(:func:`btgrid_tf.F_eff_tf` / :func:`btgrid_tf.gamma_nu_NP_tf`), driven by the
``np_model`` / ``np_model_nu`` strings, so a plotted curve is the fit's model.

A "λ set" is the pair of dicts :class:`NPLambdas` carries:

    eff = {lambda_inf, lambda2, lambda4, lambda6, delta_lambda2, np_model}   (F_eff / TMD)
    gnu = {lambda_inf_nu, lambda2_nu, lambda4_nu, np_model_nu}                (γ_ν / CS)

These map 1:1 onto the form-factor keyword arguments. A band is drawn from a
caller-supplied list of λ-set "toys" (this module takes percentiles of the
resulting curves); it never samples and never sees a covariance.

CLI (plot a raw λ set, no fit involved)::

    python -m wremnants.postprocessing.scetlib_np.np_function_plots \\
        --num-lambdas lambda2=0.4,lambda4=0.4,lambda2_nu=0.15 \\
        --num-np-model tanh_6 --num-np-model-nu tanh_2 -o /tmp/np.png

The ``--num-lambdas`` names must be λ the chosen models actually use (e.g.
``lambda6`` needs ``--num-np-model tanh_6``); naming an inert λ is a hard error
rather than a silently-ignored value. Unset λ stay at their defaults (NP-unit
point).

CLI (from a fitresults: prefit dashed, postfit solid + 68% band)::

    python -m wremnants.postprocessing.scetlib_np.np_function_plots \\
        --num-fitresult <fitresults.hdf5> -o /tmp/np.png

``--num-lambdas`` combines with ``--num-fitresult``: the λ come from the fit and the
named ones are overridden, on BOTH points and pinned in the band (so "the postfit
tune but with λ4 = 0" is ``--num-fitresult ... --num-lambdas lambda4=0``).

``--num-np-model`` also combines with ``--num-fitresult``, so a fit's λ can be
replayed under a DIFFERENT F_eff form — e.g. the postfit tune under the large-b_T
turn-off variant::

    --num-fitresult ... --num-np-model tanh_6_sigmoid \\
        --num-lambdas bT_cutoff=2,bT_cutoff_width=0.2

λ the requested form needs but the fit never stored (here the two cutoff shape
constants) fall back to the registry defaults, and every such fill is reported on
stdout by :mod:`fitresult_lambdas`; pass them in ``--num-lambdas`` to choose them.

Tune flags come from the shared ingestor :mod:`np_tune`, so they are spelled the
same here as in :mod:`sigma_gen_at_lambda` / :mod:`sigma_reco_at_lambda`:
``--num-*`` / ``--den-*``, one spelling per flag (the historical unprefixed and
``--ratio-*`` aliases were removed on 2026-08-07). Either side's fitresults may be
given as the run DIRECTORY instead of the hdf5 inside it.

CLI (TWO tunes overlaid — e.g. the same fit under two PDF sets)::

    python -m wremnants.postprocessing.scetlib_np.np_function_plots \\
        --num-fitresult <runA>/fitresults.hdf5 \\
        --den-fitresult <runB>/fitresults.hdf5 -o /tmp/np.png

Any ``--den-*`` flag turns the comparison on, exactly as in
:mod:`sigma_gen_at_lambda`: the numerator is drawn from the viridis ramp, the
denominator from plasma (see the colour contract below), each solid with its own 68%
band, and the parameter inset becomes a two-column λ table. There is no ratio panel — these are form factors, and
their ratio is not a thing anyone reads.

Each side is resolved INDEPENDENTLY (its own fitresults, its own ``--*-lambdas``,
its own forms), so a λ not set on a side falls back to that side's own source and
never to the other side's value — the same rule :mod:`np_tune` uses. Forms, unlike
λ, are inherited den-from-num when the denominator has no fitresults of its own.

**Colour contract:**

* *Colormap* = WHICH TUNE. The numerator's curves are sampled from
  :data:`NUM_CMAP` (viridis), the denominator's from :data:`DEN_CMAP` (plasma); a
  single-tune plot is the numerator, so viridis. Grey means "not one of ours", i.e.
  an external reference (the lattice CS kernel, the MAP22 F_eff).
* *Position along the ramp* = WHICH y, on the TMD panel where each tune draws one
  curve per rapidity (:data:`CMAP_SHADE_RANGE`, first y at the dark end). The CS
  panel has no y, so its curve takes the ramp's y = first colour — that is what ties
  a tune's left and right curves together visually.
* *Linestyle* = PREFIT or POSTFIT, and nothing else: a side's prefit is the DASHED
  curve of its own ramp, its postfit the SOLID one, both postfits solid. (MAP22 is
  dotted, so it cannot be mistaken for a prefit.)

A series given a plain ``color`` and no ``cmap`` — a library caller, or the raw
``--num-lambdas`` point — keeps the older behaviour of shading ONE hue by y
(:func:`_y_shades`), which is the same contract with a one-hue ramp.

The prefit (λ_central) curve is drawn ONCE, from the numerator: it is a property
of the card, not of the fit. If the two sides' prefit points actually differ — two
different cards, or a ``--den-lambdas`` / ``--den-np-model`` override, which move
that side's prefit as well as its postfit — both are drawn and the difference is
reported. ``--no-prefit`` drops them.
"""

import argparse
import os
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np

from wremnants.postprocessing.scetlib_np import btgrid_tf, np_tune
from wremnants.postprocessing.scetlib_np.params import (
    EFF_PARAMS,
    GNU_PARAMS,
    NPTune,
    active_params,
    check_lambda_overrides,
    parse_lambda_overrides,
)

# One λ point = one :class:`~params.NPTune` (forms + their λ, validated together).
# ``NPLambdas`` was this module's own dataclass with the same job, which
# :mod:`fitresult_lambdas` imported FROM the plotting module — a backwards
# dependency now gone. Kept as an alias so existing callers keep working; its
# ``.eff`` / ``.gnu`` are the same ``btgrid_tf`` kwargs, and ``from_flat`` is
# ``NPTune.create`` with the λ supplied as a base (so λ the chosen forms ignore
# are dropped, not rejected — a fitresults legitimately carries λ from a wider
# vocabulary than the plotted form uses).
NPLambdas = NPTune


@dataclass
class Series:
    """A curve to draw: a λ set, styling, and optional toys for an error band."""

    label: str
    lam: NPLambdas
    color: Optional[str] = None
    linestyle: str = "-"
    lw: float = 2.0
    toys: Optional[List[NPLambdas]] = None  # band drawn from these if present
    band_pct: Tuple[float, float] = (16.0, 84.0)
    #: Colormap this series takes its per-y TMD colours from (see the colour
    #: contract). ``None`` falls back to shades of :attr:`color`. When set and
    #: ``color`` is None, the CS-panel curve — which has no y — uses the ramp's
    #: FIRST (y = first) colour, so the two panels stay visibly the same series.
    cmap: Optional[str] = None


# ---------------------------------------------------------------------------
# Lattice reference for the CS panel
# ---------------------------------------------------------------------------
# Cridge, Marinelli & Tackmann, arXiv:2506.13874 (JHEP 12 (2025) 043) Eq. (3.34):
# a fit of THIS parametrization to the lattice CS-kernel data of refs [89-91]
# (2307.12359, 2306.06488, 2302.06502). Their Eq. (3.32) is
#
#     2 γ̃_ζ^np(b_T) = γ̃_ν^np(b_T) = −λ∞ tanh(λ2 b_T²/λ∞ + λ4 b_T⁴/λ∞)
#
# i.e. the tanh_2 form with λ6_ν = 0 — the reference is drawn in the paper's OWN
# model, NOT in whatever ``np_model_nu`` the plotted series use. This panel shows
# γ̃_ν, so the curve is used as-is; a lattice paper plotting γ̃_ζ shows HALF of it.
#
# The covariance is included only to draw a ±1σ guide band. Its provenance is
# Ref. [66] of that paper (Dehnadi, Ploessl & Tackmann, in preparation, no arXiv
# number) and the lattice authors have said it is not ready to be trusted, so it
# is a visual scale reference and must NOT be turned into a fit prior. See
# ``knowledge/30_physics_global/np_parametrization_constraints.md`` §10-§13.
LATTICE_CS_MU = np.array([1.6853, 0.0870, 0.0074])  # λ∞, λ2 [GeV²], λ4 [GeV⁴]
LATTICE_CS_SIGMA = np.array([0.5069, 0.0332, 0.0066])
LATTICE_CS_CORR = np.array(
    [
        [1.0000, 0.5212, -0.7249],
        [0.5212, 1.0000, -0.9135],
        [-0.7249, -0.9135, 1.0000],
    ]
)
LATTICE_CS_COV = np.outer(LATTICE_CS_SIGMA, LATTICE_CS_SIGMA) * LATTICE_CS_CORR


def lattice_gamma_nu(bT, lam_inf, lambda2_nu, lambda4_nu):
    """γ̃_ν^NP(b_T) in the paper's tanh_2 form, Eq. (3.32). Plain numpy: the
    reference is a fixed curve, never differentiated, so it does not go through
    ``btgrid_tf`` (which would also force the series' own ``np_model_nu``)."""
    bT = np.asarray(bT, dtype=float)
    return -lam_inf * np.tanh((lambda2_nu * bT**2 + lambda4_nu * bT**4) / lam_inf)


def lattice_cs_band(bT, nsig=1.0, n_samples=20000, seed=0):
    """Envelope of Eq. (3.32) over the ``nsig`` covariance ellipsoid of Eq. (3.34).

    Sampled on the ellipsoid SURFACE (not the interior) so the result is the
    nsig envelope of the curve, which is what a band should show. λ∞ ≤ 0 draws
    are dropped: they are not a CS kernel."""
    chol = np.linalg.cholesky(LATTICE_CS_COV)
    rng = np.random.default_rng(seed)
    u = rng.normal(size=(3, n_samples))
    u /= np.linalg.norm(u, axis=0)
    pars = LATTICE_CS_MU[:, None] + nsig * (chol @ u)
    pars = pars[:, pars[0] > 0]
    curves = lattice_gamma_nu(
        np.asarray(bT)[None, :], pars[0][:, None], pars[1][:, None], pars[2][:, None]
    )
    return curves.min(axis=0), curves.max(axis=0)


def draw_lattice_cs_reference(ax, bT, nsig=1.0, band=True):
    """Overlay the lattice-fit CS kernel on a γ̃_ν panel (drawn behind the series)."""
    if band:
        # No legend entry: the band comes from a provisional covariance (see the
        # LATTICE_CS_MU comment) and should read as a soft visual guide around the
        # central curve, not as a quoted uncertainty.
        lo, hi = lattice_cs_band(bT, nsig=nsig)
        ax.fill_between(bT, lo, hi, color="0.55", alpha=0.28, lw=0, zorder=0)
    ax.plot(
        bT,
        lattice_gamma_nu(bT, *LATTICE_CS_MU),
        color="0.15",
        lw=1.8,
        ls=(0, (5, 2)),
        zorder=1,
        label=r"lattice fit (2506.13874)",
    )


# ---------------------------------------------------------------------------
# MAP22 reference for the TMD panel
# ---------------------------------------------------------------------------
# MAP Collaboration MAP22 fit at N3LL, central replica: NangaParbat-Legacy
# ``FitResults/MAP22_N3LL/fitconfig.yaml``, parameterisation ``MAP22g52``
# (inc/NangaParbat/MAP22g52.h, the ``ifunc < 2`` TMD-PDF branch). Data: DY only
# (E605/E288/E772, CDF/D0, STAR, LHCb/CMS/ATLAS), MSHT2020 NNLO PDFs.
#
#   f_int(x,b) = [ g1 e^{-g1 b²/4} + λ² g1B² (1 - g1B b²/4) e^{-g1B b²/4}
#                  + λ2² g1C e^{-g1C b²/4} ] / [ g1 + λ² g1B² + λ2² g1C ]
#   g_{1,1B,1C}(x) = N (x/x̂)^σ ((1-x)/(1-x̂))^{α²},   x̂ = 0.1
#
# TWO THINGS TO KNOW BEFORE COMPARING:
#  * Their EVOLUTION factor exp(-g2² b² ln(ζ/Q0²)/4) is STRIPPED here. That is
#    their CS-kernel piece, which in our model is γ_ν^NP (the other panel) —
#    keeping it would double-count the CS side. What is drawn is the intrinsic
#    (boundary-condition) part only, which is the F_eff analogue.
#  * f_int is PER HADRON, while F_eff is the two-beam product. So the reference is
#    f_int(x1,b)·f_int(x2,b) with x_{1,2} = (Q/√s)e^{±Y} — which is also where the
#    Y dependence comes from, i.e. this reference PREDICTS ΔΛ₂ rather than fitting it.
#
# Scheme caveat: f_NP is defined as the remainder after a PARTICULAR perturbative
# piece. MAP22 uses its own b* (Eq. 4 of arXiv:2502.04166, with b_min ≠ 0, which we
# do not have) at N3LL, so their NP function is not identical in scheme to ours.
# Treat as an external shape reference, not as truth.
MAP22_N3LL_PARS = (
    0.2482420850967863,  # g2      (evolution; unused here, see above)
    0.31560484357030422,  # N1
    1.1349706417780698,  # alpha1
    0.51541879126481949,  # sigma1
    1.8893385924128125,  # lambda
    0.0043928548895256119,  # N3     ) FF branch,
    10.843398922979528,  # beta1    ) unused for
    0.0070355513246244906,  # delta1 ) the TMD PDF
    1.5397221152343268,  # gamma1   )
    0.061287263131035774,  # lambdaF)
    0.21739975134899683,  # N3B     )
    0.12887993135476603,  # N1B
    0.016302544233918744,  # N1C
    0.019602747712729469,  # lambda2
    4.2653349811775403,  # alpha2
    4.3239264662485306,  # alpha3
    0.41795806310999201,  # sigma2
    12.673323632665454,  # sigma3
)
MAP22_XHAT = 0.1
# Default kinematics for the Y -> x map: on-shell Z at 13 TeV.
MAP22_Q_DEFAULT = 91.1876
MAP22_SQRTS_DEFAULT = 13000.0


def map22_f_int(x, bT, pars=MAP22_N3LL_PARS):
    """MAP22 intrinsic (evolution-stripped) NP function for ONE hadron."""
    x = np.asarray(x, dtype=float)
    bT = np.asarray(bT, dtype=float)
    N1, al1, sg1, lam = pars[1], pars[2], pars[3], pars[4]
    N1B, N1C, lam2 = pars[11], pars[12], pars[13]
    al2, al3, sg2, sg3 = pars[14], pars[15], pars[16], pars[17]

    def g(N, sg, al):
        return N * (x / MAP22_XHAT) ** sg * ((1 - x) / (1 - MAP22_XHAT)) ** (al**2)

    g1, g1B, g1C = g(N1, sg1, al1), g(N1B, sg2, al2), g(N1C, sg3, al3)
    q = (bT / 2.0) ** 2
    num = (
        g1 * np.exp(-g1 * q)
        + lam**2 * g1B**2 * (1 - g1B * q) * np.exp(-g1B * q)
        + lam2**2 * g1C * np.exp(-g1C * q)
    )
    return num / (g1 + lam**2 * g1B**2 + lam2**2 * g1C)


def map22_F_eff(
    bT, Y, Q=MAP22_Q_DEFAULT, sqrt_s=MAP22_SQRTS_DEFAULT, pars=MAP22_N3LL_PARS
):
    """MAP22 analogue of F_eff(b_T, Y): the two-beam product at x = (Q/√s)e^{±Y}."""
    x1 = (Q / sqrt_s) * np.exp(+abs(Y))
    x2 = (Q / sqrt_s) * np.exp(-abs(Y))
    return map22_f_int(x1, bT, pars) * map22_f_int(x2, bT, pars)


def draw_map22_tmd_reference(ax, bT, y_values, colors, Q=None, sqrt_s=None):
    """Overlay the MAP22 F_eff analogue on a TMD panel, one curve per y."""
    Q = MAP22_Q_DEFAULT if Q is None else Q
    sqrt_s = MAP22_SQRTS_DEFAULT if sqrt_s is None else sqrt_s
    for yi, y in enumerate(y_values):
        # Label EVERY curve with its y: the reference is y-dependent (that is the
        # whole point -- it PREDICTS the rapidity dependence), so a single legend
        # entry for the whole family makes the panel unreadable.
        #
        # Colour: the caller passes the NEUTRAL GREY shade ramp (see
        # ``_y_shades`` / the colour contract in the module docstring) -- an external
        # reference must not wear a tune's colour, and grey is what the lattice
        # reference already uses on the CS panel. Like-y pairing with the tunes is by
        # shade LEVEL (both ramps go dark -> light with y), not by identical colour.
        #
        # DOTTED, not dashed, because ``plot_series_from_fitresult`` draws the PREFIT
        # series with "--"; a dashed reference would be indistinguishable from it.
        ax.plot(
            bT,
            map22_F_eff(bT, y, Q, sqrt_s),
            color=colors[yi],
            lw=1.8,
            ls=(0, (1, 1.3)),
            zorder=1,
            label=f"MAP22, y={y:g}",
        )


def gamma_nu_curve(bT, gnu):
    """γ_ν^NP(b_T) for one gnu dict (CS sector).

    ``gnu`` carries the numeric λ plus the ``np_model_nu`` key; the form reads the
    λ it needs and ignores the model key (passed explicitly as the selector)."""
    return np.asarray(
        btgrid_tf.gamma_nu_NP_tf(bT, gnu, np_model_nu=gnu["np_model_nu"]), dtype=float
    )


def f_eff_curve(bT, y, eff):
    """F_eff(y, b_T) for one eff dict at rapidity ``y`` (TMD sector)."""
    return np.asarray(
        btgrid_tf.F_eff_tf(y, bT, eff, np_model=eff["np_model"]), dtype=float
    )


def _band(curves, pct):
    """(lo, hi) percentile envelope across a stack of curves, or None if empty."""
    if not len(curves):
        return None
    stack = np.asarray(curves)
    return np.percentile(stack, pct[0], axis=0), np.percentile(stack, pct[1], axis=0)


_CORNERS = {
    "upper right": (0.97, 0.97, "top", "right"),
    "upper left": (0.03, 0.97, "top", "left"),
    "lower left": (0.03, 0.03, "bottom", "left"),
    "lower right": (0.97, 0.03, "bottom", "right"),
}


# LaTeX labels per λ name, for the parameter inset.
_EFF_LABELS = {
    "lambda2": r"\lambda_2",
    "lambda4": r"\lambda_4",
    "delta_lambda2": r"\delta\lambda_2",
    "lambda6": r"\lambda_6",
    "lambda_inf": r"\lambda_\infty",
    "bT_cutoff": r"b_T^{\rm cut}",
    "bT_cutoff_width": r"w_{\rm cut}",
    "abs_margin": r"m_{\rm fold}",
}
_GNU_LABELS = {
    "lambda2_nu": r"\lambda_2^\nu",
    "lambda4_nu": r"\lambda_4^\nu",
    "lambda6_nu": r"\lambda_6^\nu",
    "lambda_inf_nu": r"\lambda_\infty^\nu",
    "abs_margin_nu": r"m_{\rm fold}^\nu",
}


def _sector_view(lam, sector):
    """``(model, active λ names, λ dict, print order)`` for one panel's sector."""
    if sector == "gnu":
        model = lam.gnu.get("np_model_nu", "?")
        return model, active_params(np_model_nu=model), lam.gnu, GNU_PARAMS
    model = lam.eff.get("np_model", "?")
    return model, active_params(np_model=model), lam.eff, EFF_PARAMS


def _param_inset(ax, lam, sector, corner="upper right"):
    """Small text box listing only the λ that drive the panel for its model."""
    model, active, src, order = _sector_view(lam, sector)
    labels = _GNU_LABELS if sector == "gnu" else _EFF_LABELS
    lines = [rf"${labels[k]} = {src.get(k, 0):+.4f}$" for k in order if k in active]
    lines.append(rf"model: {model}")
    box = dict(boxstyle="round,pad=0.35", fc="white", ec="0.6", alpha=0.85)
    x, y, va, ha = _CORNERS[corner]
    ax.text(
        x,
        y,
        "\n".join(lines),
        transform=ax.transAxes,
        fontsize=8,
        va=va,
        ha=ha,
        bbox=box,
    )


def _inset_table_text(lams, titles, sector, title_width=26):
    """One λ per row, one λ SET per column: the two-tune version of the inset.

    Monospace and plain λ names (not the single-set box's mathtext), because what
    matters with two tunes side by side is that the columns line up. A λ a column's
    form does not use prints ``--`` rather than a value it would not read — the two
    sides may legitimately be at different forms. Column headers are the series'
    own legend labels (truncated), so the box needs no colour key to say which
    column is which curve."""
    views = [_sector_view(lam, sector) for lam in lams]
    rows = [p for p in views[0][3] if any(p in v[1] for v in views)]
    heads = [
        (t if len(t) <= title_width else "…" + t[-(title_width - 1) :]) for t in titles
    ]
    w = max([len(p) for p in rows] + [len("model")])
    cw = max([len(h) for h in heads] + [8])
    out = [" " * w + "  " + "  ".join(f"{h:>{cw}}" for h in heads)]
    for p in rows:
        cells = [f"{v[2][p]:+.4f}" if p in v[1] else "--" for v in views]
        out.append(f"{p:<{w}}  " + "  ".join(f"{c:>{cw}}" for c in cells))
    out.append(f"{'model':<{w}}  " + "  ".join(f"{v[0]:>{cw}}" for v in views))
    return "\n".join(out)


def _param_inset_compare(ax, lams, titles, sector, corner="upper right"):
    """Draw the multi-column λ table (:func:`_inset_table_text`) in a corner."""
    x, y, va, ha = _CORNERS[corner]
    ax.text(
        x,
        y,
        _inset_table_text(lams, titles, sector),
        transform=ax.transAxes,
        fontsize=6.5,
        family="monospace",
        va=va,
        ha=ha,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.6", alpha=0.88),
    )


def _fixed_color(c):
    """Resolve a ``"CN"`` property-cycle reference to a concrete colour NOW.

    :func:`plot_output.save_plot` imports ``wums.plot_tools``, which runs
    ``hep.style.use`` at IMPORT and so replaces the property cycle after every
    artist on the figure already exists. ``Line2D`` resolves ``"CN"`` at draw time
    (against the new cycle) while ``fill_between`` resolves it eagerly (the old
    one), so a banded series came out with its line and its band in DIFFERENT
    colours. Pinning the colour at creation makes the late style change a no-op."""
    import matplotlib.colors as mcolors

    return mcolors.to_hex(mcolors.to_rgba(c)) if isinstance(c, str) else c


#: One COLORMAP per side, sampled across the requested rapidities. Both are
#: perceptually uniform and colourblind-safe, and they are far apart at every shade
#: level (viridis runs blue -> teal -> yellow-green, plasma purple -> salmon ->
#: orange), so "which y" and "which tune" stay separately readable.
NUM_CMAP = "viridis"
DEN_CMAP = "plasma"

#: Where in each colormap the y ramp starts and ends. Avoids both extreme ends: the
#: near-black start and the near-white finish are unreadable as line colours.
CMAP_SHADE_RANGE = (0.25, 0.85)

#: Neutral hue for EXTERNAL references (lattice on the CS panel, MAP22 on the TMD
#: one). Part of the colour contract: a colormap identifies WHICH TUNE, grey means
#: "not one of ours".
REFERENCE_GREY = "0.35"


def _cmap_shades(name, n):
    """``n`` colours across one colormap, over :data:`CMAP_SHADE_RANGE`.

    The per-y ramp of one tune. ``n == 1`` takes the ramp's start, which is also what
    the y-less CS panel uses, so a tune is recognisably the same series in both."""
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap(name)
    lo, hi = CMAP_SHADE_RANGE
    if n <= 1:
        return [_to_hex(cmap(lo)[:3])]
    return [_to_hex(cmap(t)[:3]) for t in np.linspace(lo, hi, n)]


def _y_shades(base, n):
    """``n`` shades of ONE hue, dark (first y) -> light (last y).

    The fallback per-y ramp for a series given a plain ``color`` and no ``cmap``
    (library callers, and the raw ``--num-lambdas`` point). Same job as
    :func:`_cmap_shades`: the shade LEVEL says which rapidity."""
    rgb = np.asarray(_to_rgb(base), dtype=float)
    if n <= 1:
        return [_to_hex(rgb)]
    dark = rgb * 0.55
    # Only 45% of the way to white: enough separation between shades to read, while
    # the palest curve is still unmistakably the tune's hue on a white background.
    light = 1.0 - (1.0 - rgb) * 0.55
    return [_to_hex(dark + (light - dark) * t) for t in np.linspace(0.0, 1.0, n)]


def _series_colors(s, n_y, fallback):
    """``(CS-panel colour, [per-y TMD colours])`` for one series.

    A ``cmap`` series takes its TMD ramp from that colormap and its single CS colour
    from the ramp's y = first end (an explicit ``color`` still wins, for a caller that
    wants the CS curve somewhere else on the ramp). Otherwise the ramp is shades of
    the series' own hue."""
    if s.cmap:
        shades = _cmap_shades(s.cmap, n_y)
        return (_fixed_color(s.color) if s.color else shades[0]), shades
    base = _fixed_color(s.color or fallback)
    return base, _y_shades(base, n_y)


def _to_rgb(c):
    import matplotlib.colors as mcolors

    return mcolors.to_rgb(c)


def _to_hex(rgb):
    import matplotlib.colors as mcolors

    return mcolors.to_hex(np.clip(rgb, 0.0, 1.0))


def _legend(ax, loc):
    """Panel legend, wrapped to two columns once a two-tune overlay fills it.

    The TMD panel draws one curve per (series, y) plus the per-y reference, so a
    comparison of two fits at three rapidities is a dozen entries — a single column
    of those covers the panel.

    ``handlelength`` is set explicitly, and that is not cosmetic: linestyle carries
    the prefit-vs-postfit distinction, and matplotlib's default key is 2 font-sizes
    long (~13 pt at the small font used here) while dash patterns scale with
    linewidth — short enough that a dashed key rendered as one unbroken stub,
    indistinguishable from the solid one. Long enough for a full dash cycle, or the
    legend lies."""
    n = len(ax.get_legend_handles_labels()[1])
    ax.legend(
        loc=loc,
        fontsize=8 if n <= 8 else 6.5,
        ncol=1 if n <= 8 else 2,
        handlelength=3.6,
    )


def plot_np_functions(
    series: Sequence[Series],
    *,
    y_values: Sequence[float] = (0.0, 2.5, 5.0),
    bT_max: float = 4.0,
    n_points: int = 401,
    outpath: str,
    inset_from: Optional[Series] = None,
    insets: Optional[Sequence[Tuple[str, Series]]] = None,
    f_ymax: Optional[float] = None,
    lattice_reference: bool = True,
    lattice_nsig: float = 1.0,
    map22_reference: bool = True,
    map22_Q: Optional[float] = None,
    map22_sqrts: Optional[float] = None,
    args=None,
):
    """Draw the two NP form factors for one or more λ sets.

    Parameters
    ----------
    series
        Curves to overlay. Each series with ``toys`` draws a percentile band on
        its panels; the caller supplies the toys.
    y_values
        Rapidity values for the TMD panel (F_eff depends on y; γ_ν does not).
    bT_max, n_points
        b_T grid for the curves [GeV^-1].
    inset_from
        Series whose λ fill the per-panel parameter box (default: last series,
        typically the postfit point).
    insets
        ``[(column title, series), …]`` for a MULTI-column λ table instead of the
        single-set box — how a two-tune comparison shows both tunes. Wins over
        ``inset_from`` when given.
    lattice_reference, lattice_nsig
        Overlay the lattice-fit CS kernel (Eq. 3.34 of arXiv:2506.13874) on the
        γ̃_ν panel, on by default — it is the reference for "is this kernel
        physical / the right size". Drawn in the paper's own tanh_2 form. The
        band is a visual guide from a provisional covariance, NOT a prior; see
        the ``LATTICE_CS_MU`` comment above.
    map22_reference, map22_Q, map22_sqrts
        Overlay the MAP22 (N3LL) analogue of F_eff on the TMD panel, on by default
        — the external reference for "is this TMD boundary condition sane". Needs
        kinematics to map y -> x; defaults are on-shell Z at 13 TeV. Read the
        ``MAP22_N3LL_PARS`` comment for the two comparison caveats (evolution factor
        stripped; per-hadron vs two-beam) before drawing conclusions.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    bT = np.linspace(0.0, bT_max, n_points)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.5))

    auto_colors = [c for c in plt.rcParams["axes.prop_cycle"].by_key()["color"]]

    # NP factors evaluated at the bare b_T grid (this grid's b* prescription is
    # the identity, b_bar == b_T; see param_model / base.conf b0_over_bmax=0); a
    # b*-frozen grid would need b_T -> b_bar mapped first. F_eff runs away at large
    # b_T for λ4 < 0 toys, so the TMD panel scales to the line curves (below), not
    # a runaway band tail.
    line_fmax = 0.0

    # Reference first so it sits behind every fitted curve.
    if lattice_reference:
        draw_lattice_cs_reference(axL, bT, nsig=lattice_nsig)

    n_y = len(y_values)
    for si, s in enumerate(series):
        # One ramp per series: `color` draws the y-less CS curve, `shades` the TMD
        # ones. Same series, same colours, both panels.
        color, shades = _series_colors(s, n_y, auto_colors[si % len(auto_colors)])

        # ---- CS panel: γ_ν^NP(b_T) (no y dependence) ----
        axL.plot(
            bT,
            gamma_nu_curve(bT, s.lam.gnu),
            label=s.label,
            color=color,
            lw=s.lw,
            ls=s.linestyle,
        )
        if s.toys:
            band = _band([gamma_nu_curve(bT, t.gnu) for t in s.toys], s.band_pct)
            if band is not None:
                axL.fill_between(
                    bT,
                    band[0],
                    band[1],
                    color=color,
                    alpha=0.22,
                    label=f"{s.label} {int(s.band_pct[1]-s.band_pct[0])}% band",
                )

        # ---- TMD panel: F_eff(b_T, y) per requested y ----
        for yi, y in enumerate(y_values):
            yc = shades[yi]
            line = f_eff_curve(bT, y, s.lam.eff)
            line_fmax = max(line_fmax, float(np.nanmax(line)))
            axR.plot(
                bT, line, color=yc, lw=s.lw, ls=s.linestyle, label=f"{s.label}, y={y:g}"
            )
            if s.toys:
                band = _band([f_eff_curve(bT, y, t.eff) for t in s.toys], s.band_pct)
                if band is not None:
                    axR.fill_between(bT, band[0], band[1], color=yc, alpha=0.18)

    if map22_reference:
        # Grey ramp, shaded by y like the tunes: an external reference never wears a
        # tune's hue (same rule as the grey lattice curve on the CS panel).
        y_colors = _y_shades(REFERENCE_GREY, len(y_values))
        draw_map22_tmd_reference(axR, bT, y_values, y_colors, map22_Q, map22_sqrts)

    axL.axhline(0, color="k", lw=0.5)
    axL.set_xlabel(r"$b_T$ [GeV$^{-1}$]")
    axL.set_ylabel(r"$\tilde\gamma_\nu^{\rm NP}(b_T)$")
    axL.set_title(
        r"CS rapidity anomalous dimension $\tilde\gamma_\nu^{\rm NP}(b_T)$", fontsize=11
    )
    _legend(axL, loc="lower left")
    axL.grid(alpha=0.3)
    if args.ylim_cs is not None:
        axL.set_ylim(args.ylim_cs[0], args.ylim_cs[1])

    axR.set_xlabel(r"$b_T$ [GeV$^{-1}$]")
    axR.set_ylabel(r"$F_{\rm eff}(b_T, y)$")
    axR.set_title(r"TMD-effective NP factor $F_{\rm eff}(b_T, y)$", fontsize=11)
    _legend(axR, loc="upper right")
    axR.grid(alpha=0.3)
    if args.ylim_tmd is not None:
        axR.set_ylim(args.ylim_tmd[0], args.ylim_tmd[1])
    else:
        # Scale to the line curves so a runaway band tail (bare F_eff, λ4 < 0) can't
        # dominate the autoscale.
        top = f_ymax if f_ymax is not None else max(1.1, 1.2 * line_fmax)
        axR.set_ylim(0.0, top)

    # Param boxes diagonally opposite each panel's legend to avoid collisions:
    # CS legend lower-left -> box upper-right; TMD upper-right -> box lower-left.
    if insets:
        titles = [t for t, _ in insets]
        lams = [s.lam for _, s in insets]
        _param_inset_compare(axL, lams, titles, "gnu", corner="upper right")
        _param_inset_compare(axR, lams, titles, "eff", corner="lower left")
    else:
        inset = (
            inset_from if inset_from is not None else (series[-1] if series else None)
        )
        if inset is not None:
            _param_inset(axL, inset.lam, "gnu", corner="upper right")
            _param_inset(axR, inset.lam, "eff", corner="lower left")

    # Allow --outpath to be a directory (trailing slash or no extension): append
    # a default filename rather than erroring on a bare ".png".
    from wremnants.postprocessing.scetlib_np import plot_output

    fig.tight_layout()
    outdir, basename = plot_output.split_outpath(
        outpath, default_name="np_functions.png"
    )
    plot_output.save_plot(outdir, basename, fig=fig, args=args, dpi=140)
    plt.close(fig)
    print(f"Wrote {outdir}/{basename}.png(.pdf) + {basename}.log")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

# Baseline λ for the raw mode: SCETlib "knobs off" (NP-unit point). --num-lambdas
# overrides any subset; unset λ stay here. Carries all ALL_PARAMS so the
# form-factor kwargs are complete regardless of the chosen model.
_DEFAULT_LAMBDAS = dict(
    lambda2=0.0,
    lambda4=0.0,
    lambda6=0.0,
    delta_lambda2=0.0,
    lambda_inf=1.0,
    lambda2_nu=0.0,
    lambda4_nu=0.0,
    lambda6_nu=0.0,
    lambda_inf_nu=1.0,
    # tanh_6_sigmoid turn-off shape constants (inert for every other model)
    bT_cutoff=2.0,
    bT_cutoff_width=0.2,
    # tanh_6_abs damping-fold margins (inert for every other model)
    abs_margin=0.0,
    abs_margin_nu=0.0,
)


def make_parser():
    p = argparse.ArgumentParser(
        description="Plot SCETlib NP form factors γ_ν^NP(b_T) and F_eff(b_T,y).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # The shared tune surface (np_tune): --num-fitresult / -lambdas / -np-model /
    # -label. One spelling per flag; the unprefixed and --ratio-* aliases are gone.
    src = np_tune.add_tune_args(p, "num")
    src.add_argument(
        "--n-toys",
        type=int,
        default=500,
        help="band toys, per side (needs a fitresults carrying a covariance).",
    )
    src.add_argument("--seed", type=int, default=0, help="band RNG seed.")
    # Second tune: any --den-* flag overlays it (no ratio panel — see the module
    # docstring). Same spellings as sigma_gen_at_lambda / sigma_reco_at_lambda.
    np_tune.add_tune_args(p, "den")

    p.add_argument(
        "--y",
        type=float,
        nargs="+",
        default=[0.0, 2.5, 5.0],
        help="rapidity values for the TMD panel.",
    )
    p.add_argument("--bT-max", type=float, default=4.0)
    p.add_argument(
        "--f-ymax",
        type=float,
        default=None,
        help="fixed upper y-limit for the TMD panel (default: auto "
        "from the line curves; clips runaway negative-λ band tails).",
    )
    p.add_argument(
        "--ylim-cs",
        type=float,
        nargs=2,
        default=None,
        help="fixed y-limits for the CS panel (default: auto).",
    )
    p.add_argument(
        "--no-map22-reference",
        dest="map22_reference",
        action="store_false",
        help="drop the MAP22 (N3LL) F_eff reference from the TMD panel; it is drawn "
        "by default as the external boundary-condition reference.",
    )
    p.add_argument(
        "--map22-Q",
        type=float,
        default=None,
        help="Q [GeV] for the MAP22 reference's y->x map (default: on-shell Z).",
    )
    p.add_argument(
        "--map22-sqrts",
        type=float,
        default=None,
        help="sqrt(s) [GeV] for the MAP22 reference's y->x map (default: 13000).",
    )
    p.add_argument(
        "--no-lattice-reference",
        dest="lattice_reference",
        action="store_false",
        help="drop the lattice-fit CS kernel (2506.13874 Eq. 3.34) from the CS "
        "panel; it is drawn by default as the physicality/size reference.",
    )
    p.add_argument(
        "--lattice-nsig",
        type=float,
        default=1.0,
        help="n-sigma envelope for the lattice reference band (default 1). The "
        "covariance is provisional -- a visual guide, never a fit prior.",
    )
    p.add_argument(
        "--ylim-tmd",
        type=float,
        nargs=2,
        default=None,
        help="fixed y-limits for the TMD panel (default: auto).",
    )
    p.add_argument(
        "--no-prefit",
        dest="prefit",
        action="store_false",
        help="drop the prefit (λ_central) curve(s); it is drawn by default as the "
        "starting point the fit moved away from.",
    )
    p.add_argument("--outpath", "-o", required=True)
    return p


def _run_label(fitresult_path):
    """A fit's run name: the directory the fitresults sits in. Same convention as
    :mod:`sigma_gen_at_lambda`, so the two tools label the same run identically."""
    return os.path.basename(os.path.dirname(os.path.abspath(fitresult_path)))


def _inset_title(side, label, width=18):
    """``"num: <label>"`` for the λ-table header, keeping the side marker.

    Run-directory labels are long and differ at the END, so a too-long one is cut
    from the FRONT — after the side prefix is attached, never before it, or the
    "which column is which curve" marker is the first thing lost."""
    return f"{side}: " + (label if len(label) <= width else "…" + label[-(width - 1) :])


def _fmt_diff(v):
    """One side of a :func:`np_tune.tune_diff` entry — a form name, a λ, or absent
    (a λ only one side's form uses)."""
    if v is None:
        return "-"
    return v if isinstance(v, str) else f"{v:+.5g}"


def _series_for_side(args, side, inherit_forms=None, default_label="input"):
    """The curves for ONE tune side: ``[prefit, postfit]`` from its fitresults, or
    a single raw λ point when it has none.

    λ are resolved per side and never inherited from the other side (the
    :mod:`np_tune` rule); FORMS do inherit den-from-num when this side has no
    fitresults of its own, so overriding one form on the denominator cannot
    silently move the other. Returns ``(series, (np_model, np_model_nu))``."""
    fr = getattr(args, f"{side}_fitresult")
    overrides = parse_lambda_overrides(
        getattr(args, f"{side}_lambdas"), what=f"--{side}-lambdas"
    )
    np_model = getattr(args, f"{side}_np_model")
    np_model_nu = getattr(args, f"{side}_np_model_nu")
    if inherit_forms is not None and not fr:
        np_model = np_model or inherit_forms[0]
        np_model_nu = np_model_nu or inherit_forms[1]

    if fr:
        # Reading lives in the reader module; the plotter stays pure. The λ come
        # from the fitresults, with --*-lambdas overriding whichever it names (the
        # inert-λ check happens there, where the fit's forms are known).
        from wremnants.postprocessing.scetlib_np import fitresult_lambdas as frl

        series = frl.plot_series_from_fitresult(
            fr,
            result=getattr(args, f"{side}_result"),
            n_toys=args.n_toys,
            seed=args.seed,
            np_model=np_model,
            np_model_nu=np_model_nu,
            overrides=overrides,
        )
    else:
        np_model = np_model or "tanh_2"
        np_model_nu = np_model_nu or "tanh_2"
        check_lambda_overrides(
            overrides, np_model, np_model_nu, what=f"--{side}-lambdas"
        )
        lam = NPLambdas.from_flat(
            {**_DEFAULT_LAMBDAS, **overrides}, np_model, np_model_nu
        )
        label = getattr(args, f"{side}_label") or default_label
        series = [Series(label=label, lam=lam, color="C3")]

    # This side's colour ramp, on every one of its curves (prefit AND postfit) and in
    # both panels; `color=None` lets the CS curve take the ramp's y = first end.
    for s in series:
        s.cmap = NUM_CMAP if side == "num" else DEN_CMAP
        s.color = None

    point = series[-1]  # the curve that IS this tune (postfit, or the raw point)
    return series, (point.lam.np_model, point.lam.np_model_nu)


def _report_band(series, fitresult_path, side):
    """Say whether this side got an error band, and where to find one if not.

    The λ live in the FIT pass, the covariance in the follow-up ``cov/`` pass, so
    pointing at a run directory legitimately yields no band. Silence there would
    read as "this tune has no uncertainty"."""
    point = series[-1]
    if point.toys:
        print(f"[np_function_plots] {side}: {len(point.toys)}-toy 68% band")
        return
    hint = ""
    if fitresult_path:
        cov = os.path.join(os.path.dirname(os.path.abspath(fitresult_path)), "cov")
        if os.path.exists(os.path.join(cov, "fitresults.hdf5")):
            hint = f"; the covariance pass is at {cov} — point --{side}-fitresult there for a band"
    print(
        f"[np_function_plots] {side}: NO band — no λ covariance in "
        f"{fitresult_path or 'the raw λ point'} (or every λ frozen){hint}"
    )


def main(argv=None):
    parser = make_parser()
    args = parser.parse_args(argv)

    compare = np_tune.side_requested(args, "den")
    if args.den_result and not args.den_fitresult:
        parser.error(
            "--den-result needs --den-fitresult (no second fitresults to read)"
        )
    if args.den_label and not compare:
        parser.error("--den-label needs one of the other --den-* flags")

    try:
        # A fitresults flag may name the run directory; resolve before any read.
        np_tune.normalize_fitresult_args(args, extra=())
        num_series, num_forms = _series_for_side(args, "num")
        den_series = None
        if compare:
            den_series, _ = _series_for_side(args, "den", inherit_forms=num_forms)
    except ValueError as e:
        parser.error(str(e))

    _report_band(num_series, args.num_fitresult, "num")
    # --num-label names the tune's own curve in EVERY mode, not just the two-tune
    # one: with a single fitresults it renames "postfit" (the prefit curve keeps its
    # own name — it is the card, not the fit).
    if args.num_label:
        num_series[-1].label = args.num_label
    series, insets = num_series, None

    if compare:
        _report_band(den_series, args.den_fitresult, "den")
        num_point, den_point = num_series[-1], den_series[-1]
        # Label each curve by its run directory; a side with no fitresults (a λ
        # scan) has no run to name it with, so name it by what DIFFERS between the
        # two tunes — the shared λ say nothing about which curve is which.
        lam_num, lam_den = np_tune.tune_labels(
            num_point.lam, den_point.lam, fallback="identical tune"
        )
        num_point.label = args.num_label or (
            _run_label(args.num_fitresult) if args.num_fitresult else lam_num
        )
        den_point.label = args.den_label or (
            _run_label(args.den_fitresult) if args.den_fitresult else lam_den
        )
        # Colours are already the two sides' ramps (set in _series_for_side); all that
        # is left is that both postfits draw SOLID, with linestyle saying only
        # prefit-vs-postfit.
        num_point.linestyle = den_point.linestyle = "-"

        # The prefit is the card's λ_central, not the fit's, so ONE curve normally
        # covers both sides. When the two sides' prefit points really do differ,
        # draw both rather than hiding one behind the other.
        pre_num, pre_den = num_series[:-1], den_series[:-1]
        if pre_num and pre_den and pre_num[-1].lam == pre_den[-1].lam:
            pre_den = []
            # Say in the legend that the one dashed curve covers both sides —
            # otherwise its hue (the numerator's) reads as "the numerator's prefit"
            # and the denominator's looks missing.
            for s in pre_num:
                s.label = "prefit (λ_central, both sides)"
        elif pre_num and pre_den:
            diff = np_tune.tune_diff(pre_num[-1].lam, pre_den[-1].lam)
            print(
                "[np_function_plots] the two sides' PREFIT (λ_central) points differ ("
                + ", ".join(
                    f"{k}: {_fmt_diff(v[0])} vs {_fmt_diff(v[1])}"
                    for k, v in diff.items()
                )
                + ") — both prefit curves drawn. Expected if a --den-lambdas / "
                "--den-np-model override moved it (those apply to both fit points "
                "of their side); otherwise the two fits are on different cards."
            )
            # Both dashed: they are told apart by their side's hue, not by a third
            # dash pattern nobody can read at legend size.
            for s in pre_num:
                s.label = f"prefit ({num_point.label})"
            for s in pre_den:
                s.label = f"prefit ({den_point.label})"
        if not args.prefit:
            pre_num, pre_den = [], []
        series = [*pre_num, num_point, *pre_den, den_point]
        insets = [
            (_inset_title("num", num_point.label), num_point),
            (_inset_title("den", den_point.label), den_point),
        ]

        diff = np_tune.tune_diff(num_point.lam, den_point.lam)
        if diff:
            print("[np_function_plots] num vs den, only what differs:")
            for k, (a, b) in diff.items():
                print(f"    {k:<16} {_fmt_diff(a):>12} vs {_fmt_diff(b):>12}")
        else:
            print(
                "[np_function_plots] [warning] the two tunes are IDENTICAL — every "
                "curve will lie on top of its partner."
            )
    elif not args.prefit:
        series = num_series[-1:]

    plot_np_functions(
        series,
        lattice_reference=args.lattice_reference,
        lattice_nsig=args.lattice_nsig,
        map22_reference=args.map22_reference,
        map22_Q=args.map22_Q,
        map22_sqrts=args.map22_sqrts,
        y_values=args.y,
        bT_max=args.bT_max,
        outpath=args.outpath,
        insets=insets,
        f_ymax=args.f_ymax,
        args=args,
    )


if __name__ == "__main__":
    main()
