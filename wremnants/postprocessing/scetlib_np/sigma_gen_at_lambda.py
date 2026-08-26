"""Evaluate the SCETlib NP matched σ_gen at a given λ tune, via the core.

Builds the datacard-free core ``SigmaGenModel`` and runs its matched gen-level
prediction σ_gen(λ; g) = σ_resum(λ; g) + σ_ns on the (ptVGen, absYVGen) gen grid
— Steps 1–2 of ``param_model.py`` (the object that folds through R into the fit,
BEFORE the gen→reco fold and ratio). Prints σ_gen on the gen grid and can plot a
1-D projection (e.g. the ptZ = ptVGen distribution).

Tunes are resolved by the shared ingestor :mod:`np_tune`, which builds validated
:class:`~params.NPTune` objects (a form pair AND the λ those forms use, never
separable). Three of them:

  * **base** — ``--meta-from HDF5`` / the ``--theory-corr`` file's Nonperturbative
    runcard / else the canonical FranksVals tanh_2 default. The model is BUILT
    here (positive σ_gen, which the constructor requires). Its forms are
    **card-locked**: ``--num-np-model`` / ``--den-np-model`` never touch them,
    mirroring ``param_model``, where the R denominator must stay
    histmaker-consistent.
  * **num** (``--num-fitresult`` / ``--num-lambdas`` / ``--num-np-model``) and
    **den** (the matching ``--den-*``) — the tunes EVALUATED. Each is the base plus
    its own fitresults' postfit λ, plus its own ``--*-lambdas``, at its own forms.

λ not set on a side stay at the BASE, so the two sides agree except where you
said otherwise; forms, unlike λ, are inherited den-from-num, so overriding one
form on the denominator cannot silently move the other. Evaluating, unlike
constructing, has no positivity guard, so a weak-NP tune can be inspected;
non-positive bins are warned, not rejected.

Because forms are per-side, a fit's λ can be replayed under a DIFFERENT form —
"the same tune, but with the large-b_T turn-off" is::

    --num-fitresult <fitresults.hdf5> --den-fitresult <same fitresults.hdf5> \\
        --den-np-model tanh_6_sigmoid --den-lambdas bT_cutoff=2,bT_cutoff_width=0.1

λ the requested form needs that the source never stored (here the two cutoff
shape constants) come from the registry defaults, and every fill is reported.

Can ALSO overlay the gen distribution in an official ``TheoryCorrection`` hist
(``--theory-corr``). Those ``.pkl.lz4`` files carry the ``{generator}_hist``
object — the official SCETlib+DYTurbo prediction on a (Q, absY, qT, charge, vars)
grid — so its central (``pdf0``) entry is the same physical object the param-model
σ_gen reconstructs. Overlaying is a direct end-to-end check that the bt-grid +
on-the-fly reconstruction reproduces the official run; pick any ``vars`` label,
e.g. ``lambda21.0``, to overlay a λ-shifted official run against the model at the
matching λ.

Any ``--den-*`` flag turns on the comparison: both tunes are evaluated on the SAME
core (one bt-grid construction) and the output becomes the two σ_gen overlaid with
a ratio panel (num ÷ den). Recommended with ``--native-points``, which keeps qT
native so an unphysical tune's laundered σ(qT) ringing is visible; there the ratio
is masked where the denominator is non-positive and a difference panel is added,
since a difference still reads across those zero-crossings.

The ratio panel's y-range is picked automatically (a small window around 1 in the
binned mode, a 2–98th-percentile clip in the native one). ``--rrange LO HI`` pins
it instead — use it to zoom past a near-zero-crossing spike, or to hold the same
window across several figures being compared by eye. It applies to whichever ratio
panel is drawn (binned, native, or the theory-corr one); the difference panel and
the printed min/max always stay unclipped, so nothing is hidden silently.

``--den-lambdas`` works on its own, with no fitresults on either side: the second
tune is then the base + those λ, exactly as the first is the base +
``--num-lambdas``. That is the λ-scan comparison — e.g.

    --num-lambdas lambda_inf=2.0 --den-lambdas lambda_inf=1.0

overlays σ_gen at the two λ_∞ and ratios them.

The bT-grid is required (``--btgrid``). The gen-bin edges (ptVGen, absYVGen) are
chosen per axis, in order: explicit ``--ptv-edges`` / ``--absy-edges``, then a
``--gen-edges-from`` / ``--datacard`` hdf5, then a built-in default (1-GeV ptVGen
bins over [0, 40]; a single rapidity-inclusive absYVGen bin [0, 5]) — so the
script runs with no gen-edge input. ``--datacard`` feeds ONLY the gen edges; λ
are not sourced from it (use ``--meta-from``).

Run inside a container that binds the inputs (same as the validation scripts):

    export APPTAINER_BIND="/scratch,/cvmfs,/work,/ceph,/home"
    singularity run --cleanenv <wmassdevrolling> bash -c \\
      "source main/WRemnants/setup.sh; \\
       python3 -m wremnants.postprocessing.scetlib_np.sigma_gen_at_lambda \\
         --num-lambdas lambda2=0.4,lambda4=0.1,lambda2_nu=0.15 \\
         --theory-corr <wremnants-data>/data/TheoryCorrections/scetlib_dyturbo_..._CorrZ.pkl.lz4 \\
         --plot ~/public_html/alphaS/YYMMDD_sigmagen/ptZ.png"
"""

import argparse
import sys
import time

import numpy as np

from wremnants.postprocessing.scetlib_np import np_tune
from wremnants.postprocessing.scetlib_np.params import EFF_PARAMS, GNU_PARAMS
from wremnants.postprocessing.scetlib_np.sigma_gen import (
    _NONSING_DYTURBO_DEFAULT,
    _NONSING_FO_SING_DEFAULT,
    _default_btgrid_dir,
)

# The theory-correction reference loader/projection lives in the shared validation
# library so the CLIs share one implementation; the λ-tune resolution now lives in
# :mod:`np_tune` (one ingestor for every CLI).
from wremnants.postprocessing.scetlib_np.validation.agreement import (
    Q_HI,
    Q_LO,
    load_theory_corr_hist,
    resolve_gen_axes,
    theory_corr_projection,
)


def _lambda_box_text(eff, gnu):
    """Compact multi-line λ-tune annotation: np_model form(s) + non-zero params."""
    lines = [f"np_model = {eff.get('np_model')}"]
    if gnu.get("np_model_nu") != eff.get("np_model"):
        lines.append(f"np_model_nu = {gnu.get('np_model_nu')}")
    for p in GNU_PARAMS:
        if abs(gnu.get(p, 0.0)) > 0:
            lines.append(f"{p} = {gnu[p]:.4g}")
    for p in EFF_PARAMS:
        if abs(eff.get(p, 0.0)) > 0:
            lines.append(f"{p} = {eff[p]:.4g}")
    return "λ tune:\n" + "\n".join(lines)


def make_projection_plot(
    sigma_gen,
    gen_axes,
    axis,
    out_path,
    eff,
    gnu,
    s_corr=None,
    corr_label=None,
    model_label=None,
    eff_den=None,
    gnu_den=None,
    args=None,
    rrange=None,
):
    """Step histogram of the matched σ_gen(λ) projection onto one gen axis
    (default ptVGen = ptZ), summing over the other.

    With an overlay ``s_corr`` (bin-integrated σ on the SAME ``axis`` edges — a
    TheoryCorrection projection, or the second tune's σ_gen), it is drawn on top
    and a ratio panel is added below. Without ``s_corr`` the figure is a single
    panel. Values are plotted as DIFFERENTIAL dσ/dx (bin-integrated σ ÷ bin width)
    so the variable binning, notably the wide ptVGen overflow bin, reads correctly;
    the ratio is width-independent.

    ``rrange=(lo, hi)`` pins the ratio panel's y-range instead of the automatic
    window around 1.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = [n for n, _ in gen_axes]
    if axis not in names:
        raise SystemExit(f"--plot-axis {axis!r} not in gen axes {names}")
    ai = names.index(axis)
    other = 1 - ai  # exactly 2 gen axes (ptVGen, absYVGen)
    edges = np.asarray(gen_axes[ai][1], dtype=np.float64)
    widths = np.diff(edges)
    s = sigma_gen.sum(axis=other)
    ds = s / widths
    corr_label = corr_label or "SCETlib+DYTurbo"
    show_ratio = s_corr is not None
    if show_ratio:
        ratio = np.divide(s, s_corr, out=np.ones_like(s), where=s_corr != 0)
        fig, (ax, axr) = plt.subplots(
            2,
            1,
            sharex=True,
            figsize=(7, 6),
            gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06},
        )
    else:
        fig, ax = plt.subplots(figsize=(7, 5))
        axr = None

    lab = "p$_T^Z$ (ptVGen) [GeV]" if axis == "ptVGen" else axis
    ax.stairs(
        ds, edges, color="C3", lw=1.6, label=model_label or "σ_gen(λ) (param model)"
    )
    if s_corr is not None:
        ax.stairs(
            s_corr / widths, edges, color="C0", lw=1.6, ls=(0, (4, 2)), label=corr_label
        )
    ax.set_ylabel(r"d$\sigma_{\mathrm{gen}}$/d(" + axis + ")  [a.u.]")
    ax.margins(x=0)
    ax.legend(loc="upper right", fontsize=9)
    if eff_den is not None:
        # Two tunes: a box per tune, matching the native plot — num on the right
        # (mistyrose/C3, the model curve), den on the left (aliceblue/C0, overlay).
        ax.text(
            0.975,
            0.58,
            "num  " + _lambda_box_text(eff, gnu),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7,
            family="monospace",
            bbox=dict(boxstyle="round", fc="mistyrose", ec="C3", alpha=0.85),
        )
        ax.text(
            0.02,
            0.58,
            "den  " + _lambda_box_text(eff_den, gnu_den),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7,
            family="monospace",
            bbox=dict(boxstyle="round", fc="aliceblue", ec="C0", alpha=0.85),
        )
    else:
        ax.text(
            0.975,
            0.60,
            _lambda_box_text(eff, gnu),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7.5,
            family="monospace",
            bbox=dict(boxstyle="round", facecolor="white", edgecolor="0.7", alpha=0.9),
        )

    if show_ratio:
        # Short, fixed panel label — the curves are already named in the legend and
        # the λ boxes, and inlining the tune names here makes the y-label unreadable.
        axr.stairs(ratio, edges, color="k", lw=1.4)
        axr.axhline(1.0, color="0.5", lw=0.8, ls="--")
        axr.set_ylabel("num / den" if model_label else "model / corr")
        axr.set_xlabel(lab)
        axr.margins(x=0)
        rlo, rhi = float(np.min(ratio)), float(np.max(ratio))
        if rrange is not None:
            axr.set_ylim(*rrange)
        else:
            # Zoom around 1 to show the (often sub-%) residual, keeping 1.0 in frame;
            # small window if the ratio is flat.
            pad = max((rhi - rlo) * 0.25, 0.003)
            axr.set_ylim(min(rlo, 1.0) - pad, max(rhi, 1.0) + pad)
        rng = f"; ratio [{rlo:.4f}, {rhi:.4f}]"
        if rrange is not None:
            # The printed range is the true one; say so when the panel crops it.
            n_out = int(np.sum((ratio < rrange[0]) | (ratio > rrange[1])))
            rng += f", panel --rrange [{rrange[0]:g}, {rrange[1]:g}]" + (
                f", {n_out}/{ratio.size} bins outside" if n_out else ""
            )
    else:
        ax.set_xlabel(lab)
        rng = ""

    from wremnants.postprocessing.scetlib_np import plot_output

    outdir, basename = plot_output.split_outpath(out_path)
    plot_output.save_plot(outdir, basename, fig=fig, args=args, dpi=130)
    plt.close(fig)
    print(
        f"[plot] wrote {outdir}/{basename}.png(.pdf) + {basename}.log  "
        f"(axis={axis}, summed over {names[other]}{rng})"
    )


def make_native_projection_plot(
    sigma_YqT,
    Y_unique,
    qT_unique,
    W_absY,
    W_ptVGen,
    axis,
    out_path,
    eff,
    gnu,
    args=None,
    sigma_YqT_den=None,
    num_label=None,
    den_label=None,
    eff_den=None,
    gnu_den=None,
    rrange=None,
):
    """Line plot of σ_gen(λ) with the plotted axis kept on its NATIVE btgrid grid.

    The plotted axis (``ptVGen`` = qT, or ``absYVGen`` = signed Y) is NOT rebinned;
    the OTHER spatial axis is integrated with the same ``W_absY`` / ``W_ptVGen``
    weights the binned path uses (i.e. over the --absy-edges / --ptv-edges range),
    and Q is already folded into ``sigma_YqT`` (--q-lo/--q-hi). This keeps the
    plotted axis at full resolution so a sub-bin σ(qT)<0 dip — laundered by the
    qT->ptVGen rebin in :func:`make_projection_plot` — is visible.

    With a second tune ``sigma_YqT_den`` (same grid) the figure becomes a
    comparison: the two native σ(qT) overlaid (top), a ratio numerator ÷
    denominator (middle) MASKED where the denominator is ≤ a small positive floor
    — a native ratio is meaningless across the unphysical tune's zero-crossings —
    and a difference numerator − denominator (bottom).

    ``rrange=(lo, hi)`` pins the ratio panel's y-range instead of the robust
    percentile clip; the difference panel is left alone.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    wY = np.asarray(W_absY, dtype=np.float64).sum(axis=0)  # (NY,)
    wqT = np.asarray(W_ptVGen, dtype=np.float64).sum(axis=0)  # (NqT,)

    def _project(s):
        # s is (NY, NqT), Q already integrated; fold the OTHER axis over its edges.
        s = np.asarray(s, dtype=np.float64)
        return wY @ s if axis == "ptVGen" else s @ wqT

    if axis == "ptVGen":  # keep qT native, integrate Y over the absY edges (fold)
        x = np.asarray(qT_unique, dtype=np.float64)
        xlab = r"p$_T^Z$ = q$_T$ (native) [GeV]"
    else:  # absYVGen: keep native (signed) Y, integrate qT over the ptV edges
        x = np.asarray(Y_unique, dtype=np.float64)
        xlab = "Y (native, signed)"
    order = np.argsort(x)
    x = x[order]
    y = _project(sigma_YqT)[order]
    nneg = int(np.sum(y < 0))
    ylab = r"σ$_{\mathrm{gen}}$ (native; other axes integrated)  [a.u.]"
    has_ratio = sigma_YqT_den is not None

    if has_ratio:
        y_den = _project(sigma_YqT_den)[order]
        num_label = num_label or "σ_gen (numerator)"
        den_label = den_label or "σ_gen (denominator)"
        # Ratio masked ONLY where the denominator is non-positive — the
        # unphysical-tune zero-crossings where a ratio flips sign / diverges. The
        # native σ(qT) spectrum falls orders of magnitude, so a peak-relative floor
        # would swallow the (positive, physical) tail; the y-range is instead
        # robustly clipped below so a near-zero-crossing spike doesn't compress it.
        good = y_den > 0
        ratio = np.full_like(y, np.nan)
        np.divide(y, y_den, out=ratio, where=good)
        diff = y - y_den
        n_masked = int(np.sum(~good))

        fig, (ax, axr, axd) = plt.subplots(
            3,
            1,
            sharex=True,
            figsize=(10, 9),
            gridspec_kw={"height_ratios": [3, 1.2, 1.2], "hspace": 0.07},
        )
        ax.plot(x, y, color="C3", lw=1.5, marker=".", ms=3, label=num_label)
        ax.plot(
            x,
            y_den,
            color="C0",
            lw=1.5,
            ls=(0, (4, 2)),
            marker=".",
            ms=3,
            label=den_label,
        )
        ax.axhline(0.0, color="0.4", lw=0.8, ls="--")
        ax.set_ylabel(ylab)
        ax.margins(x=0)
        ax.legend(loc="upper right", fontsize=9)
        ax.text(
            0.975,
            0.58,
            "num  " + _lambda_box_text(eff, gnu),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7,
            family="monospace",
            bbox=dict(boxstyle="round", fc="mistyrose", ec="C3", alpha=0.85),
        )
        if eff_den is not None:
            ax.text(
                0.02,
                0.58,
                "den  " + _lambda_box_text(eff_den, gnu_den),
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=7,
                family="monospace",
                bbox=dict(boxstyle="round", fc="aliceblue", ec="C0", alpha=0.85),
            )
        axr.plot(x, ratio, color="k", lw=1.4, marker=".", ms=3)
        axr.axhline(1.0, color="0.5", lw=0.8, ls="--")
        axr.set_ylabel("num / den")
        axr.margins(x=0)
        finite = ratio[np.isfinite(ratio)]
        if rrange is not None:
            axr.set_ylim(*rrange)
        elif finite.size:
            # Robust y-range (2–98th pct) so a spike at a near-zero-crossing of the
            # denominator does not compress the bulk of the spectrum's ratio.
            rlo, rhi = (float(v) for v in np.percentile(finite, [2, 98]))
            pad = max((rhi - rlo) * 0.1, 0.02)
            axr.set_ylim(min(rlo, 1.0) - pad, max(rhi, 1.0) + pad)
        axd.plot(x, diff, color="C2", lw=1.4, marker=".", ms=3)
        axd.axhline(0.0, color="0.5", lw=0.8, ls="--")
        axd.set_ylabel("num − den")
        axd.set_xlabel(xlab)
        axd.margins(x=0)
        tail = (
            f"; {n_masked}/{y_den.size} ratio pts masked (den≤0); " f"num {nneg} pts<0"
        )
        if rrange is not None:
            # Say how much the pinned panel crops, so a zoom can't quietly hide the
            # ringing this mode exists to show.
            n_out = int(np.sum((finite < rrange[0]) | (finite > rrange[1])))
            tail += f"; panel --rrange [{rrange[0]:g}, {rrange[1]:g}]" + (
                f", {n_out}/{finite.size} pts outside" if n_out else ""
            )
    else:
        fig, ax = plt.subplots(figsize=(10, 7))
        ax.plot(x, y, color="C3", lw=1.5, marker=".", ms=3, label="σ_gen(λ)")
        ax.axhline(0.0, color="0.4", lw=0.8, ls="--")
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.margins(x=0)
        ax.legend(loc="upper right", fontsize=9)
        ax.text(
            0.975,
            0.60,
            _lambda_box_text(eff, gnu),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.6", alpha=0.85),
        )
        tail = f"; {nneg}/{y.size} points < 0"

    from wremnants.postprocessing.scetlib_np import plot_output

    outdir, basename = plot_output.split_outpath(out_path)
    plot_output.save_plot(outdir, basename, fig=fig, args=args, dpi=200)
    plt.close(fig)
    print(
        f"[plot] wrote {outdir}/{basename}.png(.pdf) + {basename}.log  "
        f"(NATIVE axis={axis}{tail})"
    )


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--btgrid", default=_default_btgrid_dir(), help="SCETlib bT-grid directory"
    )
    p.add_argument(
        "--datacard",
        default=None,
        help="hdf5 fallback for the GEN EDGES (no default). λ are NOT "
        "sourced from it — use --meta-from for that",
    )
    # λ tune: the CONSTRUCTION base (card-locked forms), then one group per side.
    p.add_argument(
        "--meta-from",
        default=None,
        help="hdf5 to read the base λ tune from (datacard/fitresults metadata); optional",
    )
    np_tune.add_tune_args(p, "num")
    # gen-edge source
    p.add_argument(
        "--ptv-edges",
        default=None,
        help="ptVGen edges 'a,b,c,...' (default: 1-GeV bins over [0,40])",
    )
    p.add_argument(
        "--absy-edges",
        default=None,
        help="absYVGen edges 'a,b,c,...' (default: single bin [0,5])",
    )
    p.add_argument(
        "--gen-edges-from",
        default=None,
        help="hdf5 whose scetlib_np auxiliary gives the gen edges "
        "(default: --datacard, else the built-in defaults)",
    )
    # TheoryCorrection overlay
    p.add_argument(
        "--theory-corr",
        default=None,
        help="TheoryCorrection .pkl.lz4 to overlay the official "
        "SCETlib+DYTurbo gen distribution (its {generator}_hist)",
    )
    p.add_argument(
        "--theory-corr-proc",
        default=None,
        help="proc key in the corr file (default: the single physics key)",
    )
    p.add_argument(
        "--theory-corr-var",
        default="pdf0",
        help="vars label to read from the corr hist (default pdf0 = central)",
    )
    p.add_argument(
        "--theory-corr-normalize",
        action="store_true",
        help="rescale the corr curve to the model σ_gen integral "
        "(shape-only comparison; default off = absolute overlay)",
    )
    # model / output
    p.add_argument("--q-lo", type=float, default=Q_LO)
    p.add_argument("--q-hi", type=float, default=Q_HI)
    p.add_argument(
        "--no-nonsingular",
        action="store_true",
        help="resum-only σ_gen (σ_ns = 0; skips the FO inputs)",
    )
    p.add_argument(
        "--nonsingular-fo-sing",
        default=_NONSING_FO_SING_DEFAULT,
        help="SCETlib singular …_nnlo_sing…combined.pkl for σ_ns = DYTurbo − singular "
        "(default: CT18Z; set the MSHT20/other-PDF path to MATCH the --btgrid PDF)",
    )
    p.add_argument(
        "--nonsingular-dyturbo",
        default=_NONSING_DYTURBO_DEFAULT,
        help="full FO prediction for σ_ns: either the DYTurbo "
        "results_…-{scale}-scetlibmatch.txt, or an NNLOjet stitched-export base "
        "name (…/nnlojet_export_stitched/ptz, the per-y __ylo__yhi.dat slices are "
        "stitched and a [q-lo, q-hi] Q bin inserted) "
        "(default: CT18Z DYTurbo; set the path MATCHING the --btgrid PDF)",
    )
    p.add_argument(
        "--plot",
        default=None,
        help="optional path (e.g. .png/.pdf) to write a 1-D projection plot "
        "of σ_gen(λ) [+ theory-corr overlay/ratio]; see --plot-axis",
    )
    p.add_argument(
        "--plot-axis",
        default="ptVGen",
        choices=["ptVGen", "absYVGen"],
        help="gen axis to project onto for --plot (default ptVGen = ptZ)",
    )
    p.add_argument(
        "--rrange",
        type=float,
        nargs=2,
        metavar=("LO", "HI"),
        default=None,
        help="y range for the ratio panel of --plot (default: automatic). Applies "
        "to whichever ratio is drawn (num/den, binned or --native-points, or the "
        "--theory-corr one); the printed min/max stay unclipped",
    )
    p.add_argument(
        "--native-points",
        action="store_true",
        help="plot --plot-axis on its NATIVE btgrid grid (skip that axis's rebin); "
        "the OTHER axis and Q are still integrated over the specified edges "
        "(--absy-edges / --ptv-edges / --q-lo/--q-hi). Reveals sub-bin negativity "
        "the gen rebin launders. No TheoryCorrection overlay in this mode.",
    )
    # Ratio of TWO tunes (numerator = --num-fitresult / --num-lambdas, denominator
    # = --den-fitresult / --den-lambdas), both evaluated on the SAME core. Either
    # --den flag turns the comparison on.
    np_tune.add_tune_args(p, "den")
    args = p.parse_args(argv)
    # The second tune is requested by EITHER flag: --den-fitresult (its postfit) or
    # --den-lambdas alone (base + those λ). Everything downstream gates on this, not
    # on --den-fitresult, so the λ-only comparison gets the same overlay/ratio/diff.
    do_ratio = np_tune.side_requested(args, "den")
    if do_ratio and args.theory_corr:
        p.error(
            "--den-* (--den-fitresult/--den-lambdas) and --theory-corr both define an "
            "overlay; pick one"
        )
    if args.den_result and not args.den_fitresult:
        p.error("--den-result needs --den-fitresult (no second fitresults to read)")
    if args.den_label and not do_ratio:
        p.error("--den-label needs one of the other --den-* flags")
    if args.rrange is not None and args.rrange[0] >= args.rrange[1]:
        p.error(f"--rrange needs LO < HI (got {args.rrange[0]}, {args.rrange[1]})")

    # ---- λ: build a PHYSICAL base tune (model CONSTRUCTED there so the
    # positive-σ_gen guard passes), then EVALUATE at base + the requested
    # overrides (--num-fitresult postfit, then --num-lambdas). Params not set stay at the
    # base, NOT at 0.
    import os

    if args.num_fitresult and not args.meta_from:
        # the fitresults carries the card λ_central metadata: construct there
        args.meta_from = args.num_fitresult

    try:
        # CONSTRUCTION tune: the card's own λ AND its own forms. --*-np-model does
        # NOT reach here — a CLI form stamped onto the card's λ would declare a
        # form whose parameters the card never carried.
        base_tune = np_tune.resolve_base_tune(args)
        num_tune, _ = np_tune.resolve_side_tune(args, "num", base_tune)
        den_tune = None
        if do_ratio:
            # inherit=num_tune: the denominator's FORMS default to the numerator's
            # (its λ still fall back to the base), so overriding one form here
            # cannot silently move the other.
            den_tune, _ = np_tune.resolve_side_tune(
                args, "den", base_tune, inherit=num_tune
            )
    except ValueError as e:
        p.error(str(e))

    base = base_tune.as_lambda_central()
    eff, gnu = num_tune.eff, num_tune.gnu
    eval_np_model, eval_np_model_nu = num_tune.np_model, num_tune.np_model_nu
    explicit = num_tune != base_tune

    eff2 = gnu2 = eval_np_model2 = eval_np_model_nu2 = None
    num_label = den_label = None
    if do_ratio:
        eff2, gnu2 = den_tune.eff, den_tune.gnu
        eval_np_model2, eval_np_model_nu2 = den_tune.np_model, den_tune.np_model_nu
        # Legend labels: the fitresults directory name where there is one, else the
        # λ/forms that DIFFER between the tunes (λ-scan mode has no run directory to
        # name it with, and the shared λ say nothing about which curve is which).
        lam_num, lam_den = np_tune.tune_labels(num_tune, den_tune)
        num_label = args.num_label or (
            os.path.basename(os.path.dirname(os.path.abspath(args.num_fitresult)))
            if args.num_fitresult
            else lam_num
        )
        den_label = args.den_label or (
            os.path.basename(os.path.dirname(os.path.abspath(args.den_fitresult)))
            if args.den_fitresult
            else lam_den
        )
        if num_tune == den_tune:
            print(
                "  [warning] the two tunes are IDENTICAL — every ratio will be 1. "
                "Each side is the base plus its OWN overrides, so give the differing "
                "λ to --num-lambdas and --den-lambdas."
            )

    gen_axes = resolve_gen_axes(args)

    from wremnants.postprocessing.scetlib_np.sigma_gen import SigmaGenModel

    print(
        "\n[core] constructing SigmaGenModel at the base tune (bt-grid integral) …",
        flush=True,
    )
    t0 = time.time()
    core = SigmaGenModel(
        btgrid_dir=args.btgrid,
        lambda_central=base,
        gen_axes=gen_axes,
        Q_lo=args.q_lo,
        Q_hi=args.q_hi,
        include_nonsingular=not args.no_nonsingular,
        nonsingular_fo_sing=args.nonsingular_fo_sing,
        nonsingular_dyturbo=args.nonsingular_dyturbo,
    )
    print(
        f"  constructed in {time.time()-t0:.1f}s; gen grid {core.gen_shape} "
        f"({[n for n, _ in core.gen_axes]})"
    )

    print(f"\n[λ] evaluating matched σ_gen at:")
    print(f"  F_eff  : {eff}")
    print(f"  γ_ν^NP : {gnu}")
    if explicit:
        diff = np_tune.tune_diff(num_tune, base_tune)
        print(
            "  (differs from the base tune in: "
            + ", ".join(f"{k}={v[0]}" for k, v in diff.items())
            + "; the rest stay at the base)"
        )
    else:
        print("  (no overrides — evaluating at the base tune itself)")

    # ---- σ_gen on the (ptVGen, absYVGen) gen grid. Reuse the construction central
    # when no overrides; else evaluate the tune (no positivity guard, so a weak-NP
    # tune can be inspected — it can dip negative at the lowest qT).
    t0 = time.time()
    if explicit:
        sigma_gen = np.asarray(
            core.sigma_gen(
                eff, gnu, np_model=eval_np_model, np_model_nu=eval_np_model_nu
            ).numpy(),
            dtype=np.float64,
        )
    else:
        sigma_gen = np.asarray(core.sigma_gen_central.numpy(), dtype=np.float64)
    print(f"  σ_gen computed in {time.time()-t0:.1f}s; shape {sigma_gen.shape}")

    n_bad = int(np.sum(sigma_gen <= 0))
    if n_bad:
        print(
            f"  [warning] {n_bad}/{sigma_gen.size} σ_gen bins are non-positive at "
            f"this tune (expected where the NP damping is weak, esp. low qT)."
        )

    print(f"\n  Σ σ_gen        = {sigma_gen.sum():.6g}")
    print(f"  per-bin σ_gen  : min {sigma_gen.min():.4g}  max {sigma_gen.max():.4g}")
    print("\n  σ_gen(λ) per (ptVGen × absY) bin:")
    with np.printoptions(precision=4, suppress=True, linewidth=140):
        print(sigma_gen)

    # ---- second tune (--den-fitresult / --den-lambdas): evaluate on the same core,
    # print the binned ratio (numerator ÷ denominator) on the gen grid, which —
    # unlike the native ratio — is laundered positive and well-behaved.
    sigma_gen2 = None
    if do_ratio:
        t0 = time.time()
        sigma_gen2 = np.asarray(
            core.sigma_gen(
                eff2, gnu2, np_model=eval_np_model2, np_model_nu=eval_np_model_nu2
            ).numpy(),
            dtype=np.float64,
        )
        print(f"\n[den] σ_gen(den) computed in {time.time()-t0:.1f}s")
        print(f"  F_eff (den)  : {eff2}")
        print(f"  γ_ν^NP (den) : {gnu2}")
        rb = np.where(sigma_gen2 != 0, sigma_gen / sigma_gen2, np.nan)
        print(
            f"  Σ num / Σ den   = {sigma_gen.sum()/sigma_gen2.sum():.5f}  "
            f"(num={num_label}, den={den_label})"
        )
        print(f"  binned num/den  : min {np.nanmin(rb):.4f}  max {np.nanmax(rb):.4f}")
        with np.printoptions(precision=4, suppress=True, linewidth=140):
            print("  binned num/den per (ptVGen × absY) bin:")
            print(rb)

    # ---- NP physical-validity detectors at THIS λ. A wrong-sign tune's pathology
    # (anti-damping NP → oscillating, negative native σ(qT)) is AVERAGED AWAY in
    # the binned σ_gen above, so check the native spectrum and form factors
    # directly. Detectors only — change nothing; fit-time enforcement is the
    # np_damping_wall.NPDampingWall regularizer.
    from wremnants.postprocessing.scetlib_np import param_model_diagnostics as ppd

    rep = ppd.np_physical_report(
        core, eff, gnu, np_model=eval_np_model, np_model_nu=eval_np_model_nu
    )
    damp, neg = rep["damp"], rep["neg"]
    print("\n  NP physical-validity:")
    print(
        f"    γ_ν^NP damping : {'OK' if not damp['gamma_nu_wrong_sign'] else 'WRONG SIGN'}"
        f"  (max γ_ν over probe bT = {damp['gamma_nu_max']:+.3g}; must be ≤ 0)"
    )
    print(
        f"    F_eff decays   : {'OK' if not damp['F_eff_growing'] else 'GROWING (bT-integral divergence sign)'}"
    )
    print(
        f"    native σ(qT)≥0 : neg_area_frac={neg['neg_area_frac']:.3g} "
        f"(λ_central {rep['central_neg_area']:.3g}), min/peak={neg['min_over_peak']:+.3g}, "
        f"n_neg_bins={neg['n_neg_bins']}"
    )

    # WHERE the negativity sits (native Y × qT). The σ(qT)<0 dip is laundered by
    # the gen-binning, so locating it is the discriminator: negative cells inside
    # the |Y|≤2.5 acceptance at accessible qT matter for interpretation; cells only
    # at |Y|>2.5 or beyond the fit's qT reach are doubly invisible. ACCEPT_ABSY is
    # the Z dilepton |yll| acceptance edge (boson-Y proxy).
    ACCEPT_ABSY = 2.5
    if neg.get("n_neg_bins") and neg.get("worst") is not None:
        w = neg["worst"]
        print(
            f"    worst neg cell : σ={w['value']:.4g} ({w['frac_of_peak']*100:+.1f}% of peak) "
            f"at Y={w['Y']:+.3g}, qT={w['qT']:.3g} GeV"
        )
        cells = neg.get("neg_bins") or []
        in_cells = [c for c in cells if abs(c["Y"]) <= ACCEPT_ABSY]
        n_in, n_out = len(in_cells), len(cells) - len(in_cells)
        qr = neg.get("neg_qT_range")
        print(
            f"    neg-cell split : {n_in} inside |Y|≤{ACCEPT_ABSY}, {n_out} outside "
            f"(|Y|≤{neg.get('neg_absY_max', float('nan')):.2g} reached); "
            f"qT∈[{qr[0]:.3g}, {qr[1]:.3g}] GeV"
            if qr
            else ""
        )
        if in_cells:
            wi = min(in_cells, key=lambda c: c["value"])
            print(
                f"    worst in-acc   : σ={wi['value']:.4g} ({wi['frac_of_peak']*100:+.1f}% of peak) "
                f"at Y={wi['Y']:+.3g}, qT={wi['qT']:.3g} GeV   "
                f"[the |Y|≤{ACCEPT_ABSY} pathology the fit region could see]"
            )
        n_show = min(10, len(cells))
        if n_show:
            print(f"    most-negative {n_show} cells (Y, qT[GeV], σ, %peak):")
            for c in cells[:n_show]:
                inacc = "in " if abs(c["Y"]) <= ACCEPT_ABSY else "out"
                print(
                    f"        [{inacc}] Y={c['Y']:+.3g}  qT={c['qT']:7.3g}  "
                    f"σ={c['value']:+.4g}  ({c['frac_of_peak']*100:+.1f}%)"
                )

    # ---- *_abs (damping-fold) forms: the damping probes above are satisfied BY
    # CONSTRUCTION, so they are no evidence about the λ. Report what the fold is
    # actually doing — an inactive fold means this really is a tanh_6 tune, an
    # active one means non-damping λ held down by the fold alone.
    fold = rep.get("fold") or {}
    if fold:
        print(
            f"    damping fold   : {'INACTIVE (≡ tanh_6 tune)' if fold['inactive'] else 'ACTIVE (λ non-damping)'}"
        )
        if "tmd_fold_frac" in fold:
            print(
                f"      TMD  folded b_T fraction {fold['tmd_fold_frac']:.3f} "
                f"{fold['tmd_fold_bT_range']}  max F_eff={fold['F_eff_max']:.4g} "
                f"(unfolded {fold['F_eff_max_bare']:.4g})"
            )
        if "cs_fold_frac" in fold:
            print(
                f"      CS   folded b_T fraction {fold['cs_fold_frac']:.3f} "
                f"{fold['cs_fold_bT_range']}  max γ_ν={fold['gamma_nu_max']:+.4g} "
                f"(unfolded {fold['gamma_nu_max_bare']:+.4g})"
            )

    if not rep["ok"]:
        print(
            "    ⚠  UNPHYSICAL NP TUNE — the differential σ(qT) is negative / the NP is "
            "anti-damping.\n       This is hidden in the binned σ_gen above; do not treat "
            "this point as a physical prediction."
        )

    # ---- optional: project an official TheoryCorrection run onto the same axis.
    s_corr = None
    corr_label = None
    if args.theory_corr:
        h_corr = load_theory_corr_hist(args.theory_corr, args.theory_corr_proc)
        s_corr = theory_corr_projection(
            h_corr,
            core.gen_axes,
            args.plot_axis,
            var=args.theory_corr_var,
            q_window=(args.q_lo, args.q_hi),
        )
        other = 1 - [n for n, _ in core.gen_axes].index(args.plot_axis)
        s_model = sigma_gen.sum(axis=other)
        sum_corr, sum_model = float(s_corr.sum()), float(s_model.sum())
        norm_label = ""
        if args.theory_corr_normalize and sum_corr != 0:
            scale = sum_model / sum_corr
            s_corr = s_corr * scale
            norm_label = f", ×{scale:.4f} (shape-normalized)"
            print(
                f"\n[theory-corr] shape-normalized to the model integral (×{scale:.5f})"
            )
        corr_label = f"SCETlib+DYTurbo ({args.theory_corr_var}){norm_label}"
        rcorr = np.divide(s_model, s_corr, out=np.ones_like(s_model), where=s_corr != 0)
        print(
            f"\n[theory-corr] projection onto {args.plot_axis} (var={args.theory_corr_var}):"
        )
        print(f"  Σ corr           = {sum_corr:.6g}")
        print(
            f"  Σ model / Σ corr = {sum_model/sum_corr:.5f}"
            + (
                "  (== 1 after --theory-corr-normalize)"
                if args.theory_corr_normalize
                else ""
            )
        )
        print(f"  model / corr per bin : min {rcorr.min():.4f}  max {rcorr.max():.4f}")
        with np.printoptions(precision=4, suppress=True, linewidth=140):
            print(f"  model/corr: {rcorr}")
        if not args.plot:
            print("[theory-corr] (pass --plot to also write the overlay figure)")

    if args.plot:
        if args.native_points:
            # Keep the plotted axis native (no rebin); integrate the other axis
            # + Q over the specified edges. Reveals the sub-bin σ(qT)<0 dip the
            # gen rebin launders. TheoryCorrection overlay (s_corr) is skipped;
            # with a second tune it is overlaid + ratio'd + differenced.
            sigma_YqT = core.sigma_YqT_native(
                eff, gnu, np_model=eval_np_model, np_model_nu=eval_np_model_nu
            )
            sigma_YqT_den = None
            if do_ratio:
                sigma_YqT_den = core.sigma_YqT_native(
                    eff2, gnu2, np_model=eval_np_model2, np_model_nu=eval_np_model_nu2
                )
            make_native_projection_plot(
                sigma_YqT,
                core.Y_unique,
                core.qT_unique,
                core.W_absY,
                core.W_ptVGen,
                args.plot_axis,
                args.plot,
                eff,
                gnu,
                args=args,
                sigma_YqT_den=sigma_YqT_den,
                num_label=num_label,
                den_label=den_label,
                eff_den=eff2,
                gnu_den=gnu2,
                rrange=args.rrange,
            )
        else:
            print("Making plot")
            # A second tune overlays its projection with a num/den ratio panel
            # (takes precedence over the mutually-exclusive theory-corr).
            if do_ratio:
                other = 1 - [n for n, _ in core.gen_axes].index(args.plot_axis)
                s_corr = sigma_gen2.sum(axis=other)
                corr_label = den_label
            make_projection_plot(
                sigma_gen,
                core.gen_axes,
                args.plot_axis,
                args.plot,
                eff,
                gnu,
                s_corr=s_corr,
                corr_label=corr_label,
                model_label=num_label if do_ratio else None,
                eff_den=eff2 if do_ratio else None,
                gnu_den=gnu2 if do_ratio else None,
                args=args,
                rrange=args.rrange,
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
