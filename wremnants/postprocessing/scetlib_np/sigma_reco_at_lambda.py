"""Evaluate the SCETlib NP matched σ_reco at postfit tune(s) and overlay data.

Reco-level companion to :mod:`sigma_gen_at_lambda`. That script stops at the
gen-level σ_gen (Steps 1–2 of :mod:`param_model`); this one carries it one step
further — folds σ_gen through the response matrix R (Step 3) to get σ_reco on the
fit's reco grid — and overlays the OBSERVED DATA:

    σ_reco(λ; b) = Σ_g  R(b, g) · σ_gen(λ; g)          (Step 3, gen → reco)

so the prediction lives on the SAME (ptll, yll) reco bins as the data and can be
compared to it directly (unlike σ_gen, whose ptVGen axis has no data counterpart).

Inputs, and why they are where they are:
  * ``--datacard`` (required) — the rabbit input tensor. Supplies EVERYTHING reco:
    the response matrix R + gen-total N_gen (its ``scetlib_np`` auxiliary), the
    λ_central BASE tune (its metadata), the reco axes, AND the observed data
    (``hdata_obs``). The data is NOT in the fitresults (those are typically written
    with ``saveHists=False``); it lives in the datacard the fit read. Two fits that
    share a datacard therefore share the data by construction.
  * the λ tunes, resolved by the shared ingestor :mod:`np_tune` exactly as in
    :mod:`sigma_gen_at_lambda`: ``--num-*`` (``--num-fitresult`` /
    ``--num-lambdas`` / ``--num-np-model`` / ``--num-label``) and the matching
    ``--den-*``. Any ``--den-*``
    adds a SECOND tune on the same core (one R, one bt-grid), for a
    wall-vs-nowall style overlay. Forms are per-side, so the same λ can be
    evaluated under two different NP forms.

Normalization: σ_reco is a normalized (weighted) cross section — its integral is
~O(400), not the ~7 M event yield of the data — and the real fit floats the signal
normalization, so ONLY THE SHAPE is physical. Each prediction is therefore scaled
so Σ σ_reco = Σ data over the plotted range (``--data-norm shape``, the default): a
pure shape comparison. ``--data-norm none`` leaves σ_reco raw (data will not
overlay). The prediction is SIGNAL-only (Z; the ~0.3 % non-signal background and
the postfit detector nuisances are not added — those are common to the tunes and
absorbed into the shape normalization); for the full postfit Data-vs-Pred use
rabbit's ``rabbit_plot_hists.py -m 'Project ch0 ptll'`` on a ``--saveHists`` run.

Run inside the container that binds the inputs (same as the validation scripts):

    export APPTAINER_BIND="/scratch,/cvmfs,/work,/ceph,/home"
    singularity run --cleanenv <wmassdevrolling> bash -c \\
      "source main/WRemnants/setup.sh; \\
       python3 -m wremnants.postprocessing.scetlib_np.sigma_reco_at_lambda \\
         --datacard <ZMassDilepton.hdf5> \\
         --num-fitresult <nowall/fitresults.hdf5> --num-label 'no wall' \\
         --den-fitresult <wall/fitresults.hdf5>   --den-label wall \\
         --plot-axis ptll \\
         --plot ~/public_html/alphaS/YYMMDD_sigmareco/ptll_data.png"
"""

import argparse
import sys
import time

import numpy as np

from wremnants.postprocessing.scetlib_np import np_tune


def _sigma_reco(model, eff, gnu):
    """σ_reco flat (N_reco,) = R · σ_gen(λ) for one tune, via the model's core."""
    import tensorflow as tf

    gen = model._sigma_gen_at(
        eff, gnu, np_model=eff["np_model"], np_model_nu=gnu["np_model_nu"]
    )
    return np.asarray(
        tf.linalg.matvec(model.R, tf.reshape(gen, [-1])).numpy(), dtype=np.float64
    )


def _project(flat, reco_shape, keep_index):
    """Sum a flat reco vector onto one reco axis, keeping ``keep_index``."""
    a = np.asarray(flat, dtype=np.float64).reshape(reco_shape)
    other = tuple(i for i in range(a.ndim) if i != keep_index)
    return a.sum(axis=other)


def make_reco_data_plot(
    edges,
    data,
    data_err,
    preds,
    axis_name,
    out_path,
    tune_boxes=None,
    norm_note="",
    args=None,
):
    """Data (points) vs one or more σ_reco predictions (steps), + a pred/data panel.

    ``preds`` is a list of ``(label, values, color)`` on the SAME ``edges`` as
    ``data``. All curves are drawn as DIFFERENTIAL dσ/dx (÷ bin width) so the
    variable ptll binning reads correctly; the ratio panel is width-independent.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    edges = np.asarray(edges, dtype=np.float64)
    widths = np.diff(edges)
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, (ax, axr) = plt.subplots(
        2,
        1,
        sharex=True,
        figsize=(7.5, 6.8),
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.06},
    )

    lab = "p$_T^{\\ell\\ell}$ (ptll) [GeV]" if axis_name == "ptll" else axis_name

    # ---- top: density (a.u.) ----
    for label, vals, color in preds:
        ax.stairs(vals / widths, edges, color=color, lw=1.7, label=label)
    ax.errorbar(
        centers,
        data / widths,
        yerr=data_err / widths,
        fmt="o",
        ms=3.5,
        color="k",
        lw=1.0,
        capsize=0,
        label="data",
        zorder=5,
    )
    ax.set_ylabel(r"d$\sigma_{\mathrm{reco}}$/d(" + axis_name + ")  [a.u.]")
    ax.margins(x=0)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper right", fontsize=9)
    if norm_note:
        ax.text(
            0.5,
            1.01,
            norm_note,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=8,
            color="0.35",
        )
    if tune_boxes:
        # num on the right (mistyrose/C3), den on the left (aliceblue/C0), matching
        # sigma_gen_at_lambda's two-tune annotation.
        placements = [
            (0.975, "right", "mistyrose", "C3"),
            (0.025, "left", "aliceblue", "C0"),
        ]
        for (prefix, text), (x, ha, fc, ec) in zip(tune_boxes, placements):
            ax.text(
                x,
                0.60,
                prefix + text,
                transform=ax.transAxes,
                ha=ha,
                va="top",
                fontsize=7,
                family="monospace",
                bbox=dict(boxstyle="round", fc=fc, ec=ec, alpha=0.85),
            )

    # ---- bottom: pred / data ----
    rel = np.divide(data_err, data, out=np.zeros_like(data), where=data > 0)
    axr.fill_between(
        edges,
        np.append(1 - rel, (1 - rel)[-1]),
        np.append(1 + rel, (1 + rel)[-1]),
        step="post",
        color="0.8",
        alpha=0.7,
        lw=0,
        label="data stat",
    )
    for label, vals, color in preds:
        ratio = np.divide(vals, data, out=np.ones_like(vals), where=data > 0)
        axr.stairs(ratio, edges, color=color, lw=1.5)
    axr.axhline(1.0, color="0.5", lw=0.8, ls="--")
    axr.set_ylabel("pred / data")
    axr.set_xlabel(lab)
    axr.margins(x=0)

    # Zoom the ratio to the spread of the predictions (keep 1.0 in frame).
    allr = np.concatenate(
        [np.divide(v, data, out=np.ones_like(v), where=data > 0) for _, v, _ in preds]
    )
    rlo, rhi = float(np.min(allr)), float(np.max(allr))
    pad = max((rhi - rlo) * 0.2, 0.01)
    axr.set_ylim(min(rlo, 1.0) - pad, max(rhi, 1.0) + pad)

    from wremnants.postprocessing.scetlib_np import plot_output

    outdir, basename = plot_output.split_outpath(out_path)
    plot_output.save_plot(outdir, basename, fig=fig, args=args, dpi=150)
    plt.close(fig)
    print(
        f"[plot] wrote {outdir}/{basename}.png(.pdf) + {basename}.log  "
        f"(axis={axis_name}; pred/data spread [{rlo:.4f}, {rhi:.4f}])"
    )


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--datacard",
        required=True,
        help="rabbit input tensor hdf5: supplies R, N_gen, λ_central base, the "
        "reco axes, AND the observed data (hdata_obs)",
    )
    # λ tunes via the shared ingestor: --num-* and --den-*. Same surface, same
    # precedence and same per-side form overrides as sigma_gen_at_lambda.
    np_tune.add_tune_args(p, "num")
    np_tune.add_tune_args(p, "den")
    p.add_argument(
        "--meta-from",
        default=None,
        help="hdf5 for the base λ tune (default: --datacard metadata)",
    )
    p.add_argument(
        "--btgrid", default=None, help="SCETlib bT-grid dir (default: model default)"
    )
    p.add_argument("--q-lo", type=float, default=60.0)
    p.add_argument("--q-hi", type=float, default=120.0)
    p.add_argument(
        "--plot-axis",
        default="ptll",
        help="reco axis to project onto (default ptll; the other reco axes are summed)",
    )
    p.add_argument(
        "--data-norm",
        default="shape",
        choices=["shape", "none"],
        help="'shape' (default): scale each prediction so Σpred=Σdata (the real fit "
        "floats the signal norm, so only shape is physical). 'none': raw σ_reco.",
    )
    p.add_argument("--plot", default=None, help="output path (.png/.pdf)")
    p.add_argument(
        "--no-central",
        action="store_true",
        help="do not overlay the λ_central (datacard) σ_reco reference curve",
    )
    args = p.parse_args(argv)

    from rabbit.debugdata import FitDebugData
    from rabbit.inputdata import FitInputData
    from wremnants.postprocessing.scetlib_np.param_model import SCETlibNPParamModel
    from wremnants.postprocessing.scetlib_np.sigma_gen_at_lambda import _lambda_box_text

    # ---- construct the param model from the datacard (R + N_gen + λ_central) ----
    print(
        "\n[core] loading the datacard + constructing SCETlibNPParamModel …", flush=True
    )
    t0 = time.time()
    indata = FitInputData(args.datacard)
    model = SCETlibNPParamModel(
        indata, btgrid_dir=args.btgrid, Q_lo=args.q_lo, Q_hi=args.q_hi
    )
    if getattr(model, "gen_level", False):
        raise SystemExit(
            "--datacard is a gen-level (σUL) card with no response matrix; there "
            "is no reco fold to make. Use sigma_gen_at_lambda.py instead."
        )
    print(f"  constructed in {time.time()-t0:.1f}s; R shape {tuple(model.R.shape)}")

    # ---- reco axes + observed data, from the SAME datacard ----
    dbg = FitDebugData(indata)
    channels = list(dbg.data_obs_hists.keys())
    if len(channels) != 1:
        raise SystemExit(
            f"expected a single fit channel, found {channels}; --plot-axis can't pick"
        )
    hdat = dbg.data_obs_hists[channels[0]]
    reco_names = [a.name for a in hdat.axes]
    reco_shape = [a.size for a in hdat.axes]
    if args.plot_axis not in reco_names:
        raise SystemExit(
            f"--plot-axis {args.plot_axis!r} not in reco axes {reco_names}"
        )
    ai = reco_names.index(args.plot_axis)
    edges = np.asarray(hdat.axes[ai].edges, dtype=np.float64)

    data = _project(hdat.values(), reco_shape, ai)
    data_var = _project(hdat.variances(), reco_shape, ai)
    data_err = np.sqrt(np.maximum(data_var, 0.0))
    print(
        f"[data] channel {channels[0]}: reco axes {list(zip(reco_names, reco_shape))}; "
        f"Σ data = {data.sum():.6g} (projected onto {args.plot_axis})"
    )

    # ---- base tune + the two tunes' σ_reco, via the shared ingestor ----
    if not args.meta_from:
        args.meta_from = args.datacard
    do_ratio = np_tune.side_requested(args, "den")
    try:
        base_tune = np_tune.resolve_base_tune(args)
        num_tune, _ = np_tune.resolve_side_tune(args, "num", base_tune)
        # forms inherit den-from-num (λ still fall back to the base), so a
        # per-side form override cannot silently move the other form
        den_tune = (
            np_tune.resolve_side_tune(args, "den", base_tune, inherit=num_tune)[0]
            if do_ratio
            else None
        )
    except ValueError as e:
        p.error(str(e))

    import os

    sides = [("num", num_tune, args.num_fitresult, args.num_label, "C3")]
    if do_ratio:
        sides.append(("den", den_tune, args.den_fitresult, args.den_label, "C0"))

    projections = []  # (label, projected-σ_reco, color)
    tune_boxes = []
    for role, tune, fr, label, color in sides:
        proj = _project(_sigma_reco(model, tune.eff, tune.gnu), reco_shape, ai)
        # default label = the fit's run directory basename (mirrors sigma_gen_at_lambda)
        auto = label or (
            os.path.basename(os.path.dirname(os.path.abspath(fr)))
            if fr
            else tune.describe()
        )
        projections.append([auto, proj, color])
        tune_boxes.append(
            (f"{role}  " if do_ratio else "", _lambda_box_text(tune.eff, tune.gnu))
        )
        print(f"[{role}] {tune.describe()}\n  Σ σ_reco (proj) = {proj.sum():.6g}")

    central_proj = None
    if not args.no_central:
        central_proj = _project(
            np.asarray(model.sigma_reco_central.numpy(), dtype=np.float64),
            reco_shape,
            ai,
        )

    # ---- normalization ----
    norm_note = ""
    if args.data_norm == "shape":
        dtot = data.sum()
        for entry in projections:
            ptot = entry[1].sum()
            if ptot != 0:
                entry[1] = entry[1] * (dtot / ptot)
        if central_proj is not None:
            ctot = central_proj.sum()
            if ctot != 0:
                central_proj = central_proj * (dtot / ctot)
        norm_note = (
            "signal σ_reco (param model), each shape-normalized to data "
            "(Σpred=Σdata); no background / postfit nuisances"
        )
    else:
        norm_note = "raw σ_reco (a.u.); NOT normalized to the data yield"

    preds = []
    if central_proj is not None:
        preds.append(("λ_central", central_proj, "0.4"))
    preds.extend((lbl, vals, col) for lbl, vals, col in projections)

    if not args.plot:
        print("\n[done] (no --plot given; nothing written)")
        return 0

    make_reco_data_plot(
        edges,
        data,
        data_err,
        preds,
        args.plot_axis,
        args.plot,
        tune_boxes=tune_boxes if do_ratio else None,
        norm_note=norm_note,
        args=args,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
