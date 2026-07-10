"""Unified entry point for the SCETlib-NP param-model agreement checks (option A).

One CLI, ``--reference {card,histmaker}``, over the two agreement references that
share the datacard-built ``SCETlibNPParamModel``:

  * ``card``      : the model's λ_central σ_reco / σ_gen vs the datacard itself
                    (``indata.norm[signal]`` and the ``N_gen`` auxiliary) — no
                    external inputs.
  * ``histmaker`` : the same vs an external histmaker ``nominal`` (+ gen MC and
                    the Corr[var] λ-variations).

The third reference — the official ``TheoryCorrection`` at an arbitrary λ tune,
gen-only, cardless — is its own CLI (``sigma_gen_at_lambda``): genuinely distinct
(no datacard/response, arbitrary λ, theory truth).

The shared reference loaders/aligners live in ``validation.agreement``; the
fit-side guard + pathology detectors and the card comparison
(``run_card_diagnostics``) live in ``param_model_diagnostics`` (numpy Layer 0).

Run inside the wmass singularity, e.g.:

    python3 -m wremnants.postprocessing.scetlib_np.validate_agreement \\
        --reference card --datacard <hdf5> [--outdir <plotdir>]

    python3 -m wremnants.postprocessing.scetlib_np.validate_agreement \\
        --reference histmaker --datacard <hdf5> --histmaker <hdf5> \\
        [--plot-out <path>] [--gen-histmaker <hdf5>] [--variation lambda21.0 ...]
"""

import argparse
import sys
import time

import numpy as np

from rabbit.inputdata import FitInputData
from wremnants.postprocessing.scetlib_np.param_model import SCETlibNPParamModel
from wremnants.postprocessing.scetlib_np.param_model_diagnostics import (
    run_card_diagnostics,
)
from wremnants.postprocessing.scetlib_np.sigma_gen import _default_btgrid_dir
from wremnants.postprocessing.scetlib_np.validation.agreement import (
    GEN_HIST,
    GEN_SAMPLE,
    NOMINAL_HIST,
    SIGNAL_PROC,
    SIGNAL_SAMPLE,
    VARIATION_HIST,
    align_gen,
    align_nominal,
    load_gen_hist,
    load_nominal,
    validate_variation,
)
from wremnants.postprocessing.scetlib_np.validation_plots import (
    plot_ptll_ratio,
    summarize,
    tf_to_hist,
)


def run_card(args):
    """``--reference card``: model σ_reco/σ_gen(λ_c) vs the datacard itself
    (``indata.norm[signal]`` + the ``N_gen`` auxiliary). No external inputs."""
    print("Loading FitInputData …", flush=True)
    t0 = time.time()
    indata = FitInputData(args.datacard)
    print(f"  loaded in {time.time()-t0:.1f}s; nproc={indata.nproc}", flush=True)

    print(
        "Constructing SCETlibNPParamModel (runs the bt integral at λ_central) …",
        flush=True,
    )
    t0 = time.time()
    kw = dict(
        signal_proc=args.signal_proc, check_agreement=False
    )  # report below, not the guard
    if args.btgrid:
        kw["btgrid_dir"] = args.btgrid
    model = SCETlibNPParamModel(indata, **kw)
    print(f"  constructed in {time.time()-t0:.1f}s", flush=True)

    run_card_diagnostics(
        model,
        indata,
        outdir=(args.outdir or None),
        do_plots=bool(args.outdir),
        args=args,
    )
    return 0


def run_histmaker(args):
    """``--reference histmaker``: model σ_reco/σ_gen(λ_c) + λ-variations vs an
    external histmaker ``nominal`` / gen MC / Corr[var]."""
    import os

    btgrid = args.btgrid or _default_btgrid_dir()

    print("=" * 70)
    print("SCETlibNPParamModel  vs  histmaker nominal  — integral validation")
    print("=" * 70)
    print(f"  datacard  : {args.datacard}")
    print(f"  histmaker : {args.histmaker}")
    print(f"  btgrid    : {btgrid}")

    print("\nLoading FitInputData …")
    t0 = time.time()
    indata = FitInputData(args.datacard)
    print(f"  loaded in {time.time()-t0:.1f}s; nproc={indata.nproc}")

    print("\nConstructing SCETlibNPParamModel (runs the bT integral at λ_central) …")
    t0 = time.time()
    model = SCETlibNPParamModel(
        indata,
        btgrid_dir=btgrid,
        signal_proc=args.signal_proc,
    )
    print(f"  constructed in {time.time()-t0:.1f}s")
    print(f"  λ_central eff : {model.eff_central}")
    print(f"  λ_central gnu : {model.gnu_central}")
    print(f"  reco axes     : {[n for n, _ in model._reco_axes_meta]}")
    print(f"  reco shape    : {model.reco_shape}")

    # Everything stays at the hist level: the nominal is reordered/cropped onto
    # the model's reco binning (a Hist), and the model's tf output is wrapped
    # into a matching Hist via tf_to_hist.
    print(f"\nLoading histmaker {args.hist!r} for {args.sample!r} …")
    h_nom = align_nominal(
        load_nominal(args.histmaker, args.sample, args.hist), model._reco_axes_meta
    )
    h_model = tf_to_hist(model.sigma_reco_central, model._reco_axes_meta)

    summarize(
        h_model.values(flow=False), h_nom.values(flow=False), model._reco_axes_meta
    )

    # Optional third curve: resum-only σ_reco, derived from the matched model
    # by subtracting the exposed σ_ns at gen level and re-folding through R
    # (σ_gen = σ_resum + σ_ns, both cached on the model — no second model).
    extra_models = None
    model_label = r"ParamModel $\sigma_{reco}$"
    if args.overlay_resum:
        import tensorflow as tf

        print("\nDeriving resum-only σ_reco (matched − σ_ns, folded through R) …")
        gen_resum = model.sigma_gen_central - model.sigma_ns
        reco_resum = tf.linalg.matvec(model.R, tf.reshape(gen_resum, [-1]))
        h_resum = tf_to_hist(reco_resum, model._reco_axes_meta)
        scale_r = float(h_nom.values().sum() / h_resum.values().sum())
        extra_models = [(h_resum, "ParamModel resum-only", "blue", scale_r)]
        model_label = r"ParamModel matched (resum$\oplus$FO)"

    if args.plot_out:
        # σ_reco_central is a theory σ (the response is normalized) while nominal
        # is a weighted event yield — density-normalize both for a pure shape
        # comparison, so there's no ad-hoc Σ/Σ scale factor (see plot_ptll_ratio).
        # Project onto each reco shape axis: ptll (-> args.plot_out) and yll
        # (-> <base>_yll<ext>); the helper sums out the remaining reco axes.
        base, ext = os.path.splitext(args.plot_out)
        for ax_name in ("ptll", "yll"):
            out = args.plot_out if ax_name == "ptll" else f"{base}_{ax_name}{ext}"
            plot_ptll_ratio(
                h_model,
                h_nom,
                axis=ax_name,
                out_path=out,
                density=True,
                model_label=model_label,
                title="",
                autorrange=0.2,
                ratio_legend=False,
                no_sci=True,
                extra_models=extra_models,
            )

    # ---- Gen-level cross-check: σ_gen(λ_central) [the bT integral, NO response]
    # vs a gen-level histmaker at the same NP tune. Tests the integral alone.
    if args.gen_histmaker:
        print(f"\nGen-level check: σ_gen vs {args.gen_hist!r} from {args.gen_sample!r}")
        gen_arr = align_gen(
            load_gen_hist(args.gen_histmaker, args.gen_sample, args.gen_hist),
            model._gen_axes_meta,
        )
        h_gen_mc = tf_to_hist(gen_arr, model._gen_axes_meta)
        h_gen_model = tf_to_hist(model.sigma_gen_central, model._gen_axes_meta)
        print("  [gen level] σ_gen(integral) vs gen MC, on the model gen grid:")
        summarize(
            h_gen_model.values(flow=False),
            h_gen_mc.values(flow=False),
            model._gen_axes_meta,
        )
        if args.plot_out:
            base, ext = os.path.splitext(args.plot_out)
            gen_out = f"{base}_gen{ext}"
            # ptVGen-projected unit-normalized (density) shape ratio model/gen MC,
            # used only to size the ratio panel — the plot itself density-normalizes.
            edges = np.asarray(model._gen_axes_meta[0][1], dtype=np.float64)
            centers = 0.5 * (edges[:-1] + edges[1:])
            mp = h_gen_model.project("ptVGen").values(flow=False)
            npj = h_gen_mc.project("ptVGen").values(flow=False)
            good = npj > 0
            rp = np.where(
                good,
                (mp / mp.sum()) / np.where(good, npj / npj.sum(), 1.0),
                np.nan,
            )
            # Focus the view on the fit region (ptll fit stops at 44); the wide
            # [44,100] gen-grid tail bin would otherwise squash the spectrum and
            # loosen the ratio range. Set the ratio range from the visible bins.
            xhi = 44.0
            vis = good & (centers < xhi)
            half = float(np.nanmax(np.abs(rp[vis] - 1.0)))
            rr = max(half * 1.3, 0.004)
            plot_ptll_ratio(
                h_gen_model,
                h_gen_mc,
                axis="ptVGen",
                out_path=gen_out,
                density=True,
                model_label=r"ParamModel $\sigma_{\mathrm{gen}}$",
                ref_label="histmaker gen MC",
                rlabel="model / gen MC",
                title="",
                rrange=(1.0 - rr, 1.0 + rr),
                xlim=(0.0, xhi),
                ratio_legend=False,
                no_sci=True,
            )

    # ---- Off-central λ-response: model rnorm vs histmaker Corr[var]/Corr[pdf0] ----
    for var_label in args.variation or []:
        out_v = None
        if args.plot_out:
            base, ext = os.path.splitext(args.plot_out)
            out_v = f"{base}_var_{var_label}{ext}"
        validate_variation(model, args, var_label, out_path=out_v)

    print("\nDone.")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--reference",
        required=True,
        choices=["card", "histmaker"],
        help="agreement reference: 'card' (the datacard itself) or 'histmaker' "
        "(an external histmaker nominal + gen MC + Corr variations)",
    )
    # ---- shared ----
    p.add_argument(
        "--datacard",
        required=True,
        help="fit-input hdf5 (FitInputData + the scetlib_np auxiliary + λ_central metadata)",
    )
    p.add_argument(
        "--btgrid",
        default=None,
        help="SCETlib bT-grid dir (default: the model's data-area copy)",
    )
    p.add_argument(
        "--signal-proc", default=SIGNAL_PROC, help="indata signal process name"
    )
    # ---- card only ----
    p.add_argument(
        "--outdir",
        default=None,
        help="[card] plot output dir ('' / unset to skip plotting)",
    )
    # ---- histmaker only ----
    p.add_argument(
        "--histmaker",
        default=None,
        help="[histmaker] histmaker hdf5 holding the 'nominal' hist",
    )
    p.add_argument(
        "--sample",
        default=SIGNAL_SAMPLE,
        help="[histmaker] histmaker sample group for the signal",
    )
    p.add_argument(
        "--hist",
        default=NOMINAL_HIST,
        help="[histmaker] histogram name to compare against",
    )
    p.add_argument(
        "--plot-out",
        default="",
        help="[histmaker] path for the ptll-projection ratio plot ('' to skip)",
    )
    p.add_argument(
        "--gen-histmaker",
        default="",
        help="[histmaker] gen-level histmaker hdf5 (same NP tune) for the gen-level "
        "cross-check (σ_gen vs gen MC, NO response matrix); '' to skip",
    )
    p.add_argument(
        "--gen-hist", default=GEN_HIST, help="[histmaker] gen-level hist name"
    )
    p.add_argument(
        "--gen-sample",
        default=GEN_SAMPLE,
        help="[histmaker] gen histmaker sample group",
    )
    p.add_argument(
        "--variation",
        nargs="*",
        default=[],
        metavar="LABEL",
        help="[histmaker] off-central λ-response check: one or more 'vars'-axis labels "
        "(e.g. lambda21.0 lambda2_nu0.25 lambda41.0). For each, compare the model's "
        "rnorm=σ_reco(λ)/σ_reco(λ_c) to the histmaker's Corr[var]/Corr[pdf0].",
    )
    p.add_argument(
        "--variation-hist",
        default=VARIATION_HIST,
        help="[histmaker] histmaker hist (with a 'vars' axis) holding the NP λ variations",
    )
    p.add_argument(
        "--overlay-resum",
        action="store_true",
        help="[histmaker] overlay the resum-only σ_reco (matched minus the exposed σ_ns, "
        "re-folded through R) as a third curve on the ptll plot.",
    )
    args = p.parse_args(argv)

    if args.reference == "card":
        return run_card(args)
    if not args.histmaker:
        p.error("--reference histmaker requires --histmaker")
    return run_histmaker(args)


if __name__ == "__main__":
    sys.exit(main())
