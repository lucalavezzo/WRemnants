"""Test 6 v2 — Random-truth 4-POI recovery test for SCETlibNPParamModel.

Bypasses rabbit entirely. For each of N random truth points drawn from the
analysis priors (rejection-sampled so each truth satisfies hard bounds +
tanh2 stability):

    1. Compute pseudo-data = model(x_truth) (all 8 params at truth).
    2. Start ALL 8 params at their prior centers (= runcard).
    3. Float only the 4 POIs (lambda2, lambda4, lambda2_nu, delta_lambda2).
       The 4 non-POI params (lambda_inf, lambda6, lambda_inf_nu, lambda4_nu)
       stay frozen at their prior centers throughout — they do NOT track
       truth, so the fit has to absorb that mismatch through the 4 POIs.
    4. Minimize ||compute(λ)[:, signal] − data[:, signal]||² with BFGS
       (autograd from TF, projected to the 4 POI dimensions).
    5. Record the 4-D residual (postfit_POI − truth_POI) per POI.

This mirrors what rabbit actually does in production: 4 POIs float, the
other 4 sit at runcard regardless of nature. The interesting question is
whether the 4 POIs can mimic shifts in the frozen 4 (degeneracy) and, if
so, by how much they get pulled away from their own truth values.

Environment knobs:
    TEST6_N      number of random truths (default 50)
    TEST6_SEED   RNG seed (default 42)
    MAX_FITS     hard cap on number of fits (debug)

Run inside singularity::

    python3 -m wremnants.postprocessing.scetlib_np.validation.truth_start_grid
"""

import csv
import os
import sys
import time
from pathlib import Path

import numpy as np

from wremnants.postprocessing.scetlib_np.sigma_gen import _default_btgrid_dir
import tensorflow as tf
from scipy.optimize import minimize

from rabbit.inputdata import FitInputData

from wremnants.postprocessing.scetlib_np.param_model import (
    ALL_PARAMS,
    SCETlibNPParamModel,
)

FIT_HDF5 = ("/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/"
            "260526_fit_NewVarsCT18ZLambda6/"
            "ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile"
            "_NewVarsCT18ZLambda6_v4_franks/ZMassDilepton.hdf5")

OUTPUT_DIR = Path(os.environ.get(
    "TEST6_OUT",
    "/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/"
    "260530_scetlib_np_truth_random_4poi",
))

# (mu, sigma) for each param's Gaussian prior. Mirror of the test 5
# generator's PRIORS so the truth distribution matches what rabbit sees.
PRIORS = {
    "lambda2":        (0.4, 0.5),
    "lambda4":        (0.4, 0.5),
    "lambda6":        (0.0, 0.05),
    "delta_lambda2":  (0.0, 0.20),
    "lambda_inf":     (1.0, 0.5),
    "lambda2_nu":     (0.15, 0.10),
    "lambda4_nu":     (0.0, 0.05),
    "lambda_inf_nu":  (2.0, 1.0),
}

# The 4 POIs that rabbit actually floats in test 5 (see
# scripts/rabbit/scetlib_np/injection_tests.sh COMMON_FREEZE).
POI_NAMES = ("lambda2", "lambda4", "lambda2_nu", "delta_lambda2")


def _hard_bounds_ok(p):
    return (p["lambda2"]        > 0.05
            and p["lambda_inf"]     > 0.05
            and p["lambda_inf_nu"]  > 0.05)


def _tanh2_stability_ok(p):
    return (p["lambda4"]
            + (1.0 / 3.0) * (p["lambda2"] ** 3) / (p["lambda_inf"] ** 2)) > 0.0


def _accept(p):
    return _hard_bounds_ok(p) and _tanh2_stability_ok(p)


def sample_truth(rng):
    """Rejection-sample a truth dict from PRIORS subject to stability."""
    rejected = 0
    while True:
        p = {k: float(rng.normal(mu, s)) for k, (mu, s) in PRIORS.items()}
        if _accept(p):
            return p, rejected
        rejected += 1
        if rejected > 10000:
            raise RuntimeError("rejection rate explosion — check priors")


def _signal_indices(indata, signal_proc="Zmumu"):
    procs = [p.decode() if isinstance(p, bytes) else str(p) for p in indata.procs]
    return procs.index(signal_proc), indata.nproc


def _build_x(values_dict, param_order):
    return np.array([values_dict[k] for k in param_order], dtype=np.float64)


def _make_loss_and_grad_partial(model, signal_proc_idx, data_ratio_tf,
                                x_template, free_idx):
    """Loss closure where only `free_idx` of x floats; others held at
    x_template. Returns (loss, grad_w.r.t._free_subset).
    """
    free_idx = np.asarray(free_idx, dtype=np.int64)

    def loss_and_grad(x_free):
        x_full = x_template.copy()
        x_full[free_idx] = np.asarray(x_free, dtype=np.float64).ravel()
        xv = tf.Variable(x_full, dtype=model.indata.dtype)
        with tf.GradientTape() as tape:
            rnorm = model.compute(xv)
            ratio = rnorm[:, signal_proc_idx]
            loss = tf.reduce_mean((ratio - data_ratio_tf) ** 2)
        grad_full = tape.gradient(loss, xv).numpy().astype(np.float64)
        return float(loss.numpy()), grad_full[free_idx]

    return loss_and_grad


def run_fit(model, signal_proc_idx, x_truth, x_start, free_idx):
    """Pseudo-data = model(x_truth). Float `free_idx`, freeze the rest at
    x_start. Returns dict with postfit (full 8-vector) and diagnostics.
    """
    x_truth_tf = tf.constant(x_truth, dtype=model.indata.dtype)
    rnorm_truth = model.compute(x_truth_tf).numpy()
    data_ratio_tf = tf.constant(
        rnorm_truth[:, signal_proc_idx], dtype=model.indata.dtype
    )

    x_template = x_start.astype(np.float64).copy()
    loss_and_grad = _make_loss_and_grad_partial(
        model, signal_proc_idx, data_ratio_tf, x_template, free_idx,
    )

    x0 = x_start[free_idx].copy()
    f0, g0 = loss_and_grad(x0)
    g0_max = float(np.max(np.abs(g0)))
    g0_norm = float(np.linalg.norm(g0))

    t0 = time.time()
    result = minimize(
        loss_and_grad, x0, jac=True, method="BFGS",
        options={"gtol": 1e-12, "maxiter": 500, "disp": False},
    )
    elapsed = time.time() - t0

    postfit = x_template.copy()
    postfit[free_idx] = result.x
    resid_full = np.abs(postfit - x_truth)
    resid_poi = np.abs(postfit[free_idx] - x_truth[free_idx])
    return {
        "postfit": postfit,
        "n_iters": int(result.nit),
        "n_fev": int(result.nfev),
        "success": bool(result.success),
        "max_resid_full": float(resid_full.max()),
        "max_resid_poi":  float(resid_poi.max()),
        "final_loss": float(result.fun),
        "start_loss": float(f0),
        "start_grad_max": g0_max,
        "start_grad_norm": g0_norm,
        "message": (result.message.decode()
                    if isinstance(result.message, bytes) else str(result.message)),
        "elapsed_s": float(elapsed),
    }


def write_csv(rows, param_order, csv_path):
    if not rows:
        return
    fieldnames = (
        ["truth_idx", "n_iters", "n_fev", "success",
         "max_resid_full", "max_resid_poi", "final_loss", "start_loss",
         "start_grad_max", "start_grad_norm", "elapsed_s", "message"]
        + [f"truth_{n}"  for n in param_order]
        + [f"start_{n}"  for n in param_order]
        + [f"postfit_{n}" for n in param_order]
    )
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def make_plots(rows, param_order, output_dir, run_meta=None):
    """For each POI: scatter of (truth, postfit) and histogram of residuals
    in units of σ_prior. Also a per-POI bias/pull summary plot.

    ``run_meta`` (optional) is recorded in each plot's ``.log`` alongside the
    command line (this script is env-var driven, so the config lives there).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from wremnants.postprocessing.scetlib_np import plot_output

    plots = []
    n_ok = [r for r in rows if np.isfinite(r["final_loss"])]
    if not n_ok:
        return plots

    for poi in POI_NAMES:
        truths   = np.array([r[f"truth_{poi}"]   for r in n_ok])
        postfits = np.array([r[f"postfit_{poi}"] for r in n_ok])
        resids   = postfits - truths
        sigma = PRIORS[poi][1]
        pulls = resids / sigma

        fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
        # (1) scatter truth vs postfit
        ax = axes[0]
        lo = min(truths.min(), postfits.min()) - 0.05 * sigma
        hi = max(truths.max(), postfits.max()) + 0.05 * sigma
        ax.plot([lo, hi], [lo, hi], "k--", lw=1, alpha=0.5, label="postfit = truth")
        ax.scatter(truths, postfits, s=18, alpha=0.7)
        ax.set_xlabel(f"truth {poi}")
        ax.set_ylabel(f"postfit {poi}")
        ax.set_title(f"{poi}: truth vs postfit")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

        # (2) residual vs truth
        ax = axes[1]
        ax.axhline(0, color="k", lw=1, alpha=0.5)
        ax.scatter(truths, resids, s=18, alpha=0.7)
        ax.set_xlabel(f"truth {poi}")
        ax.set_ylabel(f"postfit − truth ({poi})")
        ax.set_title(f"{poi}: residual vs truth")
        ax.grid(alpha=0.3)

        # (3) pull histogram
        ax = axes[2]
        ax.hist(pulls, bins=20, alpha=0.75, edgecolor="k")
        ax.axvline(0, color="k", lw=1)
        ax.axvline(pulls.mean(), color="C3", ls="--",
                   label=f"mean = {pulls.mean():+.2f}")
        ax.set_xlabel(f"(postfit − truth) / σ_prior  [σ_prior = {sigma}]")
        ax.set_ylabel("count")
        ax.set_title(f"{poi}: pull")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

        fig.suptitle(f"{poi}: N={len(n_ok)} random truths, 4-POI BFGS")
        fig.tight_layout()
        plot_output.save_plot(
            str(output_dir), f"recover_{poi}", fig=fig, meta_info=run_meta, dpi=110
        )
        plt.close(fig)
        plots.append((poi, f"recover_{poi}.pdf", f"recover_{poi}.png"))

    # summary box: bias (mean residual / sigma) and RMS per POI
    fig, ax = plt.subplots(figsize=(7, 4))
    bias = []
    rms = []
    for poi in POI_NAMES:
        truths   = np.array([r[f"truth_{poi}"]   for r in n_ok])
        postfits = np.array([r[f"postfit_{poi}"] for r in n_ok])
        resids = postfits - truths
        sigma = PRIORS[poi][1]
        bias.append(resids.mean() / sigma)
        rms.append(np.sqrt(np.mean(resids**2)) / sigma)
    x = np.arange(len(POI_NAMES))
    ax.bar(x - 0.18, bias, width=0.36, label="bias / σ_prior")
    ax.bar(x + 0.18, rms,  width=0.36, label="RMS / σ_prior")
    ax.axhline(0, color="k", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(POI_NAMES, rotation=20)
    ax.set_ylabel("residual (units of σ_prior)")
    ax.set_title(f"4-POI fit summary, N={len(n_ok)}")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    plot_output.save_plot(
        str(output_dir), "summary_bias_rms", fig=fig, meta_info=run_meta, dpi=110
    )
    plt.close(fig)
    plots.append(("summary", "summary_bias_rms.pdf", "summary_bias_rms.png"))

    return plots


def write_readme(rows, plots, output_dir):
    n = len(rows)
    n_ok = sum(1 for r in rows if r["success"] and np.isfinite(r["final_loss"]))
    lines = [
        "# Test 6 v2 — random-truth 4-POI recovery",
        "",
        f"- N truths: {n}",
        f"- successful fits: {n_ok}",
        f"- POIs floated: {', '.join(POI_NAMES)}",
        f"- frozen at prior centers: " + ", ".join(
            p for p in ALL_PARAMS if p not in POI_NAMES
        ),
        "",
        "## Plots",
        "",
    ]
    for tag, pdf, png in plots:
        lines.append(f"- **{tag}**: [{pdf}]({pdf}) · ![{tag}]({png})")
    (output_dir / "README.md").write_text("\n".join(lines) + "\n")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    seed = int(os.environ.get("TEST6_SEED", "42"))
    N = int(os.environ.get("TEST6_N", "50"))
    max_fits_env = os.environ.get("MAX_FITS", "").strip()
    if max_fits_env:
        N = min(N, int(max_fits_env))

    rng = np.random.default_rng(seed)

    print("Loading FitInputData …")
    t0 = time.time()
    indata = FitInputData(FIT_HDF5)
    print(f"  loaded in {time.time()-t0:.1f}s; nproc={indata.nproc}")

    print("Constructing SCETlibNPParamModel …")
    t0 = time.time()
    model = SCETlibNPParamModel(
        indata,
        btgrid_dir=_default_btgrid_dir(),
        signal_proc="Zmumu",
    )
    print(f"  constructed in {time.time()-t0:.1f}s")
    signal_proc_idx, _ = _signal_indices(indata)
    param_order = list(model._param_order)
    print(f"  param order: {param_order}")

    free_idx = np.array(
        [param_order.index(n) for n in POI_NAMES], dtype=np.int64
    )
    print(f"  POI indices: {dict(zip(POI_NAMES, free_idx.tolist()))}")

    # x_start: all 8 at prior centers (= runcard / rabbit starting point).
    prior_center = {k: PRIORS[k][0] for k in param_order}
    x_start = _build_x(prior_center, param_order)
    print(f"  x_start (prior centers): {dict(zip(param_order, x_start.tolist()))}")

    print(f"\nN truths to draw: {N}  (seed={seed})\n")

    rows = []
    total_rejected = 0
    for i in range(N):
        truth_dict, rej = sample_truth(rng)
        total_rejected += rej
        x_truth = _build_x(truth_dict, param_order)
        print(
            f"[{i+1}/{N}] truth: "
            + ", ".join(f"{p}={truth_dict[p]:+.3f}" for p in POI_NAMES)
            + f"  (rej={rej})",
            flush=True,
        )

        try:
            res = run_fit(model, signal_proc_idx, x_truth, x_start, free_idx)
        except Exception as e:
            print(f"  FIT RAISED: {e!r}", flush=True)
            res = {
                "postfit": np.full_like(x_truth, np.nan),
                "n_iters": 0, "n_fev": 0, "success": False,
                "max_resid_full": float("nan"),
                "max_resid_poi":  float("nan"),
                "final_loss": float("nan"),
                "start_loss": float("nan"),
                "start_grad_max": float("nan"),
                "start_grad_norm": float("nan"),
                "message": f"EXC: {e!r}",
                "elapsed_s": 0.0,
            }

        resid_poi_per = res["postfit"][free_idx] - x_truth[free_idx]
        print(
            f"  iters={res['n_iters']:3d} nfev={res['n_fev']:3d} "
            f"success={res['success']} "
            f"max|Δpoi|={res['max_resid_poi']:.3g} loss={res['final_loss']:.3g} "
            f"(start loss={res['start_loss']:.3g}) "
            f"t={res['elapsed_s']:.1f}s  "
            + "  ".join(
                f"Δ{n}={resid_poi_per[k]:+.3f}"
                for k, n in enumerate(POI_NAMES)
            ),
            flush=True,
        )

        row = {
            "truth_idx": i,
            "n_iters": res["n_iters"],
            "n_fev": res["n_fev"],
            "success": int(res["success"]),
            "max_resid_full": res["max_resid_full"],
            "max_resid_poi":  res["max_resid_poi"],
            "final_loss": res["final_loss"],
            "start_loss": res["start_loss"],
            "start_grad_max": res["start_grad_max"],
            "start_grad_norm": res["start_grad_norm"],
            "elapsed_s": res["elapsed_s"],
            "message": res["message"],
        }
        for j, name in enumerate(param_order):
            row[f"truth_{name}"]   = float(x_truth[j])
            row[f"start_{name}"]   = float(x_start[j])
            row[f"postfit_{name}"] = float(res["postfit"][j])
        rows.append(row)

    print(f"\nTotal samples rejected by stability: {total_rejected}")

    csv_path = OUTPUT_DIR / "summary.csv"
    write_csv(rows, param_order, csv_path)
    print(f"Wrote {csv_path}")

    run_meta = {
        "run_config": {
            "seed": seed,
            "N": N,
            "POIs": POI_NAMES,
            "fit_hdf5": str(FIT_HDF5),
            "btgrid_dir": str(_default_btgrid_dir()),
        }
    }
    plots = make_plots(rows, param_order, OUTPUT_DIR, run_meta=run_meta)
    for tag, pdf, png in plots:
        print(f"Wrote {OUTPUT_DIR / pdf}")
        print(f"Wrote {OUTPUT_DIR / png}")

    write_readme(rows, plots, OUTPUT_DIR)
    print(f"Wrote {OUTPUT_DIR / 'README.md'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
