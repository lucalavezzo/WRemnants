"""Analyze SCETlib NP injection-test fit results.

Walks ``$OUT_BASE`` (the directory produced by
``injection_tests.sh``) and produces:

- Block (1) / (2): one pull histogram per (block, param). x-axis
  ``(λ_post − λ_truth) / σ_prior``, 10 entries each (one per toy).
- Block (3) / (4): one scatter plot per (block, param). x = λ_start,
  y = λ_postfit. Reference lines: horizontal at truth (data-only
  recovery), diagonal y=x (prior-only / no motion). Block (4) gets a
  linear-fit overlay whose slope is the data-vs-prior tug-of-war
  number.

Also writes ``summary.csv`` with one row per (block, fit_or_toy, param)
containing postfit / prefit values, pulls, and per-block aggregates.

Run inside the wmass singularity with setup.sh sourced::

    python3 scripts/rabbit/scetlib_np/injection_analyze.py \\
        /path/to/scetlib_np_injection_<timestamp>/
"""

import argparse
import glob
import os
import re
import sys
import csv

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wums import ioutils

from wremnants.postprocessing.scetlib_np import plot_output


# ---- inputs ----------------------------------------------------------------

# The four free ParamModel POUs we care about, in the order to plot.
PARAMS = ["lambda2", "lambda4", "lambda2_nu", "delta_lambda2"]

# Theorist-recommended prior widths (must match DEFAULT_PRIOR_SIGMAS in
# scetlib_np.param_model).
SIGMA_PRIOR = {
    "lambda2_nu":    0.10,
    "lambda2":       0.50,
    "lambda4":       0.50,
    "delta_lambda2": 0.20,
}

# Truth = runcard λ_central (the data-generating point).
# These are the FranksVals defaults that the histmaker bakes into pseudoData.
TRUTH = {
    "lambda2":       0.40,
    "lambda4":       0.40,
    "lambda2_nu":    0.15,
    "delta_lambda2": 0.00,
}

DEFAULT_OUTDIR = ("/home/submit/lavezzo/public_html/alphaS/"
                  "260530_scetlib_np_injection_tests")


# ---- block (5) constants --------------------------------------------------

# Block (5) lets all 8 ParamModel POUs float, so we need to track all of
# them when reading per-fit postfits and building pull histograms.
PARAMS_B5 = [
    "lambda2", "lambda4", "lambda6", "delta_lambda2",
    "lambda_inf", "lambda2_nu", "lambda4_nu", "lambda_inf_nu",
]

TRUTH_B5 = {
    "lambda2":       0.40,
    "lambda4":       0.40,
    "lambda6":       0.00,
    "delta_lambda2": 0.00,
    "lambda_inf":    1.00,
    "lambda2_nu":    0.15,
    "lambda4_nu":    0.00,
    "lambda_inf_nu": 2.00,
}

# σ widths for normalizing pulls. The 4 priored POUs reuse SIGMA_PRIOR;
# the 4 currently-unpriored POUs get wide defaults matching the test (5)
# Gaussian draws.
SIGMA_REF_B5 = {
    "lambda2":       0.50,
    "lambda4":       0.50,
    "lambda6":       0.05,
    "delta_lambda2": 0.20,
    "lambda_inf":    0.50,
    "lambda2_nu":    0.10,
    "lambda4_nu":    0.05,
    "lambda_inf_nu": 1.00,
}


# ---- HDF5 readers ---------------------------------------------------------

def _values_by_name(parms_hist, names):
    """Pull values out of a 1D hist whose axis bins are parameter names."""
    axis = parms_hist.axes[0]
    labels = [b.decode() if isinstance(b, bytes) else str(b) for b in axis]
    values = np.asarray(parms_hist.values())
    out = {}
    for n in names:
        if n in labels:
            out[n] = float(values[labels.index(n)])
    return out


def _resolve(maybe_proxy):
    """Materialize an H5PickleProxy if needed."""
    return maybe_proxy.get() if hasattr(maybe_proxy, "get") else maybe_proxy


def load_postfits(path, params):
    """Open ``path`` and eagerly extract {postfit, prefit} dicts for every
    ``results_*`` group, materializing H5PickleProxies before the file is
    closed. Returns ``{group_name: {"postfit": {...}, "prefit": {...}}}``.
    """
    out = {}
    with h5py.File(path, "r") as f:
        for k in f.keys():
            if not k.startswith("results_"):
                continue
            res = ioutils.pickle_load_h5py(f[k])
            try:
                parms     = _resolve(res["parms"])
                preparms  = _resolve(res["parms_prefit"])
            except (KeyError, ValueError) as exc:
                print(f"  [warn] {path}:{k}: cannot read parms: {exc}")
                continue
            out[k] = {
                "postfit": _values_by_name(parms, params),
                "prefit":  _values_by_name(preparms, params),
            }
    return out


# ---- block collectors -----------------------------------------------------

def collect_toys(block_dir):
    """Block (1)/(2): single fitresults.hdf5 with N ``results_toy<i>`` groups."""
    fit_files = sorted(glob.glob(os.path.join(block_dir, "fitresults*.hdf5")))
    if not fit_files:
        return []
    groups = load_postfits(fit_files[0], PARAMS)
    toy_groups = {k: v for k, v in groups.items() if "toy" in k}

    def _toy_idx(name):
        m = re.search(r"toy(\d+)", name)
        return int(m.group(1)) if m else 1 << 30

    out = []
    for k in sorted(toy_groups, key=_toy_idx):
        out.append({
            "toy":     _toy_idx(k),
            "postfit": toy_groups[k]["postfit"],
            "prefit":  toy_groups[k]["prefit"],
        })
    return out


def collect_shifts(out_base, block_prefix, params=None):
    """Block (3)/(4)/(5): N subdirs ``<prefix>__<slug>``, each one fit.

    ``params`` lets a caller (e.g. block 5) request a different set of
    POUs than the default ``PARAMS``.
    """
    if params is None:
        params = PARAMS
    pattern = os.path.join(out_base, f"{block_prefix}__*")
    out = []
    for subdir in sorted(glob.glob(pattern)):
        slug = os.path.basename(subdir).split("__", 1)[1]
        fit_files = sorted(glob.glob(os.path.join(subdir, "fitresults*.hdf5")))
        if not fit_files:
            print(f"  [warn] no fitresults under {subdir!r}; skipping")
            continue
        groups = load_postfits(fit_files[0], params)
        if not groups:
            print(f"  [warn] no result groups in {fit_files[0]!r}; skipping")
            continue
        # single result group for pseudoData runs (e.g. results_nominal)
        res = next(iter(groups.values()))
        out.append({
            "slug":    slug,
            "postfit": res["postfit"],
            "prefit":  res["prefit"],
        })
    return out


# ---- plotting -------------------------------------------------------------

def plot_pull_histogram(values, param, block_label, outdir, args=None):
    truth = TRUTH[param]
    sigma = SIGMA_PRIOR[param]
    pulls = (np.asarray(values) - truth) / sigma

    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    # 24 bins of width 0.25σ over [-3, 3] (finer than the default 12 × 0.5σ).
    bins = np.linspace(-3.0, 3.0, 25)
    ax.hist(pulls, bins=bins, edgecolor="black", color="C0", alpha=0.7)
    ax.axvline(0, color="black", lw=0.8)
    ax.axvline(-1, color="gray", lw=0.5, ls="--")
    ax.axvline(+1, color="gray", lw=0.5, ls="--")
    ax.set_xlabel(
        r"$(\lambda_\mathrm{post} - \lambda_\mathrm{truth}) "
        r"\,/\, \sigma_\mathrm{prior}$"
    )
    ax.set_ylabel("count")
    ax.set_title(f"{block_label}: {param}")

    median = float(np.median(pulls))
    rms    = float(np.sqrt(np.mean(pulls**2)))
    info = (
        f"N = {len(pulls)}\n"
        f"median pull = {median:+.3f}\n"
        f"RMS pull = {rms:.3f}\n"
        rf"$\sigma_\mathrm{{prior}}$ = {sigma}"
        + f"\n" + rf"$\lambda_\mathrm{{truth}}$ = {truth}"
    )
    ax.text(0.02, 0.98, info, transform=ax.transAxes, va="top", ha="left",
            fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
    ax.set_xlim(-3, 3)
    plt.tight_layout()
    plot_output.save_plot(
        outdir, f"{block_label}__{param}__pull", fig=fig, args=args, dpi=150
    )
    plt.close(fig)
    return median, rms


def plot_pull_histogram_b5(values, param, block_label, outdir,
                           truth_map, sigma_map, args=None):
    """Block-(5) variant of plot_pull_histogram that takes explicit truth
    and sigma maps (so we don't have to mutate the module-level TRUTH /
    SIGMA_PRIOR which blocks (1)-(4) depend on).
    """
    truth = truth_map[param]
    sigma = sigma_map[param]
    pulls = (np.asarray(values) - truth) / sigma

    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    bins = np.linspace(-3.0, 3.0, 25)  # 24 bins of width 0.25σ
    ax.hist(pulls, bins=bins, edgecolor="black", color="C0", alpha=0.7)
    ax.axvline(0, color="black", lw=0.8)
    ax.axvline(-1, color="gray", lw=0.5, ls="--")
    ax.axvline(+1, color="gray", lw=0.5, ls="--")
    ax.set_xlabel(
        r"$(\lambda_\mathrm{post} - \lambda_\mathrm{truth}) "
        r"\,/\, \sigma_\mathrm{ref}$"
    )
    ax.set_ylabel("count")
    ax.set_title(f"{block_label}: {param}")

    median = float(np.median(pulls))
    rms    = float(np.sqrt(np.mean(pulls ** 2)))
    info = (
        f"N = {len(pulls)}\n"
        f"median pull = {median:+.3f}\n"
        f"RMS pull = {rms:.3f}\n"
        rf"$\sigma_\mathrm{{ref}}$ = {sigma}"
        + f"\n" + rf"$\lambda_\mathrm{{truth}}$ = {truth}"
    )
    ax.text(0.02, 0.98, info, transform=ax.transAxes, va="top", ha="left",
            fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
    ax.set_xlim(-3, 3)
    plt.tight_layout()
    plot_output.save_plot(
        outdir, f"{block_label}__{param}__pull", fig=fig, args=args, dpi=150
    )
    plt.close(fig)
    return median, rms, pulls


def plot_shift_scatter(starts, posts, param, block_label, outdir,
                       fit_slope=False, args=None):
    """x = λ_start, y = λ_postfit. Truth + diagonal as reference."""
    starts = np.asarray(starts)
    posts  = np.asarray(posts)
    truth  = TRUTH[param]
    sigma  = SIGMA_PRIOR[param]

    # Plot window
    all_vals = np.concatenate([starts, posts, [truth]])
    lo, hi = float(all_vals.min()), float(all_vals.max())
    margin = max(0.1 * (hi - lo), 0.1 * sigma)
    xs = np.array([lo - margin, hi + margin])

    fig, ax = plt.subplots(figsize=(6.0, 5.5))
    # diagonal y=x → "no motion" / prior dominates
    ax.plot(xs, xs, "--", color="gray", lw=1, label=r"$y = x$ (no motion)")
    # truth horizontal → "data fully recovers"
    ax.axhline(truth, color="C2", lw=1.5, ls="-",
               label=rf"$y = \lambda_\mathrm{{truth}}$ = {truth}")
    # points
    ax.scatter(starts, posts, color="C0", s=55, zorder=10,
               edgecolor="black", linewidths=0.6, label=r"postfit")

    annot_lines = [
        f"N = {len(starts)}",
        rf"$\lambda_\mathrm{{truth}}$ = {truth}",
        rf"$\sigma_\mathrm{{prior}}$ = {sigma}",
    ]
    slope = None
    sigma_data_eff = None
    if fit_slope and len(starts) >= 2 and np.ptp(starts) > 0:
        slope_, intercept = np.polyfit(starts, posts, 1)
        slope = float(slope_)
        xfit = np.array([starts.min(), starts.max()])
        ax.plot(xfit, slope * xfit + intercept, "-", color="C3", lw=1.5,
                label=rf"linear fit (slope = {slope:.3f})")
        annot_lines.append(f"slope = {slope:.3f}")
        if 0.0 < slope < 1.0:
            sigma_data_eff = float(sigma * np.sqrt(slope / (1.0 - slope)))
            annot_lines.append(
                rf"$\sigma_\mathrm{{data}}^\mathrm{{eff}}$ = {sigma_data_eff:.3g}"
            )

    ax.set_xlim(lo - margin, hi + margin)
    ax.set_ylim(lo - margin, hi + margin)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(rf"$\lambda_\mathrm{{start}}$ ({param})")
    ax.set_ylabel(rf"$\lambda_\mathrm{{post}}$ ({param})")
    ax.set_title(f"{block_label}: {param}")
    ax.legend(loc="best", fontsize=9, framealpha=0.85)

    ax.text(0.02, 0.98, "\n".join(annot_lines),
            transform=ax.transAxes, va="top", ha="left", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    plt.tight_layout()
    plot_output.save_plot(
        outdir, f"{block_label}__{param}__scatter", fig=fig, args=args, dpi=150
    )
    plt.close(fig)
    return slope, sigma_data_eff


# ---- summary --------------------------------------------------------------

def write_summary_csv(rows, outdir):
    """One row per (block, fit_id, param) with prefit / postfit / pull."""
    path = os.path.join(outdir, "summary.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["block", "fit_id", "param", "prefit", "postfit",
                    "truth", "sigma_prior", "pull"])
        for r in rows:
            w.writerow(r)
    print(f"  wrote {path}")


def write_summary_md(stats, outdir, block_params=None):
    """A small markdown summary of per-param aggregates.

    ``block_params`` (optional) maps block-label -> list of params to
    iterate. Defaults to ``PARAMS`` for any unlisted block. Used so that
    block (5) can show all 8 POUs while blocks (1)-(4) keep the original
    4-column table.
    """
    block_params = block_params or {}
    path = os.path.join(outdir, "README.md")
    lines = ["# SCETlib NP injection-test summary", ""]
    for block, params in stats.items():
        lines.append(f"## {block}")
        lines.append("")
        lines.append("| Param | N | median pull | RMS pull | slope | σ_data_eff |")
        lines.append("|---|---|---|---|---|---|")
        plist = block_params.get(block, PARAMS)
        for p in plist:
            d = params.get(p, {})
            lines.append(
                f"| {p} | {d.get('N', '-')} | "
                f"{d.get('median_pull', '-')} | "
                f"{d.get('rms_pull', '-')} | "
                f"{d.get('slope', '-')} | "
                f"{d.get('sigma_data_eff', '-')} |"
            )
        lines.append("")
    with open(path, "w") as f:
        f.write("\n".join(lines))
    print(f"  wrote {path}")


# ---- main -----------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("out_base",
                        help="$OUT_BASE — the dir written by the injection script")
    parser.add_argument("-o", "--outdir", default=DEFAULT_OUTDIR,
                        help=f"output dir (default: {DEFAULT_OUTDIR})")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print(f"Reading from: {args.out_base}")
    print(f"Writing to:   {args.outdir}")

    csv_rows = []
    stats = {}

    # ---- blocks (1) and (2): pull histograms -----------------------------
    for name, label in [
        ("01_asimov_priors",   "Block 1 — Asimov toys, priors ON"),
        ("02_asimov_nopriors", "Block 2 — Asimov toys, priors OFF"),
    ]:
        block_dir = os.path.join(args.out_base, name)
        toys = collect_toys(block_dir)
        if not toys:
            print(f"\n[skip] {name}: no toy results found")
            continue
        print(f"\n[{name}] {len(toys)} toys")
        stats[label] = {}
        for p in PARAMS:
            values = [t["postfit"].get(p) for t in toys
                      if p in t.get("postfit", {})]
            if not values:
                continue
            median, rms = plot_pull_histogram(values, p, name, args.outdir, args=args)
            stats[label][p] = {
                "N":            len(values),
                "median_pull":  f"{median:+.3f}",
                "rms_pull":     f"{rms:.3f}",
            }
            print(f"   {p:>14s}  median={median:+.3f}  RMS={rms:.3f}")
            for t in toys:
                if p not in t.get("postfit", {}):
                    continue
                truth = TRUTH[p]; sigma = SIGMA_PRIOR[p]
                pull = (t["postfit"][p] - truth) / sigma
                csv_rows.append([name, f"toy{t['toy']}", p,
                                 t["prefit"].get(p, ""), t["postfit"][p],
                                 truth, sigma, f"{pull:.6f}"])

    # ---- blocks (3) and (4): shift scatter -------------------------------
    for name, label, fit_slope in [
        ("03_inject_nopriors", "Block 3 — pseudoData shifts, priors OFF", False),
        ("04_inject_priors",   "Block 4 — pseudoData shifts, priors ON",  True),
    ]:
        shifts = collect_shifts(args.out_base, name)
        if not shifts:
            print(f"\n[skip] {name}: no shift results found")
            continue
        print(f"\n[{name}] {len(shifts)} shifts")
        stats[label] = {}
        for p in PARAMS:
            starts = [s["prefit"].get(p)  for s in shifts
                      if p in s.get("prefit", {}) and p in s.get("postfit", {})]
            posts  = [s["postfit"].get(p) for s in shifts
                      if p in s.get("prefit", {}) and p in s.get("postfit", {})]
            if not starts:
                continue
            slope, sde = plot_shift_scatter(starts, posts, p, name,
                                            args.outdir, fit_slope=fit_slope, args=args)
            entry = {"N": len(starts)}
            if slope is not None:
                entry["slope"] = f"{slope:.3f}"
            if sde is not None:
                entry["sigma_data_eff"] = f"{sde:.3g}"
            stats[label][p] = entry
            extra = ""
            if slope is not None:
                extra = f"  slope={slope:.3f}"
                if sde is not None:
                    extra += f"  σ_data_eff={sde:.3g}"
            print(f"   {p:>14s}  N={len(starts)}{extra}")
            for s in shifts:
                if p not in s.get("postfit", {}):
                    continue
                truth = TRUTH[p]; sigma = SIGMA_PRIOR[p]
                pull = (s["postfit"][p] - truth) / sigma
                csv_rows.append([name, s["slug"], p,
                                 s["prefit"].get(p, ""), s["postfit"][p],
                                 truth, sigma, f"{pull:.6f}"])

    # ---- block (5): multi-param random + corner shifts -------------------
    block5_label = (
        "Block 5 — pseudoData multi-shift, all 8 POUs floating"
    )
    block5_name = "05_test5"
    shifts5 = collect_shifts(args.out_base, block5_name, params=PARAMS_B5)
    block_params = {}
    if not shifts5:
        print(f"\n[skip] {block5_name}: no shift results found")
    else:
        print(f"\n[{block5_name}] {len(shifts5)} shifts")
        stats[block5_label] = {}
        block_params[block5_label] = PARAMS_B5
        for p in PARAMS_B5:
            values = [s["postfit"].get(p) for s in shifts5
                      if p in s.get("postfit", {})]
            if not values:
                continue
            median, rms, pulls = plot_pull_histogram_b5(
                values, p, block5_name, args.outdir,
                truth_map=TRUTH_B5, sigma_map=SIGMA_REF_B5, args=args,
            )
            stats[block5_label][p] = {
                "N":           len(values),
                "median_pull": f"{median:+.3f}",
                "rms_pull":    f"{rms:.3f}",
            }
            # Per-param console log: median, RMS, and any |pull| > 3
            # entries (potential optimizer failures).
            print(f"   {p:>14s}  N={len(values)}  "
                  f"median={median:+.3f}  RMS={rms:.3f}")
            for s, pull in zip(shifts5, pulls):
                if p not in s.get("postfit", {}):
                    continue
                if abs(pull) > 3.0:
                    print(f"      [|pull|>3] {s['slug']}: "
                          f"postfit={s['postfit'][p]:.4g}  pull={pull:+.3f}")
            # CSV rows
            for s in shifts5:
                if p not in s.get("postfit", {}):
                    continue
                truth = TRUTH_B5[p]
                sigma = SIGMA_REF_B5[p]
                pull = (s["postfit"][p] - truth) / sigma
                csv_rows.append([block5_name, s["slug"], p,
                                 s["prefit"].get(p, ""), s["postfit"][p],
                                 truth, sigma, f"{pull:.6f}"])

    # ---- write summaries -------------------------------------------------
    print()
    write_summary_csv(csv_rows, args.outdir)
    write_summary_md(stats, args.outdir, block_params=block_params)

    print(f"\nDone. {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
