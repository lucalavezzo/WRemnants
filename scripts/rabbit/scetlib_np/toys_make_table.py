#!/usr/bin/env python3
"""Tabulate SCETlib NP toy postfits as a PDF table.

One row per toy (postfit per lambda + the pdfAlphaS POI), with `truth` and
`prefit` reference rows on top and `mean`/`RMS` summary rows on the bottom.
Reads each <indir>/toy*/fitresults.hdf5 (skips unreadable/partial ones), pulls
postfit (`parms`) and prefit (`parms_prefit`) for the floating params, and
writes a PDF table (+ markdown + CSV) into <outdir>. Run inside the singularity
with setup.sh sourced.

    python3 scripts/rabbit/scetlib_np/toys_make_table.py \\
        --indir /ceph/.../scetlib_np_toys_<ts> --outdir /path/to/outdir

NOTE on pdfAlphaS: it is the single *unconstrained* parameter of interest
(`hsystsnoconstraint == ['pdfAlphaS']`), NOT a constrained nuisance. In these
toys its truth is fixed at nominal (0) by `--pseudoData nominal`; it is never
thrown -- the fit floats it freely to absorb each toy's Poisson stats, so the
spread of its postfit across toys IS the statistical sigma on alpha_S.

Because the toys run with `--toysDataMode observed`, pdfAlphaS is BLINDED: the
stored postfit = (true fit value) + a fixed per-name N(0, 5) offset (identical
in every toy). So the absolute pdfAlphaS values below are offset by an unknown
constant; only their *spread* (the RMS row) is meaningful without unblinding.
A blind-safe bias check needs an Asimov reference carrying the same offset, or
`--unblind pdfAlphaS`. With `--noHessian` there is also no per-toy sigma.
"""
import argparse
import csv
import glob
import os

import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# pdfAlphaS first (the headline POI), then the floating SCETlib NP lambdas.
PARAMS = ("pdfAlphaS", "lambda2", "lambda4", "lambda2_nu", "delta_lambda2")
PARAM_TEX = (
    r"pdfAlphaS$^{*}$",
    r"$\lambda_2$",
    r"$\lambda_4$",
    r"$\lambda_2^{\nu}$",
    r"$\delta\lambda_2$",
)
# Injected truth: lambdas = FranksVals lambda_central; pdfAlphaS = nominal (0).
TRUTH = {
    "pdfAlphaS": 0.00,
    "lambda2": 0.40,
    "lambda4": 0.40,
    "lambda2_nu": 0.15,
    "delta_lambda2": 0.00,
}
# pdfAlphaS prefit is the per-toy frequentist throw of its (inert) constraint
# centre, not a fixed reference -- blank it in the single prefit row.
PER_TOY_PREFIT = ("pdfAlphaS",)
# Rows that are summaries / references rather than individual toys.
REF_ROWS = ("truth", "prefit", "mean", "RMS")


def _fmt(v):
    return "—" if (v is None or not np.isfinite(v)) else f"{v:+.5f}"


def _result_key(path):
    """Single 'results_*' group key in a (toy) fitresults file."""
    with h5py.File(path, "r") as f:
        keys = [k for k in f.keys() if k.startswith("results_")]
    toy = [k for k in keys if "toy" in k]
    return (toy or keys)[0][len("results_"):]


def _vals(hist):
    names = [n.decode() if isinstance(n, bytes) else str(n) for n in list(hist.axes[0])]
    v = hist.values()
    return {p: (float(v[names.index(p)]) if p in names else float("nan")) for p in PARAMS}


def read_toy(outdir):
    """Return (prefit_dict, postfit_dict) for one toy, or raise if unreadable."""
    from rabbit import io_tools

    f = os.path.join(outdir, "fitresults.hdf5")
    res = _result_key(f)
    fr = io_tools.get_fitresult(f, result=res)
    return _vals(fr["parms_prefit"].get()), _vals(fr["parms"].get())


def summarize(rows):
    """mean / RMS (ddof=1) of the postfit per param over the toy rows only."""
    mean, rms = {}, {}
    for p in PARAMS:
        vals = np.array([v[p] for _, v in rows if np.isfinite(v[p])])
        mean[p] = float(np.mean(vals)) if len(vals) else float("nan")
        rms[p] = float(np.std(vals, ddof=1)) if len(vals) > 1 else float("nan")
    return mean, rms


def render_pdf(table, out_pdf, title, caption=None, args=None):
    """table: list of (label, {param: value}). Writes a matplotlib PDF table."""
    labels = [lab for lab, _ in table]
    cell = [[_fmt(vals[p]) for p in PARAMS] for _, vals in table]

    nrows = len(table)
    fig, ax = plt.subplots(figsize=(8.5, 0.5 * nrows + 1.3))
    ax.axis("off")
    ax.set_title(title, fontsize=11, pad=14)
    tbl = ax.table(
        cellText=cell,
        rowLabels=labels,
        colLabels=list(PARAM_TEX),
        loc="center",
        cellLoc="center",
        rowLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(11)
    tbl.scale(1.0, 1.6)
    # style header + reference/summary rows
    for (r, c), cellobj in tbl.get_celld().items():
        if r == 0:  # column header
            cellobj.set_facecolor("#404a5c")
            cellobj.set_text_props(color="white", fontweight="bold")
        if c == -1:  # row labels
            cellobj.set_text_props(fontweight="bold")
        lab = labels[r - 1] if r >= 1 else None
        if lab in ("truth", "prefit"):  # shade reference rows
            cellobj.set_facecolor("#e8edf4")
        if lab in ("mean", "RMS"):  # shade summary rows
            cellobj.set_facecolor("#fbeede")
    if caption:
        fig.text(0.5, 0.01, caption, ha="center", va="bottom", fontsize=8, wrap=True)
    from wremnants.postprocessing.scetlib_np import plot_output

    outdir, basename = plot_output.split_outpath(out_pdf)
    plot_output.save_plot(outdir, basename, fig=fig, args=args)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, help="toy suite dir to read (contains toy*/)")
    ap.add_argument("--outdir", required=True, help="dir to write the PDF/MD/CSV table into")
    ap.add_argument("--name", default="toys_postfit_table", help="output basename")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    toydirs = sorted(
        d for d in glob.glob(os.path.join(args.indir, "toy[0-9]*")) if os.path.isdir(d)
    )
    # detect whether the toys were run with --unblind (changes how to read pdfAlphaS)
    unblinded = False
    for d in toydirs:
        lg = os.path.join(d, "log.txt")
        if os.path.exists(lg):
            with open(lg) as fh:
                unblinded = "--unblind" in fh.read()
            break
    prefit = None
    rows = []
    for d in toydirs:
        slug = os.path.basename(d)
        try:
            pre, post = read_toy(d)
        except Exception as e:  # noqa: BLE001
            print(f"# skip {slug}: {e}")
            continue
        if prefit is None:
            prefit = dict(pre)
            for p in PER_TOY_PREFIT:  # not a fixed reference -> blank it
                prefit[p] = float("nan")
        rows.append((slug, post))

    if not rows:
        raise SystemExit(f"no readable toys under {args.indir}")

    mean, rms = summarize(rows)

    table = [("truth", TRUTH)]
    if prefit is not None:
        table.append(("prefit", prefit))
    table += rows
    table.append(("mean", mean))
    table.append(("RMS", rms))

    if unblinded:
        caption = (
            f"pdfAlphaS$^{{*}}$: unconstrained POI, truth=0 (nominal), UNBLINDED "
            f"-> mean-0 is the bias, RMS is the stat sigma on alpha_S. "
            f"mean/RMS over {len(rows)} toy(s); --noHessian (no per-toy sigma)."
        )
    else:
        caption = (
            f"pdfAlphaS$^{{*}}$: unconstrained POI, truth=0 (nominal), BLINDED "
            f"(postfit offset by a fixed constant) -> only its RMS is meaningful here. "
            f"mean/RMS over {len(rows)} toy(s); --noHessian (no per-toy sigma)."
        )

    base = os.path.join(args.outdir, args.name)
    # markdown
    md = ["| | " + " | ".join(PARAMS) + " |", "|---|" + "---:|" * len(PARAMS)]
    for label, vals in table:
        md.append(f"| {label} | " + " | ".join(_fmt(vals[p]) for p in PARAMS) + " |")
    md.append("")
    md.append(f"> {caption.replace('$', '').replace('^{{*}}', '*').replace('{', '').replace('}', '')}")
    with open(base + ".md", "w") as fh:
        fh.write("\n".join(md) + "\n")
    # csv
    with open(base + ".csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["row", *PARAMS])
        for label, vals in table:
            w.writerow([label, *[vals[p] for p in PARAMS]])
    # pdf
    render_pdf(
        table,
        base + ".pdf",
        title=f"SCETlib NP toy postfits — {os.path.basename(args.indir.rstrip('/'))}",
        caption=caption,
        args=args,
    )

    print("\n".join(md))
    print(f"\n# wrote {base}.pdf(.png) + {os.path.basename(base)}.log")
    print(f"# wrote {base}.md")
    print(f"# wrote {base}.csv")


if __name__ == "__main__":
    main()
