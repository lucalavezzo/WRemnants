#!/usr/bin/env python3
"""SCETlib NP ParamModel TOY ensemble: statistical-spread / bias test.

Sibling to closure_suite.py. Instead of one Asimov fit, this throws
an ensemble of statistically-fluctuated toys and looks at the distribution of
the recovered floating lambdas across toys:

  * each toy is a Poisson throw around the *nominal-at-truth* pseudodata
    (``-t 1 --pseudoData nominal --toysDataMode observed``), so the truth the
    toys scatter around is the datacard's FranksVals lambda_central -- NOT the
    (shifted) optimizer start. With --toysDataMode expected the central would
    instead be expected_yield() at the xparam_default start shift, i.e. the
    shifted point, which is NOT what we want here.
  * constrained nuisances are randomized the standard frequentist way
    (--toysSystRandomize frequentist randomizes the constraint minima / global
    observables); the data is Poisson-fluctuated (--toysDataRandomize poisson).
    The lambdas are unconstrained POU/POIs and are NEVER thrown -- they are
    what we measure.
  * the optimizer START is shifted to the t5_all4 point by default (all four
    lambdas off truth) via the xparam_default=... spec token, matching the
    5th closure test. This doubles as a "converges from a bad start even with noise" check;
    pass --shift '' for a truth start (faster, isolates the pure stat spread).
  * --noHessian --noEDM (values only; the full Hessian OOMs for this model).
    With no per-toy sigma, the ENSEMBLE itself gives the uncertainty: the RMS
    of the fitted lambdas across toys is the statistical sigma, and the mean
    minus truth is the bias.

Each toy runs as its OWN rabbit_fit process with a distinct --seed (so they can
run concurrently and are individually reproducible). NB: a single process with
``-t N`` would run the N toys SERIALLY off one seed stream -- fine but slow.

After all toys finish, reads every fitresults.hdf5 and writes toys_table.{md,csv}
(per-toy raw values + per-param mean / RMS / bias).

Run inside the singularity with setup.sh sourced, e.g.:

    singularity exec --cleanenv <wmassdevrolling.img> bash -c \\
      'source WRemnants/setup.sh && \\
       python3 scripts/rabbit/scetlib_np/toys_suite.py --outbase /ceph/.../scetlib_np_toys_<ts>'
"""

import argparse
import csv
import os
import subprocess
import sys
import time

import numpy as np
import h5py

# ---- datacard truth (FranksVals lambda_central) ---------------------------
TRUTH = {
    "lambda2": 0.40,
    "lambda4": 0.40,
    "lambda2_nu": 0.15,
    "delta_lambda2": 0.00,
}
FLOAT_PARAMS = ("lambda2", "lambda4", "lambda2_nu", "delta_lambda2")
FREEZE = ["lambda_inf", "lambda6", "lambda_inf_nu", "lambda4_nu", "^scetlibNP.*"]

# Default optimizer-start shift = the t5_all4 closure test (all four off truth).
T5_SHIFT = {"lambda2": 0.60, "lambda4": 0.55, "lambda2_nu": 0.18, "delta_lambda2": 0.10}

# ---- inputs (override via env vars; same defaults as the closure suite) ----
FIT_HDF5 = os.environ.get(
    "FIT_HDF5",
    "/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/260528_debug_SCETlibPOIModel/"
    "ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile_excludeSCETlibNP/"
    "ZMassDilepton.hdf5",
)
BTGRID_DIR = os.environ.get(
    "BTGRID_DIR",
    "/ceph/submit/data/user/l/lavezzo/zstuff/Z_COM13_CT18Z_N3p0LL_btgrid_fineall/",
)
PARAMMODEL = "wremnants.postprocessing.scetlib_np.SCETlibNPParamModel"


def parse_shift(s):
    """'lambda2=0.6,lambda4=0.55' -> {'lambda2':0.6,...}; '' -> {}."""
    s = (s or "").strip()
    if not s:
        return {}
    out = {}
    for tok in s.split(","):
        k, v = tok.split("=")
        out[k.strip()] = float(v)
    return out


def shift_str(shift):
    return ",".join(f"{k}={v}" for k, v in shift.items())


def build_cmd(outdir, maxiter, seed, no_hessian, shift, unblind=False, priors=False):
    cmd = [
        "rabbit_fit.py",
        FIT_HDF5,
        "-v", "4",
        # All model knobs are spec tokens (the model does not read env vars):
        # inputs as key=value, the start shift via xparam_default=..., and
        # Gaussian priors on the floating lambdas (DEFAULT_PRIOR_SIGMAS) via
        # priors=1 (rabbit#133 dropped the --paramModelPriors CLA). The priors
        # are centred on the *start* (xparamdefault == the shift) -- i.e.
        # anchored at the offset theorist-central, NOT at the (unknown) truth
        # the toys are thrown around.
        "--paramModel", PARAMMODEL,
        f"btgrid_dir={BTGRID_DIR}",
        *([f"xparam_default={shift_str(shift)}"] if shift else []),
        *(["priors=1"] if priors else []),
        "--freezeParameters", *FREEZE,
        "--minimizerMaxiter", str(maxiter),
        # ---- toy throw: 1 toy, around nominal-at-truth pseudodata ----
        "-t", "1",
        "--pseudoData", "nominal",
        "--toysDataMode", "observed",
        "--toysSystRandomize", "frequentist",
        "--toysDataRandomize", "poisson",
        "--seed", str(seed),
        "-o", outdir + "/",
    ]
    if no_hessian:
        cmd += ["--noHessian", "--noEDM"]
    # observed-mode toys are auto-blinded (blinds the unconstrained pdfAlphaS POI);
    # this is pseudodata, so unblind to read pdfAlphaS absolutely (bias/closure).
    if unblind:
        cmd += ["--unblind"]
    return cmd


def launch(slug, outbase, maxiter, threads, seed, shift, no_hessian, unblind=False, priors=False):
    outdir = os.path.join(outbase, slug)
    os.makedirs(outdir, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env["TF_NUM_INTRAOP_THREADS"] = str(threads)
    env["TF_NUM_INTEROP_THREADS"] = "2"
    logf = open(os.path.join(outdir, "log.txt"), "w")
    cmd = build_cmd(outdir, maxiter, seed, no_hessian, shift, unblind, priors)
    logf.write(f"# slug={slug}  seed={seed}\n")
    logf.write(f"# xparam_default={shift_str(shift)}\n")
    logf.write("# " + " ".join(cmd) + "\n")
    logf.flush()
    proc = subprocess.Popen(cmd, env=env, stdout=logf, stderr=subprocess.STDOUT)
    return {"slug": slug, "seed": seed, "outdir": outdir, "proc": proc, "logf": logf}


def _toy_result_key(path):
    """Discover the single 'results_*' group key in a toy fitresults file."""
    with h5py.File(path, "r") as f:
        keys = [k for k in f.keys() if k.startswith("results_")]
    if not keys:
        raise ValueError(f"no results_* group in {path}")
    # prefer the toy group if multiple are present
    toy = [k for k in keys if "toy" in k]
    key = (toy or keys)[0]
    return key[len("results_"):]


def read_postfit(outdir):
    """Return {param: postfit_value} for the floating params (one toy)."""
    from rabbit import io_tools

    f = os.path.join(outdir, "fitresults.hdf5")
    result = _toy_result_key(f)
    fr = io_tools.get_fitresult(f, result=result)
    post = fr["parms"].get()
    names = [n.decode() if isinstance(n, bytes) else str(n) for n in list(post.axes[0])]
    pv = post.values()
    return {p: (float(pv[names.index(p)]) if p in names else np.nan) for p in FLOAT_PARAMS}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outbase", required=True, help="output base dir for the toy suite")
    ap.add_argument("--ntoys", type=int, default=10)
    ap.add_argument("--threads", type=int, default=32, help="CPU threads per fit")
    ap.add_argument("--max-concurrent", type=int, default=5,
                    help="max toys running at once (throttle for a busy node)")
    ap.add_argument("--maxiter", type=int, default=50)
    ap.add_argument("--seed-base", type=int, default=123456789)
    ap.add_argument("--shift", default=shift_str(T5_SHIFT),
                    help="optimizer-start shift 'p=v,...' (default: t5_all4; '' for truth start)")
    ap.add_argument("--no-hessian", action="store_true", default=True)
    ap.add_argument("--hessian", dest="no_hessian", action="store_false")
    ap.add_argument("--unblind", action="store_true",
                    help="pass --unblind to rabbit_fit (pseudodata: read pdfAlphaS absolutely). "
                         "observed-mode toys are auto-blinded otherwise.")
    ap.add_argument("--priors", action="store_true",
                    help="enable Gaussian priors on the floating lambdas (priors=1 spec token) "
                         "(DEFAULT_PRIOR_SIGMAS), centred on the t5 start.")
    ap.add_argument("--tabulate-only", action="store_true")
    ap.add_argument(
        "--resume",
        action="store_true",
        help="skip toys whose fitresults.hdf5 is already complete/readable; "
        "(re)launch only the missing/partial ones (removing any partial file first)",
    )
    args = ap.parse_args()

    shift = parse_shift(args.shift)
    slugs = [f"toy{i:02d}" for i in range(1, args.ntoys + 1)]
    seeds = {slug: args.seed_base + i for i, slug in enumerate(slugs, start=1)}

    os.makedirs(args.outbase, exist_ok=True)
    print(f"[toys] outbase   = {args.outbase}", flush=True)
    print(f"[toys] ntoys     = {args.ntoys}  threads/fit = {args.threads}  max_concurrent = {args.max_concurrent}", flush=True)
    print(f"[toys] truth     = {TRUTH}", flush=True)
    print(f"[toys] start shift = {shift_str(shift) or '(truth start)'}", flush=True)
    print(f"[toys] no_hessian = {args.no_hessian}", flush=True)
    print(f"[toys] unblind   = {args.unblind}", flush=True)
    print(f"[toys] priors    = {args.priors}", flush=True)

    if not args.tabulate_only:
        t0 = time.time()
        pending = list(slugs)
        if args.resume:
            keep = []
            for slug in pending:
                outdir = os.path.join(args.outbase, slug)
                try:
                    read_postfit(outdir)  # readable + has the toy result group
                    print(f"[toys] resume: {slug} already complete, skipping", flush=True)
                    continue
                except Exception:  # noqa: BLE001
                    fr = os.path.join(outdir, "fitresults.hdf5")
                    if os.path.exists(fr):
                        os.remove(fr)  # drop the partial/unreadable file before re-running
                        print(f"[toys] resume: removed partial {slug}/fitresults.hdf5", flush=True)
                    keep.append(slug)
            pending = keep
            print(f"[toys] resume: {len(pending)} toy(s) to (re)run: {pending}", flush=True)
        running = []
        while pending or running:
            while pending and len(running) < args.max_concurrent:
                slug = pending.pop(0)
                print(f"[toys] launching {slug} (seed={seeds[slug]}) ...", flush=True)
                running.append(launch(slug, args.outbase, args.maxiter, args.threads,
                                       seeds[slug], shift, args.no_hessian, args.unblind,
                                       args.priors))
            # wait for at least one running job to finish
            done = []
            for job in running:
                rc = job["proc"].poll()
                if rc is not None:
                    job["logf"].close()
                    print(f"[toys]   {job['slug']} done rc={rc}", flush=True)
                    done.append(job)
            running = [j for j in running if j not in done]
            if running and not done:
                # block on the first running job to avoid a busy spin
                job = running[0]
                job["proc"].wait()
            time.sleep(2)
        print(f"[toys] all toys finished in {time.time()-t0:.0f}s", flush=True)

    # ---- tabulate -----------------------------------------------------------
    per_toy = {}  # slug -> {param: value}
    for slug in slugs:
        outdir = os.path.join(args.outbase, slug)
        try:
            per_toy[slug] = read_postfit(outdir)
        except Exception as e:  # noqa: BLE001
            print(f"[toys] WARNING: could not read {slug}: {e}", flush=True)

    # per-toy raw CSV
    raw_path = os.path.join(args.outbase, "toys_raw.csv")
    with open(raw_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["toy", "seed", *FLOAT_PARAMS])
        for slug in slugs:
            if slug in per_toy:
                w.writerow([slug, seeds[slug], *[per_toy[slug][p] for p in FLOAT_PARAMS]])

    # per-param summary
    summary = {}
    for p in FLOAT_PARAMS:
        vals = np.array([per_toy[s][p] for s in slugs if s in per_toy and np.isfinite(per_toy[s][p])])
        n = len(vals)
        mean = float(np.mean(vals)) if n else np.nan
        std = float(np.std(vals, ddof=1)) if n > 1 else np.nan  # ensemble RMS = stat sigma
        bias = mean - TRUTH[p]
        err_mean = std / np.sqrt(n) if (n and np.isfinite(std)) else np.nan
        bias_sig = bias / err_mean if (err_mean and np.isfinite(err_mean) and err_mean > 0) else np.nan
        summary[p] = dict(n=n, mean=mean, std=std, bias=bias, err_mean=err_mean, bias_sig=bias_sig)

    sum_path = os.path.join(args.outbase, "toys_table.csv")
    with open(sum_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["param", "truth", "start", "n", "mean", "rms_statsigma", "mean_minus_truth", "err_on_mean", "bias_over_errmean"])
        for p in FLOAT_PARAMS:
            s = summary[p]
            w.writerow([p, TRUTH[p], shift.get(p, TRUTH[p]), s["n"], s["mean"], s["std"], s["bias"], s["err_mean"], s["bias_sig"]])

    md = []
    md.append("# SCETlib NP toy ensemble\n")
    md.append(f"- datacard: `{os.path.basename(FIT_HDF5)}`")
    md.append(f"- truth (FranksVals): {TRUTH}")
    md.append(f"- start shift: {shift_str(shift) or '(truth start)'}")
    md.append(f"- toys: {len(per_toy)}/{args.ntoys} read; `-t 1 --pseudoData nominal --toysDataMode observed`, --noHessian\n")
    md.append("| param | truth | start | N | mean | RMS (stat σ) | mean−truth | err(mean) | bias/err |")
    md.append("|---|---:|---:|:--:|---:|---:|---:|---:|---:|")
    for p in FLOAT_PARAMS:
        s = summary[p]
        md.append(
            "| {p} | {t:.3f} | {st:+.3f} | {n} | {m:+.5f} | {sd:.5f} | {b:+.5f} | {em:.5f} | {bs:+.2f} |".format(
                p=p, t=TRUTH[p], st=shift.get(p, TRUTH[p]), n=s["n"],
                m=s["mean"], sd=s["std"], b=s["bias"], em=s["err_mean"], bs=s["bias_sig"],
            )
        )
    md.append("\n## per-toy fitted values\n")
    md.append("| toy | seed | " + " | ".join(FLOAT_PARAMS) + " |")
    md.append("|---|---:|" + "---:|" * len(FLOAT_PARAMS))
    for slug in slugs:
        if slug in per_toy:
            md.append(f"| {slug} | {seeds[slug]} | " + " | ".join(f"{per_toy[slug][p]:+.5f}" for p in FLOAT_PARAMS) + " |")
    md_text = "\n".join(md) + "\n"
    md_path = os.path.join(args.outbase, "toys_table.md")
    with open(md_path, "w") as fh:
        fh.write(md_text)

    print("\n" + md_text, flush=True)
    print(f"[toys] wrote {raw_path}", flush=True)
    print(f"[toys] wrote {sum_path}", flush=True)
    print(f"[toys] wrote {md_path}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
