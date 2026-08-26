#!/usr/bin/env python3
"""SCETlib NP ParamModel closure suite: inject shifted starts, recover truth.

Five closure tests. Each one:
  * injects a *shifted* start point for the floating lambda via the
    xparam_default=... token in the --paramModel spec (this also moves the
    prior mean, but we run priors OFF here so it is purely the optimizer
    start),
  * fits Asimov pseudo-data generated at the datacard truth
    (``-t 0 --pseudoData nominal``); the truth is the datacard's lambda_central
    (FranksVals), auto-detected from its metadata,
  * runs WITH the Hessian (no --noHessian/--noEDM) so postfit uncertainties
    are available,
and we then check that the postfit lambda recover the truth within uncertainty.

The four floating params are lambda2, lambda4, lambda2_nu, delta_lambda2; the
rest are frozen. delta_lambda2 is shifted in several tests (it floats but its
truth is 0).

Fits run in PARALLEL (the node has hundreds of cores), each capped to
--threads CPU threads so they don't thrash. After all fits finish, the script
reads every fitresults.hdf5 and prints + writes a recovery/uncertainty table
(CSV + markdown).

Run inside the singularity with setup.sh sourced, e.g.:

    singularity exec --cleanenv <wmassdevrolling.img> bash -c \\
      'source WRemnants/setup.sh && \\
       python3 scripts/rabbit/scetlib_np/closure_suite.py --outbase /ceph/.../scetlib_np_closure_<ts>'
"""

import argparse
import csv
import os
import subprocess
import sys
import time

import numpy as np

# ---- datacard truth (FranksVals lambda_central) ---------------------------
# These are the values baked into the datacard (tag
# scetlib_dyturbo_LatticeNPLambda4Bugfix_FranksVals_...); the ParamModel
# auto-detects the same from the fit hdf5 metadata. Only the floating params
# are listed (the frozen ones can't move from truth anyway).
TRUTH = {
    "lambda2": 0.40,
    "lambda4": 0.40,
    "lambda2_nu": 0.15,
    "delta_lambda2": 0.00,
}

FLOAT_PARAMS = ("lambda2", "lambda4", "lambda2_nu", "delta_lambda2")

# Freeze the 4 unconstrained POUs + the old discrete-template NP systs (the
# latter would double-count with our continuous lambda).
FREEZE = ["lambda_inf", "lambda6", "lambda_inf_nu", "lambda4_nu", "^scetlibNP.*"]

# ---- the 5 tests: (slug, {param: shifted_start}) --------------------------
# Each shift is the *start* value; truth is in TRUTH above. Magnitudes are kept
# within ~1 prior sigma so the integral stays physical.
TESTS = [
    ("t1_l2_up", {"lambda2": 0.55}),  # single, +0.15
    ("t2_l4_dn", {"lambda4": 0.25}),  # single, -0.15
    ("t3_delta_up", {"delta_lambda2": 0.15}),  # delta only, +0.15 from 0
    ("t4_l2nu_delta", {"lambda2_nu": 0.20, "delta_lambda2": -0.10}),  # mixed incl delta
    (
        "t5_all4",
        {"lambda2": 0.60, "lambda4": 0.55, "lambda2_nu": 0.18, "delta_lambda2": 0.10},
    ),  # all four off truth
]

# ---- inputs (override via env vars) ---------------------------------------
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


def shift_str(shift):
    return ",".join(f"{k}={v}" for k, v in shift.items())


def build_cmd(outdir, maxiter, priors, no_hessian, shift):
    cmd = [
        "rabbit_fit.py",
        FIT_HDF5,
        "-v",
        "4",
        # All model knobs are spec tokens (env vars are not supported by the
        # model anymore): inputs as key=value, the start shift via
        # xparam_default=..., priors via priors=1.
        "--paramModel",
        PARAMMODEL,
        f"btgrid_dir={BTGRID_DIR}",
        *([f"xparam_default={shift_str(shift)}"] if shift else []),
        *(["priors=1"] if priors else []),
        "--freezeParameters",
        *FREEZE,
        "--minimizerMaxiter",
        str(maxiter),
        "-t",
        "0",
        "--pseudoData",
        "nominal",
        "-o",
        outdir + "/",
    ]
    # The full Hessian OOMs for this param model (rabbit pfor-tiles the 2000-pt
    # bT intermediate over all ~3754 params -> ~33 TB tensor). With --noHessian
    # --noEDM we get postfit VALUES only (variances = nan); uncertainties then
    # need a likelihood scan or a chunked Hessian (--hessianParallelIterations).
    if no_hessian:
        cmd += ["--noHessian", "--noEDM"]
    return cmd


def launch(slug, shift, outbase, maxiter, threads, priors, no_hessian):
    outdir = os.path.join(outbase, slug)
    os.makedirs(outdir, exist_ok=True)
    env = os.environ.copy()
    # Thread caps so parallel fits don't oversubscribe.
    env["OMP_NUM_THREADS"] = str(threads)
    env["TF_NUM_INTRAOP_THREADS"] = str(threads)
    env["TF_NUM_INTEROP_THREADS"] = "2"
    logf = open(os.path.join(outdir, "log.txt"), "w")
    cmd = build_cmd(outdir, maxiter, priors, no_hessian, shift)
    logf.write(f"# slug={slug}\n# xparam_default={shift_str(shift)}\n")
    logf.write("# " + " ".join(cmd) + "\n")
    logf.flush()
    proc = subprocess.Popen(cmd, env=env, stdout=logf, stderr=subprocess.STDOUT)
    return {"slug": slug, "shift": shift, "outdir": outdir, "proc": proc, "logf": logf}


def read_postfit(outdir):
    """Return {param: (start, postfit, err)} for the floating params."""
    from wremnants.postprocessing.scetlib_np import fitresult_lambdas as frl

    f = os.path.join(outdir, "fitresults.hdf5")
    readout = frl.read_lambdas(f, result="nominal", params=FLOAT_PARAMS)
    out = {}
    for p in FLOAT_PARAMS:
        d = readout["params"].get(p, {})
        if d.get("present"):
            out[p] = (d["prefit"], d["postfit"], d["postfit_sigma"])
        else:
            out[p] = (np.nan, np.nan, np.nan)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outbase", required=True, help="output base dir for the suite")
    ap.add_argument("--threads", type=int, default=96, help="CPU threads per fit")
    ap.add_argument("--maxiter", type=int, default=50)
    ap.add_argument("--priors", action="store_true", help="enable Gaussian priors on the lambdas (priors=1 spec token)")
    ap.add_argument(
        "--no-hessian",
        action="store_true",
        help="pass --noHessian --noEDM (values only; avoids the full-Hessian OOM)",
    )
    ap.add_argument("--serial", action="store_true", help="run fits one at a time")
    ap.add_argument(
        "--tabulate-only",
        action="store_true",
        help="skip running; just read existing fitresults under --outbase and tabulate",
    )
    args = ap.parse_args()

    os.makedirs(args.outbase, exist_ok=True)
    print(f"[suite] outbase   = {args.outbase}", flush=True)
    print(f"[suite] FIT_HDF5  = {FIT_HDF5}", flush=True)
    print(f"[suite] truth     = {TRUTH}", flush=True)
    print(f"[suite] floating  = {FLOAT_PARAMS}", flush=True)
    print(f"[suite] tests     = {[s for s, _ in TESTS]}", flush=True)
    print(f"[suite] priors    = {args.priors}  threads/fit = {args.threads}", flush=True)

    if not args.tabulate_only:
        t0 = time.time()
        if args.serial:
            for slug, shift in TESTS:
                print(f"[suite] running {slug} ({shift_str(shift)}) ...", flush=True)
                job = launch(slug, shift, args.outbase, args.maxiter, args.threads, args.priors, args.no_hessian)
                rc = job["proc"].wait()
                job["logf"].close()
                print(f"[suite]   {slug} done rc={rc}", flush=True)
        else:
            jobs = []
            for slug, shift in TESTS:
                print(f"[suite] launching {slug} ({shift_str(shift)}) ...", flush=True)
                jobs.append(launch(slug, shift, args.outbase, args.maxiter, args.threads, args.priors, args.no_hessian))
            for job in jobs:
                rc = job["proc"].wait()
                job["logf"].close()
                print(f"[suite]   {job['slug']} done rc={rc}", flush=True)
        print(f"[suite] all fits finished in {time.time()-t0:.0f}s", flush=True)

    # ---- tabulate -----------------------------------------------------------
    rows = []
    for slug, shift in TESTS:
        outdir = os.path.join(args.outbase, slug)
        try:
            res = read_postfit(outdir)
        except Exception as e:  # noqa: BLE001
            print(f"[suite] WARNING: could not read {slug}: {e}", flush=True)
            continue
        for p in FLOAT_PARAMS:
            start, post, err = res[p]
            truth = TRUTH[p]
            bias = post - truth
            biassig = bias / err if (err and np.isfinite(err) and err > 0) else np.nan
            shifted = p in shift
            rows.append(
                {
                    "test": slug,
                    "param": p,
                    "shifted": "yes" if shifted else "no",
                    "truth": truth,
                    "start": start,
                    "postfit": post,
                    "err": err,
                    "post_minus_truth": bias,
                    "bias_over_sigma": biassig,
                }
            )

    # write CSV
    csv_path = os.path.join(args.outbase, "closure_table.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else [])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # write + print markdown
    md = []
    md.append(f"# SCETlib NP closure suite\n")
    md.append(f"- datacard: `{os.path.basename(FIT_HDF5)}`")
    md.append(f"- truth (FranksVals): {TRUTH}")
    md.append(f"- floating: {', '.join(FLOAT_PARAMS)}; priors: {args.priors}\n")
    md.append(
        "| test | param | shifted | truth | start | postfit | ±σ | post−truth | bias/σ |"
    )
    md.append("|---|---|:--:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        md.append(
            "| {test} | {param} | {shifted} | {truth:.3f} | {start:+.3f} | "
            "{postfit:+.5f} | {err:.5f} | {pmt:+.5f} | {bos:+.2f} |".format(
                test=r["test"],
                param=r["param"],
                shifted=r["shifted"],
                truth=r["truth"],
                start=r["start"],
                postfit=r["postfit"],
                err=r["err"],
                pmt=r["post_minus_truth"],
                bos=r["bias_over_sigma"],
            )
        )
    md_text = "\n".join(md) + "\n"
    md_path = os.path.join(args.outbase, "closure_table.md")
    with open(md_path, "w") as fh:
        fh.write(md_text)

    print("\n" + md_text, flush=True)
    print(f"[suite] wrote {csv_path}", flush=True)
    print(f"[suite] wrote {md_path}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
