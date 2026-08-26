"""Measure the bare-metal cost of one σ_gen evaluation, alone, offline.

This script bypasses rabbit entirely. It builds a SCETlibNPParamModel,
then times:

  (1) forward only `_sigma_gen_at(eff, gnu)` — the bT integral itself.
  (2) forward only `compute(param)` — adds R · σ_gen and the (N_reco, N_proc)
      ratio shaping.
  (3) forward + backward of `compute(param)` under tf.GradientTape — what
      rabbit's loss_val_grad invokes per L-BFGS step.

For each, we run N+1 calls and report (first-call wall, mean / std of the
rest). First call wall ≈ XLA compile time; the rest are steady-state.

Run inside the wmass singularity:

    singularity exec --bind /scratch,/work,/home,/ceph \
        /cvmfs/.../wmassdevrolling:testing \
        python3 -m wremnants.postprocessing.scetlib_np.validation.timing
"""

import json
import os
import sys
import time

import numpy as np

from wremnants.postprocessing.scetlib_np.sigma_gen import _default_btgrid_dir
import tensorflow as tf

# FranksVals λ_central fallback: this lets us
# build the model even if the upstream correction pkl isn't locally available.
# Passed to the constructor below via the ``lambda_central`` kwarg.
LAMBDA_CENTRAL = {
    "eff_params": {
        "np_model": "tanh_2", "lambda_inf": 1.0, "lambda2": 0.4,
        "lambda4": 0.4, "lambda6": 0.0, "delta_lambda2": 0.0,
    },
    "gnu_params": {
        "np_model_nu": "tanh_2", "lambda_inf_nu": 2.0,
        "lambda2_nu": 0.15, "lambda4_nu": 0.0,
    },
}
from rabbit.inputdata import FitInputData

from wremnants.postprocessing.scetlib_np.param_model import SCETlibNPParamModel


FIT_HDF5 = ("/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/"
            "260528_debug_SCETlibPOIModel/"
            "ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile"
            "_LatticeNPLambda4Bugfix_FranksVals/ZMassDilepton.hdf5")


N_WARMUP = 1   # first call(s) discarded; XLA-compile bucket
N_TIMED = 10   # number of timed steady-state calls


def _stats(times):
    a = np.array(times)
    return float(a.mean()), float(a.std(ddof=0)), float(a.min()), float(a.max())


def time_block(label, fn, n_warmup, n_timed):
    """Run fn() once per call, time wall clock. Returns dict of stats."""
    warm = []
    print(f"\n=== {label} ===", flush=True)
    for k in range(n_warmup):
        t0 = time.perf_counter()
        fn()
        dt = time.perf_counter() - t0
        warm.append(dt)
        print(f"  warmup {k+1}/{n_warmup}: {dt:.3f} s", flush=True)

    timed = []
    for k in range(n_timed):
        t0 = time.perf_counter()
        fn()
        dt = time.perf_counter() - t0
        timed.append(dt)
        print(f"  call   {k+1}/{n_timed}: {dt:.3f} s", flush=True)

    mean, std, lo, hi = _stats(timed)
    print(f"  --> steady-state: mean {mean:.3f} ± {std:.3f} s "
          f"(min {lo:.3f}, max {hi:.3f}) over {n_timed} calls", flush=True)
    return {
        "label": label,
        "warmup_times": warm,
        "timed_times": timed,
        "mean": mean, "std": std, "min": lo, "max": hi,
    }


def main():
    print("Loading FitInputData …", flush=True)
    indata = FitInputData(FIT_HDF5)

    print("Building SCETlibNPParamModel …", flush=True)
    model = SCETlibNPParamModel(
        indata,
        lambda_central=LAMBDA_CENTRAL,
        btgrid_dir=_default_btgrid_dir(),
        signal_proc="Zmumu",
    )

    lambda_central_np = model.xparamdefault.numpy().astype(np.float64)
    eff_c, gnu_c = model._eff_gnu_from_array(lambda_central_np)

    # As tf.constant (matches what _sigma_gen_at receives from _unpack_params).
    eff_c_tf = {k: (tf.constant(v, dtype=model.indata.dtype)
                    if isinstance(v, float) else v) for k, v in eff_c.items()}
    gnu_c_tf = {k: (tf.constant(v, dtype=model.indata.dtype)
                    if isinstance(v, float) else v) for k, v in gnu_c.items()}

    param_tf = tf.constant(lambda_central_np, dtype=model.indata.dtype)
    param_var = tf.Variable(lambda_central_np, dtype=model.indata.dtype)

    # ---- (1) forward only: _sigma_gen_at -----------------------------------
    def f_sigma_gen():
        out = model._sigma_gen_at(eff_c_tf, gnu_c_tf)
        # Force materialization so we time the actual compute, not graph build.
        _ = out.numpy()

    r1 = time_block("(1) forward  _sigma_gen_at(eff, gnu)",
                    f_sigma_gen, N_WARMUP, N_TIMED)

    # ---- (2) forward only: compute(param) -----------------------------------
    def f_compute():
        out = model.compute(param_tf)
        _ = out.numpy()

    r2 = time_block("(2) forward  compute(param)",
                    f_compute, N_WARMUP, N_TIMED)

    # ---- (3) forward + backward: gradient of sum(compute) w.r.t. param ------
    def f_grad():
        with tf.GradientTape() as tape:
            out = model.compute(param_var)
            loss = tf.reduce_sum(out)
        g = tape.gradient(loss, param_var)
        _ = g.numpy()

    r3 = time_block("(3) fwd+bwd  d/dparam sum(compute(param))",
                    f_grad, N_WARMUP, N_TIMED)

    # ---- summary -----------------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY (mean steady-state per call, full Hankel, no rabbit)")
    print("=" * 70)
    for r in (r1, r2, r3):
        print(f"  {r['label']:50s} {r['mean']*1000:8.1f} ms ± {r['std']*1000:.1f}")
    print("\nfirst-call (compile) wall:")
    for r in (r1, r2, r3):
        wt = r["warmup_times"][0] if r["warmup_times"] else float("nan")
        print(f"  {r['label']:50s} {wt:7.2f} s")
    print()
    print(f"ratio fwd+bwd / fwd_compute  : {r3['mean']/r2['mean']:.2f}x")
    print(f"ratio compute / _sigma_gen_at: {r2['mean']/r1['mean']:.2f}x")

    return 0


if __name__ == "__main__":
    sys.exit(main())
