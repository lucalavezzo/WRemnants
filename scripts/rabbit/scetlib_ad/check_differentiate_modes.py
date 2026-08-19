#!/usr/bin/env python3
"""Verify that the straight-through path gives the SAME derivatives as
differentiating through the SCETlib bridge.

The param model's default (`differentiate=straightthrough`) evaluates SCETlib
once per parameter point and hands autodiff an exact local quadratic:

    d     = p - stop_gradient(p)                # value 0, unit gradient
    sigma = stop_gradient(val) + J.d + 0.5 d^T K d

That looks like it throws derivative information away, so this script checks that
it does not. It builds the model both ways on the same card and compares, at the
same parameter point:

    value          compute(p)
    gradient       d/dp  sum(compute(p))
    Hessian        d2/dp2 sum(compute(p))    (GradientTape.jacobian of the gradient)

`through` is the reference: it calls ScetlibCachedXsecTF inside the graph and lets
TF drive the C++ callbacks via its nested custom_gradient wrappers. Timings are
reported too, since the cost -- not the accuracy -- is the reason for the default.

    python scripts/rabbit/scetlib_ad/check_differentiate_modes.py \
        --card <card>.hdf5 --conf <runcard>.conf --cache <cache>.npz
"""

import argparse
import os
import sys
import time

import numpy as np
import tensorflow as tf

sys.path.insert(
    0,
    os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ),
)

from rabbit.inputdata import FitInputData  # noqa: E402

from wremnants.postprocessing.scetlib_ad.param_model import (  # noqa: E402
    SCETlibADParamModel,
)


def derivatives(model, p):
    """value / gradient / HVP / Hessian of sum(compute(p)), each independently.

    Layered on purpose, because rabbit uses different autodiff machinery at each
    level and they do not all survive a py_function:

    * gradient          -- GradientTape.gradient        (Fitter._loss_val_grad)
    * HVP               -- nested tapes, reverse-over-reverse
                           (Fitter._loss_val_grad_hessp_revrev, the minimizer's
                           default)
    * Hessian           -- GradientTape.jacobian of the gradient
                           (Fitter.loss_val_grad_hess, the postfit covariance)

    Each is attempted separately and a failure is reported rather than raised, so
    the table shows WHICH level a mode stops working at.
    """
    x = tf.Variable(p, dtype=model.indata.dtype)
    out = {}

    out["value"] = float(tf.reduce_sum(model.compute(x)))

    try:
        with tf.GradientTape() as t:
            loss = tf.reduce_sum(model.compute(x))
        out["gradient"] = np.asarray(t.gradient(loss, x).numpy())
    except Exception as e:  # noqa: BLE001 - reporting, not handling
        out["gradient"] = e

    v = np.ones(int(p.size), dtype=np.float64)
    try:
        with tf.GradientTape() as t2:
            with tf.GradientTape() as t1:
                loss = tf.reduce_sum(model.compute(x))
            g = t1.gradient(loss, x)
        hvp = t2.gradient(g, x, output_gradients=tf.constant(v, dtype=x.dtype))
        out["hvp"] = None if hvp is None else np.asarray(hvp.numpy())
    except Exception as e:  # noqa: BLE001
        out["hvp"] = e

    try:
        with tf.GradientTape() as t2:
            with tf.GradientTape() as t1:
                loss = tf.reduce_sum(model.compute(x))
            g = t1.gradient(loss, x)
        h = t2.jacobian(g, x)
        out["hessian"] = None if h is None else np.asarray(h.numpy())
    except Exception as e:  # noqa: BLE001
        out["hessian"] = e

    return out


def rel(a, b):
    """max|a - b| / max|b| -- norm-relative, so legitimately zero entries do not
    blow it up."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--card", required=True)
    ap.add_argument("--conf", required=True)
    ap.add_argument("--cache", required=True)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument(
        "--displace",
        type=float,
        default=0.05,
        help="evaluate the derivatives this far (relatively) off the anchor, so "
        "the comparison is not done at a special point",
    )
    args = ap.parse_args()

    indata = FitInputData(args.card)
    common = dict(
        cache=args.cache,
        conf=args.conf,
        gen_level=True,
        threads=args.threads,
        jitCompile="off",
    )

    modes = (
        ("through", dict(differentiate="through")),
        ("straightthrough", dict(differentiate="straightthrough")),
    )
    results, names = {}, None
    for label, extra in modes:
        m = SCETlibADParamModel(indata, **common, **extra)
        names = list(m._param_order)
        p = np.asarray(m.xparamdefault.numpy(), dtype=np.float64) * (
            1.0 + args.displace
        )
        t0 = time.time()
        results[label] = derivatives(m, p)
        results[label]["_time"] = time.time() - t0

    def status(v):
        if isinstance(v, Exception):
            return f"FAILED: {type(v).__name__}"
        if v is None:
            return "None (no derivative recorded)"
        return "ok"

    print("\nwhich autodiff levels each mode supports:")
    print(
        f"   {'mode':<30}" + "".join(f"{k:>34}" for k in ("gradient", "hvp", "hessian"))
    )
    for label in results:
        print(
            f"   {label:<30}"
            + "".join(
                f"{status(results[label][k]):>34}"
                for k in ("gradient", "hvp", "hessian")
            )
        )

    ref = results["through"]
    print("\nagreement with differentiate=through, where through works:")
    for label, r in results.items():
        if label == "through":
            continue
        bits = [f"value {abs(r['value'] / ref['value'] - 1):.2e}"]
        for k in ("gradient", "hvp", "hessian"):
            a, b = r[k], ref[k]
            if (
                isinstance(a, Exception)
                or isinstance(b, Exception)
                or a is None
                or b is None
            ):
                continue
            bits.append(f"{k} {rel(a, b):.2e}")
        print(f"   {label:<30} " + "   ".join(bits))

    print("\nper-parameter gradient:")
    print(f"   {'parameter':<18}" + "".join(f"{label_[:22]:>24}" for label_ in results))
    for i, n in enumerate(names):
        row = ""
        for label_ in results:
            g = results[label_]["gradient"]
            row += f"{'n/a':>24}" if isinstance(g, Exception) else f"{g[i]:24.10g}"
        print(f"   {n:<18}" + row)

    print("\ncost of one (value, gradient, HVP, Hessian) sweep:")
    for label_, r in results.items():
        print(f"   {label_:<30} {r['_time']:7.2f} s")

    g_ref = ref["gradient"]
    if not isinstance(g_ref, Exception):
        tol = 1e-8
        bad = [
            label_
            for label_, r in results.items()
            if label_ != "through"
            and not isinstance(r["gradient"], Exception)
            and rel(r["gradient"], g_ref) > tol
        ]
        if bad:
            raise SystemExit(f"gradient mismatch beyond {tol:g} for: {bad}")
        print(f"\nGradients agree with the reference to better than {tol:g}.")


if __name__ == "__main__":
    main()
