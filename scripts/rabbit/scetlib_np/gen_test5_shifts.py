#!/usr/bin/env python3
"""Generate the shift list for SCETlib NP injection test (5).

Test (5) probes the multi-param parameter space, with particular
emphasis on the negative-lambda4 / small-lambda_inf corner where the
ParamModel's tanh-based factorisation can become unstable.

For each shift we draw all 8 lambdas from independent Gaussians, then
reject draws violating either hard physical bounds (denominator
positivity) or the tanh_2 stability constraint::

    lambda4 + (1/3) * lambda2**3 / lambda_inf**2 > 0

We collect ``--N`` (default 60) accepted random draws plus 5 hand-crafted
corner cases that sit deep inside the (negative lambda4, small
lambda_inf) region. Output is one entry per line in the format consumed
by the injection-test shell driver::

    lambda2=...,lambda4=...,...  t5_<slug>

Slugs: ``t5_rand00..t5_rand<N-1>`` for random draws,
``t5_corner_neg_0..t5_corner_neg_4`` for the corner cases.

Stdlib + numpy only. Does not import any analysis machinery.
"""

import argparse
import sys

import numpy as np


# Order matters: we emit keys in this order on every line.
PARAM_ORDER = [
    "lambda2",
    "lambda4",
    "lambda6",
    "delta_lambda2",
    "lambda_inf",
    "lambda2_nu",
    "lambda4_nu",
    "lambda_inf_nu",
]

# (mu, sigma) for each lambda's Gaussian draw.
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

# Hard physical bounds: denominators in the bT integral must stay > 0.
def _hard_bounds_ok(p):
    return (p["lambda2"]       > 0.05
            and p["lambda_inf"]    > 0.05
            and p["lambda_inf_nu"] > 0.05)


def _tanh2_stability_ok(p):
    """tanh_2 coefficient at large bT must be positive."""
    return (p["lambda4"]
            + (1.0 / 3.0) * (p["lambda2"] ** 3) / (p["lambda_inf"] ** 2)) > 0.0


def _accept(p):
    return _hard_bounds_ok(p) and _tanh2_stability_ok(p)


def _draw_one(rng):
    return {k: float(rng.normal(mu, sigma)) for k, (mu, sigma) in PRIORS.items()}


def _format_entry(p, slug):
    parts = [f"{k}={p[k]:.6g}" for k in PARAM_ORDER]
    return ",".join(parts) + "  " + slug


# Corner cases: all 8 lambdas explicit; only lambda4, lambda_inf differ
# from the runcard centrals.
CORNER_DEFAULTS = {
    "lambda2":        0.4,
    "lambda6":        0.0,
    "delta_lambda2":  0.0,
    "lambda2_nu":     0.15,
    "lambda4_nu":     0.0,
    "lambda_inf_nu":  2.0,
}

CORNERS = [
    # (lambda4, lambda_inf, descriptor)
    (-0.50, 0.20, "deep negative, well inside stability"),
    (-0.20, 0.30, "moderate negative"),
    (-0.10, 0.30, "mild"),
    (-0.80, 0.15, "extreme"),
    (-0.10, 0.40, "near boundary"),
]


def _build_corner(lambda4, lambda_inf):
    p = dict(CORNER_DEFAULTS)
    p["lambda4"]    = lambda4
    p["lambda_inf"] = lambda_inf
    return p


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--output", default="-",
                    help="output path (default: stdout)")
    ap.add_argument("--seed", type=int, default=42,
                    help="RNG seed (default: 42)")
    ap.add_argument("--N", type=int, default=60,
                    help="number of accepted random draws (default: 60)")
    args = ap.parse_args(argv)

    rng = np.random.default_rng(args.seed)

    # ---- random draws ------------------------------------------------------
    accepted = []
    rejected = 0
    while len(accepted) < args.N:
        cand = _draw_one(rng)
        if _accept(cand):
            accepted.append(cand)
        else:
            rejected += 1
        if rejected > 100 * args.N:
            raise RuntimeError(
                f"rejection rate too high: {rejected} rejected after "
                f"{len(accepted)} accepted; check priors / bounds"
            )

    # ---- corner cases ------------------------------------------------------
    corners = []
    for i, (lam4, lam_inf, descr) in enumerate(CORNERS):
        p = _build_corner(lam4, lam_inf)
        if not _accept(p):
            raise RuntimeError(
                f"corner case {i} ({descr}) violates stability: "
                f"lambda4={lam4}, lambda_inf={lam_inf}"
            )
        corners.append(p)

    # ---- emit --------------------------------------------------------------
    lines = []
    for i, p in enumerate(accepted):
        lines.append(_format_entry(p, f"t5_rand{i:02d}"))
    for i, p in enumerate(corners):
        lines.append(_format_entry(p, f"t5_corner_neg_{i}"))

    out_text = "\n".join(lines) + "\n"
    if args.output == "-" or args.output is None:
        sys.stdout.write(out_text)
    else:
        with open(args.output, "w") as f:
            f.write(out_text)

    print(
        f"[scetlib_np_gen_test5_shifts] accepted={len(accepted)} "
        f"rejected={rejected} corners={len(corners)} seed={args.seed}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
