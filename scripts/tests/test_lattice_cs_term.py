#!/usr/bin/env python3
"""Checks of wremnants/postprocessing/scetlib_ad/lattice_cs_term.py (the lattice CS-kernel term whose kernel SCETlib
computes at every step). Needs a SCETlib build with DrellYan.gamma_nu_points (scetlib-cms branch gamma-nu-points);
skips (exit 0) otherwise. No big cache: the calculation is configured live from a runcard, and a small cache can be
given to exercise the bitwise snapshot check against real rules.

    python scripts/tests/test_lattice_cs_term.py [--conf <runcard>] [--cache <small cache.npz>]

  1. the data file holds exactly the per-ensemble ASWZ points and block covariance of the phase-3 inputs file;
  2. the stat-only lattice-only fit on SCETlib's kernel converges (gradient ~ 0) and the J-mapped systematics
     leave its minimum unchanged (Woodbury);
  3. the term's chi2 == the explicit k1 fit of the same residuals (k1 floated by least squares);
  4. TF gradient of the penalty == central finite differences; the parameter the kernel does not depend on has an
     exactly zero gradient;
  5. exp(2 tau) compensation (tau= on the -r line).
Exits 1 on any failure."""

import argparse
import os
import sys
import types

import numpy as np

CONF = (
    "/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/ad_scetlib_caches/"
    "pdf62_y35_260921_y25/cache.conf"
)
FAIL = []


def check(name, ok, msg=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {msg}", flush=True)
    if not ok:
        FAIL.append(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conf", default=CONF)
    ap.add_argument("--cache", default=None)
    args = ap.parse_args()
    try:
        import scetlib_tf

        scetlib_tf.ScetlibGammaNuTF  # noqa: B018
    except (ImportError, AttributeError) as e:
        print(f"[SKIP] no ScetlibGammaNuTF in this SCETlib build ({e})")
        return 0
    import tensorflow as tf

    from wremnants.postprocessing.scetlib_ad import lattice_cs_term as LT
    from wremnants.postprocessing.scetlib_ad import xsec_backend as xb

    conf = (
        args.conf
        if args.cache is None
        else args.cache.replace("cache.npz", "cache.conf")
    )
    _, sigma = xb.configure(conf, threads=4)
    sing, nons = sigma.sub_pieces()
    if args.cache:
        n_eig = sum(
            1
            for n in xb.ScetlibADXsec.cache_param_names(args.cache)
            if n.startswith("pdf_eig")
        )
        sing.set_pdf_eig_params(n_eig)
        nons.set_pdf_eig_params(n_eig)
        fn = scetlib_tf.ScetlibCachedXsecTF.load(args.cache, sing, nons)
    else:
        fn = types.SimpleNamespace(_sing=sing)
    names = list(sing.gradient_param_names())
    anchor = np.array(sing.gradient_central(), float)

    fit = [
        ("alphaS", "alphas", 0.118, 0.002),
        ("lambda2_nu", "np_gnu_lambda2", 0.15, 0.1),
        ("lambda4_nu", "np_gnu_lambda4", 0.0, 0.5),
        ("resumTNP_gamma_nu", "tnp_gamma_nu", 0.0, 1.0),
        (
            "lambda2",
            "np_eff_lambda2",
            float(anchor[names.index("np_eff_lambda2")]),
            0.5,
        ),
    ]
    idx = np.array([names.index(f[1]) for f in fit])
    c0 = np.array([f[2] for f in fit])
    c1 = np.array([f[3] for f in fit])
    held = anchor.copy()
    held[idx] = 0.0
    S = np.zeros((len(names), len(fit)))
    S[idx, np.arange(len(fit))] = 1.0

    class PM:
        params = np.array([f[0].encode() for f in fit])
        nparams = len(fit)
        scetlib_names = names
        core = types.SimpleNamespace(tf_fn=fn)
        _p_base_anchor = anchor.copy()

        @staticmethod
        def scetlib_full_vector_tf(param):
            p = tf.constant(c0) + tf.constant(c1) * tf.cast(
                param[: len(fit)], tf.float64
            )
            return tf.constant(held) + tf.linalg.matvec(tf.constant(S), p)

    def term(**kw):
        m = LT.LatticeCSTermMapping.parse_args(
            types.SimpleNamespace(channel_info={}, procs=[]),
            *[f"{k}={v}" for k, v in kw.items()],
        )
        t = LT.LatticeCSTerm(m, tf.float64)
        t.param_model_override = PM()
        t.set_expectations(None, None, parms=PM.params)
        return t

    rules = "require" if args.cache else "any"
    T = term(syst="Jnf+Jbt", tau=8.0, rules=rules)
    c = T.core

    old = np.load(
        os.path.join(os.path.dirname(LT.__file__), "data", "lattice_aswz_inputs.npz")
    )
    check(
        "1. data == phase-3 inputs (points, covariance)",
        all(
            np.array_equal(getattr(c, a), old[b])
            for a, b in (
                ("b_fm", "b_fm"),
                ("a_fm", "a_fm"),
                ("y", "y"),
                ("cov_stat", "cov_stat"),
            )
        ),
    )
    fnom = c.fit_nominal
    check(
        "2. stat-only lattice fit converged",
        bool(np.max(np.abs(fnom["grad"])) < 1e-7),
        f"|grad| {np.max(np.abs(fnom['grad'])):.1e}, chi2 {fnom['chi2']:.4f}, lambda {fnom['lam']}",
    )
    dmin = float(np.max(np.abs(c.fit_final["lam"] - fnom["lam"])))
    check("2. J-mapped systematics keep the minimum", dmin < 1e-7, f"{dmin:.1e}")

    th = np.array([0.5, -0.8, 0.012, 0.6, 0.3])
    p = held + S @ (c0 + c1 * th)
    r0 = c.zeta(p) - c.y
    W = np.linalg.inv(c.cov)
    k1 = -(c.v @ W @ r0) / (c.v @ W @ c.v)
    r = r0 + k1 * c.v
    chi2_fit = float(r @ W @ r)
    chi2_term = float(T.chi2_tf(tf.constant(th)))
    check(
        "3. chi2 == explicit k1 fit",
        abs(chi2_term - chi2_fit) < 1e-10 * max(1, chi2_fit),
        f"{chi2_term - chi2_fit:+.1e}",
    )

    x = tf.Variable(th)
    with tf.GradientTape() as tape:
        pen = T.compute_nll_penalty(x, None)
    g = tape.gradient(pen, x).numpy()
    steps = np.array([1e-3, 1e-3, 1e-5, 1e-2, 1e-2])
    fd = np.zeros(len(th))

    def cd(i, h):
        e = np.eye(len(th))[i] * h
        up = float(T.chi2_tf(tf.constant(th + e)))
        dn = float(T.chi2_tf(tf.constant(th - e)))
        return (up - dn) / (2 * h) * 0.5 * np.exp(-16.0)

    for i, h in enumerate(steps):
        fd[i] = (4 * cd(i, h / 2) - cd(i, h)) / 3  # Richardson
    rel = float(np.max(np.abs(g[:4] - fd[:4]) / np.max(np.abs(fd[:4]))))
    check("4. TF gradient == FD", rel < 1e-8, f"rel {rel:.1e}")
    check("4. no kernel dependence -> exactly zero gradient", g[4] == 0.0)
    half = 0.5 * (chi2_term - T.offset)
    check(
        "5. penalty == 1/2 (chi2 - offset) exp(-2 tau)",
        abs(float(pen) - half * np.exp(-16.0)) <= 1e-14 * abs(half * np.exp(-16.0)),
    )
    print("ALL PASS" if not FAIL else f"FAILED: {FAIL}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
