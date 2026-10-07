#!/usr/bin/env python3
"""Self-contained checks of wremnants/postprocessing/scetlib_ad/lattice_cs_chi2.py (no SCETlib cache needed).

    python scripts/tests/test_lattice_cs_chi2.py [--card <scetlib_ad datacard>]

Numpy core (always):
  1. the analytic k1 profile == an explicit least-squares fit with k1 floated, at several tunes;
  2. the lattice-only minimum (syst=J, lambda_inf_nu = 2) == the reference fit (chi2 6.7127 / 18 dof,
     lambda2_nu 0.18443, lambda4_nu -0.005920; WRemnantsHelpers lattice-cs-kernel/260923-scetlib-kernel-fit);
  3. Woodbury: the J-mapped systematics leave the minimum unchanged and give the Gauss-Newton covariance
     stat + sum d d^T (sigma 0.05666 / 0.00400, rho -0.911 for syst=J);
  4. syst components: J == Jnf+Jk+Jbt; the default is Jnf+Jbt;
  5. d chi2 / d alpha_s (analytic) == finite differences of the linear model.
TF regularizer (needs --card, a scetlib_ad card carrying the theory-correction meta):
  6. TF chi2 == numpy, gradient == finite differences (lambdas, alphaS, CS TNPs; pert=alphas and pert=live);
  7. exp(2 tau) compensation through a fake fitter (live variable), tau= cross-check, mismatch raises;
  8. pert=live takes the alphaS and TNP maps from the fitter's param model, refuses an alphaS anchor != the table;
  9. ydata=asimov: the term vanishes at the truth.
Exits 1 on any failure."""

import argparse
import os
import sys

import numpy as np

from wremnants.postprocessing.scetlib_ad import lattice_cs_chi2 as L

CARD_A = (
    "/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/260916_Z_2D_card_adcorr/"
    "ZMassDilepton_ptll_yll_adexclpdf/ZMassDilepton.hdf5"
)
TUNES = [(0.18, -0.006), (0.06427, 0.0), (-1e-6, 0.0444), (0.0351, 0.0037), (0.15, 0.0)]
FAIL = []


def check(name, ok, msg=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {msg}", flush=True)
    if not ok:
        FAIL.append(name)


def lam(l2, l4):
    return dict(lambda_inf_nu=2.0, lambda2_nu=l2, lambda4_nu=l4)


def gn_cov(c, l2, l4):
    s2 = 1.0 / np.cosh((l2 * c.u + l4 * c.u**2) / 2.0) ** 2
    X = np.column_stack([-0.5 * s2 * c.u, -0.5 * s2 * c.u**2, c.v])
    return np.linalg.inv(X.T @ c.W @ X)[:2, :2]


def numpy_checks():
    c = L.LatticeCSCore(syst="none")
    md = 0.0
    for l2, l4 in TUNES:
        r0 = c.r0(lam(l2, l4))
        k = -(c.v @ c.W @ r0) / (c.v @ c.W @ c.v)  # explicit 1-parameter GLS
        r = r0 + k * c.v
        md = max(md, abs(r @ c.W @ r - c.chi2(lam(l2, l4))))
    check("1 analytic k1 profile == explicit GLS", md < 1e-10, f"max |d| {md:.1e}")
    cj = L.LatticeCSCore(syst="J")
    res, lm = cj.fit()
    ok = abs(res.fun - 6.712728) < 1e-5 and abs(lm["lambda2_nu"] - 0.184433) < 1e-5
    ok = ok and abs(lm["lambda4_nu"] + 0.0059205) < 1e-6
    check(
        "2 lattice-only minimum == 260923 reference",
        ok,
        f"chi2 {res.fun:.6f}, l2 {lm['lambda2_nu']:.6f}, l4 {lm['lambda4_nu']:.7f}",
    )
    rs, ls = c.fit()
    C = gn_cov(cj, lm["lambda2_nu"], lm["lambda4_nu"])
    s = np.sqrt(np.diag(C))
    ok = (
        abs(rs.fun - res.fun) < 1e-8 and abs(ls["lambda2_nu"] - lm["lambda2_nu"]) < 1e-6
    )
    ok = ok and abs(s[0] - 0.056663) < 1e-5 and abs(s[1] - 0.0039976) < 1e-6
    check(
        "3 Woodbury: J leaves the minimum, GN cov = stat+syst",
        ok,
        f"sigma {s}, rho {C[0, 1] / s[0] / s[1]:+.4f}",
    )
    c3 = L.LatticeCSCore(syst="Jnf+Jk+Jbt")
    ok = np.allclose(cj.cov, c3.cov, rtol=0, atol=1e-15) and L.DEFAULT_SYST == "Jnf+Jbt"
    check("4 syst components: J == Jnf+Jk+Jbt; default Jnf+Jbt", ok)
    h, md = 1e-5, 0.0
    for l2, l4 in TUNES:
        an = cj.dchi2_dalphas(lam(l2, l4))
        fd = (cj.chi2(lam(l2, l4), h) - cj.chi2(lam(l2, l4), -h)) / (2 * h)
        md = max(md, abs(an - fd) / max(1.0, abs(fd)))
    check("5 d chi2 / d alpha_s analytic == FD", md < 1e-8, f"max rel {md:.1e}")


NAMES = (
    "alphaS",
    "lambda2",
    "lambda2_nu",
    "lambda4_nu",
    "pdfEig0",
    "resumTNP_gamma_nu",
    "resumTNP_gamma_cusp",
)


def tf_checks(card):
    import tensorflow as tf

    from rabbit import inputdata
    from rabbit.mappings import helpers as mh
    from rabbit.regularization import helpers as rh

    indata = inputdata.FitInputData(card)
    names = np.array(NAMES)
    n = len(names)
    mod = "wremnants.postprocessing.scetlib_ad.lattice_cs_chi2"

    def build(*args):
        m = mh.load_mapping(f"{mod}.LatticeCSMapping", indata, *args)
        return rh.load_regularizer(f"{mod}.LatticeCSChi2", m, dtype=tf.float64)

    def theta(t_as, l2, l4, tn=0.0, tc=0.0):
        return np.array([t_as, 0.0, (l2 - 0.15) / 0.1, l4 / 0.5, 0.0, tn, tc])

    def arm(r):
        r.set_expectations(tf.zeros(n, tf.float64), None, parms=names)
        return r

    def ad_vs_fd(reg, x0, idx):
        x = tf.Variable(x0)
        with tf.GradientTape() as tp:
            p = reg.compute_nll_penalty(x, None)
        g = tp.gradient(p, x).numpy()
        worst = 0.0
        for i in idx:
            e = 1e-6
            xp, xm = x0.copy(), x0.copy()
            xp[i] += e
            xm[i] -= e
            fp = float(reg.compute_nll_penalty(tf.constant(xp), None))
            fm = float(reg.compute_nll_penalty(tf.constant(xm), None))
            worst = max(worst, abs(g[i] - (fp - fm) / (2 * e)) / max(1.0, abs(g[i])))
        return float(p), worst

    # 6a pert=alphas (the phase-2 LATLIVE mode, spelled both ways)
    for spell in ("alphas=live", "pert=alphas"):
        reg = arm(build("syst=Jnf+Jbt+mu0scale", "offset=0", spell))
        core = L.LatticeCSCore(syst="Jnf+Jbt+mu0scale")
        mv, mg = 0.0, 0.0
        for t_as in (-0.7, 0.0, 0.5):
            for l2, l4 in TUNES:
                p, w = ad_vs_fd(reg, theta(t_as, l2, l4), (0, 2, 3))
                mv = max(mv, abs(2 * p - core.chi2(lam(l2, l4), 0.002 * t_as)))
                mg = max(mg, w)
        check(
            f"6a TF == numpy, AD == FD ({spell})",
            mv < 1e-9 and mg < 1e-5,
            f"{mv:.1e}, {mg:.1e}",
        )
    # 6b pert=live: TNPs and alpha_s second order
    reg = arm(build("offset=0", "pert=live"))
    core = L.LatticeCSCore()
    mv, mg = 0.0, 0.0
    for t_as in (-0.7, 0.5):
        for tn, tc in ((0.0, 0.0), (1.5, -2.0), (-2.0, 1.0)):
            for l2, l4 in TUNES[:3]:
                p, w = ad_vs_fd(reg, theta(t_as, l2, l4, tn, tc), (0, 2, 3, 5, 6))
                tnp = {"resumTNP_gamma_nu": tn, "resumTNP_gamma_cusp": tc}
                mv = max(
                    mv, abs(2 * p - core.chi2(lam(l2, l4), 0.002 * t_as, tnp, "live"))
                )
                mg = max(mg, w)
    check(
        "6b pert=live: TF == numpy, AD == FD incl. TNPs",
        mv < 1e-9 and mg < 1e-5,
        f"{mv:.1e}, {mg:.1e}",
    )
    r_al = arm(build("offset=0", "pert=alphas"))
    x = tf.constant(theta(0.4, 0.06, 0.003, 1.0, 1.0))
    check(
        "6c pert=alphas ignores the TNPs",
        r_al.compute_nll_penalty(x, None)
        == r_al.compute_nll_penalty(tf.constant(theta(0.4, 0.06, 0.003)), None),
    )

    class PM:
        _scetlib_order = tuple(NAMES)
        _rp_quad = np.array([True, True, True, True, False, False, False])
        _rp_id = ~_rp_quad
        _rp_scale = np.ones(7)
        npoi, npou = 1, 6
        params = np.array([nm.encode() for nm in _scetlib_order])

        def __init__(self, a0=0.118, tscale=1.0):
            self._rp_c = np.array(
                [
                    [a0, 0.4, 0.15, 0.0, 0.0, 0.0, 0.0],
                    [0.002, 0.5, 0.1, 0.5, 0.0, 0.0, 0.0],
                    [0.0] * 7,
                ]
            )
            self._rp_scale = np.ones(7)
            self._rp_scale[5] = tscale

    class Fitter:
        def __init__(self, tau, reg, pm=None):
            self.tau = tf.Variable(tau, dtype=tf.float64)
            self.regularizers = [reg]
            self.param_model = pm

        def arm_regularizers(self):
            for r in self.regularizers:
                r.set_expectations(tf.zeros(n, tf.float64), None, parms=names)

    x = tf.constant(theta(0.3, 0.06427, 0.0))
    r1 = arm(build("offset=0"))
    base = float(r1.compute_nll_penalty(x, None))
    r2 = build("offset=0")
    f = Fitter(8.0, r2)
    f.arm_regularizers()
    v8 = float(r2.compute_nll_penalty(x, None)) * np.exp(16.0)
    f.tau.assign(5.0)
    v5 = float(r2.compute_nll_penalty(x, None)) * np.exp(10.0)
    check(
        "7a exp(2 tau) compensation with the live fitter.tau",
        abs(v8 / base - 1) < 1e-12 and abs(v5 / base - 1) < 1e-12,
    )
    try:
        Fitter(5.0, build("tau=8")).arm_regularizers()
        check("7b tau= mismatch raises", False)
    except ValueError:
        check("7b tau= mismatch raises", True)
    r3 = build("pert=live")
    Fitter(8.0, r3, PM(tscale=2.0)).arm_regularizers()
    ok = r3.as_map == (0.118, 0.002) and r3._tnp_map["resumTNP_gamma_nu"] == (
        0.0,
        2.0,
        0.0,
    )
    check(
        "8a pert=live maps (alphaS and TNPs) taken from the param model",
        ok,
        f"{r3.as_map} {r3._tnp_map}",
    )
    try:
        Fitter(8.0, build("pert=live"), PM(0.119)).arm_regularizers()
        check("8b anchor != pert table raises", False)
    except ValueError:
        check("8b anchor != pert table raises", True)
    for mode in ("pert=alphas", "pert=live"):
        ra = arm(build("offset=0", mode, "ydata=asimov"))
        v0 = float(ra.compute_nll_penalty(tf.zeros(n, tf.float64), None))
        v1 = float(ra.compute_nll_penalty(tf.constant(theta(0.5, 0.15, 0.0)), None))
        check(
            f"9 ydata=asimov vanishes at the truth ({mode})",
            abs(v0) < 1e-12 and v1 > 0,
            f"{v0:.1e}, {v1:.2e}",
        )


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--card", default=CARD_A)
    args = ap.parse_args()
    numpy_checks()
    if os.path.exists(args.card):
        tf_checks(args.card)
    else:
        print(f"[SKIP] TF checks: card {args.card} not found")
    print("ALL PASS" if not FAIL else f"FAILED: {FAIL}")
    sys.exit(1 if FAIL else 0)


if __name__ == "__main__":
    main()
