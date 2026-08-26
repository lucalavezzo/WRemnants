"""Validation of the ``tanh_6_abs`` damping-fold NP forms (``btgrid_tf``).

``tanh_6_abs`` is ``tanh_6`` with its damping exponent pushed through
``abs_fold_tf`` (|exponent + margin| − margin). What that is supposed to buy, and
what this module asserts:

  1. PARITY — on the physical (damping) region the fold is the identity, so
     tanh_6_abs IS tanh_6, to the last bit, for F_eff and γ_ν and hence for σ_gen.
     This is what lets the card closure / λ-response validated under tanh_6 carry
     over to the folded form unchanged (contrast tanh_6_sigmoid, which only
     reduces to tanh_6 in the bT_cutoff → ∞ limit).
  2. BOUNDEDNESS — for EVERY λ, of either sign, λ∞ included:
     F_eff ≤ 1 + ``abs_margin`` and γ_ν ≤ ``abs_margin_nu``, uniformly in b_T (the
     margin IS the allowed excursion of the function — the ``F_eff ≤ 1.2`` /
     ``γ_ν ≤ 0.2`` cap, as an identity rather than a penalty). This is the wall-free
     property: no infeasible λ is left, so no NPDampingWall term is needed.
  3. DIFFERENTIABILITY — J and K stay finite across the fold, including the
     forward-over-forward nesting the full-K Hessian path uses under
     ``@tf.function``, where per-op JVP support is not a given (the reason
     ``btgrid_tf._frozen_eq_zero`` exists, and why the fold is a frozen-condition
     ``tf.where`` and not ``tf.abs``).
  4. MIRROR DEGENERACY — the honest cost. The tanh argument is ODD under flipping
     (λ2, λ4, λ6, δλ2) together, so at margin = 0 the folded form is EXACTLY
     invariant under λ → −λ: two mirror-image minima, and the sign of a postfit λ
     is a branch choice, not a measurement. A margin > 0 breaks it (the fold's
     offset is not producible by the polynomial, which vanishes at b_T = 0).
  5. NOT MONOTONE, NOT A tanh_6 SHAPE — a folded tune damps but has a kink where
     the exponent crosses −margin, a shape no physical-λ tanh_6 can produce. The
     fold makes the PREDICTION physical, not the λ. Quantified here so nobody
     reads "no wall needed" as "any λ is now a fine NP tune".

Pure form-factor level (no bT-grid, no datacard, no rabbit), so it runs anywhere
TF does:

    python -m wremnants.postprocessing.scetlib_np.validation.abs_fold
"""

import numpy as np
import tensorflow as tf

from wremnants.postprocessing.scetlib_np import btgrid_tf as fz_tf
from wremnants.postprocessing.scetlib_np.params import (
    EFF_MODEL_PARAMS,
    GNU_MODEL_PARAMS,
    active_params,
    const_params,
)

# The fit's own b_T grid (combined_btgrid: 2000 points, 1e-3 … 50 GeV⁻¹, b_bar ==
# b_T since b0_over_bmax = 0), so the bounds are asserted where they are used.
BT = np.geomspace(1e-3, 50.0, 2000)
Y_PROBE = (0.0, 1.5, 2.5, 5.0)

# A physical (damping) tune: every coefficient ≥ 0 ⇒ the fold is inactive.
PHYS = dict(
    lambda2=0.4,
    lambda4=0.3,
    lambda6=0.016,
    delta_lambda2=0.0,
    lambda_inf=1.0,
    lambda2_nu=0.15,
    lambda4_nu=0.05,
    lambda6_nu=0.0007,
    lambda_inf_nu=2.0,
)

# The np-wall-local-minima study's unphysical optima, where plain tanh_6 blows up:
# "C" (2D, best NLL) and the 1D tanh_6_sigmoid escape (δλ2 → −3.35).
UNPHYS = {
    "C (λ2=1.10, λ4=−0.61)": dict(PHYS, lambda2=1.10, lambda4=-0.61, lambda6=0.01),
    "C + δλ2=−3.35": dict(
        PHYS, lambda2=1.10, lambda4=-0.61, lambda6=0.01, delta_lambda2=-3.35
    ),
    "λ6 < 0": dict(PHYS, lambda6=-0.5),
    "λ2 < 0 (no small-b turn-on)": dict(PHYS, lambda2=-0.5),
    "λ∞ < 0 (both sides)": dict(PHYS, lambda_inf=-1.0, lambda_inf_nu=-2.0),
    "CS λ2ν, λ4ν, λ6ν < 0": dict(
        PHYS, lambda2_nu=-0.5, lambda4_nu=-0.3, lambda6_nu=-0.01
    ),
}


def _eff(values, margin=0.0):
    d = {k: values[k] for k in values if k in active_params(np_model="tanh_6_abs")}
    d["abs_margin"] = margin
    return d


def _gnu(values, margin=0.0):
    d = {k: values[k] for k in values if k in active_params(np_model_nu="tanh_6_abs")}
    d["abs_margin_nu"] = margin
    return d


def _F(values, y, model="tanh_6_abs", margin=0.0):
    return fz_tf.F_eff_tf(y, BT, _eff(values, margin), np_model=model).numpy()


def _G(values, model="tanh_6_abs", margin=0.0):
    return fz_tf.gamma_nu_NP_tf(BT, _gnu(values, margin), np_model_nu=model).numpy()


def main():
    fails = []

    def chk(name, cond, detail=""):
        print(("PASS" if cond else "FAIL"), name, detail)
        if not cond:
            fails.append(name)

    print("=" * 78)
    print("1. registry / vocabulary")
    print("=" * 78)
    chk(
        "tanh_6_abs adds exactly abs_margin to tanh_6 (TMD)",
        set(EFF_MODEL_PARAMS["tanh_6_abs"]) - set(EFF_MODEL_PARAMS["tanh_6"])
        == {"abs_margin"},
    )
    chk(
        "tanh_6_abs adds exactly abs_margin_nu to tanh_6 (CS)",
        set(GNU_MODEL_PARAMS["tanh_6_abs"]) - set(GNU_MODEL_PARAMS["tanh_6"])
        == {"abs_margin_nu"},
    )
    chk(
        "both margins are declared shape CONSTANTS",
        const_params(np_model="tanh_6_abs", np_model_nu="tanh_6_abs")
        == {"abs_margin", "abs_margin_nu"},
    )

    print()
    print("=" * 78)
    print("2. parity: on the physical region tanh_6_abs IS tanh_6 (bit-exact)")
    print("=" * 78)
    for y in Y_PROBE:
        a = _F(PHYS, y, "tanh_6_abs")
        b = _F(PHYS, y, "tanh_6")
        chk(
            f"F_eff parity at y={y}",
            np.array_equal(a, b),
            f"max|Δ|={np.max(np.abs(a - b)):.3g}",
        )
    ga, gb = _G(PHYS, "tanh_6_abs"), _G(PHYS, "tanh_6")
    chk("γ_ν parity", np.array_equal(ga, gb), f"max|Δ|={np.max(np.abs(ga - gb)):.3g}")
    # a nonzero margin must not disturb parity either: fold(x)=x for x ≥ −margin
    for margin in (0.0, 0.01, 0.1):
        a = _F(PHYS, 2.5, "tanh_6_abs", margin=margin)
        chk(
            f"F_eff parity holds at margin={margin}",
            np.array_equal(a, _F(PHYS, 2.5, "tanh_6")),
        )
        g = _G(PHYS, "tanh_6_abs", margin=margin)
        chk(
            f"γ_ν parity holds at margin={margin}",
            np.array_equal(g, _G(PHYS, "tanh_6")),
        )
    # NP off (λ∞ = 0) is untouched: fold(0) = 0 for any margin ≥ 0
    off = dict(PHYS, lambda_inf=0.0, lambda_inf_nu=0.0)
    chk("λ∞=0 → F_eff ≡ 1", np.array_equal(_F(off, 2.5, margin=0.05), np.ones_like(BT)))
    chk(
        "λ∞_ν=0 → γ_ν ≡ 0",
        np.array_equal(_G(off, margin=0.05), np.zeros_like(BT)),
    )

    print()
    print("=" * 78)
    print("3. boundedness: F_eff ≤ 1 + margin, γ_ν ≤ margin_nu, for ANY λ")
    print("=" * 78)
    for margin in (0.0, 0.02, 0.2):
        tol = 1e-12 * (1.0 + margin)
        for label, vals in UNPHYS.items():
            Fmax = max(np.max(_F(vals, y, margin=margin)) for y in Y_PROBE)
            Fbare = max(np.max(_F(vals, y, "tanh_6")) for y in Y_PROBE)
            gmax = float(np.max(_G(vals, margin=margin)))
            gbare = float(np.max(_G(vals, "tanh_6")))
            chk(
                f"[margin={margin}] {label}: F_eff ≤ 1 + margin",
                Fmax <= 1.0 + margin + tol,
                f"max F_eff={Fmax:.6g} (unfolded {Fbare:.4g})",
            )
            chk(
                f"[margin={margin}] {label}: γ_ν ≤ margin",
                gmax <= margin + tol,
                f"max γ_ν={gmax:+.6g} (unfolded {gbare:+.4g})",
            )
            chk(
                f"[margin={margin}] {label}: finite everywhere",
                np.all(np.isfinite(_F(vals, 2.5, margin=margin))) and np.isfinite(gmax),
            )
    # the fold must actually be DOING something in those cases (else the test is
    # vacuous): at least one probe tune anti-damps under the unfolded parent
    chk(
        "unfolded tanh_6 does blow up on these tunes",
        max(max(np.max(_F(v, y, "tanh_6")) for y in Y_PROBE) for v in UNPHYS.values())
        > 10.0,
    )
    # and the cap is TIGHT, not just an upper bound nothing approaches: a tune that
    # wants to anti-damp saturates it. (What the user sets IS what F_eff reaches.)
    for margin in (0.02, 0.2):
        reach = max(
            np.max(_F(v, y, margin=margin)) for v in UNPHYS.values() for y in Y_PROBE
        )
        chk(
            f"[margin={margin}] the cap is reached (tight, not vacuous)",
            reach > 1.0 + 0.9 * margin,
            f"max F_eff over tunes = {reach:.6g} vs cap {1 + margin:.6g}",
        )

    print()
    print("=" * 78)
    print("4. derivatives: J and K finite across the fold (incl. @tf.function)")
    print("=" * 78)
    ORDER = ("lambda2", "lambda4", "lambda6", "delta_lambda2", "lambda_inf")

    def _f_of_vec(v, margin=0.0):
        vals = {n: v[i] for i, n in enumerate(ORDER)}
        vals["abs_margin"] = tf.constant(margin, tf.float64)
        return fz_tf.F_eff_tf(2.5, BT, vals, np_model="tanh_6_abs")

    @tf.function
    def _jvp_hvp(v):
        """One forward-over-forward pass, exactly as _ratio_compact_hess nests."""
        t = tf.one_hot(1, len(ORDER), dtype=tf.float64)  # λ4 direction
        with tf.autodiff.ForwardAccumulator(v, t) as acc_j:
            with tf.autodiff.ForwardAccumulator(v, t) as acc_i:
                out = _f_of_vec(v)
            di = acc_i.jvp(out)
        return di, acc_j.jvp(di)

    for label, vals in list(UNPHYS.items()) + [("physical", PHYS)]:
        v = tf.constant([vals[n] for n in ORDER], dtype=tf.float64)
        J, K = _jvp_hvp(v)
        chk(
            f"J, K finite ({label})",
            bool(np.all(np.isfinite(J.numpy())) and np.all(np.isfinite(K.numpy()))),
            f"|J|max={np.max(np.abs(J.numpy())):.3g} |K|max={np.max(np.abs(K.numpy())):.3g}",
        )
    # exactly on the fold: λ ≡ 0 puts the exponent at 0 for every b_T
    v0 = tf.constant([0.0, 0.0, 0.0, 0.0, 1.0], dtype=tf.float64)
    J0, K0 = _jvp_hvp(v0)
    chk(
        "J, K finite ON the fold (all λ = 0)",
        bool(np.all(np.isfinite(J0.numpy())) and np.all(np.isfinite(K0.numpy()))),
    )

    print()
    print("=" * 78)
    print("5. mirror degeneracy (margin=0) and its breaking (margin>0)")
    print("=" * 78)
    flip = {k: (-v if k != "lambda_inf" else v) for k, v in PHYS.items()}
    flip["lambda_inf_nu"] = PHYS["lambda_inf_nu"]
    for y in Y_PROBE:
        chk(
            f"F_eff(λ) == F_eff(−λ) at margin=0, y={y}",
            np.array_equal(_F(PHYS, y), _F(flip, y)),
        )
    chk("γ_ν(λ) == γ_ν(−λ) at margin=0", np.array_equal(_G(PHYS), _G(flip)))
    d = np.max(np.abs(_F(PHYS, 2.5, margin=0.05) - _F(flip, 2.5, margin=0.05)))
    chk("margin=0.05 BREAKS the mirror", d > 0, f"max|Δ|={d:.3g}")

    print()
    print("=" * 78)
    print("6. what the fold does NOT buy (reported, not asserted)")
    print("=" * 78)
    for label, vals in UNPHYS.items():
        F = _F(vals, 2.5)
        # a folded F_eff still decays overall, but not monotonically: the kink
        # where the exponent crosses −margin turns it around
        rises = int(np.sum(np.diff(F) > 1e-15))
        Fb = _F(vals, 2.5, "tanh_6")
        folded = float(np.mean(Fb > 1.0))
        print(
            f"    {label:32s} folded b_T fraction {folded:5.3f} | "
            f"{rises:4d}/{BT.size - 1} rising steps | F_eff(b_max)={F[-1]:.3g}"
        )
    print(
        "    → a folded tune damps but is NOT a tanh_6 shape (kinked, locally\n"
        "      rising). Judge such a fit with param_model_diagnostics.fold_activity,\n"
        "      not with the wall's silence."
    )

    print()
    print("=" * 78)
    print("7. the margin is a CONSTANT: SCETlibNPParamModel refuses to float it")
    print("=" * 78)
    from wremnants.postprocessing.scetlib_np.param_model import SCETlibNPParamModel

    check = SCETlibNPParamModel._check_shape_constants

    class _Stub:  # just the attributes the guard reads
        _np_model_fit = "tanh_6_abs"
        _np_model_nu_fit = "tanh_6_abs"
        _param_order = tuple(
            active_params(np_model="tanh_6_abs", np_model_nu="tanh_6_abs")
        )

    stub = _Stub()
    start = np.array(
        [0.0 if not n.startswith("abs_margin") else 0.02 for n in stub._param_order]
    )

    def _raises(freeze, values=start):
        try:
            check(stub, freeze, values)
            return False
        except ValueError:
            return True

    chk("unfrozen margins RAISE", _raises(["lambda_inf", "lambda_inf_nu"]))
    chk(
        "explicitly frozen margins pass",
        not _raises(["lambda_inf", "lambda_inf_nu", "abs_margin", "abs_margin_nu"]),
    )
    chk("regex-frozen margins pass", not _raises(["abs_margin.*", "lambda_inf.*"]))
    chk("no freeze list at all (standalone build) passes", not _raises(None))
    neg = np.where(
        [n.startswith("abs_margin") for n in stub._param_order], -0.01, start
    )
    chk("negative margin RAISES even when frozen", _raises(["abs_margin.*"], neg))

    print()
    if fails:
        print(f"{len(fails)} FAILURES:", fails)
        raise SystemExit(1)
    print("ALL PASS")


if __name__ == "__main__":
    main()
