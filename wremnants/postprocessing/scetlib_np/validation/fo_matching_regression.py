"""Regression test for the fixed-order matching onto R's gen bins.

The matching grid M (the FO inputs' own binning, where σ_ns is data) and the gen
grid G (R's binning, where σ_gen must be delivered) are independent. Two regimes,
picked by :meth:`SigmaGenModel._setup_matching` from the grids themselves:

  * M refines G   -> additive:       σ_gen(g) = σ_resum(g) + Σ_{B⊂g} σ_ns(B)
  * otherwise     -> multiplicative: σ_gen(g) = Σ_{c⊂g} σ_resum(c)·f(B(c)),
                                     f(B) = 1 + σ_ns(B)/σ_resum(B)

Every DYTurbo correction is in the first regime (those FO grids are finer than the
analysis binning) and MUST stay bit-identical — that is what most of this test
pins. The coarse NNLOjet stitched exports are in the second; there the test checks
positivity and the conservation property (the matching-grid total is preserved at
any λ), which is what the old centre-in-bin projection violated: on the analysis
gen binning it double-counted 11 FO qT bins, left three |Y| bins with identically
zero σ_ns, and produced 7 non-positive gen bins.

Run inside the container (needs the btgrids on /scratch and the FO inputs):

    python -m wremnants.postprocessing.scetlib_np.validation.fo_matching_regression
"""

import sys

import numpy as np
import tensorflow as tf

from wremnants.postprocessing.scetlib_np import sigma_gen as sg
from wremnants.postprocessing.scetlib_np import theory_correction as tc

# λ_central of the FranksVals tanh_2 tune (the value is irrelevant to the grid
# algebra under test; it only has to be a physical point the model can build at).
LAM = {
    "eff_params": {
        "np_model": "tanh_2",
        "lambda2": 0.4,
        "lambda4": 0.4,
        "delta_lambda2": 0.0,
        "lambda_inf": 1.0,
    },
    "gnu_params": {
        "np_model_nu": "tanh_2",
        "lambda2_nu": 0.15,
        "lambda4_nu": 0.0,
        "lambda_inf_nu": 2.0,
    },
}
# The analysis gen binning: ptZ_binning (0.5-wide from 1 to 12 GeV) x the folded
# 20-quantile |Y| edges. Hardcoded rather than imported so the expected numbers
# below cannot silently change with the binning module.
PTZ = np.array(
    [
        0,
        1,
        1.5,
        2,
        2.5,
        3,
        3.5,
        4,
        4.5,
        5,
        5.5,
        6,
        6.5,
        7,
        7.5,
        8,
        8.5,
        9,
        9.5,
        10,
        10.5,
        11,
        11.5,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        22,
        24,
        26,
        28,
        30,
        33,
        37,
        44,
    ],
    dtype=np.float64,
)
QUANTILE_Y = np.array([0, 0.15, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.8, 2.5])
UNIFORM_PTV = np.linspace(0, 40, 41)  # the sigma_gen_at_lambda CLI default

NNLOJET_TAG = (
    "scetlib_nnlojet_LatticeNPLambda4Bugfix_FranksValsVars_MSHT20aN3LO_N4p0LL_N3LO"
)


def _build(inputs, ptv_edges, absy_edges):
    return sg.SigmaGenModel(
        lambda_central=LAM,
        gen_axes=[("ptVGen", ptv_edges), ("absYVGen", absy_edges)],
        **inputs,
    )


def main():
    failures = []

    def check(name, ok, detail=""):
        print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}", flush=True)
        if not ok:
            failures.append(name)

    ct18 = dict(
        btgrid_dir=sg._default_btgrid_dir(),
        nonsingular_fo_sing=sg._NONSING_FO_SING_DEFAULT,
        nonsingular_dyturbo=sg._NONSING_DYTURBO_DEFAULT,
    )
    print("\n=== CT18Z / DYTurbo: additive regime, values frozen")
    c_fit = _build(ct18, PTZ, QUANTILE_Y)
    check("fit binning -> additive", c_fit._matching["mode"] == "additive")
    c_uni = _build(ct18, UNIFORM_PTV, np.array([0.0, 2.5]))
    check("uniform binning -> additive", c_uni._matching["mode"] == "additive")
    tot = float(c_uni.sigma_gen_central.numpy().sum())
    ns = float(c_uni.sigma_ns.numpy().sum())
    check(
        "Σσ_gen == 1178.7845840518", abs(tot - 1178.7845840517757) < 1e-8, f"({tot!r})"
    )
    check("Σσ_ns == -31.9875527325", abs(ns - (-31.98755273251006)) < 1e-8, f"({ns!r})")

    # The NNLOjet inputs live in a collaborator's area, recorded in the corr's own
    # build metadata; skip rather than fail when they are not reachable.
    try:
        inp = tc.read_corr_inputs(NNLOJET_TAG, proc="Z")
        nnlojet = dict(
            btgrid_dir=tc.btgrid_for(inp["pdf_set"]),
            nonsingular_fo_sing=inp["nnlo_sing"],
            nonsingular_dyturbo=inp["dyturbo_fo"],
        )
        reachable = all(
            sg.fo_input_exists(p) for p in (inp["nnlo_sing"], inp["dyturbo_fo"])
        )
    except Exception as exc:  # noqa: BLE001 - any resolution problem = skip
        print(f"\n=== NNLOjet part SKIPPED (cannot resolve {NNLOJET_TAG}: {exc})")
        reachable = False

    if reachable:
        print("\n=== MSHT20aN3LO / NNLOjet: multiplicative regime")
        a_fit = _build(nnlojet, PTZ, QUANTILE_Y)
        check(
            "fit binning -> multiplicative", a_fit._matching["mode"] == "multiplicative"
        )
        gen = a_fit.sigma_gen_central.numpy()
        check(
            "no non-positive gen bins", bool((gen > 0).all()), f"(min {gen.min():.4f})"
        )
        resum = a_fit._resum_gen(a_fit.sigma_YqT_central).numpy()
        check(
            "Σσ_gen == Σresum + Σσ_ns(covered)",
            abs(gen.sum() - (1202.2235 - 55.5492)) < 2e-3,
            f"({gen.sum():.4f})",
        )
        check(
            "Σresum through the refinement == direct fit-grid resum",
            abs(resum.sum() - 1202.2235) < 2e-3,
            f"({resum.sum():.4f})",
        )

        # Conservation at a SHIFTED λ: the whole point of keeping the factor on the
        # matching grid is that Σ_{c⊂B} σ_gen(c; λ) = σ_resum(B; λ) + σ_ns(B) for
        # every λ, not just λ_central.
        eff = dict(LAM["eff_params"], lambda2=0.55)
        gnu = dict(LAM["gnu_params"], lambda2_nu=0.22)
        s_native = a_fit.sigma_YqT_native(eff, gnu)
        res_ref = a_fit._resum_ref(s_native)
        res_M = tf.matmul(
            a_fit._S_ptV, tf.matmul(res_ref, a_fit._S_absY, transpose_b=True)
        ).numpy()
        shifted = a_fit.sigma_gen(eff, gnu, sigma_YqT=s_native).numpy()
        ns_M = a_fit._ns_M.numpy()
        own_q = a_fit._own_q.numpy()
        own_y = a_fit._own_y.numpy()
        covered = np.zeros(ns_M.shape, dtype=bool)
        covered[
            np.ix_(
                np.unique(own_q[own_q < ns_M.shape[0]]),
                np.unique(own_y[own_y < ns_M.shape[1]]),
            )
        ] = True
        target = float((res_M + ns_M)[covered].sum())
        check(
            "Σσ_gen(λ') == Σ_B[σ_resum(B;λ') + σ_ns(B)]",
            abs(shifted.sum() - target) < 1e-6 * abs(target),
            f"({shifted.sum():.6f} vs {target:.6f})",
        )
        check(
            "shifted-λ gen bins all positive",
            bool((shifted > 0).all()),
            f"(min {shifted.min():.4f})",
        )

        print("\n=== MSHT20aN3LO on the CLI's uniform grid: aligned -> additive")
        a_uni = _build(nnlojet, UNIFORM_PTV, np.array([0.0, 2.5]))
        check("uniform binning -> additive", a_uni._matching["mode"] == "additive")
        check(
            "Σσ_gen == 1119.5999",
            abs(float(a_uni.sigma_gen_central.numpy().sum()) - 1119.5999) < 1e-3,
            f"({float(a_uni.sigma_gen_central.numpy().sum()):.4f})",
        )

    print("\n=== sigma_gen_central - sigma_ns is the resum-only gen spectrum")
    models = [("CT18Z", c_fit)] + ([("aN3LO", a_fit)] if reachable else [])
    for name, mdl in models:
        got = (mdl.sigma_gen_central - mdl.sigma_ns).numpy()
        want = mdl._resum_gen(mdl.sigma_YqT_central).numpy()
        check(
            f"{name}: identity holds",
            bool(np.allclose(got, want, rtol=1e-12, atol=1e-9)),
            f"(max |Δ| {np.abs(got - want).max():.2e})",
        )

    print("\n" + ("ALL PASS" if not failures else f"FAILURES: {failures}"))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
