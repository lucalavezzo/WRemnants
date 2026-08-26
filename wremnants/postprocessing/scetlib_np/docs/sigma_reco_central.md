# How `sigma_reco_central` is built

Reference notes for `SCETlibNPParamModel.sigma_reco_central` — the per-reco-bin
reference cross section at λ_central. It is the denominator of the per-step
`rnorm` ratio and the quantity the histmaker closure
(`validation/histmaker_validation.py`) compares against the nominal
reco spectrum. Defined in `param_model.py` (`tf.linalg.matvec(self.R, gen_flat)`).

> **The formula lives in one place.** The full response fold — `P = R_raw/N_gen`,
> `σ_reco = Σ_g P·σ_gen`, the factor definitions, and the helicitySig SUM-vs-UL
> subtlety — is written out in the `param_model.py` module docstring (single
> source of truth). This file is supplementary: the closure logic, measured
> numbers, and verification pointers that don't belong in the docstring.

## Helicity: SUM for R, UL for N_gen (NOT both UL)

R and N_gen reduce the **same** `helicitySig` axis in **opposite** ways, because
they are filled with different tensors:

- **R** is filled with `nominal_weight_helicity` (= nominal_weight ×
  helWeight_tensor): a PARTITION of the event weight into 8 helicity pieces that
  ADD BACK UP. `project(*reco_axes, *gen_axes)` SUMS helicitySig → the physical
  angular-resolved reco×gen yield (the angular dependence lives in the
  cosThetaStar*/phiStar* reco bins).
- **N_gen** is filled with `csAngularMoments`: a moment expansion whose A_i bins
  (0..7) do NOT sum to σ (some are negative). Only UL (`helicitySig = -1`) is the
  angular-integrated total, so N_gen takes UL.

Taking UL of R (treating it like N_gen) would discard the angular partition and
inflate the closure ~15×. (An earlier version of this file claimed "UL on both
sides" — that was wrong; the code sums helicitySig for R. See the inline comment
in `response_matrix.load_R`.)

## Why N_gen is a separate hist, not the reco-marginal of R

`P(b|g)` is the fraction of gen-bin-g events reconstructed into reco bin b
(efficiency × migration). The full gen total is NOT recoverable by summing R
over reco bins — that gives `ε(g)·N_gen(g)`, about 36% of the gen total here,
and reco over/underflow doesn't change it (reco failures are simply absent from
R, not in its flow). ε ranges **0.07–0.54** across gen bins (median 0.47) —
strongly gen-dependent, which is why dividing by the reco-marginal
(migration-only) closes far worse. That gen-dependence is the whole reason
N_gen is a second histogram.

## Units / why the closure works

`sigma_reco_central` carries σ_gen's **cross-section** units; the histmaker
nominal is in **event counts**, so the comparison is meaningful only up to one
overall scale `k`. It closes because the Z MC was generated with SCETlib at
λ_central, so `σ_gen(λ_c; g) ≈ k·N_gen(g)` (same shape). Substituting:

```
σ_reco_central(b) ≈ Σ_g [R_raw(b,g)/N_gen(g)] · k·N_gen(g)
                  =  k · Σ_g R_raw(b,g)
                  =  k · (nominal reco yield)
```

So the closure is really a test that **σ_gen(λ_c; g) has the right gen shape** —
i.e. that `σ_gen(λ_c; g)/N_gen(g)` is ≈ flat in g. In the fit only the ratio
`rnorm = σ_reco(λ)/σ_reco_central` survives, so `k` and any overall constant
cancel; `sigma_reco_central` is purely the reference the λ-varied prediction is
normalized against.

## Verification scripts
(env: set PYTHONPATH to include the WRemnants base and `wums/`)

- reco-marginal vs N_gen (efficiency, over/underflow): `/tmp/verify_reco_marginal.py`
- helicitySig sum vs UL in R: `/tmp/verify_helicity.py`
