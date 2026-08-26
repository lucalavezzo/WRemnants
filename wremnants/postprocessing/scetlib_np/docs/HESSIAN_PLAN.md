# Exact compact-derivative path for the SCETlib NP ParamModel Hessian

Status: IMPLEMENTED — the straight-through param model + two-pass GN Hessian
described here are live (see the UPDATE sections below and the
`hessian_straightthrough` / `hessian_gn` spec tokens). Author: closure-fit work,
2026-06.

## TL;DR

rabbit's postfit covariance step OOMs for this ParamModel: it tries to allocate a
single tensor of shape `[3754, 546840, 2000]` (~32.8 TB). The fix lives entirely
in `param_model.py`: replace the bT fold in `compute()` with a *straight-through*
wrap that keeps the exact forward value but hands autodiff the **compact** local
derivatives of the response w.r.t. the `<= 8` lambda. Then the differentiated
graph never contains the `[546840, 2000]` slab, and the Hessian is feasible.
rabbit is not modified.

---

## 1. The objects and where the cost is

The bT fold (`reconstruct_batch_tf`, `btgrid_tf.py`) builds, per btgrid bin and
per bT point, the integrand

```
integrand[g, t] = bT_J0_kernel[g,t] * I_pert[g,t] * exp(C_nu[g,t] * gNP[t]) * Feff[g,t]
sigma[g]        = sum_t  w[t] * integrand[g, t]            # Simpson over bT
```

with shapes (verified from `combined_btgrid.pkl`):

```
g  = btgrid bin index,  Ng  = 546840      # (Q, Y, qT) grid points
t  = bT index,          Nbt = 2000
```

So `integrand` is `[546840, 2000]` and bT is summed out immediately. Everything
downstream is bT-free:

```
sigma[g]  --(Q-integral, |Y|-fold, qT->ptVGen rebin)-->  sigma_gen[gen bins]
sigma_reco[b] = sum_gen R[b, gen] * sigma_gen[gen]        # response matrix R
ratio[b]      = sigma_reco[b] / sigma_reco_central[b]      # b = reco bin, Nreco
rnorm[b,p]    = ratio[b] on the signal column, 1 elsewhere
```

`Nreco` (reco/fit bins) is a separate, much smaller number, and `546840 x 2000`
is purely *internal* to the fold. Call the differentiable map

```
r(lambda) := ratio(lambda)      # shape [Nreco], depends on the <= 8 lambda only
```

`r` is the only place the lambda (hence the bT fold) enter the likelihood.

---

## 2. The likelihood, and exactly which derivatives need the fold

Let `theta` be all the ordinary nuisances (~3750 of them), `lambda` the NP params
(`nlam <= 8`, 4 floating here). Poisson NLL:

```
L(lambda, theta) = sum_b [ mu_b - n_b * ln(mu_b) ] + constraints(theta)
```

The prediction couples lambda in only one place -- multiplicatively on the signal
through `r_b(lambda)`:

```
mu_b(lambda, theta) = r_b(lambda) * s_b * nu_b(theta)  +  (background, no lambda)
```

(`s_b` = signal nominal, `nu_b(theta)` = the usual nuisance response.)

### Gradient

```
dL/dlambda_i = sum_b (1 - n_b/mu_b) * dmu_b/dlambda_i ,   dmu_b/dlambda_i = (dr_b/dlambda_i) * s_b * nu_b
dL/dtheta_a  = sum_b (1 - n_b/mu_b) * dmu_b/dtheta_a + ... ,  dmu_b/dtheta_a = r_b * s_b * (dnu_b/dtheta_a)
```

Key point already visible at first order:
- `dL/dlambda` needs **`dr/dlambda`** (a fold derivative).
- `dL/dtheta`  needs only **`r`** (the fold *value*, already computed in the forward pass) -- NOT a fold derivative.

That is why the *fit* (gradient only) never OOMs: the fold is evaluated once and
differentiated only in the `<= 8` lambda directions.

### Hessian blocks

```
H_{lambda_i lambda_j} = sum_b [ (n_b/mu_b^2)(dmu_b/dlambda_i)(dmu_b/dlambda_j)
                                + (1 - n_b/mu_b) * d2mu_b/dlambda_i dlambda_j ]
      with  d2mu_b/dlambda_i dlambda_j = (d2 r_b / dlambda_i dlambda_j) * s_b * nu_b

H_{lambda_i theta_a}  = sum_b [ (n_b/mu_b^2)(dmu_b/dlambda_i)(dmu_b/dtheta_a)
                                + (1 - n_b/mu_b) * d2mu_b/dlambda_i dtheta_a ]
      with  d2mu_b/dlambda_i dtheta_a = (d r_b / dlambda_i) * s_b * (dnu_b/dtheta_a)

H_{theta_a theta_b}   = sum_b [ (n_b/mu_b^2)(dmu_b/dtheta_a)(dmu_b/dtheta_b)
                                + (1 - n_b/mu_b) * d2mu_b/dtheta_a dtheta_b ]
      with  d2mu_b/dtheta_a dtheta_b  built from  r_b  and  nu-derivatives
```

So, reading off what each block needs from the fold:

| block            | fold quantity needed                | size                       |
|------------------|-------------------------------------|----------------------------|
| `H_{lambda lambda}` | `dr/dlambda` **and** `d2r/dlambda^2` | `[Nreco, nlam]`, `[Nreco, nlam, nlam]` |
| `H_{lambda theta}`  | `dr/dlambda` only                   | `[Nreco, nlam]` (the *same* object)    |
| `H_{theta theta}`   | `r` (value) only -- **no** fold derivative | `[Nreco]`              |

**This corrects my earlier hand-wave.** It is *not* true that "every nuisance row
carries bT dependence" in the sense of needing its own fold computation. The
nuisances DO correlate with lambda -- the cross block `H_{lambda theta}` is
generically nonzero, which is exactly why the profiled uncertainty needs the full
matrix (Schur complement) -- but the only fold derivative that block needs is the
*shared* `dr/dlambda`. The nuisance-nuisance block needs no fold derivative at all.

In short, the entire Hessian's dependence on the bT fold is captured by just two
compact objects:

```
J = dr/dlambda      shape [Nreco, nlam]          (nlam <= 8)
K = d2r/dlambda^2   shape [Nreco, nlam, nlam]
```

both indexed by the `<= 8` lambda -- never by the 3750 nuisances, never by the
`[546840, 2000]` slab.

---

## 3. Why rabbit OOMs anyway (the real inefficiency)

rabbit forms the covariance from

```
hess = t2.jacobian(grad, self.x)        # fitter.py, experimental_use_pfor=True (default)
```

`jacobian` with `pfor` vectorizes the backward pass over **all 3754 components of
`grad` at once**. The forward value `r` appears in the loss graph for *every*
parameter's gradient (it multiplies the signal in every signal bin), so when the
backward of `grad` is traced, it re-enters the fold graph -- which contains the
`[546840, 2000]` integrand -- and pfor stamps a copy of that intermediate for each
of the 3754 rows:

```
[ 3754 , 546840 , 2000 ]  =  (params)  x  (the fold's internal slab)  ~= 32.8 TB
```

The waste is **not** that the nuisances "need the grid." It is that pfor cannot see
that the fold's derivative content is *low rank* -- fully described by `J` and `K`
over `<= 8` directions -- so it mechanically re-materializes the full internal slab
for all 3754 rows instead of reusing the compact `J`/`K`. (Note: for **Asimov** data
`n_b ~= mu_b`, so the `(1 - n_b/mu_b)` terms vanish and even `K` drops out; the
covariance reduces to the Fisher / Gauss-Newton form `J^T A J`, needing only `J`.)

---

## 4. The fix: straight-through exact value + compact derivatives (in `param_model.py`)

Wrap the `lambda -> ratio` map so the exact fold value is preserved but autodiff
sees a compact quadratic carrying `J` and `K`:

```python
val = fold(lam)                          # exact ratio [Nreco]; value only
J   = jac_compact(lam)                    # dr/dlambda     [Nreco, nlam]   (stop_gradient)
K   = hess_compact(lam)                   # d2r/dlambda^2  [Nreco, nlam, nlam] (stop_gradient)
d   = lam - tf.stop_gradient(lam)         # identically 0 in value, unit gradient
ratio = tf.stop_gradient(val) \
        + tf.linalg.matvec(J, d) \
        + 0.5 * tf.einsum('rij,i,j->r', K, d, d)
```

At the evaluation point (`d = 0`):
- value      = `val`  -> exact (forward bit-identical to current `compute()`),
- 1st deriv  = `J`    -> exact local Jacobian,
- 2nd deriv  = `K`    -> exact local curvature.

Because `fold(lam)` is `stop_gradient`'d, rabbit's `jacobian(grad, x)` never
differentiates through it -> **no slab, no OOM**. What pfor now tiles is
`[3754, Nreco, nlam]` ~ tens of GB. Plugging `ratio` into the blocks of section 2
reproduces every block exactly:
- `H_{lambda lambda}` from the `J*d` and `0.5 d^T K d` terms,
- `H_{lambda theta}`  from `J` (via `dmu/dlambda`),
- `H_{theta theta}`   from the preserved value `r` (no fold).

### Computing `J` and `K` compactly

Use **forward-mode AD over the `<= 8` lambda inputs** (not the big outputs):
- `J`: `nlam` JVP passes (`tf.autodiff.ForwardAccumulator`), each peaking at one
  `[546840, 2000]` slab (~9 GB), done sequentially.
- `K`: `~nlam^2` forward-over-forward (or forward-over-reverse) directions, same
  per-direction peak.

Memory peak is one slab; cost is a few dozen fold passes -- fine as a one-time
operation.

### Avoiding any fit slowdown

`J`/`K` must NOT be recomputed every minimizer step. Do the covariance as a
separate **Hessian-only pass at the fixed postfit point**:
- Phase A: fit with `--noHessian` (exact model, cheap) -> postfit lambda. (Done.)
- Phase B: load the postfit and evaluate the covariance once, with the
  straight-through path active. Per-call `J`/`K` cost is then one-time.

Phase B options (decide in P0): (b1) rabbit `--externalPostfit` / `--noFit` if it
drives `jacobian(grad, x)`; (b2) a small standalone covariance script that rebuilds
the loss at the postfit -- fully under our control.

---

## 5. Gauss-Newton / Fisher simplification (start here)

For Asimov closure fits the residual `(1 - n_b/mu_b) -> 0`, so the covariance is the
Fisher information `J^T A J` and needs **only `J`, not `K`**. Therefore:
1. Implement the **linear** straight-through (`J` only) first -- exact for Asimov,
   much simpler.
2. Add the `0.5 d^T K d` term afterward for full correctness on real data.

---

## 6. Phased plan & validation gates

- **P0 - scoping (no code).** Confirm the Phase-B hook (b1 vs b2); decide GN-first
  vs full-K; fix the wrap boundary (`lambda -> ratio`, leaving the one-hot
  signal-column placement outside); write out exact shapes/contractions.
- **P1 - oracle.** Downsample the btgrid (and/or reduce reco bins) so the *full*
  pfor Hessian fits in memory. This is the exact reference to validate against.
- **P2 - implement** the straight-through wrap in `compute()`; assert the forward
  is **bit-identical** to current `compute()` (preserves the validated 0.14% closure).
- **P3 - validate derivatives.** `J` vs finite-difference and vs a tape gradient on
  the small grid; `K` vs finite-difference of the gradient; then the **full
  covariance vs the P1 oracle** (must match to numerical precision).
- **P4 - real problem.** Run Phase B on the actual fit: no OOM, sane sigma;
  cross-check sigma on 1-2 params against an independent **likelihood scan**.
- **P5 - integrate.** Clean flag/entry point, docstring + this note, memory update.

---

## 7. Risks & mitigations

- Second-order correctness is the hard part -> P1 exact small-grid Hessian oracle +
  finite-difference checks before trusting anything.
- `--externalPostfit` may not drive the Hessian as needed -> fallback (b2)
  standalone covariance script.
- Forward must stay exact -> P2 bit-identity assertion.
- `K` cost/accuracy on real data -> GN-first (`J` only) is already exact for Asimov;
  add `K` only when real-data fits need it.

**UPDATE (full-K now works under @tf.function).** `_ratio_compact_hess` (nested
forward-mode `K`) used to raise `IndexError: list index out of range` under
rabbit's `@tf.function`. Root cause (pinned by isolation micro-tests): it is NOT
the `SelectV2`/`tf.where` op — it is the `Equal` op inside the `lambda_inf==0` /
`den==0` masks receiving a *tangent-carrying* input. Building that `Equal` JVP in
nested forward-over-forward mode under `@tf.function` is the TF bug. Fix:
`btgrid_tf._frozen_eq_zero(x) = tf.equal(tf.stop_gradient(x), 0)` freezes the
comparison input in all three masks (`_safe_div`, `F_eff_tf`, `gamma_nu_NP_tf`).
The comparison is a measure-zero boundary `tf.where` never differentiates, so
value and first/second derivatives are unchanged; only the spurious tangent into
`Equal` is removed. Validation (isolation, float64): forward value + 1st-grad
bit-identical to the old `tf.where` code (Δ=0 exactly); full-K under `@tf.function`
matches the exact reverse-mode (jacobian-of-jacobian) Hessian oracle to ≤3e-16
relative, on both a synthetic ratio and the real `reconstruct_batch_tf` path. So
`hessian_straightthrough=1` *without* `hessian_gn=1` is now a
working full-K covariance path; GN stays the Asimov default (exact there, cheaper).

## 8. Sequencing

Hold until the `--noHessian` recovery suite finishes (it provides the postfit lambda
Phase B loads, and frees the node). Then P0 -> P1 are quick; P3 (oracle match) is the
gate that confirms the approach before touching the full problem.

## 9. Full-K for DATA: precomputed chunked K (design, 2026-06-11)

### Why this section exists

Two facts established 2026-06-11:

1. **The in-rabbit full-K pass has never run at full grid scale.** The §7
   validation used a 400-bin grid SLICE (`/tmp/realgrid_fullk.py`, "~1 min
   smoke test, not the full cov pass") — it proved correctness of nested-FA K
   under `@tf.function`, not feasibility. The first full-scale attempt
   (cov pass with `STRAIGHTTHROUGH=1`, no `_GN=1`) was OOM-KILLED at graph
   build: `_ratio_compact_hess`'s Python double loop unrolls all 64 nested-FA
   bT-fold passes INTO rabbit's traced `loss_val_grad_hess` graph. At
   (Nu=284605, Nbt=2000) fp64 that is a ~TB-scale compile/execution peak.

2. **GN (J-only) is an approximation on real/toy data.** Exact for Asimov
   (the K term enters as Σ_b (1 − n_b/μ_b)·K_bij and residuals vanish there),
   but on data the residuals are nonzero. The error is statistically
   suppressed — residuals are O(1/√n), sign-fluctuating, and partially cancel
   over ~50k reco bins, and the omission affects only the λ–λ block (J and
   all nuisance/cross terms stay exact) — but it is not zero. Data
   covariances should use full-K, or at minimum quantify GN-vs-K once on
   toys (gate G4 below).

### Key insight: K never needs to be inside rabbit's graph

In the two-pass recipe the cov pass evaluates everything at the FIXED loaded
postfit λ̂ (`--externalPostfit ... --noFit`). J and K enter
`_ratio_straightthrough` under `tf.stop_gradient` — they are *values at λ̂*,
not differentiated subgraphs. So K can be computed ONCE, eagerly, in small
chunks, at model construction — and enter rabbit's graph as a constant.
The nested-FA passes then never coexist in one graph.

### Implementation plan (all in param_model.py; no rabbit changes)

1. **λ̂ source**: the ParamModel already receives all CLI args as kwargs
   (`ph.load_models` passes `vars(args)`). When `STRAIGHTTHROUGH=1`, `_GN`
   unset, and `externalPostfit` is set: read λ̂ for the 8 compact params by
   name from the externalPostfit fitresults hdf5 (same file rabbit loads
   into `fitter.x`).

2. **Chunked eager K** — `_compact_hess_precompute(lambda_hat)`:
   - one jitted pair kernel, compiled ONCE (tangents are *inputs*, not
     trace-time constants):
     ```python
     @tf.function(jit_compile=True)
     def _pair_jvp(param, ti, tj):
         with tf.autodiff.ForwardAccumulator(param, tj) as acc_j:
             with tf.autodiff.ForwardAccumulator(param, ti) as acc_i:
                 r = self._ratio_from_param(param)
             di = acc_i.jvp(r)
         return acc_j.jvp(di)          # (N_reco,)
     ```
   - Python loop over the 36 (i ≤ j) pairs (K symmetric), eager between
     calls so each pass's intermediates are freed; fill K[r,i,j] = K[r,j,i].
   - cache `self._K_cached = tf.constant(K)` keyed on
     `lambda_hat.tobytes()`.

3. **Use in the graph**: `_ratio_straightthrough(..., use_curvature=True)`
   uses `tf.stop_gradient(self._K_cached)` instead of calling
   `_ratio_compact_hess(param)` inline. Validity condition: the graph is
   only evaluated at λ̂ — true by construction in the `--noFit` cov pass.
   Keep the existing docstring WARNING (never enable during a FIT); add the
   stronger note that cached-K is *only* valid for the one-point cov pass.

4. **Optional**: restrict pairs to unfrozen λ (rabbit masks frozen rows
   anyway): 4 floating → 10 pairs → ~30 s. Default: all 36 (~1–3 min).

### Cost & memory

- Memory: one pair pass ≈ one forward + two tangent streams ≈ 15–20 GB
  (vs ~TB unrolled). Fits any production node.
- Time: 36 × ~2–3 s ≈ 1–3 min one-off per cov pass, + one ~30 s compile.
  ~2–3× the arithmetic of the (infeasible) unrolled graph because each call
  recomputes the primal — irrelevant at this absolute scale.

### Validation gates

- **G1 (oracle)**: chunked-eager K == exact reverse-mode (jac-of-jac)
  Hessian on the 400-bin slice (reuse `validate_k_fix.py` machinery);
  require ≤1e-13 rel, matching §7.
- **G2 (scale)**: full-grid chunked K completes on a production node;
  record peak RSS and wall time.
- **G3 (composition, Asimov)**: in-rabbit cov pass with cached-K on Asimov
  must reproduce the GN reference (`260609_scetlib_np_hessian_GN_mainport_asimov`)
  — the K term multiplies zero residuals there, so any difference is a bug.
- **G4 (physics, data)**: on toys: ST-full-K vs ST-GN vs an exact
  frozen-λ rabbit Hessian + likelihood scan (the §7 open item). This
  quantifies the actual GN error on data BEFORE unblinding — if GN-vs-K is
  negligible against the λ uncertainties, GN may be acceptable for speed;
  that becomes a measured statement instead of an assumption.
