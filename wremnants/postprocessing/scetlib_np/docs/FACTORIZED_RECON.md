# Factorized bT-grid reconstruction (GPU memory fix)

**Status:** implemented 2026-06-11; the only reconstruction path. The physics
now lives in `SigmaGenModel` (`sigma_gen.py`); `SCETlibNPParamModel` holds one as
`self.core`. The legacy `(Nbins, Nbt)` path / `legacy_recon` token were removed
2026-06-26 (the dedup is bit-exact-verified, so there is nothing to fall back to).
The derived layout is now memoized in `combined_btgrid.factorized.npz` next to the
combined pickle — see "Derived cache" below.

## Problem

The original reconstruction (`btgrid_tf.reconstruct_batch_tf`) lays every
λ-dependent tensor out as `(Nbins, Nbt)` = (546840, 2000) float64 = **8.75 GB
each**. On GPU the fit dies at ParamModel construction (`build_bT_J0_kernel`,
`Op:BesselJ0` OOM on a V100-32GB): `I_pert`, `C_nu` and the kernel's `arg`
intermediate are already resident (3 × 8.75 GB) when the J0 output asks for a
4th, on top of rabbit's 5.98 GB `hlogk`. Even past construction, each fit
step would materialize several more `(Nbins, Nbt)` intermediates
(`exp_g_factor`, gathered `F_eff`, `integrand`, binary-op temps, tape
activations) — a ~50–70 GB peak. The layout fundamentally cannot fit 32 GB.

## Two exact observations

The per-bin integrand is (see the `param_model` module docstring for the full
formula):

    sigma_i = qT_i * SUM_b  w_b * bT_b * J0(qT_i * bT_b)
                          * I_pert[i, b] * exp(C_nu[i, b] * g_nu(b; lambda))
                          * F_eff(Y_i, b; lambda)

1. **qT enters the λ-dependent part only through the J0 kernel** (and the row
   index). The kernel therefore only needs the `NqT = 141` *unique* qT
   values, not 546840 per-bin rows: `(141, 2000)` ≈ 2.3 MB instead of
   8.75 GB.

2. **`I_pert` / `C_nu` rows are bit-identical across qT below the
   profile-scale transition.** SCETlib's qT scale provider
   (`src/qT/Scale_provider.cpp`, `_f_run` / transition function) is piecewise
   in `x = qT/Q`: for `x < x1` the scales are exactly canonical, independent
   of qT. Our grid was produced with `transition_points = [0.2, 0.6, 1.0]`
   and `scale_setting = spectrum` (stored in `grid["config"]`), so for
   `qT < 0.2·Q` (≈ 15–24 GeV across the Q grid) SCETlib evaluates the *same*
   b-space function for every qT bin and the cached rows are byte-for-byte
   copies. Measured on the fineall grid: 546840 bins → **284605 unique
   (I_pert, C_nu) row pairs (1.9× dedup)** — ~6.8k shared low-qT rows used by
   many bins + ~278k qT-unique rows above the transition. (At large qT,
   `C_nu → 0` as resummation switches off — also visible in the grid.)

## Implementation

Three changes, all in `wremnants/postprocessing/scetlib_np/` (rabbit
untouched):

1. **`btgrid_tf.dedup_grid_rows`** (construction-time, numpy): hashes the raw
   bytes of each `(I_pert row, C_nu row, F_eff-Y-index)` triple, builds the
   `bin → unique row` map, then **verifies the grouping by direct
   `np.array_equal` comparison** — the dedup is guaranteed bit-exact and
   makes *no assumption* about profiles, transition points or qT ranges. A
   regenerated grid with different profiles (or no degeneracy) self-adapts;
   the achieved dedup factor is logged at construction so a regression is
   visible.

2. **J0 kernel on `qT_unique`** with the Simpson weights and the per-bin qT
   prefactor folded in: `KwqT = qT_u ⊗ (bT · J0(qT_u · bT)) · w_simpson`,
   shape `(141, 2000)`.

3. **`btgrid_tf.reconstruct_batch_factorized_tf`**: evaluates the λ-dependent
   block on the unique rows only,

       M = I_u * exp(C_u * g_nu) * F_eff_u          # (Nu, Nbt)
       S = M @ KwqT^T                               # (Nu, 141)
       sigma_flat = gather_nd(S, [u(i), qT_idx(i)]) # (Nbins,)

   No `(Nbins, Nbt)` tensor is ever materialized. Downstream (sparse→dense,
   Q-integration, Y/qT rebins, R contraction) is unchanged.

`SCETlibNPParamModel.__init__` builds the factorized constants by default;
`_sigma_YqT_native_at` dispatches on `self.factorized`.

**C_ν sub-dedup (added same day):** C_ν is the ν-evolution coefficient and
depends only on (Q, profile-qT) — not Y — so its standalone unique-row count
is just **1888** (vs 284605 joint rows, 150×). `dedup_grid_rows` therefore
deduplicates the C rows a second time; `exp(C·g)` runs on the (1888, 2000)
table (30 MB) and is gathered back to (Nu, Nbt). Bit-identical (exp of
identical rows + row-replicating gather), ~150× fewer fp64 transcendentals,
and the (Nu, Nbt) C constant never exists on device (−4.55 GB GPU).

## Measured timings (model branch alone, CPU 64 threads, XLA, fp64)

`/tmp` bench, fineall grid, reduce_sum(σ_flat) proxy; legacy layout not
benchmarked (it OOMs GPUs and was never the GPU question):

| pass | factorized | + C_ν sub-dedup |
|---|---|---|
| forward | 684 ms | 815 ms |
| forward + grad (4 λ) | 2944 ms | 2922 ms |
| HVP (rev-over-rev) | 3935 ms | **3072 ms** |

Reading: the branch is memory-bandwidth-bound on CPU (the forward's ~20 GB of
(Nu, Nbt) streams ÷ ~30 GB/s ≈ its runtime), so killing transcendentals only
pays in the second-order pass (−22% HVP — the dominant call in rabbit's
trust-krylov iterations). The (Nu, Nbt) elementwise stream itself is the
exact-arithmetic floor: shrinking it further would mean fewer bT points or
fp32 — both approximations, both ruled out. Fit-iteration share: with
grad ≈ 3 s + a few HVPs ≈ 3 s each vs observed 8–30 s/iteration, the model
branch dominates the CPU fit; rabbit core is the minority cost.

## Exactness

- The dedup + gathers are **bit-exact** (verified at construction; a gather
  replicates rows, and its backward scatter-adds cotangents, so λ-gradients
  are mathematically unchanged).
- The Simpson-as-matmul evaluates the same integrand with the same weights on
  the same sampling; only the floating-point multiplication grouping and
  summation order change (elementwise-then-`reduce_sum` vs `matmul`).
  Expected agreement ≲1e-14 relative — far below the fp64 quadrature error
  already in the grid, and five orders below the 0.14 % reco-closure
  validation. This is a rounding-order effect, **not an approximation**: no
  term is dropped, no precision is reduced (everything stays fp64).

Parity is asserted by `wremnants/postprocessing/scetlib_np/validation/factorized_parity.py`
(legacy vs factorized; σ values and dσ/dλ for the floating λ; λ_central, an
off-central `delta_lambda2 ≠ 0` point, and the static fast path). Result of
the 2026-06-11 run on the fineall grid: see the "Parity results" section
below.

## Memory budget (V100-32GB, fp64)

| | legacy | factorized |
|---|---|---|
| constants | I_pert + C_nu + kernel = 26.3 GB | I_u + C_u = 9.1 GB, KwqT 2.3 MB |
| per-step λ-intermediates | several × 8.75 GB | several × 4.55 GB (Nu, Nbt), S 0.3 GB |
| construction peak | > 35 GB (OOM) | ~14 GB |

With rabbit's own ~6–8 GB (hlogk 5.98 GB + working set) the fit-step peak is
expected to be borderline-OK under XLA fusion. **If it still OOMs**, the next
exact lever is gradient checkpointing (`tf.recompute_grad` around the
`M`/matmul block): tape memory drops to ~constants + one transient, at ~2×
the model-branch FLOPs (tens of ms/step). fp32 is deliberately NOT used
anywhere.

## Derived cache

Building the factorized layout from the raw `combined_btgrid.pkl` (sanitize →
`dense_index_map` → `dedup_grid_rows` with its bit-exact verification → weighted
J0 kernel) is the slow part of construction (~18 GB load + the dedup), and is a
**pure function of the combined grid**. `SigmaGenModel` memoizes the derived
arrays in `combined_btgrid.factorized.npz` next to the pickle, so repeat
constructions skip the raw load and the dedup. Invalidated when the combined
pickle changes (mtime+size, recorded in the npz header) **or** when the
derivation code changes (bump `_FACTORIZED_SCHEMA_VERSION` in `sigma_gen.py`).
Freshness also requires the combined pickle itself to be up to date vs its
shards (`btgrid_cache.is_combined_fresh`). Writes are PID-temp + atomic
`os.replace`, so concurrent fit jobs don't clobber each other.

## Knobs

None. The reconstruction is unconditionally the factorized path (the legacy
`(Nbins, Nbt)` path and the `legacy_recon` token were removed — the dedup is
bit-exact-verified, so there is no parity fallback to keep), and the derived
cache is always on with automatic staleness detection (no flag). The combined
pickle must be present alongside the `.npz`: it is what the cache is verified
against (mtime+size + fresh-vs-shards), and if it is absent the layout is rebuilt
rather than the `.npz` trusted blindly.

## Parity results

2026-06-11, fineall grid, CPU, fp64; tolerance 1e-12 (rtol + atol, see below):

- **dedup:** 546840 bins → 284605 unique rows (1.92×), verified bit-exact at
  construction.
- **dσ/dλ** (λ2, λ4, δλ2, λ2_ν; λ_central + off-central δλ2≠0 point): all
  agree at **≤ 2e-14** relative.
- **σ_flat:** 542413 of 542960 non-zero bins agree at < 1e-12 relative. The
  remaining 547 bins (all at qT ∈ [48, 100], i.e. the fixed-order tail) show
  up to 3.6e-10 *relative-to-themselves* — these are cancellation-dominated
  Hankel sums: the worst bin (Q=96, qT=100) has σ ≈ 3e-7 (≈1e-8 of the
  spectrum peak) and a Simpson condition number κ = Σ|wᵢyᵢ|/|Σwᵢyᵢ| ≈ 6.5e6,
  so summation-order noise is κ·ε_fp64 ≈ 7e-10 — exactly what is observed.
  In **absolute** terms the worst difference is ~7e-17 of the spectrum peak:
  below the fp64 resolution of the peak bins, and below what fp64 determines
  in the LEGACY path itself (reordering its own sum moves those bins as
  much). The Q-integrated (NY, NqT) object shows the same picture (max
  1.1e-10 relative, 25 entries > 1e-12, all in the same high-qT tail).

The parity script therefore uses an rtol+atol metric: a bin fails only if it
disagrees relatively AND the absolute difference exceeds tol × the array's
peak scale. Gradient agreement at 1e-14 is the load-bearing check — a real
indexing/mapping bug would corrupt gradients O(1), not reshuffle rounding in
ill-conditioned tail bins.

Rerun with the C_ν sub-dedup (151×, verified bit-exact): ALL OK; the σ_flat
diffs are bit-identical to the pre-sub-dedup run (forward is bit-identical
by construction), λ_ν-gradient diffs shift within 1e-14 (the gather's
backward segment-sums cotangents in a different order — rounding only).
