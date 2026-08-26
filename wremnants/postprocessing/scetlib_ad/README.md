# `scetlib_ad` — a fully differentiable SCETlib prediction for rabbit

A rabbit `ParamModel` in which **every theory parameter SCETlib exposes is a
continuous fit parameter with exact derivatives**, rather than a discrete template
morph whose joint response with the others is an outer product:

| parameter group | count (Z) | notes |
|---|---|---|
| `alphaS` | 1 | physical units; with a PDF α_s member pair it is the PDF-consistent coupling, not α_s at fixed PDF |
| nonperturbative λ | 8 | Collins–Soper and TMD form factors, `tanh_2`/`tanh_6` |
| theory nuisance parameters | 10 | `gamma_cusp`, `gamma_mu_q`, `gamma_nu`, `s`, `h_qqV` + 5 beam-function TNPs |
| PDF eigenvector coefficients | `n_eig` | extra differentiable columns, exact at `c_e = 0, ±1` |

Which of these exist is a property of the **cache**, not of this code: the model
reads `gradient_param_names()` and registers what it finds. Only the
profile-scale parameters — `kappaFO`, `kappaf`, `muf`, the transition points —
are outside SCETlib's autodiff (they need d/dμ of the PDF convolution grids) and
still need template nuisances.

The physics comes from the SCETlib `autodiff-sigmaul` branch: `ScetlibCachedXsecTF`
(`scetlib-cms/py/scetlib_tf.py`) replays a prepared cache — compressed bin rules
for the resummed piece plus a frozen fixed-order grid for the nonsingular one —
and returns exact first and second derivatives from clad.

| file | role |
|---|---|
| `params.py` | SCETlib ↔ rabbit name map, prior σ, POI/POU defaults, impact groups |
| `xsec_backend.py` | `ScetlibADXsec` (configure + cache load + value/J/K) and `GenFold` (cache bins → the card's gen grid) |
| `response.py` | the reco fold: `P(b|g) = R_raw/N_gen`, the datacard response auxiliary, the positivity floor |
| `param_model.py` | `SCETlibADParamModel` — the rabbit adapter |

Scripts live in `scripts/rabbit/scetlib_ad/`. The BUILD itself is SCETlib's --
`scetlib-cms/examples/matched_ad/prepare_cache.py` holds the steps
(`plan_variations`, `build_prologue`, `build_variations`, `write_cache`,
`fork_member_build`) and `scetlib-cms/py/scetlib_cache.py` the cache file
format and both merges. The two scripts here are the rabbit side: the gen axes
off a card, the runcard, `--subset`, and the shard scheduling.


| script | role |
|---|---|
| `backend_check.py` | standalone cache sanity: anchor round trip, FD-checked Jacobian, Hessian symmetry, fold sum rule |
| `prepare_cache_for_card.py` | build a cache for a card's gen binning, or an explicit `--grid-json` |
| `build_cache_parallel.py` | split the BINS across processes (`--bin-groups`) and merge the shards (`--merge-bins`, `--merge-only`) |
| `make_debug_card.py` | a self-contained gen-level card built from a cache, for closure tests |
| `compare_to_scetlib_run.py` | validate the resummed piece against a native SCETlib production run |
| `conf/Z_CT18Z_N3p0LL_FranksVals.conf` | runcard reproducing the current analysis central (see below) |

## Running

Inside the WRemnants singularity, `source setup.sh` then the SCETlib one:

```bash
source $WREM_BASE/scetlib-cms/setup.sh     # PYTHONPATH, LD_LIBRARY_PATH, ulimit -s

# ONE job: minimize, postfit Hessian, impacts.
rabbit_fit.py <gen_card>.hdf5 -v 3 \
  --paramModel wremnants.postprocessing.scetlib_ad.SCETlibADParamModel \
      cache=<cache>.npz conf=<runcard>.conf gen_level=1 threads=32 \
  --jitCompile off -t 0 --doImpacts -o <dir>
```

`--jitCompile off` is **mandatory** — the model reaches SCETlib through
`tf.py_function`, which XLA cannot compile. The model refuses to construct
otherwise, with that message.

**One pass, not two.** The fit and the postfit covariance are a single job. The
model always includes the exact second-derivative term, so the composite Hessian
is exact and there is nothing to configure — no `--noHessian` fit followed by an
`--externalPostfit --noFit` pass. Measured on a 30-bin debug card, one pass
reproduces the two-pass numbers exactly (α_s = 0.1195 ± 0.00045, identical λ
uncertainties) with a much harder-converged EDM (6.1e-22 vs 3.5e-17).

A Gauss-Newton variant dropping that term was measured ~5× faster and, on Asimov,
numerically identical — but it is an approximation on real data, so it is
deliberately **not** offered. The exact Hessian costs ~1 s/bin of serial work,
scaling as `1 + P(P+1)/2` in the parameter count: 3–80 s per minimiser iteration
on 64 threads at the few-hundred-to-1200-bin gen binnings we fit on. Reintroduce a
switch only if a fit on the full correction grid ever needs it.

Two things keep that affordable: the `(value, J, K)` triple is cached on the
parameter vector, so one C++ Hessian is computed per *distinct* point rather than
per HVP (the minimiser takes many HVPs at fixed `x`); and the exact Hessian makes
the minimiser converge harder, so it needs fewer iterations.

## How the derivatives get into TensorFlow

`ScetlibCachedXsecTF` is an ordinary TF-differentiable function — its backward pass
is itself a `custom_gradient` whose own gradient contracts Hessian-vector
products — so nested `GradientTape`s work and TF drives every C++ call. The model
just calls it inside the graph, exactly as `examples/matched_ad/tf_gradients.py`
does. **There is no surrogate anywhere**: autodiff differentiates the real
prediction, and rabbit's postfit Hessian is the real Hessian.

One requirement that imposes, and it is easy to break by accident:

> Map rabbit's fit vector into SCETlib's layout with a **constant 0/1 matrix
> multiply**, never `tensor_scatter_nd_update`.

Some mapping is unavoidable — rabbit's vector holds only the fitted parameters,
POIs first, while SCETlib's holds every registered parameter in registry order. But
a scatter's backward pass contains a gather, whose gradient TF represents as
`tf.IndexedSlices`, and the bridge's second-order py_function payloads call
`.numpy()` on the incoming cotangent and fail on it — so **everything past first
order breaks**, while first order keeps working, which makes it a nasty way to
fail. Isolated: a nested-tape HVP works on a bare `Variable`, and fails with a
scatter in front even when the scatter covers the whole vector. A constant matmul
is bit-identical (the entries are exactly 0 and 1) and free at these sizes.

Fixing it upstream — densifying the incoming cotangent before `.numpy()` in
`_uhvp_py` — would remove the trap rather than leave us relying on a comment.

## Validating a cache

Two independent checks, both cheap relative to a fit:

- `compare_to_scetlib_run.py` — against a native SCETlib production run with
  `calculation_piece = sing`. That pkl **is** the resummed cross section,
  bin-integrated, so replaying only our cache's rules gives the same object: no
  matching, no fixed-order generator, no MC, no correction file in between. Any
  disagreement is therefore ours — runcard, quadrature, Q integration, rule
  compression. The reference's finer bins are summed onto ours, which is exact
  because both sides are bin-integrated, and the script refuses to run unless our
  edges really are a subset of the reference's.
- `backend_check.py` — the cache's own consistency, plus the fold sum rule.

**Use an analysis runcard, not `examples/matched_ad/matched.conf`.** The example
runcard is plain N3LL with SCETlib's default profile scales; the analysis is
N³⁺⁰LL — every TNP at `theta = 0` with `'level0'` — with `lambda = 0`, transition
points `[0.2, 0.6, 1.0]`, scale floors, `compensate_fo` and `collins_soper4`.
Measured on a 30-bin grid, that difference moves the λ response by 7–35% of the
response itself. `conf/Z_CT18Z_N3p0LL_FranksVals.conf` transcribes the config of
the production run the current analysis correction is built from, with the
provenance in comments.

Note the consequence: the analysis order is *defined* by the `[TNPs]` block, and a
non-`off` TNP scheme is exactly what registers a TNP as a gradient parameter. So
an analysis-faithful cache has 19 parameters, not 9. They are not fitted by
default — pass them in `fit_params`, and then `priors=1` is required (the model
refuses to float a TNP free, and refuses any parameter whose Jacobian column is
identically zero, which `resumTNP_b_qqDS` is for the Z).

Study logbook: `WRemnantsHelpers/studies/scetlib-ad-param-model/LOGBOOK.md`.
