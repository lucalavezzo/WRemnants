"""Smoke test for SCETlibNPParamModel gen_level=1 (no response-matrix fold).

Verifies the gen-level σUL mode end to end:
  1. gen_level construction reads the gen binning from the (single) fit channel
     and needs NO scetlib_np auxiliary / R / N_gen: R is None,
     sigma_reco_central absent, sigma_gen_central_flat present, gen_shape ==
     the channel binning.
  2. compute() at λ_central → ratio is exactly 1 in the signal column (to machine
     precision) and exactly 1 in every other process column; shape (N_gen, nproc).
  3. compute() with one λ perturbed → finite, non-trivial per-gen-bin ratio, and
     it equals the raw σ_gen(λ)/σ_gen(λ_central) (i.e. NO reco fold was applied).

The gen binning (ptVGen, absYVGen edges) is read from the unfolding histmaker
output; λ_central is passed explicitly (the shipped example) since older cards
predate the λ_central-into-metadata propagation. The model strings only select
the F_eff/γ_ν functional form — the btgrid caches the NP-off integrands, so the
exact central values do not affect this code-path test (the ratio is taken vs
the same λ_central, so it is 1 at central by construction).

Run inside the wmass singularity with setup.sh sourced, e.g.:

    singularity exec --cleanenv --bind /scratch,/work,/home,/ceph,/cvmfs \\
      /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/bendavid/cmswmassdocker/wmassdevrolling:latest \\
      bash -c 'source setup.sh && \\
        python3 -m wremnants.postprocessing.scetlib_np.validation.gen_level_smoke'
"""

import sys
import time
import types

import numpy as np

from wremnants.postprocessing.scetlib_np.sigma_gen import _default_btgrid_dir
import tensorflow as tf

from rabbit.inputdata import FitInputData

from wremnants.postprocessing.scetlib_np import response_matrix as fz_R
from wremnants.postprocessing.scetlib_np.param_model import SCETlibNPParamModel

# Any Z-dilepton fit input works — gen_level only reads procs/dtype/systs and the
# (swapped-in) gen channel from it; the scetlib_np auxiliary is NOT needed.
FIT_HDF5 = (
    "/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/260528_debug_SCETlibPOIModel/"
    "ZMassDilepton_ptll_yll_cosThetaStarll_quantile_phiStarll_quantile_excludeSCETlibNP/"
    "ZMassDilepton.hdf5"
)
# Unfolding histmaker output — used here only to read the real (ptVGen, absYVGen)
# gen-bin edges for the synthetic gen channel.
UNFOLD_HDF5 = (
    "/ceph/submit/data/group/cms/store/user/lavezzo/alphaS/"
    "260411_histmaker_dilepton_unfolding/mz_dilepton.hdf5"
)
# Explicit λ_central: the real production central (read from the metadata of a
# recent 260623 histmaker output), passed here because the older FitInputData
# card predates the λ_central-into-metadata propagation. Using the physical
# central keeps σ_gen(λ_central) > 0 on the gen binning (the positivity guard).
LAMBDA_CENTRAL = {
    "eff_params": {
        "np_model": "tanh_2",
        "lambda_inf": 1.0,
        "lambda2": 0.4,
        "lambda4": 0.4,
        "lambda6": 0.0,
        "delta_lambda2": 0.0,
    },
    "gnu_params": {
        "np_model_nu": "tanh_2",
        "lambda_inf_nu": 2.0,
        "lambda2_nu": 0.15,
        "lambda4_nu": 0.0,
    },
}


def main():
    print("Loading FitInputData …", flush=True)
    t0 = time.time()
    indata = FitInputData(FIT_HDF5)
    procs = [p.decode() if isinstance(p, bytes) else str(p) for p in indata.procs]
    signal_proc = "Zmumu" if "Zmumu" in procs else procs[0]
    print(f"  loaded in {time.time()-t0:.1f}s; nproc={indata.nproc}; "
          f"procs={procs}; signal={signal_proc}", flush=True)

    # Real gen-bin edges (ptVGen, absYVGen) from the unfolding output.
    print("\nReading gen-bin edges from the unfolding output …", flush=True)
    gen_axes_meta = fz_R.load_R(UNFOLD_HDF5)["gen_axes"]  # [(name, edges), ...]
    print(f"  gen axes: {[ (n, len(e)-1) for n, e in gen_axes_meta ]}", flush=True)

    # Build a synthetic single 2-axis GEN channel from those edges; gen_level
    # reads its binning from the (single, non-masked) fit channel via
    # _fit_reco_axes, which only needs objects with .name / .edges.
    gen_axis_objs = [
        types.SimpleNamespace(name=name, edges=np.asarray(edges))
        for name, edges in gen_axes_meta
    ]
    indata.channel_info = {"gen_sigmaUL": {"axes": gen_axis_objs}}  # one, non-masked
    n_gen_expected = int(np.prod([len(e) - 1 for _, e in gen_axes_meta]))

    # === [1/3] gen_level construction: no R / N_gen needed =====================
    print("\n[1/3] gen_level=1 construction (no response fold) …", flush=True)
    t0 = time.time()
    gen = SCETlibNPParamModel(
        indata,
        btgrid_dir=_default_btgrid_dir(),
        signal_proc=signal_proc,
        lambda_central=LAMBDA_CENTRAL,
        gen_level=True,
    )
    print(f"  constructed in {time.time()-t0:.1f}s", flush=True)
    ok = True
    if gen.gen_level is not True:
        print("  FAIL: gen_level flag not set"); ok = False
    if gen.R is not None:
        print("  FAIL: R should be None in gen_level"); ok = False
    if getattr(gen, "sigma_reco_central", None) is not None:
        print("  FAIL: sigma_reco_central should be absent in gen_level"); ok = False
    if not hasattr(gen, "sigma_gen_central_flat"):
        print("  FAIL: sigma_gen_central_flat missing"); ok = False
    if int(np.prod(gen.gen_shape)) != n_gen_expected:
        print(f"  FAIL: gen_shape {gen.gen_shape} != expected {n_gen_expected}"); ok = False
    if not ok:
        return 1
    print(
        f"  OK: R is None, sigma_gen_central_flat set "
        f"({int(gen.sigma_gen_central_flat.shape[0])} gen bins), gen_shape={gen.gen_shape}",
        flush=True,
    )

    # === [2/3] compute() at λ_central → ratio == 1 ============================
    print("\n[2/3] compute() at λ_central → 1 in signal col, 1 elsewhere", flush=True)
    rnorm_c = gen.compute(gen.xparamdefault).numpy()
    if rnorm_c.shape != (n_gen_expected, gen.nproc):
        print(f"  FAIL: shape {rnorm_c.shape} != {(n_gen_expected, gen.nproc)}")
        return 1
    sig = rnorm_c[:, gen.signal_proc_idx]
    others = np.delete(rnorm_c, gen.signal_proc_idx, axis=1)
    sig_err = float(np.max(np.abs(sig - 1.0)))
    oth_err = float(np.max(np.abs(others - 1.0)))
    print(f"  signal col max|ratio-1|={sig_err:.3g}; other cols max|val-1|={oth_err:.3g}",
          flush=True)
    if sig_err > 1e-10 or oth_err > 0:
        print("  FAIL"); return 1
    print("  OK", flush=True)

    # === [3/3] perturbed λ → finite, non-trivial, and == raw gen ratio ========
    print("\n[3/3] perturb lambda2_nu by +10% → non-trivial gen-level ratio", flush=True)
    pert = gen.xparamdefault.numpy().copy()
    idx = gen._param_order.index("lambda2_nu")
    pert[idx] *= 1.10
    rnorm_p = gen.compute(tf.constant(pert, dtype=indata.dtype)).numpy()
    sig_p = rnorm_p[:, gen.signal_proc_idx]
    if not np.all(np.isfinite(sig_p)):
        print("  FAIL: non-finite ratio"); return 1
    if float(np.max(np.abs(sig_p - 1.0))) < 1e-6:
        print("  FAIL: ratio is trivially 1 everywhere (λ had no effect)"); return 1
    # Cross-check: gen_level ratio must be the RAW gen ratio (no R fold).
    eff_p, gnu_p = gen._unpack_params(tf.constant(pert, dtype=indata.dtype))
    sigma_gen_p = tf.reshape(gen._sigma_gen_at(eff_p, gnu_p), [-1])
    raw_gen_ratio = (sigma_gen_p / gen.sigma_gen_central_flat).numpy()
    mismatch = float(np.max(np.abs(sig_p - raw_gen_ratio)))
    print(f"  ratio range [{sig_p.min():.4f}, {sig_p.max():.4f}]; "
          f"max|signal - raw σ_gen ratio|={mismatch:.3g}", flush=True)
    if mismatch > 1e-9:
        print("  FAIL: signal column != raw gen ratio (an unexpected fold crept in)")
        return 1
    print("  OK", flush=True)

    print("\nALL CHECKS PASSED ✓", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
