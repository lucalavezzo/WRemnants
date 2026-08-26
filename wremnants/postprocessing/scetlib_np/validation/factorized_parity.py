"""Parity check: legacy (Nbins, Nbt) vs factorized bT-grid reconstruction.

Compares :func:`scetlib_np.btgrid_tf.reconstruct_batch_tf` (legacy layout)
against :func:`scetlib_np.btgrid_tf.reconstruct_batch_factorized_tf`
(deduplicated rows + unique-qT J0 kernel + Simpson-as-matmul; the default
ParamModel path since the GPU-memory restructuring — see
scetlib_np/FACTORIZED_RECON.md) on the real combined btgrid.

The two paths evaluate the SAME integrand with the SAME Simpson weights on
the SAME sampling; they may differ only in floating-point multiplication
grouping and summation order (elementwise-then-reduce_sum vs matmul), i.e.
at the few-1e-15 relative level. This script asserts that, for σ values AND
for dσ/dλ (every unfrozen λ), at λ_central and at an off-central point that
exercises the per-bin (delta_lambda2 != 0) F_eff branch.

Run on CPU inside the wmassdev container (the legacy path needs ~45 GB host
RAM for its (Nbins, Nbt) tensors):

    singularity exec --cleanenv $IMG bash -c \
      "source setup.sh && python3 wremnants/postprocessing/scetlib_np/validation/factorized_parity.py"

Options:
    --btgrid-dir DIR   override the btgrid directory
    --tolerance X      max allowed |rel diff| (default 1e-12, generous vs
                       the expected ~1e-15; catches real bugs, not noise)
"""

import argparse
import time

import numpy as np
import tensorflow as tf

from wremnants.postprocessing.scetlib_np import btgrid_cache
from wremnants.postprocessing.scetlib_np import btgrid_integrate as fz_int
from wremnants.postprocessing.scetlib_np import btgrid_tf as fz_tf

# λ test points: central-ish values (tanh_6 model, as in production), plus an
# off-central point with delta_lambda2 != 0 to exercise the unique-Y F_eff
# gather branch in BOTH implementations.
LAMBDA_POINTS = {
    "central": dict(
        eff=dict(lambda_inf=1.0, lambda2=0.4, lambda4=0.4, lambda6=0.4, delta_lambda2=0.0),
        gnu=dict(lambda2_nu=0.15, lambda4_nu=0.0, lambda_inf_nu=2.0),
    ),
    "off_central_dl2": dict(
        eff=dict(lambda_inf=1.0, lambda2=0.6, lambda4=0.55, lambda6=0.4, delta_lambda2=0.1),
        gnu=dict(lambda2_nu=0.18, lambda4_nu=0.0, lambda_inf_nu=2.0),
    ),
}
NP_MODEL = "tanh_6"
NP_MODEL_NU = "tanh_6"
# λ that float in production fits (gradients checked for these).
GRAD_PARAMS = ["lambda2", "lambda4", "delta_lambda2", "lambda2_nu"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--btgrid-dir", default=None)
    ap.add_argument("--tolerance", type=float, default=1e-12)
    args = ap.parse_args()

    btgrid_dir = args.btgrid_dir
    if btgrid_dir is None:
        from wremnants.postprocessing.scetlib_np.param_model import _default_btgrid_dir

        btgrid_dir = _default_btgrid_dir()

    grid = btgrid_cache.load(btgrid_dir)
    idx_map = fz_int.dense_index_map(grid["bins"])
    bins = grid["bins"]
    qT_pb = np.array([b[2] for b in bins])
    Y_pb = np.array([b[1] for b in bins])
    Y_u_np, Y_inv_np = np.unique(Y_pb, return_inverse=True)
    Y_inv_np = Y_inv_np.reshape(-1).astype(np.int32)
    w_simpson = fz_tf.simpson_weights(grid["bT"])

    # ---- shared constants
    bT = tf.constant(grid["bT"], dtype=fz_tf.DTYPE)
    b_bar = tf.constant(grid["b_bar"], dtype=fz_tf.DTYPE)
    Y_unique = tf.constant(Y_u_np, dtype=fz_tf.DTYPE)

    # ---- legacy constants ((Nbins, Nbt) — host RAM heavy)
    t0 = time.time()
    qT_per_bin = tf.constant(qT_pb, dtype=fz_tf.DTYPE)
    Y_per_bin = tf.constant(Y_pb, dtype=fz_tf.DTYPE)
    I_pert = tf.constant(grid["I_pert"][0], dtype=fz_tf.DTYPE)
    C_nu = tf.constant(grid["C_nu"][0], dtype=fz_tf.DTYPE)
    kernel = fz_tf.build_bT_J0_kernel(qT_per_bin, bT)
    Y_inverse_idx = tf.constant(Y_inv_np, dtype=tf.int32)
    print(f"[parity] legacy constants built in {time.time()-t0:.0f}s", flush=True)

    # ---- factorized constants (mirrors param_model.__init__)
    t0 = time.time()
    dd = fz_tf.dedup_grid_rows(grid["I_pert"][0], grid["C_nu"][0], Y_inv_np)
    I_pert_u = tf.constant(dd["I_u"], dtype=fz_tf.DTYPE)
    C_nu_uu = tf.constant(dd["C_uu"], dtype=fz_tf.DTYPE)
    c_of_u = tf.constant(dd["c_of_u"], dtype=tf.int32)
    feff_idx_u = tf.constant(dd["feff_idx_u"], dtype=tf.int32)
    qT_idx = np.searchsorted(idx_map["qT_unique"], qT_pb)
    assert np.array_equal(idx_map["qT_unique"][qT_idx], qT_pb)
    gather_idx = tf.constant(
        np.stack([dd["row_uid"].astype(np.int64), qT_idx], axis=1), dtype=tf.int32
    )
    qT_u = tf.constant(idx_map["qT_unique"], dtype=fz_tf.DTYPE)
    KwqT = (
        qT_u[:, tf.newaxis]
        * fz_tf.build_bT_J0_kernel(qT_u, bT)
        * tf.constant(w_simpson, dtype=fz_tf.DTYPE)[tf.newaxis, :]
    )
    print(f"[parity] factorized constants built in {time.time()-t0:.0f}s", flush=True)

    def run_legacy(eff_vars, gnu_vars):
        return fz_tf.reconstruct_batch_tf(
            qT_per_bin=qT_per_bin,
            bT=bT,
            I_pert=I_pert,
            C_nu=C_nu,
            b_bar=b_bar,
            Y_per_bin=Y_per_bin,
            eff_params=eff_vars,
            gnu_params=gnu_vars,
            np_model=NP_MODEL,
            np_model_nu=NP_MODEL_NU,
            bT_J0_kernel=kernel,
            bT_simpson_weights=w_simpson,
            Y_unique=Y_unique,
            Y_inverse_idx=Y_inverse_idx,
        )

    def run_factorized(eff_vars, gnu_vars):
        return fz_tf.reconstruct_batch_factorized_tf(
            b_bar=b_bar,
            I_pert_u=I_pert_u,
            C_nu_uu=C_nu_uu,
            c_of_u=c_of_u,
            eff_params=eff_vars,
            gnu_params=gnu_vars,
            np_model=NP_MODEL,
            np_model_nu=NP_MODEL_NU,
            KwqT=KwqT,
            gather_idx=gather_idx,
            Y_unique=Y_unique,
            feff_idx_u=feff_idx_u,
        )

    def check_np(a, b, tol):
        """rtol+atol comparison. A bin fails only if it disagrees relatively
        AND the absolute difference is significant against the array's own
        peak scale. The atol floor is required because the high-qT Hankel
        sums are cancellation-dominated (Simpson condition number up to
        ~6e6 measured on the fineall grid): summation-order noise there is
        κ·ε ≈ 1e-10 *relative to the tiny bin value*, i.e. ≲1e-16 of the
        spectrum peak — below what fp64 determines in EITHER path. See
        FACTORIZED_RECON.md ("Parity results")."""
        a = np.atleast_1d(np.asarray(a))
        b = np.atleast_1d(np.asarray(b))
        diff = np.abs(a - b)
        scale = np.maximum(np.abs(a), np.abs(b))
        peak = scale.max()
        mask = scale > 0
        max_rel = np.max(diff[mask] / scale[mask]) if mask.any() else 0.0
        max_vs_peak = diff.max() / peak if peak > 0 else 0.0
        failed = bool(np.any((diff > tol * scale) & (diff > tol * peak)))
        return max_rel, max_vs_peak, failed

    def _report(tag, what, a, b, tol):
        max_rel, max_vs_peak, failed = check_np(a, b, tol)
        status = "FAIL" if failed else "OK"
        print(
            f"[parity] {tag:18s} {what:16s} max rel = {max_rel:.3e}  "
            f"max |diff|/peak = {max_vs_peak:.3e}  {status}",
            flush=True,
        )
        return failed

    n_fail = 0

    # Static fast path (python-float λ, delta_lambda2 == 0): value-only check.
    pt = LAMBDA_POINTS["central"]
    s_legacy = run_legacy(dict(pt["eff"]), dict(pt["gnu"]))
    s_fact = run_factorized(dict(pt["eff"]), dict(pt["gnu"]))
    n_fail += _report(
        "central_staticdl2", "sigma_flat", s_legacy.numpy(), s_fact.numpy(), args.tolerance
    )

    for tag, pt in LAMBDA_POINTS.items():
        # tf.Variables so a single tape pass gives every λ-gradient
        eff_v = {k: tf.Variable(v, dtype=fz_tf.DTYPE, name=k) for k, v in pt["eff"].items()}
        gnu_v = {k: tf.Variable(v, dtype=fz_tf.DTYPE, name=k) for k, v in pt["gnu"].items()}
        grad_vars = [v for d in (eff_v, gnu_v) for k, v in d.items() if k in GRAD_PARAMS]

        with tf.GradientTape() as tape_l:
            s_legacy = run_legacy(eff_v, gnu_v)
            sum_l = tf.reduce_sum(s_legacy)  # scalar proxy: d(Σσ)/dλ
        g_legacy = tape_l.gradient(sum_l, grad_vars)

        with tf.GradientTape() as tape_f:
            s_fact = run_factorized(eff_v, gnu_v)
            sum_f = tf.reduce_sum(s_fact)
        g_fact = tape_f.gradient(sum_f, grad_vars)

        n_fail += _report(tag, "sigma_flat", s_legacy.numpy(), s_fact.numpy(), args.tolerance)
        for v, gl, gf in zip(grad_vars, g_legacy, g_fact):
            n_fail += _report(
                tag,
                f"d/d{v.name.split(':')[0]}",
                gl.numpy(),
                gf.numpy(),
                args.tolerance,
            )

    if n_fail:
        raise SystemExit(f"[parity] {n_fail} comparison(s) FAILED")
    print("[parity] ALL OK", flush=True)


if __name__ == "__main__":
    main()
