"""Export the SCETlib bt-grid ParamModel resummed σ onto a SCETlib run's
(Q, Y, qT) binning, as a pickle that ``plot_precision_compare.py`` can overlay.

This finishes the SCETlib precision-triangulation: it turns the *production*
``SCETlibNPParamModel`` reconstruction (the bt-grid Hankel + arctan-Q² Simpson
integration at λ_central — i.e. the on-the-fly σ_gen path used in the fit) into
a curve directly comparable to the spectrum-mode SCETlib runs that differ only
in ``target_precision_rel``. Overlaying it answers the question "does the
bt-grid + on-the-fly reconstruction throw away precision relative to running
SCETlib directly at high precision?".

It uses ONLY committed WRemnants code for the param-model side:

  * ``SCETlibNPParamModel.sigma_YqT_central`` — the native, **resum-only**,
    Q-integrated σ on the signed-(Y, qT) bt-grid grid, BEFORE the |Y|-fold and
    gen rebin. This is the object the precision runs (which are
    ``calculation_piece = sing``, i.e. resummed singular only) are directly
    comparable to. We deliberately do NOT use ``sigma_gen_central``: that is the
    *matched* spectrum (resum ⊕ DYTurbo nonsingular) rebinned to the coarse gen
    grid, which is neither the same physical object as the singular runs nor on
    a binning fine enough to read off the precision differences.

  * ``btgrid_integrate.rebin_weights`` + ``rebin_axis_tf`` — the same committed
    Simpson rebin the model uses internally to reach the gen grid, here pointed
    at the SCETlib run's own Y/qT bin edges so the exported curve lands on
    EXACTLY those bins. The bt-grid grid is the union of edges+centres of that
    experimental grid, so each target bin gets a 3-point (4th-order) Simpson
    integration — the documented <0.05 % shape residual.

The output pickle matches the schema ``plot_precision_compare.py``'s ``load_run``
expects: ``{"hist": <hist.Hist over (Q, Y, qT, vars)>, "config", "meta_data"}``,
with a single Q bin ``[Q_lo, Q_hi]`` (``sigma_YqT_central`` is already
Q-integrated over that window) and a one-entry ``"paramodel"`` vars axis. The
stored values are BIN-INTEGRATED σ (matching the SCETlib runs), so the plot
script's ``/width`` recovers the differential and the ratio panel reads ~1.

Inputs live under /work and /ceph — run inside a container that binds both:

    export APPTAINER_BIND="/scratch,/cvmfs,/work,/ceph,/home"
    singularity run <wmassdevrolling> bash -c \\
      "source main/WRemnants/setup.sh; \\
       python3 -m wremnants.postprocessing.scetlib_np.validation.export_spectrum \\
         --out /work/submit/lavezzo/alphaS/scetlib-cms-newnp-lambda4fix/prod/scetlib_run/paramodel_btgrid_resum"

Then, in the scetlib_run env (no container/WRemnants needed):

    ./plot_precision_compare.py \\
        com13_ct18z_newnps_n3+0ll_lattice_lambda4bugfix_franksvals_fine \\
        com13_ct18z_newnps_n3+0ll_lattice_lambda4bugfix_franksvals_fine_rel1em3 \\
        com13_ct18z_newnps_n3+0ll_lattice_lambda4bugfix_franksvals_fine_rel1em4 \\
        paramodel_btgrid_resum \\
        --ref com13_..._franksvals_fine --qmin 60 --qmax 120
"""

import argparse
import glob
import os
import pickle
import sys
import time

import numpy as np

from wremnants.postprocessing.scetlib_np.sigma_gen import _default_btgrid_dir

# ---- Defaults (FranksVals CT18Z N3p0LL setup; mirror the validation scripts) ----
# A spectrum-mode SCETlib run whose (Q, Y, qT) binning the curve is exported onto
# (the same fine grid the precision runs share). Only its axis edges + NP config
# are read; the variation values are not used.
REF_RUN = (
    "/work/submit/lavezzo/alphaS/scetlib-cms-newnp-lambda4fix/prod/scetlib_run/"
    "com13_ct18z_newnps_n3+0ll_lattice_lambda4bugfix_franksvals_fine"
)
SIGNAL_PROC = "Zmumu"
Q_LO, Q_HI = 60.0, 120.0
NOMINAL_VARS = ("pdf0", "central", "nominal")


def _resolve_pkl(path):
    """A run dir or pkl -> a single pkl path (prefer *_combined.pkl)."""
    if os.path.isfile(path):
        return path
    if not os.path.isdir(path):
        raise FileNotFoundError(path)
    cand = sorted(glob.glob(os.path.join(path, "*_combined.pkl")))
    cand += sorted(glob.glob(os.path.join(path, "*.pkl")))
    if not cand:
        raise FileNotFoundError(f"no .pkl in {path}")
    return cand[0]


def read_ref_binning(ref_run):
    """Read (Q_edges, Y_edges, qT_edges, ref_total_sigma) from a SCETlib run.

    ``ref_total_sigma`` is the central-variation σ summed over the Z Q-window
    (a coarse normalisation sanity check vs the exported curve); ``None`` if it
    can't be computed.
    """
    pkl = _resolve_pkl(ref_run)
    with open(pkl, "rb") as fh:
        h = pickle.load(fh)["hist"]
    Q_edges = np.asarray(h.axes["Q"].edges, dtype=np.float64)
    Y_edges = np.asarray(h.axes["Y"].edges, dtype=np.float64)
    qT_edges = np.asarray(h.axes["qT"].edges, dtype=np.float64)

    ref_total = None
    try:
        names = list(h.axes["vars"])
        sel = next((v for v in NOMINAL_VARS if v in names), names[0])
        hv = h[{"vars": sel}]
        qi = int(np.argmin(np.abs(Q_edges - Q_LO)))
        qj = int(np.argmin(np.abs(Q_edges - Q_HI)))
        hv = hv[{"Q": slice(qi, qj, sum)}]
        ref_total = float(hv.values().sum())
    except Exception as exc:  # noqa: BLE001 - sanity print only
        print(f"[ref] could not compute reference total σ: {exc}", flush=True)
    return Q_edges, Y_edges, qT_edges, ref_total, pkl


def rebin_sigma_to_edges(model, sigma_YqT, Y_edges, qT_edges, tol):
    """Simpson-rebin a native (Y, qT) σ array onto (Y_edges, qT_edges).

    Uses only committed WRemnants primitives. The input ``sigma_YqT`` (e.g.
    ``model.sigma_YqT_central`` or ``model._sigma_YqT_native_at(eff, gnu)`` at an
    arbitrary λ) is already Q-integrated over [Q_lo, Q_hi], so the result is the
    bin-integrated σ over (Q-window × Y-bin × qT-bin). Returns an (NY, NqT) array.
    """
    from wremnants.postprocessing.scetlib_np import btgrid_integrate as fz_int

    Wy = fz_int.rebin_weights(model.Y_unique, Y_edges, name="Y->Y", tol=tol)
    Wq = fz_int.rebin_weights(model.qT_unique, qT_edges, name="qT->qT", tol=tol)
    sig = fz_int.rebin_axis_tf(sigma_YqT, axis=0, weights=Wy)
    sig = fz_int.rebin_axis_tf(sig, axis=1, weights=Wq)
    return np.asarray(sig.numpy(), dtype=np.float64)


def _apply_lambda_overrides(eff, gnu, spec):
    """Apply a ``name=value,...`` string onto copies of the eff/gnu λ dicts.

    Routed by :data:`params.EFF_PARAMS` / :data:`params.GNU_PARAMS` membership
    (NOT by presence in the card dicts), so a λ the card form does not carry —
    e.g. ``lambda6`` for a tanh_6 ``--np-model`` evaluation — can be added."""
    from wremnants.postprocessing.scetlib_np.params import EFF_PARAMS, GNU_PARAMS

    eff, gnu = dict(eff), dict(gnu)
    if spec:
        for tok in spec.split(","):
            if not tok.strip():
                continue
            k, v = tok.split("=", 1)
            k, v = k.strip(), float(v.strip())
            if k in EFF_PARAMS:
                eff[k] = v
            elif k in GNU_PARAMS:
                gnu[k] = v
            else:
                raise KeyError(f"--lambdas: unknown NP parameter {k!r} "
                               f"(known: {list(EFF_PARAMS) + list(GNU_PARAMS)})")
    return eff, gnu


def build_hist(sig, Q_lo, Q_hi, Y_edges, qT_edges):
    """Pack the (NY, NqT) bin-integrated σ into a (Q, Y, qT, vars) hist.Hist.

    Single Q bin = [Q_lo, Q_hi] makes explicit that the curve only covers the Z
    window; a one-entry ``paramodel`` vars axis so load_run finds a slot.
    """
    import hist

    qax = hist.axis.Variable([Q_lo, Q_hi], name="Q", flow=False)
    yax = hist.axis.Variable(Y_edges, name="Y", flow=False)
    qtax = hist.axis.Variable(qT_edges, name="qT", flow=False)
    vax = hist.axis.StrCategory(["paramodel"], name="vars")
    h = hist.Hist(qax, yax, qtax, vax, storage=hist.storage.Double())
    h.view(flow=False)[...] = sig[np.newaxis, :, :, np.newaxis]
    return h


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--datacard", required=True, help="fit-input hdf5 (λ_central + R auxiliary)")
    p.add_argument("--btgrid", default=_default_btgrid_dir(), help="SCETlib bT-grid directory")
    p.add_argument("--ref-run", default=REF_RUN,
                   help="spectrum-mode SCETlib run (dir or pkl) whose (Q,Y,qT) edges "
                        "the curve is exported onto")
    p.add_argument("--signal-proc", default=SIGNAL_PROC)
    p.add_argument("--q-lo", type=float, default=Q_LO)
    p.add_argument("--q-hi", type=float, default=Q_HI)
    p.add_argument("--tol", type=float, default=1e-6,
                   help="edge-matching tolerance for the Simpson rebin (default 1e-6)")
    p.add_argument("--out", required=True,
                   help="output RUN DIR; the pkl is written inside it and its basename "
                        "becomes the curve label in plot_precision_compare.py")
    p.add_argument("--name", default=None,
                   help="pkl basename (default: <out dir name>)")
    p.add_argument("--cache", default=None,
                   help="npz cache of the rebinned σ (default <out>/rebinned_sig.npz); "
                        "reused unless --rebuild or the ref binning changed")
    p.add_argument("--rebuild", action="store_true",
                   help="force model construction even if the cache is present")
    p.add_argument("--lambdas", default=None,
                   help="comma-separated name=value overrides of the NP lambdas vs "
                        "λ_central, e.g. 'lambda2=0.4793,delta_lambda2=-0.0119,"
                        "lambda4=-0.0207'. Unspecified lambdas keep their λ_central "
                        "(= param_model prior-mean) values. When set, the cache is "
                        "bypassed (always rebuilt at these lambdas).")
    p.add_argument("--np-model", default=None,
                   help="F_eff form for the EVALUATION (default: the card form; "
                        "e.g. tanh_6 to export a numerator-form-override point — "
                        "its extra λ, e.g. lambda6, must come via --lambdas)")
    p.add_argument("--np-model-nu", default=None,
                   help="γ_ν^NP form for the EVALUATION (default: the card form)")
    args = p.parse_args(argv)

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)
    name = args.name or os.path.basename(os.path.normpath(out_dir))
    cache = args.cache or os.path.join(out_dir, "rebinned_sig.npz")

    print("=" * 74)
    print("EXPORT ParamModel resummed σ (bt-grid, on-the-fly) onto SCETlib run binning")
    print("=" * 74)
    print(f"  datacard : {args.datacard}")
    print(f"  btgrid   : {args.btgrid}")
    print(f"  ref-run  : {args.ref_run}")
    print(f"  out      : {out_dir}  (label '{name}')")

    Q_edges, Y_edges, qT_edges, ref_total, ref_pkl = read_ref_binning(args.ref_run)
    print(f"  ref pkl  : {ref_pkl}")
    print(f"  ref bins : NQ={Q_edges.size-1} NY={Y_edges.size-1} NqT={qT_edges.size-1}")

    eff_central = gnu_central = None
    sig = None
    if (os.path.exists(cache) and not args.rebuild and not args.lambdas
            and not args.np_model and not args.np_model_nu):
        z = np.load(cache, allow_pickle=True)
        if (z["Y_edges"].shape == Y_edges.shape and np.allclose(z["Y_edges"], Y_edges)
                and z["qT_edges"].shape == qT_edges.shape
                and np.allclose(z["qT_edges"], qT_edges)
                and float(z["q_lo"]) == args.q_lo and float(z["q_hi"]) == args.q_hi):
            sig = z["sig"]
            eff_central = z["eff_central"].item() if "eff_central" in z else None
            gnu_central = z["gnu_central"].item() if "gnu_central" in z else None
            print(f"\n[cache] reusing rebinned σ from {cache}  (shape {sig.shape})")
        else:
            print(f"\n[cache] {cache} present but binning/Q-window changed -> rebuilding")

    if sig is None:
        from rabbit.inputdata import FitInputData
        from wremnants.postprocessing.scetlib_np.param_model import SCETlibNPParamModel

        print("\n[model] constructing SCETlibNPParamModel (runs the bt-grid integral "
              "at λ_central) …", flush=True)
        t0 = time.time()
        indata = FitInputData(args.datacard)
        model = SCETlibNPParamModel(
            indata,
            btgrid_dir=args.btgrid,
            signal_proc=args.signal_proc,
            Q_lo=args.q_lo,
            Q_hi=args.q_hi,
        )
        print(f"  constructed in {time.time()-t0:.1f}s; "
              f"native σ(Y,qT) shape {tuple(model.sigma_YqT_central.shape)}")
        print(f"  λ_central F_eff   : {dict(model.eff_central)}")
        print(f"  λ_central γ_ν^NP  : {dict(model.gnu_central)}")

        # λ to evaluate at: λ_central, optionally overridden by --lambdas. The
        # unspecified ones keep λ_central (= param_model prior means).
        eff_central, gnu_central = _apply_lambda_overrides(
            model.eff_central, model.gnu_central, args.lambdas
        )
        if args.lambdas or args.np_model or args.np_model_nu:
            if args.np_model or args.np_model_nu:
                print(f"  EVAL forms        : np_model="
                      f"{args.np_model or model.core.np_model}, np_model_nu="
                      f"{args.np_model_nu or model.core.np_model_nu}")
            print(f"  λ OVERRIDE F_eff  : {eff_central}")
            print(f"  λ OVERRIDE γ_ν^NP : {gnu_central}")
            sigma_YqT = model._sigma_YqT_native_at(
                eff_central,
                gnu_central,
                np_model=args.np_model,
                np_model_nu=args.np_model_nu,
            )
        else:
            sigma_YqT = model.sigma_YqT_central

        print("\n[rebin] Simpson-rebinning native σ(Y,qT) onto the ref (Y, qT) edges "
              "via btgrid_integrate.rebin_weights …", flush=True)
        sig = rebin_sigma_to_edges(model, sigma_YqT, Y_edges, qT_edges, args.tol)
        np.savez(
            cache,
            sig=sig,
            Y_edges=Y_edges,
            qT_edges=qT_edges,
            q_lo=args.q_lo,
            q_hi=args.q_hi,
            eff_central=np.array(eff_central, dtype=object),
            gnu_central=np.array(gnu_central, dtype=object),
        )
        print(f"  wrote cache {cache}")

    print(f"\n[norm] Σ exported σ (Q∈[{args.q_lo:g},{args.q_hi:g}]) = {sig.sum():.6g}")
    if ref_total is not None:
        print(f"       Σ ref-run central σ (same Q-window)    = {ref_total:.6g}")
        print(f"       ratio exported/ref                     = {sig.sum()/ref_total:.4f}  "
              "(coarse total; the per-bin ratio is what the plot shows)")

    h = build_hist(sig, args.q_lo, args.q_hi, Y_edges, qT_edges)
    out_pkl = os.path.join(out_dir, f"{name}_combined.pkl")
    payload = {
        "hist": h,
        "config": {
            "source": "SCETlibNPParamModel.sigma_YqT_central (resum-only, bt-grid)",
            "datacard": args.datacard,
            "btgrid": args.btgrid,
            "ref_run": ref_pkl,
            "q_window": [args.q_lo, args.q_hi],
            "eff_central": eff_central,
            "gnu_central": gnu_central,
        },
        "meta_data": {
            "made_by": "export_spectrum.py",
            "note": "BIN-INTEGRATED resum-only σ on the ref (Q,Y,qT) binning; "
                    "single Q bin = [q_lo,q_hi]; valid only for that Q-window.",
        },
    }
    with open(out_pkl, "wb") as fh:
        pickle.dump(payload, fh)
    print(f"\n[write] {out_pkl}")
    print("Done. Overlay it in plot_precision_compare.py with --qmin "
          f"{args.q_lo:g} --qmax {args.q_hi:g}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
