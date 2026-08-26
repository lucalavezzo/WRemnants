"""Shared reference-loading + comparison library for the SCETlib-NP validation CLIs.

This is the Layer-1 library for the validation reorg: the per-reference loaders and
aligners that were previously duplicated-in-spirit across the standalone scripts now
live here, so the thin CLIs (``validate_agreement`` for card/histmaker, and
``sigma_gen_at_lambda`` for theory-corr) share ONE implementation.

Layering (do not invert):
  * Layer 0 — ``param_model_diagnostics`` : the numpy-only fit-side core (the in-fit
    ``run_reco_guard``, the pathology detectors, ``_shape_residual``, the card
    references). Stays put; the fit imports it. This module may import FROM it.
  * Layer 1 — THIS module : reference providers + comparison helpers.
  * Layer 2 — the CLIs.

Every reference here aligns onto the MODEL's grid and returns plain ndarrays / hists,
so the same ``_shape_residual`` / ``summarize`` compare the model to any reference.

Heavy deps (``hist``, ``wums``, ``h5py``, ``tensorflow``) are imported lazily inside
the functions that need them, so importing this module stays cheap.
"""

import numpy as np

# Shared fine->coarse rebin + plotting/array helpers (single implementation).
from wremnants.postprocessing.scetlib_np.validation_plots import (
    AXIS_LABELS,
    _merge_matrix,
    plot_ptll_ratio,
    project_inrange,
    summarize,
    tf_to_hist,
)
from wremnants.postprocessing.scetlib_np.params import (
    EFF_PARAMS,
    GNU_PARAMS,
)

# =============================================================================
# histmaker references (moved verbatim from validation.histmaker_validation)
# =============================================================================
NOMINAL_HIST = "nominal"
SIGNAL_SAMPLE = "Zmumu_2016PostVFP"
SIGNAL_PROC = "Zmumu"

# Off-central λ-response check (--variation): the histmaker stores, in this hist,
# the nominal reco yield reweighted by each NP variation of the scetlib_dyturbo
# correction, on a 'vars' axis. vars='pdf0' is the central (== 'nominal'); labels
# like 'lambda21.0', 'lambda2_nu0.25', 'lambda40.0', 'delta_lambda20.02' are
# single-λ variations. We compare Corr[var]/Corr[pdf0] (the histmaker's reco-level
# λ-reweighting) to the model's rnorm = σ_reco(λ)/σ_reco(λ_c) for the same λ —
# both are variation/central ratios, so the comparison needs no normalization and
# is non-circular (it tests the model's λ-sensitivity, exactly what the fit uses).
VARIATION_HIST = (
    "nominal_scetlib_dyturbo_LatticeNPLambda4Bugfix_FranksVals_CT18Z_N3p0LL_N2LO_Corr"
)

# ---- Gen-level cross-check: compares σ_gen(λ_central) (the bT integral, NO
# response matrix) directly to a gen-level histmaker at the same NP tune. This
# isolates the integral from R. The gen hist is rebinned onto the model's gen
# grid (the model's gen edges must be a sub-binning of the gen hist's).
# NOTE: must match the datacard's NP tune (λ_central). The datacard here is
# plain FranksVals (lambda6=0.0); use the matching gen file, NOT the
# "_Lambda60p4" (lambda6=0.4) variant.
GEN_HIST = "nominal_gen"          # (massVgen, absYVgen, ptVgen, chargeVgen, helicity)
GEN_SAMPLE = "Zmumu"
GEN_MASS_WINDOW = (60.0, 120.0)   # match the model's Q-integration window


def load_nominal(histmaker_path, sample_key, hist_name=NOMINAL_HIST):
    """Load the histmaker 'nominal' Hist for the signal sample."""
    import h5py
    from wums import ioutils as wums_io

    with h5py.File(histmaker_path, "r") as f:
        if sample_key not in f:
            raise KeyError(
                f"{histmaker_path}: no {sample_key!r} group. "
                f"top-level: {[k for k in f.keys()][:12]}"
            )
        sample = wums_io.pickle_load_h5py(f[sample_key])
        out = sample["output"]
        if hist_name not in out:
            raise KeyError(
                f"{sample_key}: no {hist_name!r} hist. available: "
                f"{list(out.keys())[:15]}"
            )
        proxy = out[hist_name]
        h = proxy.get() if hasattr(proxy, "get") else proxy
    return h


def align_nominal(h, reco_axes_meta, tol=1e-6):
    """Reorder + crop the nominal Hist onto the model's reco binning.

    Stays at the hist level: ``project(*names)`` reorders (and sums out any
    unlisted axes), then an integer-slice per axis handles the case where the
    histmaker axis is a *superset* of the fit axis (e.g. ptll has a trailing
    [44, 100] bin while the fit stops at 44).

    reco_axes_meta : list of (name, edges) in the model's reco-axis order
        (i.e. model._reco_axes_meta — the fit edges).

    Returns a hist.Hist whose axes match the model (same names/order/edges).
    """
    names = [n for n, _ in reco_axes_meta]
    have = [a.name for a in h.axes]
    missing = [n for n in names if n not in have]
    if missing:
        raise ValueError(f"nominal hist missing axes {missing}; has {have}")
    h = h.project(*names)  # reorder

    crop = {}
    for name, medges in reco_axes_meta:
        medges = np.asarray(medges, dtype=np.float64)
        hedges = np.asarray(h.axes[name].edges, dtype=np.float64)
        nb = medges.size - 1
        hits = np.where(np.isclose(hedges, medges[0], atol=tol))[0]
        if hits.size == 0:
            raise ValueError(
                f"axis {name}: model low edge {medges[0]} not found in hist "
                f"edges [{hedges[0]} .. {hedges[-1]}]"
            )
        i0 = int(hits[0])
        if i0 + nb + 1 > hedges.size or not np.allclose(
            hedges[i0 : i0 + nb + 1], medges, atol=tol
        ):
            raise ValueError(
                f"axis {name}: hist edges from index {i0} don't match model "
                f"edges. hist={hedges[i0 : i0 + nb + 1]} model={medges}"
            )
        if i0 != 0 or nb != h.axes[name].size:
            crop[name] = slice(i0, i0 + nb)  # integer (bin-index) slice
    return h[crop] if crop else h


def load_gen_hist(gen_path, sample_key, hist_name=GEN_HIST):
    """Load the gen-level Hist from a gen-distributions histmaker output."""
    import h5py
    from wums import ioutils as wums_io

    with h5py.File(gen_path, "r") as f:
        if sample_key not in f:
            raise KeyError(
                f"{gen_path}: no {sample_key!r} group. top-level: "
                f"{[k for k in f.keys()][:12]}"
            )
        out = wums_io.pickle_load_h5py(f[sample_key])["output"]
        if hist_name not in out:
            raise KeyError(
                f"{sample_key}: no {hist_name!r} hist. available: "
                f"{[k for k in out.keys() if 'gen' in k.lower()][:15]}"
            )
        proxy = out[hist_name]
        h = proxy.get() if hasattr(proxy, "get") else proxy
    return h


def align_gen(h, gen_axes_meta, mass_window=GEN_MASS_WINDOW, tol=1e-6):
    """Reduce a gen-level Hist to (NptVGen, NabsYVGen) on the model's gen grid.

    Selects the UL helicity component (angular-integrated total), the mass
    window bin, sums charge, then rebins (ptVgen, absYvgen) onto the model's
    gen edges (which must be a sub-binning of the gen hist's). Returns an
    ndarray with shape == model.gen_shape, axes (ptVGen, absYVGen).
    """
    axnames = [a.name for a in h.axes]
    # 1. UL helicity component (value -1): the angular-integrated cross section.
    if "helicity" in axnames:
        h = h[{"helicity": h.axes["helicity"].index(-1)}]
    elif "helicitySig" in axnames:
        h = h[{"helicitySig": h.axes["helicitySig"].index(-1)}]
    # 2. Mass-window bin (the gen hist's coarse mass axis brackets [60,120]).
    for mname in ("massVgen", "massVGen"):
        if mname in [a.name for a in h.axes]:
            medges = np.asarray(h.axes[mname].edges, dtype=np.float64)
            lo = np.where(np.isclose(medges, mass_window[0], atol=1e-3))[0]
            hi = np.where(np.isclose(medges, mass_window[1], atol=1e-3))[0]
            if lo.size and hi.size and hi[0] == lo[0] + 1:
                h = h[{mname: int(lo[0])}]  # the single [60,120] bin
            else:
                raise ValueError(
                    f"align_gen: mass window {mass_window} is not a single bin "
                    f"of {mname} (edges {medges[:6]}…)."
                )
            break
    # 3. Project to the gen kinematic axes (sums charge + anything left).
    gen_names = [n for n, _ in gen_axes_meta]  # ("ptVGen", "absYVGen")
    src_names = []
    for want in gen_names:
        # model uses ptVGen/absYVGen (capital G); gen hist uses ptVgen/absYvgen
        cand = next(
            (a.name for a in h.axes if a.name.lower() == want.lower()), None
        )
        if cand is None:
            raise ValueError(f"align_gen: gen hist has no axis matching {want!r}")
        src_names.append(cand)
    h = h.project(*src_names)
    vals = h.values(flow=False)  # (ptVgen_fine, absYvgen_fine), model-axis order

    # 4. Rebin each axis onto the model grid (model edges ⊆ gen-hist edges).
    Ws = []
    for src, (mname, medges) in zip(src_names, gen_axes_meta):
        Ws.append(_merge_matrix(h.axes[src].edges, medges, name=mname))
    # out[i,j] = Σ_{p,q} W0[i,p] vals[p,q] W1[j,q]
    out = Ws[0] @ vals @ Ws[1].T
    return out


def _parse_variation_label(label, param_order):
    """Parse a vars-axis label into ``(param_name, value)``.

    Labels are ``<param><value>`` e.g. ``lambda21.0`` → (lambda2, 1.0),
    ``lambda2_nu0.25`` → (lambda2_nu, 0.25), ``delta_lambda2-0.02`` →
    (delta_lambda2, -0.02). Longest-prefix match against the model's parameter
    names so ``lambda2_nu`` wins over ``lambda2``.
    """
    cands = sorted(
        (p for p in param_order if label.startswith(p)), key=len, reverse=True
    )
    if not cands:
        raise ValueError(
            f"variation label {label!r} doesn't start with any model param "
            f"{tuple(param_order)} — is it a λ variation (not a scale/FO one)?"
        )
    name = cands[0]
    try:
        value = float(label[len(name) :])
    except ValueError as exc:
        raise ValueError(
            f"can't parse a value from {label!r} after param {name!r}: {exc}"
        ) from exc
    return name, value


def validate_variation(model, args, var_label, out_path=None):
    """Off-central λ-response check for one ``vars`` label.

    Compares the model's ``rnorm = σ_reco(λ)/σ_reco(λ_c)`` (via ``model.compute``
    with that single λ overridden) to the histmaker's ``Corr[var]/Corr[pdf0]``.
    Both are variation/central ratios → no normalization, non-circular: it tests
    whether the model's λ-sensitivity (what the fit applies) matches the MC's.
    """
    import os

    import numpy as np
    import tensorflow as tf

    # ---- histmaker λ-reweighting ratio, aligned onto the model reco binning ----
    corr = load_nominal(args.histmaker, args.sample, args.variation_hist)
    va = corr.axes["vars"]
    labels = list(va)
    for need in ("pdf0", var_label):
        if need not in labels:
            raise KeyError(f"vars label {need!r} not in {args.variation_hist}: {labels}")
    h_cen = align_nominal(corr[{"vars": va.index("pdf0")}], model._reco_axes_meta)
    h_var = align_nominal(corr[{"vars": va.index(var_label)}], model._reco_axes_meta)

    # ---- model rnorm for the same single-λ override ----
    pname, pval = _parse_variation_label(var_label, list(model._param_order))
    idx = list(model._param_order).index(pname)
    param = model.xparamdefault.numpy().copy()
    central = float(param[idx])
    param[idx] = pval
    rnorm = model.compute(tf.constant(param, dtype=model.indata.dtype))
    ratio_model = (
        rnorm[:, model.signal_proc_idx].numpy().reshape(model.reco_shape)
    )

    cen = h_cen.values(flow=False).astype(np.float64)
    var = h_var.values(flow=False).astype(np.float64)
    good = cen > 0
    ratio_hist = np.where(good, var / np.where(good, cen, 1.0), np.nan)
    dr = ratio_model / ratio_hist  # double ratio model/hist, per reco bin
    finite = np.isfinite(dr)
    drf, wf = dr[finite], var[finite]

    print(f"\n=== variation {var_label!r}:  {pname}  {central:g} -> {pval:g} ===")
    print(f"  model rnorm range : {np.nanmin(ratio_model):.4f} .. {np.nanmax(ratio_model):.4f}")
    print(f"  hist  rnorm range : {np.nanmin(ratio_hist):.4f} .. {np.nanmax(ratio_hist):.4f}")
    print(
        f"  double-ratio model/hist : mean={drf.mean():.5f} median={np.median(drf):.5f} "
        f"std={drf.std():.5f} min/max={drf.min():.5f}/{drf.max():.5f}"
    )
    print(
        f"  yield-weighted mean|model/hist - 1| : "
        f"{np.average(np.abs(drf - 1.0), weights=wf):.5f}"
    )

    if out_path:
        import hist
        from wums import plot_tools

        from wremnants.postprocessing.scetlib_np import plot_output

        # Project YIELDS then ratio (ratio of sums, not sum of per-bin ratios).
        h_sig_c = tf_to_hist(model.sigma_reco_central, model._reco_axes_meta)
        h_sig_l = tf_to_hist(
            model.sigma_reco_central.numpy().reshape(model.reco_shape) * ratio_model,
            model._reco_axes_meta,
        )
        reco_edges = {
            n: np.asarray(e, dtype=np.float64) for n, e in model._reco_axes_meta
        }

        def _rnorm_1d(hnum, hden, axis):
            # flow-safe projection (see project_inrange): align_nominal parks the
            # cropped out-of-range ptll content in overflow, which a plain
            # project() would fold into the yll spectrum.
            num = project_inrange(hnum, axis).values(flow=False)
            den = project_inrange(hden, axis).values(flow=False)
            ax = hist.axis.Variable(
                reco_edges[axis], name=axis, underflow=False, overflow=False
            )
            h = hist.Hist(ax, storage=hist.storage.Double())
            h.view(flow=False)[...] = np.where(den > 0, num / den, np.nan)
            return h

        outdir = os.path.dirname(out_path) or "."
        os.makedirs(outdir, exist_ok=True)
        base0 = os.path.splitext(os.path.basename(out_path))[0]
        # Project the λ-response onto each reco shape axis: ptll (the qT response)
        # and yll (the rapidity response — should be reproduced, not spurious).
        for ax_name in ("ptll", "yll"):
            h_rn_hist = _rnorm_1d(h_var, h_cen, ax_name)
            h_rn_model = _rnorm_1d(h_sig_l, h_sig_c, ax_name)
            # Main-panel y-range: centered on 1, padded to the actual rnorm span so
            # the variation is visible (auto-range from 0 squashes it at the top);
            # floor it so a near-flat yll response isn't zoomed onto numerical noise.
            rn_vals = np.concatenate(
                [h_rn_hist.values(flow=False), h_rn_model.values(flow=False)]
            )
            m = float(np.nanmax(np.abs(rn_vals - 1.0)))
            pad = max(1.15 * m, 0.005)
            ylim = [1.0 - pad, 1.0 + pad]
            fig = plot_tools.makePlotWithRatioToRef(
                [h_rn_hist, h_rn_model],
                labels=[
                    f"histmaker  ({pname}: {central:g}→{pval:g})",
                    f"ParamModel  ({pname}: {central:g}→{pval:g})",
                ],
                colors=["black", "#e42536"],  # CMS red for the model
                xlabel=AXIS_LABELS.get(ax_name, ax_name),
                ylabel=r"$\sigma_{reco}(\lambda)\,/\,\sigma_{reco}(\lambda_c)$",
                rlabel=["model / histmaker"],
                rrange=[[0.99, 1.01]],
                autorrange=0.3,  # tighten onto the actual ~0.02-0.05% agreement
                ylim=ylim,
                binwnorm=None,
                nlegcols=1,
                grid=True,
                yerr=False,
                ratio_legend=False,
                plot_title=None,
            )
            base = base0 if ax_name == "ptll" else f"{base0}_{ax_name}"
            plot_output.save_plot(outdir, base, fig=fig, args=args)
            print(f"  rnorm plot ({ax_name}) written to: {outdir}/{base}.png(.pdf) + {base}.log")


# =============================================================================
# theory-correction reference + λ-tune resolvers
# (moved verbatim from sigma_gen_at_lambda)
# =============================================================================
Q_LO, Q_HI = 60.0, 120.0
# Canonical FranksVals (CT18Z N3+0LL lattice λ4-bugfix) tanh_2 runcard — the
# production λ_central. Construction BASE when no base λ is sourced: the model
# must be built at a PHYSICAL tune (positive σ_gen, so the constructor's response
# guard passes), and the requested λ evaluated on top. Source of truth: a
# correction file's Nonperturbative section (file_meta_data → config →
# Nonperturbative); the LatticeNPLambda4Bugfix_FranksVals_CT18Z values.
CANONICAL_BASE = {
    "eff_params": {
        "np_model": "tanh_2", "lambda2": 0.4, "lambda4": 0.4,
        "lambda6": 0.0, "delta_lambda2": 0.0, "lambda_inf": 1.0,
    },
    "gnu_params": {
        "np_model_nu": "tanh_2", "lambda2_nu": 0.15,
        "lambda4_nu": 0.0, "lambda6_nu": 0.0, "lambda_inf_nu": 2.0,
    },
}

# Built-in gen grid used when neither explicit edges nor a datacard/hdf5 source is
# given: 1-GeV ptVGen bins over [0, 40], rapidity-inclusive in a single absYVGen
# bin [0, 5] (5.0 is a TheoryCorrection absY edge, so the overlay still aligns).
DEFAULT_PTV_EDGES = np.arange(0.0, 41.0, 1.0)
DEFAULT_ABSY_EDGES = np.array([0.0, 5.0])

# corr-hist axis ↔ model gen axis (the TheoryCorrection _hist uses SCETlib names).
_CORR_AXIS = {"ptVGen": "qT", "absYVGen": "absY"}


def _parse_edges(s):
    """``a,b,c,...`` -> float ndarray of bin edges."""
    return np.array([float(x) for x in s.split(",") if x.strip()], dtype=np.float64)


def resolve_base_lambda(args):
    """Physical BASE λ tune (eff_params/gnu_params) the model is CONSTRUCTED at.

    Priority: ``--meta-from HDF5`` > the ``--theory-corr`` file's embedded
    Nonperturbative runcard > the canonical FranksVals tanh_2 default. Always a
    complete physical tune (never None), so construction lands on a positive-σ_gen
    point (the constructor's response guard); the requested λ are evaluated on top.
    """
    import os

    from wremnants.postprocessing.scetlib_np import lambda_central as lc

    if args.meta_from:
        print(f"[λ base] from hdf5 metadata {args.meta_from}")
        return lc.read_lambda_central(args.meta_from)
    if args.theory_corr:
        import pickle

        import lz4.frame

        with lz4.frame.open(args.theory_corr) as fh:
            corr = pickle.load(fh)
        base = lc.extract_lambda_central(
            corr, tag=os.path.basename(args.theory_corr),
            proc=args.theory_corr_proc or "Z",
        )
        print(f"[λ base] from the --theory-corr Nonperturbative runcard "
              f"({base.get('basename')})")
        return {"eff_params": base["eff_params"], "gnu_params": base["gnu_params"]}
    print("[λ base] none given -> canonical FranksVals tanh_2 default")
    return CANONICAL_BASE


def assemble_tune(base, overrides):
    """Full (eff_params, gnu_params) for the EVAL point = base tune + overrides,
    plus the explicitly-set names.

    ``base`` is a physical lambda_central dict (with the np_model form strings);
    params not in ``overrides`` stay at the base value (NOT 0). Each override is
    routed to eff or gnu by membership."""
    eff = dict(base["eff_params"])
    gnu = dict(base["gnu_params"])
    explicit = {}
    for name, val in overrides.items():
        if name in EFF_PARAMS:
            eff[name] = val
        elif name in GNU_PARAMS:
            gnu[name] = val
        else:
            raise SystemExit(
                f"unknown λ {name!r}; valid: {list(GNU_PARAMS) + list(EFF_PARAMS)}"
            )
        explicit[name] = val
    return eff, gnu, explicit


def resolve_gen_axes(args):
    """gen_axes = [(ptVGen, edges), (absYVGen, edges)], chosen per axis in order:
    explicit --ptv-edges/--absy-edges, then a --gen-edges-from/--datacard hdf5,
    then the built-in defaults (1-GeV ptVGen [0,40]; single absYVGen [0,5])."""
    ptv = _parse_edges(args.ptv_edges) if args.ptv_edges else None
    absy = _parse_edges(args.absy_edges) if args.absy_edges else None
    src = args.gen_edges_from or args.datacard

    src_axes = None
    if (ptv is None or absy is None) and src:
        print(f"[gen-axes] reading the scetlib_np auxiliary of {src}")
        from rabbit.inputdata import FitInputData

        from wremnants.postprocessing.scetlib_np.param_model import (
            _R_info_from_auxiliary,
        )

        indata = FitInputData(src)
        src_axes = {
            n: np.asarray(e, dtype=np.float64)
            for n, e in _R_info_from_auxiliary(indata)["gen_axes"]
        }

    def pick(name, explicit, default):
        if explicit is not None:
            print(f"[gen-axes] {name}: explicit ({explicit.size - 1} bins)")
            return explicit
        if src_axes is not None and name in src_axes:
            print(f"[gen-axes] {name}: from {src} ({src_axes[name].size - 1} bins)")
            return src_axes[name]
        print(f"[gen-axes] {name}: built-in default ({default.size - 1} bins, "
              f"[{default[0]:g}, {default[-1]:g}])")
        return default

    return [
        ("ptVGen", pick("ptVGen", ptv, DEFAULT_PTV_EDGES)),
        ("absYVGen", pick("absYVGen", absy, DEFAULT_ABSY_EDGES)),
    ]


def load_theory_corr_hist(path, proc=None):
    """Load the ``{generator}_hist`` (SCETlib+DYTurbo) Hist from a TheoryCorrection
    ``.pkl.lz4``.

    The file maps ``corr[proc][histname]``. ``proc`` defaults to the single
    physics key (``meta_data`` / ``file_meta_data`` excluded); the hist is the
    lone ``*_hist`` that is not ``minnlo_ref_hist`` (the prediction, not the
    MiNNLO reference or the ratio).
    """
    import os
    import pickle

    import lz4.frame

    with lz4.frame.open(path) as fh:
        corr = pickle.load(fh)

    meta_keys = {"meta_data", "file_meta_data"}
    procs = [k for k in corr.keys() if k not in meta_keys]
    if proc is None:
        if len(procs) != 1:
            raise SystemExit(
                f"--theory-corr-proc needed: {os.path.basename(path)} has procs {procs}"
            )
        proc = procs[0]
    elif proc not in corr:
        raise SystemExit(f"proc {proc!r} not in {list(corr.keys())}")

    entry = corr[proc]
    cands = [k for k in entry if k.endswith("_hist") and k != "minnlo_ref_hist"]
    if len(cands) != 1:
        raise SystemExit(
            f"expected one {{generator}}_hist in {proc}; found {list(entry.keys())}"
        )
    histname = cands[0]

    print(f"[theory-corr] {os.path.basename(path)} :: {proc} / {histname}")
    return entry[histname]


def theory_corr_projection(h, gen_axes, plot_axis, var="pdf0", q_window=(Q_LO, Q_HI),
                           tol=1e-6):
    """Project a TheoryCorrection ``_hist`` onto the model's ``plot_axis`` gen bins.

    Reduces the (Q, absY, qT, charge, vars) Hist to a 1-D bin-integrated σ on the
    model's ``plot_axis`` edges, restricted to the model's gen-grid extent on the
    OTHER axis so it covers the same phase space the model σ_gen projection does:

      1. select the ``vars`` entry (default ``pdf0`` = central tune);
      2. sum the Q bins whose centre falls in ``q_window`` (in-range only);
      3. sum the charge axis (in-range), if present;
      4. sum the OTHER gen axis over the model's extent [0, other_max];
      5. rebin the projection axis onto the model's ``plot_axis`` edges (model
         edges must be a sub-binning of the corr hist's: qT is fine enough that
         ptVGen always aligns; absY uses SCETlib's binning so absYVGen may not).

    Returns an ndarray of length ``len(plot_axis edges) - 1`` (bin-integrated σ).
    """
    names = [n for n, _ in gen_axes]
    if plot_axis not in names or plot_axis not in _CORR_AXIS:
        raise SystemExit(f"--plot-axis {plot_axis!r} not a model gen axis {names}")
    edges_by_name = {n: np.asarray(e, dtype=np.float64) for n, e in gen_axes}
    other_model = names[1] if plot_axis == names[0] else names[0]
    proj_corr = _CORR_AXIS[plot_axis]
    other_corr = _CORR_AXIS[other_model]

    have = [a.name for a in h.axes]
    for need in ("Q", proj_corr, other_corr, "vars"):
        if need not in have:
            raise SystemExit(
                f"theory-corr hist missing {need!r} axis; has {have}"
            )

    # 1. vars selection.
    vlist = list(h.axes["vars"])
    if var not in vlist:
        raise SystemExit(
            f"--theory-corr-var {var!r} not in corr hist vars; have {vlist}"
        )
    h = h[{"vars": vlist.index(var)}]

    # 2. Q window (sum in-range bins whose centre is inside the window).
    qe = np.asarray(h.axes["Q"].edges, dtype=np.float64)
    qc = 0.5 * (qe[:-1] + qe[1:])
    qsel = np.where((qc >= q_window[0] - tol) & (qc <= q_window[1] + tol))[0]
    if not qsel.size:
        raise SystemExit(f"no Q bins in window {q_window}; corr Q edges {qe}")
    h = h[{"Q": slice(int(qsel[0]), int(qsel[-1]) + 1, sum)}]

    # 3. charge sum (in-range), if a charge axis is present.
    if "charge" in [a.name for a in h.axes]:
        h = h[{"charge": slice(0, h.axes["charge"].size, sum)}]

    # 4. sum the OTHER axis over the model's extent [0, other_max]. Non-coinciding
    #    upper edge (SCETlib's absY binning): cut at the nearest corr edge, warn.
    other_max = edges_by_name[other_model][-1]
    oe = np.asarray(h.axes[other_corr].edges, dtype=np.float64)
    oc = 0.5 * (oe[:-1] + oe[1:])
    osel = np.where(oc <= other_max + tol)[0]
    if not osel.size:
        raise SystemExit(
            f"theory-corr {other_corr} has no bins below the model {other_model} "
            f"max {other_max}; corr edges {oe}"
        )
    cut_idx = int(osel[-1]) + 1
    actual_edge = oe[cut_idx]
    if abs(actual_edge - other_max) > tol:
        print(
            f"[theory-corr] WARNING: model {other_model} max {other_max} does not "
            f"coincide with a corr {other_corr} edge; summing corr up to "
            f"{actual_edge} ({abs(actual_edge - other_max):.3g} off)."
        )
    h = h[{other_corr: slice(0, cut_idx, sum)}]

    # 5. rebin the projection axis onto the model's plot_axis edges.
    W = _merge_matrix(
        np.asarray(h.axes[proj_corr].edges, dtype=np.float64),
        edges_by_name[plot_axis],
        name=plot_axis,
        tol=tol,
    )
    return W @ np.asarray(h.values(flow=False), dtype=np.float64)
