"""Shared plotting/array helpers for the SCETlib-NP validation scripts.

Extracted here so the packaged ``param_model_diagnostics`` (the default-on
in-fit reco guard) and the validation scripts share one implementation without a
packaged module importing an untracked script. numpy at import; ``hist`` and
``wums.plot_tools`` are imported lazily so importing this module stays cheap.
"""

import numpy as np

from wremnants.postprocessing.scetlib_np.params import bin_sum_matrix

AXIS_LABELS = {
    "ptll": r"$p_{T}^{\ell\ell}$ [GeV]",
    "yll": r"$y^{\ell\ell}$",
    "ptVGen": r"$p_{T}^{V,\,gen}$ [GeV]",
    "absYVGen": r"$|y^{V,\,gen}|$",
}


def _merge_matrix(fine_edges, coarse_edges, name="axis", tol=1e-6):
    """(N_coarse, N_fine) 0/1 matrix summing fine bins into coarse bins.

    Every coarse edge must coincide with a fine edge (the coarse binning is a
    sub-binning of the fine one) — exact merge, no interpolation. The matrix
    itself is :func:`params.bin_sum_matrix` on the fine-bin centres.
    """
    fine_edges = np.asarray(fine_edges, dtype=np.float64)
    coarse_edges = np.asarray(coarse_edges, dtype=np.float64)
    for e in coarse_edges:
        if not np.any(np.isclose(fine_edges, e, atol=tol)):
            raise ValueError(
                f"_merge_matrix[{name}]: coarse edge {e} not a fine edge "
                f"(gen hist binning is not a refinement of the model grid)."
            )
    centers = 0.5 * (fine_edges[:-1] + fine_edges[1:])
    return bin_sum_matrix(centers, coarse_edges, tol)


def project_inrange(h, axis):
    """1D Hist on ``axis``, summing the OTHER axes over their in-range bins only.

    Unlike ``hist.project()``, which folds the summed axes' flow bins back in:
    ``align_nominal`` crops out-of-fit-range ptll content into the ptll overflow,
    so projecting with flow would re-add it and tilt the spectrum. The fit and
    ``summarize`` use only the in-range (flow=False) bins, so the shape plots
    must too. Identical to ``project`` for a no-flow hist (e.g. the model).
    """
    import hist as _hist

    names = [a.name for a in h.axes]
    ai = names.index(axis)
    vals = h.values(flow=False)
    other = tuple(i for i in range(vals.ndim) if i != ai)
    out = _hist.Hist(
        _hist.axis.Variable(
            h.axes[axis].edges, name=axis, underflow=False, overflow=False
        ),
        storage=_hist.storage.Double(),
    )
    out.view(flow=False)[...] = vals.sum(axis=other)
    return out


def plot_ptll_ratio(
    h_model,
    h_nominal,
    axis="ptll",
    out_path=None,
    scale=1.0,
    ref_label="histmaker nominal",
    model_label=r"ParamModel $\sigma_{reco}$",
    rlabel="model / nominal",
    title=None,
    rrange=(0.9, 1.1),
    autorrange=None,
    xlim=None,
    density=False,
    extra_text=None,
    extra_text_loc=None,
    ratio_legend=True,
    no_sci=False,
    extra_models=None,
    args=None,
):
    """Plot the ``axis``-projected model vs reference spectra + ratio via wums.

    Both inputs are ``hist.Hist`` on the same axes (see ``tf_to_hist``); each is
    projected onto ``axis`` and passed to ``makePlotWithRatioToRef`` (handles the
    main/ratio panels and bin-width normalization). ``h_nominal`` is the
    denominator. ``scale`` multiplies the model (e.g. Σref/Σmodel to compare
    shape rather than absolute normalization).

    ``density`` normalizes each curve to unit integral (a pure shape comparison,
    removing the model↔histmaker unit convention and any total-σ difference);
    ``scale`` is then ignored. ``extra_models`` overlays additional curves as
    ``(hist, label, color, scale)`` tuples. No title is drawn by default; pass a
    non-empty ``title`` to add one.
    """
    import os

    from wums import plot_tools

    from wremnants.postprocessing.scetlib_np import plot_output

    def _unit(h):
        return h * (1.0 / float(h.values(flow=False).sum()))

    h_n = project_inrange(h_nominal, axis)
    h_m = project_inrange(h_model, axis)
    if density:
        h_n, h_m = _unit(h_n), _unit(h_m)
        suffix = ""
    else:
        if scale != 1.0:
            h_m = h_m * scale
        suffix = rf" (×{scale:.4g})" if scale != 1.0 else ""

    hists = [h_n, h_m]
    labels = [ref_label, model_label + suffix]
    colors = ["black", "#e42536"]  # CMS red for the model curve
    for h_e, lab_e, col_e, sc_e in extra_models or []:
        h_ep = project_inrange(h_e, axis)
        if density:
            h_ep = _unit(h_ep)
            labels.append(lab_e)
        else:
            if sc_e != 1.0:
                h_ep = h_ep * sc_e
            labels.append(lab_e + (rf" (×{sc_e:.4g})" if sc_e != 1.0 else ""))
        hists.append(h_ep)
        colors.append(col_e)

    # No title by default; pass an explicit non-empty ``title`` to force one.
    plot_title = title or None

    fig = plot_tools.makePlotWithRatioToRef(
        hists,
        labels=labels,
        colors=colors,
        xlabel=AXIS_LABELS.get(axis, axis),
        ylabel=(
            (
                r"$\frac{1}{\sigma}\,\frac{d\sigma}{dy}$"
                if axis in ("yll", "absYVGen")
                else r"$\frac{1}{\sigma}\,\frac{d\sigma}{dp_{T}}$ [1/GeV]"
            )
            if density else r"$\sigma$/bin"
        ),
        rlabel=[rlabel],
        rrange=[list(rrange)],
        autorrange=autorrange,
        xlim=xlim,
        binwnorm=1.0,
        nlegcols=1,
        grid=True,
        yerr=False,
        extra_text=extra_text,
        extra_text_loc=extra_text_loc,
        ratio_legend=ratio_legend,
        no_sci=no_sci,
        plot_title=plot_title,
    )

    outdir, basename = plot_output.split_outpath(out_path)
    plot_output.save_plot(outdir, basename, fig=fig, args=args)
    print(f"\n  {axis}-projection ratio plot written to: "
          f"{outdir}/{basename}.png(.pdf) + {basename}.log")


def summarize(model_sigma, nominal, reco_axes_meta):
    """Print shape-comparison diagnostics between the two reco tensors."""
    names = [n for n, _ in reco_axes_meta]
    assert model_sigma.shape == nominal.shape, (model_sigma.shape, nominal.shape)

    m = model_sigma.astype(np.float64)
    n = nominal.astype(np.float64)
    msum, nsum = m.sum(), n.sum()
    print(f"\n  Σ model σ_reco : {msum:.6g}")
    print(f"  Σ nominal      : {nsum:.6g}")
    print(f"  global scale Σmodel/Σnominal : {msum / nsum:.6g}")

    # Single global scale so we compare SHAPES, not absolute normalization.
    scale = nsum / msum
    m_scaled = m * scale

    good = n > 0
    n_bad = int((~good).sum())
    if n_bad:
        print(f"  ({n_bad} of {n.size} bins have nominal<=0; excluded from ratio stats)")
    ratio = np.full(m.shape, np.nan)
    ratio[good] = m_scaled[good] / n[good]
    r = ratio[good]

    print("\n  per-bin shape ratio  (scale·model / nominal), should ~ 1:")
    print(f"    bins compared : {r.size}")
    print(f"    mean          : {r.mean():.5f}")
    print(f"    median        : {np.median(r):.5f}")
    print(f"    std           : {r.std():.5f}")
    print(f"    min / max     : {r.min():.5f} / {r.max():.5f}")
    for q in (1, 5, 50, 95, 99):
        print(f"    p{q:<2d}           : {np.percentile(r, q):.5f}")
    # weighted (by yield) mean abs deviation — what matters for the fit
    w = n[good]
    wmad = np.average(np.abs(r - 1.0), weights=w)
    print(f"    yield-weighted mean|ratio-1| : {wmad:.5f}")

    # 1D projections (normalized) per axis
    ndim = m.ndim
    for ax in range(ndim):
        other = tuple(i for i in range(ndim) if i != ax)
        mp = m_scaled.sum(axis=other)
        npj = n.sum(axis=other)
        with np.errstate(divide="ignore", invalid="ignore"):
            rp = np.where(npj > 0, mp / npj, np.nan)
        print(f"\n  --- 1D projection on {names[ax]} (n={mp.size} bins) ---")
        print("    idx     model(scaled)        nominal        ratio")
        for i in range(mp.size):
            rr = f"{rp[i]:.4f}" if np.isfinite(rp[i]) else "   nan"
            print(f"    {i:3d}  {mp[i]:15.6g}  {npj[i]:15.6g}    {rr}")


def tf_to_hist(tensor, reco_axes_meta):
    """Build a ``hist.Hist`` on the model's reco axes from a model output.

    ``tensor`` is a model output for one process — a ``tf.Tensor``/ndarray flat
    over the reco bins (e.g. ``model.sigma_reco_central``) or already reco-shaped.
    ``reco_axes_meta`` is ``model._reco_axes_meta`` (the ``(name, edges)`` list).
    Returns a Double-storage Hist on the model's axes so ``.project(...)`` and
    ``divideHists`` work.
    """
    import hist

    vals = tensor.numpy() if hasattr(tensor, "numpy") else np.asarray(tensor)
    vals = np.asarray(vals, dtype=np.float64)

    axes = [
        hist.axis.Variable(
            np.asarray(edges, dtype=np.float64),
            name=name,
            underflow=False,
            overflow=False,
        )
        for name, edges in reco_axes_meta
    ]
    shape = tuple(ax.size for ax in axes)
    if vals.size != int(np.prod(shape)):
        raise ValueError(
            f"tf_to_hist: tensor has {vals.size} entries but reco axes imply "
            f"{int(np.prod(shape))} ({shape})."
        )
    h = hist.Hist(*axes, storage=hist.storage.Double())
    h.view(flow=False)[...] = vals.reshape(shape)
    return h
