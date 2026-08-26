"""Native-binning validation of the SCETlib bT-grid factorization against a
SCETlib spectrum-mode reference run — resummed only, NO projection.

This is the WRemnants-side reproduction of ``scetlib-np-factorize.py validate``
(``cmd_validate``), the comparison the HANDOFF documents at the 0.06-0.08 %
numerical floor. It deliberately stays in the reference's *native* fine binning
(per-(Y, qT) bin in the Z Q-bin) — NO |Y|-fold, NO rebin to the coarse gen grid,
NO ptVGen projection — so it isolates the bT-integral + Q-integration from the
projection layer that the ParamModel adds on top.

Curves (all **resummed only**, same NP, native binning):
  (1) reference : the SCETlib spectrum-mode combined.pkl, read with
                  ``input_tools.read_scetlib_hist`` (the WRemnants reader),
                  verified bit-identical to ``factorize.load_spectrum_reference``
                  (the reader cmd_validate uses) — so read_scetlib_hist IS the
                  #1 reference.
  (2) factorize : the bT-grid Hankel reconstruction via the numpy ``factorize``
                  library (reconstruct_grid_QYqT + integrate_over_Q + per-bin
                  Simpson), exactly as cmd_validate.
  (3) ParamModel: the SAME reconstruction through the TF model, read from
                  ``model.sigma_YqT_central`` (native (NY, NqT), BEFORE the
                  fold/rebin). Only built when --paramodel is given.

The reference is bin-INTEGRATED; the bT-grid is point-sampled on the edges+centres
of those bins, so the factorized/model sides are Simpson-integrated over each
reference bin to compare like-for-like. Curve 3 is also compared element-wise to
curve 2 on the native grid (a parity check of the TF port at λ_central).

Inputs live under /work and /ceph — run inside a container that binds both, e.g.
    export APPTAINER_BIND="/scratch,/cvmfs,/work,/ceph,/home"
    singularity run <wmassdevrolling:v46> bash -c \\
      "source main/WRemnants/setup.sh; \\
       python3 -m wremnants.postprocessing.scetlib_np.validation.native_validation \\
         --paramodel"
"""

import argparse
import configparser
import importlib.util
import json
import os
import pickle
import sys
import time

import numpy as np

from wremnants.postprocessing.scetlib_np.sigma_gen import _default_btgrid_dir

# The factorize library (pure numpy; no SCETlib runtime needed).
SCETLIB_RUN = "/work/submit/lavezzo/alphaS/scetlib-cms-newnp-lambda4fix/prod/scetlib_run"
if SCETLIB_RUN not in sys.path:
    sys.path.insert(0, SCETLIB_RUN)
from scetlib_run import factorize  # noqa: E402

from wremnants.utilities.io_tools import input_tools  # noqa: E402

# ---- Defaults (the documented native baseline) ----
LAMBDA6_REF = (
    "/ceph/submit/data/user/l/lavezzo/zstuff/"
    "Z_COM13_CT18Z_N3p0LL_NewNPs_Lattice_Newvals_Lambda6_Fine/combined_out/"
    "inclusive_Z_COM13_CT18Z_N3+0LL_newvals_lambda6_fine_combined.pkl"
)
FRANKS_REF = (
    "/work/submit/areimers/wmass/TheoryCorrections/SCETlib/"
    "com13_ct18z_newnps_n3+0ll_lattice_lambda4bugfix_franksvals_fine/"
    "inclusive_Z_COM13_CT18Z_N3+0LL_lattice_lambda4bugfix_franksvals_fine_combined.pkl"
)
# The ParamModel fit inputs (datacard carries λ_central + binning).
SIGNAL_PROC = "Zmumu"
Q_LO, Q_HI = 60.0, 120.0
# The scetlib+DYTurbo *matched* spectrum (resum + nonsingular FO) — read directly
# from the final theory-correction product. d["Z"]["..._hist"] is the absolute
# matched cross section on Q[60,120] × absY × qT (native 70-bin qT, |Y| folded);
# d["Z"]["..._minnlo_ratio"] is the correction ratio; "minnlo_ref_hist" is MiNNLO.
CORRZ = (
    "/home/submit/lavezzo/alphaS/main/WRemnants/wremnants-data/data/TheoryCorrections/"
    "scetlib_dyturbo_LatticeNPLambda4Bugfix_FranksVals_CT18Z_N3p0LL_N2LO_CorrZ.pkl.lz4"
)
# Fixed-order inputs for the nonsingular σ_ns = DYTurbo − SCETlib_singular
# (NP-independent), shipped in wremnants-data (same defaults as the ParamModel).
from wremnants.postprocessing.scetlib_np.param_model import (  # noqa: E402
    _NONSING_DYTURBO_DEFAULT as DYTURBO,
    _NONSING_FO_SING_DEFAULT as FO_SING,
)


def _import_cli_helpers():
    """Import the (hyphenated-filename) CLI module to reuse its reference
    extraction helpers, so we match cmd_validate's bin-handling exactly."""
    cli_path = os.path.join(SCETLIB_RUN, "scetlib-np-factorize.py")
    spec = importlib.util.spec_from_file_location("scetlib_np_factorize_cli", cli_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _np_params_from_config(config_dict):
    """eff/gnu NP param dicts from the reference pickle's config (cmd_validate fallback)."""
    cp = configparser.ConfigParser(inline_comment_prefixes="#")
    for sec_name, sec_dict in config_dict.items():
        cp[sec_name] = {k: str(v) for k, v in sec_dict.items()}
    return factorize.eff_params_from_conf(cp), factorize.gnu_params_from_conf(cp)


def load_btgrid_dense(btgrid_dir):
    """Load btgrid shards, return (grid, Q_grid, Y_grid, qT_grid) with grid['I_pert']/
    ['C_nu'] NaN-filled to a full (1, NQ*NY*NqT, NbT) product."""
    grid = factorize.load_btgrid_shards(btgrid_dir)
    Q_grid = np.array(sorted({b[0] for b in grid["bins"]}))
    Y_grid = np.array(sorted({b[1] for b in grid["bins"]}))
    qT_grid = np.array(sorted({b[2] for b in grid["bins"]}))
    NQ, NY, NqT = Q_grid.size, Y_grid.size, qT_grid.size
    n_present, n_expected = len(grid["bins"]), NQ * NY * NqT
    if n_present != n_expected:
        NbT = grid["bT"].size
        I_full = np.full((NQ * NY * NqT, NbT), np.nan)
        C_full = np.full_like(I_full, np.nan)
        Qi = {q: i for i, q in enumerate(Q_grid)}
        Yi = {y: i for i, y in enumerate(Y_grid)}
        Ti = {q: i for i, q in enumerate(qT_grid)}
        Isrc, Csrc = grid["I_pert"][0], grid["C_nu"][0]
        for r, b in enumerate(grid["bins"]):
            flat = ((Qi[b[0]] * NY) + Yi[b[1]]) * NqT + Ti[b[2]]
            I_full[flat] = Isrc[r]
            C_full[flat] = Csrc[r]
        grid["I_pert"] = I_full[np.newaxis]
        grid["C_nu"] = C_full[np.newaxis]
    return grid, Q_grid, Y_grid, qT_grid


def native_ratios(sigma_YqT, ref_spectra, ref_widths, Y_grid, qT_grid, Q_ref):
    """Match a factorized/model σ(Y,qT) (native grid) to the bin-integrated
    reference per (Y,qT) bin in the Z Q-bin, Simpson-integrating the factorized
    side over each bin (cmd_validate's logic). Returns (ratios, Y_arr, qT_arr,
    ref_arr, fac_arr, n_simpson, n_point) — ref_arr/fac_arr are the per-cell
    bin-integrated reference and factorized values."""
    Y_idx = {y: i for i, y in enumerate(Y_grid)}
    qT_idx = {q: i for i, q in enumerate(qT_grid)}
    tol = 1e-9
    ratios, keys, ref_vals, fac_vals = [], [], [], []
    n_simpson = n_point = 0
    for Y_c in sorted({k[1] for k in ref_spectra if k[0] == Q_ref}):
        for qT_c in sorted({k[2] for k in ref_spectra if k[0] == Q_ref}):
            key = (Q_ref, Y_c, qT_c)
            if key not in ref_spectra or Y_c not in Y_idx or qT_c not in qT_idx:
                continue
            ref_val = ref_spectra[key]
            wQ, wY, wqT = ref_widths[key]
            Y_lo, Y_hi = Y_c - 0.5 * wY, Y_c + 0.5 * wY
            qT_lo, qT_hi = qT_c - 0.5 * wqT, qT_c + 0.5 * wqT
            nY_in = int(((Y_grid >= Y_lo - tol) & (Y_grid <= Y_hi + tol)).sum())
            nqT_in = int(((qT_grid >= qT_lo - tol) & (qT_grid <= qT_hi + tol)).sum())
            if nY_in >= 2 and nqT_in >= 2:
                try:
                    s_Y = factorize.integrate_over_qT_bin(sigma_YqT, qT_grid, qT_lo, qT_hi)
                    fac_val = float(factorize.integrate_over_Y_bin(
                        s_Y[:, np.newaxis], Y_grid, Y_lo, Y_hi)[0])
                    n_simpson += 1
                except Exception:
                    fac_val = sigma_YqT[Y_idx[Y_c], qT_idx[qT_c]] * wY * wqT
                    n_point += 1
            else:
                fac_val = sigma_YqT[Y_idx[Y_c], qT_idx[qT_c]] * wY * wqT
                n_point += 1
            if fac_val != 0 and np.isfinite(ref_val) and np.isfinite(fac_val):
                ratios.append(ref_val / fac_val)
                keys.append((Y_c, qT_c))
                ref_vals.append(ref_val)
                fac_vals.append(fac_val)
    ratios = np.array(ratios)
    Y_arr = np.array([k[0] for k in keys])
    qT_arr = np.array([k[1] for k in keys])
    ref_arr = np.array(ref_vals)
    fac_arr = np.array(fac_vals)
    return ratios, Y_arr, qT_arr, ref_arr, fac_arr, n_simpson, n_point


def _qT_spectrum(qT_arr, val_arr, Y_arr, Y_cut):
    """Sum bin-integrated values over |Y|<Y_cut at each qT centre -> (qT, sigma(qT))."""
    m = np.abs(Y_arr) < Y_cut
    qTs = np.array(sorted(set(qT_arr[m].tolist())))
    spec = np.array([val_arr[m & (qT_arr == q)].sum() for q in qTs])
    return qTs, spec


def _spectrum_hist(qT_edges, qT_arr, val_arr, Y_arr, Y_cut, qT_max):
    """Build a 1D qT hist.Hist (native qT binning, qT<=qT_max) filled with the
    bin-integrated value summed over |Y|<Y_cut at each qT centre."""
    import hist

    edges = np.asarray(qT_edges, float)
    edges = edges[edges <= qT_max + 1e-9]
    centres = 0.5 * (edges[:-1] + edges[1:])
    qTs, spec = _qT_spectrum(qT_arr, val_arr, Y_arr, Y_cut)
    lookup = {round(q, 6): s for q, s in zip(qTs, spec)}
    h = hist.Hist(hist.axis.Variable(edges, name="qT"), storage=hist.storage.Double())
    h.view()[...] = np.array([lookup.get(round(c, 6), 0.0) for c in centres])
    return h


def _absY_qT_hist_to_gen(h, gen_axes_meta):
    """Bin-sum a (absY, qT) bin-integrated hist onto the gen grid (NptVGen, NabsYVGen)."""
    from wremnants.postprocessing.scetlib_np.param_model import _bin_sum_matrix

    qT_c = np.asarray(h.axes["qT"].centers, float)
    absY_c = np.asarray(h.axes["absY"].centers, float)
    v = h.project("qT", "absY").values(flow=False)  # (qT, absY)
    Wp = _bin_sum_matrix(qT_c, np.asarray(gen_axes_meta[0][1], float))
    Wa = _bin_sum_matrix(absY_c, np.asarray(gen_axes_meta[1][1], float))
    return Wp @ v @ Wa.T  # (NptVGen, NabsYVGen)


def corrz_matched_to_gen(corrz_path, gen_axes_meta):
    """CorrZ matched absolute spectrum projected to the gen grid (bin-sum)."""
    return _absY_qT_hist_to_gen(load_matched_corrz(corrz_path), gen_axes_meta)


def resum_ref_to_gen(reference, gen_axes_meta, q_lo, q_hi):
    """SCETlib resum reference projected to the gen grid (bin-sum), |Y|-folded —
    the same bin-sum projection as corrz_matched_to_gen, for a like-for-like
    resum vs matched decomposition on the gen grid."""
    return _absY_qT_hist_to_gen(resum_absY_qT(reference, q_lo, q_hi), gen_axes_meta)


def _central_var(h):
    if "vars" in h.axes.name:
        names = list(h.axes["vars"])
        idx = 0
        for cand in ("central", "pdf0", "nominal"):
            if cand in names:
                idx = names.index(cand)
                break
        h = h[{"vars": idx}]
    return h


def load_matched_corrz(corrz_path):
    """Read the scetlib+DYTurbo matched absolute spectrum from a CorrZ.pkl.lz4.
    Returns the central hist on (absY, qT) (Q already [60,120], charge summed)."""
    import lz4.frame

    with lz4.frame.open(corrz_path, "rb") as f:
        d = pickle.load(f)
    Z = d["Z"]
    key = next(k for k in Z if k.endswith("_hist") and k != "minnlo_ref_hist")
    h = _central_var(Z[key])
    if "charge" in h.axes.name:
        h = h[{"charge": slice(0, h.axes["charge"].size, sum)}]
    if "Q" in h.axes.name:
        # slice(...,sum) excludes flow; plain {"Q": sum} would add any Q under/overflow.
        h = h[{"Q": slice(0, h.axes["Q"].size, sum)}]
    return h  # (absY, qT)


def resum_absY_qT(reference, q_lo, q_hi):
    """Resum reference projected to (absY, qT): central var, Q-window, |Y|-folded."""
    from wums import boostHistHelpers as hh

    h = _central_var(input_tools.read_scetlib_hist(reference, charge=0))
    if "charge" in h.axes.name:
        h = h[{"charge": sum}]
    Qe = np.asarray(h.axes["Q"].edges, float)
    qi = int(np.argmin(np.abs(Qe - q_lo)))
    qj = int(np.argmin(np.abs(Qe - q_hi)))
    # slice(...,sum) sums ONLY the in-range Q bins; a plain slice + {"Q":sum}
    # re-adds the sliced-out bins via underflow (the low-mass [10,60] bin → 4× leak).
    h = h[{"Q": slice(qi, qj, sum)}]
    return hh.makeAbsHist(h, "Y")  # signed Y -> absY


def matched_vs_resum_plot(out_path, corrz_path, reference, q_lo, q_hi, Y_cut, qT_max,
                          args=None):
    """Overlay σ(qT) of the resummed reference and the scetlib+DYTurbo matched
    spectrum (summed over |Y|<~Y_cut), on the shared native qT grid, via wums.
    Ratio panel = matched/resummed = the effect of adding the fixed-order/DYTurbo
    matching. Both hists are put on a common absY binning so the |Y| cut is
    edge-aligned and identical for the two."""
    import hist
    from wums import boostHistHelpers as hh
    from wums import plot_tools

    h_m = load_matched_corrz(corrz_path)          # (absY, qT)
    h_r = resum_absY_qT(reference, q_lo, q_hi)     # (absY, qT)
    h_m, h_r = hh.rebinHistsToCommon([h_m, h_r], "absY")
    h_m, h_r = hh.rebinHistsToCommon([h_m, h_r], "qT")

    absY_edges = np.asarray(h_m.axes["absY"].edges, float)
    cut_edge = absY_edges[absY_edges <= Y_cut + 1e-9].max()  # largest edge <= Y_cut
    absYc = h_m.axes["absY"].centers
    ymask = absYc < cut_edge
    qT_edges = np.asarray(h_m.axes["qT"].edges, float)
    keep = qT_edges <= qT_max + 1e-9
    qe = qT_edges[keep]

    def to_hist(h):
        v = h.project("qT", "absY").values(flow=False)  # (qT, absY)
        sig = v[:, ymask].sum(axis=1)[: qe.size - 1]
        hh1 = hist.Hist(hist.axis.Variable(qe, name="qT"), storage=hist.storage.Double())
        hh1.view()[...] = sig
        return hh1

    hr, hm = to_hist(h_r), to_hist(h_m)
    fig = plot_tools.makePlotWithRatioToRef(
        [hr, hm],
        labels=["SCETlib resummed (our validated reco)", "scetlib + DYTurbo matched"],
        colors=["black", "#2ca02c"],
        xlabel=r"$q_T$ [GeV]", ylabel=r"$\sigma$/bin",
        rlabel=["matched / resummed"], rrange=[[0.5, 1.5]],
        binwnorm=1.0, nlegcols=1, grid=True, yerr=False, logy=True,
        plot_title=(f"Effect of DYTurbo matching: resummed vs matched, |Y|<{cut_edge:g}, "
                    f"Z Q-bin (native qT)"),
    )
    from wremnants.postprocessing.scetlib_np import plot_output

    outdir, basename = plot_output.split_outpath(out_path)
    plot_output.save_plot(outdir, basename, fig=fig, args=args)
    print(f"[plot] wrote {outdir}/{basename}.png(.pdf) + {basename}.log  (|Y|<{cut_edge:g})")


def native_matched_model_plot(out_path, corrz_path, reference, q_lo, q_hi, Y_cut,
                              qt_cutoff, qT_max, args=None):
    """Native validation that resummed + DYTurbo = CorrZ matched.

    Direct-hist approach (robust, on a common absY binning): resummed reference
    (= our reconstruction, validated to 0.08% in native_qT.png) vs CorrZ matched,
    plus (resummed + DYTurbo nonsingular). The green resummed-only curve shows the
    fixed-order gap; the red (resummed + DYTurbo) recovers the matched. Summed over
    |Y|<~Y_cut (edge-aligned), Z Q-bin, native qT."""
    import hist
    from wums import boostHistHelpers as hh
    from wums import plot_tools

    h_m = load_matched_corrz(corrz_path)            # (absY, qT) — CorrZ matched
    h_r = resum_absY_qT(reference, q_lo, q_hi)      # (absY, qT) — resummed reference
    h_m, h_r = hh.rebinHistsToCommon([h_m, h_r], "absY")
    h_m, h_r = hh.rebinHistsToCommon([h_m, h_r], "qT")

    absY_edges = np.asarray(h_m.axes["absY"].edges, float)
    cut_edge = absY_edges[absY_edges <= Y_cut + 1e-9].max()
    ym = np.asarray(h_m.axes["absY"].centers, float) < cut_edge
    qT_edges = np.asarray(h_m.axes["qT"].edges, float)
    qT_cen = np.asarray(h_m.axes["qT"].centers, float)
    vmatched = h_m.project("qT", "absY").values(flow=False)[:, ym].sum(axis=1)
    vresum = h_r.project("qT", "absY").values(flow=False)[:, ym].sum(axis=1)
    nonsing = vmatched - vresum
    nonsing[qT_cen < qt_cutoff] = 0.0
    vresum_plus = vresum + nonsing                  # resummed + DYTurbo

    keep = qT_edges <= qT_max + 1e-9
    qe = qT_edges[keep]
    n = qe.size - 1

    def H(v):
        h1 = hist.Hist(hist.axis.Variable(qe, name="qT"), storage=hist.storage.Double())
        h1.view()[...] = v[:n]
        return h1

    fig = plot_tools.makePlotWithRatioToRef(
        [H(vmatched), H(vresum), H(vresum_plus)],
        labels=["CorrZ matched (truth)", "resummed (our reco, validated)",
                "resummed + DYTurbo"],
        colors=["black", "#2ca02c", "red"],
        xlabel=r"$q_T$ [GeV]", ylabel=r"$\sigma$/bin",
        rlabel=["/ CorrZ matched"], rrange=[[0.5, 1.5]],
        binwnorm=1.0, nlegcols=1, grid=True, yerr=False, logy=True,
        plot_title=f"Native: resummed (+DYTurbo) vs CorrZ matched, |Y|<{cut_edge:g}, Z Q-bin",
    )
    from wremnants.postprocessing.scetlib_np import plot_output

    outdir, basename = plot_output.split_outpath(out_path)
    plot_output.save_plot(outdir, basename, fig=fig, args=args)
    m = qT_cen[:n] <= qT_max
    with np.errstate(divide="ignore", invalid="ignore"):
        rg = (vresum / vmatched)[:n][m]
        rr = (vresum_plus / vmatched)[:n][m]
    print(f"[native-matched] |Y|<{cut_edge:g}: resummed-only/CorrZ median {np.nanmedian(rg):.3f} "
          f"(range {np.nanmin(rg):.3f}–{np.nanmax(rg):.3f}); "
          f"resummed+DYTurbo/CorrZ median {np.nanmedian(rr):.4f}, "
          f"max|r-1| {np.nanmax(np.abs(rr-1))*100:.3f}%")
    print(f"[plot] wrote {outdir}/{basename}.png(.pdf) + {basename}.log")


def plot_native_qT(out_path, qT_edges, curves, Y_cut, qT_max, title, args=None):
    """curves: list of dicts {label, qT_arr, ref_arr, fac_arr, Y_arr, color}.
    Reference (curve[0]'s ref_arr) + each curve's factorized side as σ(qT)
    summed over |Y|<Y_cut, via wums.plot_tools.makePlotWithRatioToRef (ratio to
    reference, binwnorm=1.0)."""
    from wums import plot_tools

    h_ref = _spectrum_hist(qT_edges, curves[0]["qT_arr"], curves[0]["ref_arr"],
                           curves[0]["Y_arr"], Y_cut, qT_max)
    hists = [h_ref]
    labels = ["(1) SCETlib reference"]
    colors = ["black"]
    for c in curves:
        hists.append(_spectrum_hist(qT_edges, c["qT_arr"], c["fac_arr"],
                                    c["Y_arr"], Y_cut, qT_max))
        labels.append(c["label"])
        colors.append(c["color"])

    fig = plot_tools.makePlotWithRatioToRef(
        hists, labels=labels, colors=colors,
        xlabel=r"$q_T$ [GeV]", ylabel=r"$\sigma$/bin",
        rlabel=["factorized / ref"], rrange=[[0.97, 1.03]],
        binwnorm=1.0, nlegcols=1, grid=True, yerr=False, logy=True,
        plot_title=title,
    )
    from wremnants.postprocessing.scetlib_np import plot_output

    outdir, basename = plot_output.split_outpath(out_path)
    plot_output.save_plot(outdir, basename, fig=fig, args=args)
    print(f"[plot] wrote {outdir}/{basename}.png(.pdf) + {basename}.log")


def summarize_ratios(label, ratios, Y_arr, qT_arr, qT_cut, Y_cut):
    print(f"\n[{label}] N matched (Y,qT) cells = {len(ratios)}")
    print(f"  full set:              median = {np.nanmedian(ratios):.6f}   "
          f"std = {np.nanstd(ratios):.3e}")
    core = (np.abs(Y_arr) < Y_cut) & (qT_arr < qT_cut)
    rr = ratios[core]
    rr = rr[np.isfinite(rr)]
    print(f"  resum core |Y|<{Y_cut:g}, qT<{qT_cut:g}:  N = {core.sum()}")
    if rr.size:
        print(f"    median ratio       = {np.median(rr):.6f}   (target 1.0)")
        print(f"    |ratio-1| median   = {np.median(np.abs(rr-1))*100:.4f} %")
        print(f"    |ratio-1| 95%ile   = {np.percentile(np.abs(rr-1), 95)*100:.4f} %")
        print(f"    |ratio-1| max      = {np.max(np.abs(rr-1))*100:.4f} %")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--btgrid", default=_default_btgrid_dir(), help="SCETlib bT-grid directory")
    p.add_argument("--reference", default=FRANKS_REF,
                   help=f"SCETlib spectrum combined.pkl (#1 ref). lambda6-fine: {LAMBDA6_REF}")
    p.add_argument("--np-json", default=None,
                   help="JSON {eff_params, gnu_params} NP override (default: reference config; "
                        "with --paramodel the model's λ_central takes priority)")
    p.add_argument("--q-lo", type=float, default=Q_LO)
    p.add_argument("--q-hi", type=float, default=Q_HI)
    p.add_argument("--q-method", choices=("arctan_Q2", "simpson", "trapz"), default="arctan_Q2")
    p.add_argument("--qT-cut", type=float, default=30.0)
    p.add_argument("--Y-cut", type=float, default=1.0)
    # ParamModel (curve 3)
    p.add_argument("--paramodel", action="store_true",
                   help="also build the ParamModel and add its native σ(Y,qT) as curve 3")
    p.add_argument("--datacard", required=True, help="fit-input hdf5 (λ_central + binning)")
    p.add_argument("--signal-proc", default=SIGNAL_PROC)
    p.add_argument("--plot-out",
                   default="/home/submit/lavezzo/public_html/alphaS/260603_native_validation/native_qT.png",
                   help="output path for the σ(qT) overlay+ratio plot ('' to skip)")
    p.add_argument("--plot-qT-max", type=float, default=50.0,
                   help="max qT (GeV) shown in the plot (default 50)")
    p.add_argument("--matched", action="store_true",
                   help="also make a resummed-vs-(scetlib+DYTurbo matched) σ(qT) plot "
                        "to show the effect of the fixed-order matching")
    p.add_argument("--corrz", default=CORRZ,
                   help="CorrZ.pkl.lz4 holding the matched absolute spectrum")
    p.add_argument("--fo-sing", default=FO_SING,
                   help="SCETlib singular fixed-order combined.pkl (for σ_ns)")
    p.add_argument("--dyturbo", default=DYTURBO,
                   help="DYTurbo fixed-order results txt ('{scale}' -> mur1-muf1)")
    args = p.parse_args(argv)

    print("=" * 74)
    print("NATIVE-binning validation: (1) SCETlib reference  vs  (2) numpy factorize"
          + ("  vs  (3) ParamModel" if args.paramodel else ""))
    print("  resummed only, reference's native fine binning, Z Q-bin, NO projection")
    print("=" * 74)
    print(f"  btgrid    : {args.btgrid}")
    print(f"  reference : {args.reference}")

    cli = _import_cli_helpers()

    # ---- (1) reference, read with the WRemnants reader, verified == #1 ----
    print("\n[ref] reading reference two ways and asserting equality …")
    ref = factorize.load_spectrum_reference(args.reference)  # cmd_validate reader
    h_load = ref["hist"]
    h_read = input_tools.read_scetlib_hist(args.reference, nonsing="none")  # WRemnants reader
    v_load, v_read = h_load.values(flow=False), h_read.values(flow=False)
    if v_load.shape != v_read.shape:
        raise RuntimeError(f"reference shape mismatch: load {v_load.shape} vs read {v_read.shape}")
    dmax = float(np.nanmax(np.abs(v_load - v_read)))
    print(f"      load_spectrum_reference vs read_scetlib_hist: max|diff| = {dmax:.3e}  "
          f"({'IDENTICAL' if dmax == 0.0 else 'DIFFER'})")
    if cli._is_pointwise_reference(h_read):
        raise RuntimeError("pointwise reference — this script handles bins-mode only")
    ref_spectra, ref_widths = cli._extract_reference_spectrum(h_read)
    print(f"      {len(ref_spectra)} reference (Q,Y,qT) bin centres (bin-integrated)")
    ref_NP = _np_params_from_config(ref["config"])
    print(f"      reference NP: F_eff={ref_NP[0]}")

    # ---- ParamModel (curve 3) built first if requested: its λ_central sets the NP ----
    model = None
    if args.paramodel:
        from rabbit.inputdata import FitInputData
        from wremnants.postprocessing.scetlib_np.param_model import SCETlibNPParamModel
        print("\n[model] constructing SCETlibNPParamModel (runs btgrid integral at λ_central) …")
        t0 = time.time()
        indata = FitInputData(args.datacard)
        model = SCETlibNPParamModel(
            indata, btgrid_dir=args.btgrid,
            signal_proc=args.signal_proc, Q_lo=args.q_lo, Q_hi=args.q_hi,
            nonsingular_fo_sing=args.fo_sing, nonsingular_dyturbo=args.dyturbo,
        )
        ns = model.sigma_ns.numpy()
        sg = model.sigma_gen_central.numpy()
        print(f"      matched σ_gen (nonsingular always included); "
              f"Σ|σ_ns|/Σσ_gen = {np.abs(ns).sum()/sg.sum():.3f}")
        print(f"      constructed in {time.time()-t0:.1f}s; "
              f"native σ(Y,qT) shape {tuple(model.sigma_YqT_central.shape)}")
        print(f"      model λ_central: F_eff={model.eff_central}")
        print(f"                       γ_ν^NP={model.gnu_central}")

    # ---- NP for the numpy factorize (curve 2): model λ_central > json > ref config ----
    if model is not None:
        eff_params, gnu_params = dict(model.eff_central), dict(model.gnu_central)
        np_src = "model λ_central"
    elif args.np_json:
        with open(args.np_json) as f:
            j = json.load(f)
        eff_params, gnu_params = j["eff_params"], j["gnu_params"]
        np_src = f"json {args.np_json}"
    else:
        eff_params, gnu_params = ref_NP
        np_src = "reference pickle config"
    print(f"\n[NP] numpy factorize reconstructed at: {np_src}")
    print(f"     F_eff       : {eff_params}")
    print(f"     gamma_nu^NP : {gnu_params}")

    # ---- (2) numpy factorize reconstruction + Q-integration ----
    print(f"\n[load] btgrid shards from {args.btgrid}")
    grid, Q_grid, Y_grid, qT_grid = load_btgrid_dense(args.btgrid)
    print(f"[axes] NQ={Q_grid.size} NY={Y_grid.size} NqT={qT_grid.size}")
    t0 = time.time()
    sigma_QYqT = factorize.reconstruct_grid_QYqT(
        Q_grid, Y_grid, qT_grid, bT=grid["bT"], I_pert=grid["I_pert"][0],
        C_nu=grid["C_nu"][0], b_bar=grid["b_bar"],
        eff_params=eff_params, gnu_params=gnu_params,
    )
    sigma_YqT_np = factorize.integrate_over_Q(sigma_QYqT, Q_grid, args.q_lo, args.q_hi,
                                              method=args.q_method)
    print(f"[recon] numpy factorize reconstructed + Q-integrated in {time.time()-t0:.1f}s")

    Q_ref = min(sorted({k[0] for k in ref_spectra}),
                key=lambda Q: abs(Q - 0.5 * (args.q_lo + args.q_hi)))
    print(f"[match] Z Q-bin centre = {Q_ref}")

    r2, Y2, qT2, ref2, fac2, ns2, npn2 = native_ratios(
        sigma_YqT_np, ref_spectra, ref_widths, Y_grid, qT_grid, Q_ref)
    print(f"[bin-int] (2) Simpson/(Y,qT): {ns2}  | point×width: {npn2}")
    summarize_ratios("2 numpy factorize / 1 reference", r2, Y2, qT2, args.qT_cut, args.Y_cut)

    curves = [dict(label="(2) numpy factorize", qT_arr=qT2, ref_arr=ref2, fac_arr=fac2,
                   Y_arr=Y2, color="red")]

    # ---- (3) ParamModel native σ(Y,qT) ----
    if model is not None:
        sigma_YqT_tf = model.sigma_YqT_central.numpy()  # (NY, NqT), 0-filled missing cells
        # grids must align with the numpy side (both sorted-unique from same btgrid)
        Yg_m = np.asarray(model.Y_unique, float)
        qTg_m = np.asarray(model.qT_unique, float)
        same_grid = (Yg_m.shape == Y_grid.shape and np.allclose(Yg_m, Y_grid)
                     and qTg_m.shape == qT_grid.shape and np.allclose(qTg_m, qT_grid))
        print(f"\n[model] native grid matches numpy grid: {same_grid}")

        # (3) vs (1): Simpson per bin against the reference
        r3, Y3, qT3, ref3, fac3, ns3, npn3 = native_ratios(
            sigma_YqT_tf, ref_spectra, ref_widths, Yg_m, qTg_m, Q_ref)
        print(f"[bin-int] (3) Simpson/(Y,qT): {ns3}  | point×width: {npn3}")
        summarize_ratios("3 ParamModel / 1 reference", r3, Y3, qT3, args.qT_cut, args.Y_cut)
        curves.append(dict(label="(3) ParamModel", qT_arr=qT3, ref_arr=ref3, fac_arr=fac3,
                           Y_arr=Y3, color="#1f77b4"))

        # (3) vs (2): element-wise parity on the native grid (resum core, finite & nonzero)
        if same_grid:
            both = np.isfinite(sigma_YqT_np) & (sigma_YqT_tf != 0.0)
            core = (np.abs(Yg_m)[:, None] < args.Y_cut) & (qTg_m[None, :] < args.qT_cut)
            m = both & core
            rel = np.abs(sigma_YqT_tf[m] / sigma_YqT_np[m] - 1.0)
            print(f"\n[3 vs 2] element-wise parity (ParamModel TF vs numpy factorize, "
                  f"native, resum core, N={int(m.sum())}):")
            print(f"    max |TF/np - 1|    = {np.nanmax(rel):.3e}")
            print(f"    median |TF/np - 1| = {np.nanmedian(rel):.3e}")

    # ---- plot: σ(qT) overlay + ratio-to-reference (|Y|<Y_cut), via wums ----
    if args.plot_out:
        qT_edges = np.asarray(h_read.axes["qT"].edges, dtype=float)
        title = (f"SCETlib resummed σ(qT), native binning, |Y|<{args.Y_cut:g}, "
                 f"Z Q-bin — reference vs factorization")
        plot_native_qT(args.plot_out, qT_edges, curves, args.Y_cut, args.plot_qT_max, title,
                       args=args)

    # ---- matched (resum + DYTurbo) overlay: effect of the fixed-order matching ----
    if args.matched and args.plot_out:
        base, ext = os.path.splitext(args.plot_out)
        matched_out = base + "_matched" + (ext or ".png")
        matched_vs_resum_plot(matched_out, args.corrz, args.reference,
                              args.q_lo, args.q_hi, args.Y_cut, args.plot_qT_max,
                              args=args)
        # native validation: model's own resummed + DYTurbo vs CorrZ matched
        mm_out = base + "_model_matched" + (ext or ".png")
        native_matched_model_plot(mm_out, args.corrz, args.reference, args.q_lo,
                                   args.q_hi, args.Y_cut, qt_cutoff=1.0,
                                   qT_max=args.plot_qT_max, args=args)

    # ---- GEN-grid decomposition: model σ_gen vs resum-ref vs CorrZ-matched ----
    # Rides on the gen-grid projection (the known rebinning issue, deferred). Printed
    # to localize whether any discrepancy is the model's W-projection of the resum
    # (model vs resum-ref) or the nonsingular (resum-ref vs CorrZ).
    if model is not None:
        model_gen = model.sigma_gen_central.numpy()              # matched (always)
        # resum-only derived from the matched model: σ_resum = σ_gen − σ_ns
        model_resum_gen = model_gen - model.sigma_ns.numpy()
        resum_gen = resum_ref_to_gen(args.reference, model._gen_axes_meta, args.q_lo, args.q_hi)
        corrz_gen = corrz_matched_to_gen(args.corrz, model._gen_axes_meta)
        m_ptv, mr_ptv, r_ptv, c_ptv = (
            a.sum(axis=1) for a in (model_gen, model_resum_gen, resum_gen, corrz_gen)
        )
        tag = "matched"
        print(f"\n[gen] σ_gen on the gen grid (ptVGen-summed); model is matched, "
              f"resum-only derived as matched − σ_ns:")
        print(f"    Σ model(matched) = {m_ptv.sum():.5g}   Σ model(resum-only) = {mr_ptv.sum():.5g}   "
              f"Σ resum-ref = {r_ptv.sum():.5g}   Σ CorrZ-matched = {c_ptv.sum():.5g}")
        print(f"    Σ ratios:  model(resum)/resum-ref = {mr_ptv.sum()/r_ptv.sum():.4f}   "
              f"resum-ref/CorrZ = {r_ptv.sum()/c_ptv.sum():.4f}   "
              f"model(matched)/CorrZ = {m_ptv.sum()/c_ptv.sum():.4f}")
        with np.errstate(divide="ignore", invalid="ignore"):
            print("    per-ptVGen model(resum)/resum-ref:",
                  np.array2string(mr_ptv / r_ptv, precision=3, max_line_width=200))
            print("    per-ptVGen resum-ref/CorrZ:",
                  np.array2string(r_ptv / c_ptv, precision=3, max_line_width=200))
        if args.plot_out:
            import hist
            from wums import plot_tools

            ptV_edges = np.asarray(model._gen_axes_meta[0][1], float)

            def _ptv_hist(v):
                h1 = hist.Hist(hist.axis.Variable(ptV_edges, name="ptVGen"),
                               storage=hist.storage.Double())
                h1.view()[...] = v
                return h1
            base, ext = os.path.splitext(args.plot_out)
            gen_out = base + "_gen_matched" + (ext or ".png")
            fig = plot_tools.makePlotWithRatioToRef(
                [_ptv_hist(c_ptv), _ptv_hist(r_ptv), _ptv_hist(m_ptv)],
                labels=["CorrZ matched (truth)", "resum-ref (bin-sum)",
                        f"ParamModel σ_gen ({tag})"],
                colors=["black", "#2ca02c", "red"],
                xlabel=r"$p_T^{V,gen}$ [GeV]", ylabel=r"$\sigma$/bin",
                rlabel=["/ CorrZ"], rrange=[[0.9, 1.1]],
                binwnorm=1.0, nlegcols=1, grid=True, yerr=False, logy=True,
                plot_title="GEN grid: ParamModel matched σ_gen (resum⊕DYTurbo) vs CorrZ matched",
            )
            from wremnants.postprocessing.scetlib_np import plot_output

            outdir, basename = plot_output.split_outpath(gen_out)
            plot_output.save_plot(outdir, basename, fig=fig, args=args)
            print(f"[plot] wrote {outdir}/{basename}.png(.pdf) + {basename}.log")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
