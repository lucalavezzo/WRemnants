#!/usr/bin/env python3
"""Compare a CORRECTED and an UNCORRECTED datacard, to justify moving the
MiNNLO->SCETlib theory correction out of the histmaker and into the param model.

The correction is applied by the histmaker as a per-event weight looked up from
the CorrZ histogram -- a BIN LOOKUP, not an interpolation
(``correctionsTensor_helper.py``: "takes a histogram and returns what is in the
bin"). So the event-weight route is itself piecewise-constant and

    Sum_events w(g_event)  ==  Sum_g w(g) * R_raw(b, g)

is an IDENTITY at matched binning: the two routes are the same object written
differently. What differs in practice is grid coarseness -- the card's response
gen grid (ptVGen 21 x absYVGen 10) against the correction's own grid
(qT 70 x absY 17) -- so the residual is the WITHIN-GEN-BIN variation of the
correction, expected worst at low qT where it moves fastest.

Two tests, neither needing a cache, SCETlib, or the corr file:

T1  row-sum audit, per card
      R_rowsum(b) / norm_signal(b),   R_rowsum = Sum_g R_raw(b,g) == R @ N_gen
    The theory cancels identically, so this tests ONLY the plumbing: reco
    marginalization, cropping, axis order, channel slicing, and how much gen
    truth leaks outside the response grid. Must give the same profile on both
    cards, since the leakage is a property of R's coverage, not of the weights.

T2  granularity -- THE DECISIVE TEST
      cbar(g) = N_gen_corr(g) / N_gen_unc(g)          the coarse-grained correction
      M(b)    = [Sum_g R_raw(b,g) cbar(g)] / [Sum_g R_raw(b,g)]   matrix route
      E(b)    = norm_corr(b) / norm_unc(b)                        event route
      M/E - 1 IS the event-vs-matrix residual, with the theory, the pb/fb
      factor, the Y convention and sigma_SCETlib all cancelled out.
    Decision rule: adopt the matrix route if the yield-weighted |M/E - 1| is
    below ~0.1% and no sensitivity-carrying bin exceeds ~1%. Otherwise the
    response's gen binning needs refining toward the correction's grid.

NB the two cards must differ ONLY by the correction -- same reco axes, bin
count, processes and zero pattern -- or E(b) is not the event route and T2 is
void. That is asserted, not assumed.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
)

from wremnants.postprocessing.scetlib_ad.response import (  # noqa: E402
    R_info_from_auxiliary,
    crop_R_to_fit,
    marginalize_R_reco,
)

AXIS_LABELS = {"ptll": r"$p_{T}^{\ell\ell}$ [GeV]", "yll": r"$y^{\ell\ell}$"}


def load_card(path, signal_proc):
    """R (cropped to the fit channel), N_gen, the signal norm column, and axes."""
    from rabbit.inputdata import FitInputData

    ind = FitInputData(path)
    info = R_info_from_auxiliary(ind)

    non_masked = [
        (n, i) for n, i in ind.channel_info.items() if not i.get("masked", False)
    ]
    if len(non_masked) != 1:
        raise SystemExit(
            f"{path}: expected one non-masked channel, got {len(non_masked)}"
        )
    chan, cinfo = non_masked[0]
    fit_axes = [(a.name, np.asarray(a.edges, float)) for a in cinfo["axes"]]

    R_full, R_axes = marginalize_R_reco(
        info["R"], info["reco_axes"], [n for n, _ in fit_axes]
    )
    R = crop_R_to_fit(R_full, R_axes, fit_axes)
    n_reco = int(np.prod(R.shape[: len(fit_axes)]))
    R = R.reshape(n_reco, -1)

    procs = [p.decode() if isinstance(p, bytes) else str(p) for p in ind.procs]
    if signal_proc not in procs:
        raise SystemExit(f"{path}: signal {signal_proc!r} not in {procs}")
    norm = np.asarray(
        ind.norm.numpy() if hasattr(ind.norm, "numpy") else ind.norm, dtype=float
    )
    # Slice the channel's own rows rather than assuming it starts at 0.
    lo = int(cinfo.get("start", 0) or 0)
    hi = int(cinfo.get("stop", lo + n_reco) or lo + n_reco)
    norm_sig = norm[lo:hi, procs.index(signal_proc)]

    return dict(
        path=path,
        R=R,
        N_gen=np.asarray(info["N_gen"], float).reshape(-1),
        norm_sig=norm_sig,
        fit_axes=fit_axes,
        gen_axes=info["gen_axes"],
        procs=procs,
        channel=chan,
        lumi=cinfo.get("lumi"),
    )


def summarize(name, ratio, weights, axes, shape, top_n=6):
    """Yield-weighted stats plus the worst bins, with their coordinates."""
    good = np.isfinite(ratio) & (weights > 0)
    r, w = ratio[good], weights[good]
    dev = np.abs(r - 1.0)
    ywm = float(np.average(dev, weights=w))
    print(
        f"  {name}: bins {r.size}  mean {r.mean():.6f}  "
        f"min {r.min():.6f}  max {r.max():.6f}\n"
        f"      YIELD-WEIGHTED mean|ratio-1| = {ywm:.6f}   "
        f"({100*ywm:.4f}%)   worst |dev| = {dev.max():.3e}"
    )
    idx = np.argsort(-dev)[:top_n]
    flat = np.where(good)[0]
    print(f"      worst {top_n} bins:")
    for j in idx:
        b = flat[j]
        coord = np.unravel_index(b, shape)
        loc = ", ".join(f"{nm}[{e[c]:g},{e[c+1]:g}]" for (nm, e), c in zip(axes, coord))
        print(f"         {loc:<44} ratio {r[j]:.6f}  yield {w[j]:.4g}")
    return ywm


def profile(name, num, den, axis_idx, axes, shape):
    """1-D profile: sum numerator and denominator over the OTHER axes, then divide.

    Not a mean of per-bin ratios -- that would weight the low-yield forward bins
    like the peak.
    """
    nm, edges = axes[axis_idx]
    other = tuple(i for i in range(len(shape)) if i != axis_idx)
    n = num.reshape(shape).sum(axis=other)
    d = den.reshape(shape).sum(axis=other)
    r = np.where(d > 0, n / np.where(d == 0, np.nan, d), np.nan)
    print(f"\n  {name} profiled in {nm}:")
    print(f"      {'bin':>16} {'ratio':>12} {'dev':>12}")
    for k in range(r.size):
        print(
            f"      [{edges[k]:6g},{edges[k+1]:6g}] {r[k]:12.6f} " f"{r[k]-1.0:+12.3e}"
        )
    return r


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--card-corrected", required=True)
    ap.add_argument("--card-uncorrected", required=True)
    ap.add_argument("--signal-proc", default="Zmumu")
    ap.add_argument("--profile-axis", default="ptll")
    args = ap.parse_args()

    print("loading cards ...", flush=True)
    C = load_card(args.card_corrected, args.signal_proc)
    U = load_card(args.card_uncorrected, args.signal_proc)

    # ---- preconditions: the cards must differ ONLY by the correction --------
    print("\n=== preconditions ===")
    errs = []
    if [n for n, _ in C["fit_axes"]] != [n for n, _ in U["fit_axes"]]:
        errs.append("fit axis NAMES differ")
    for (nc, ec), (nu, eu) in zip(C["fit_axes"], U["fit_axes"]):
        if ec.shape != eu.shape or not np.allclose(ec, eu):
            errs.append(f"fit axis {nc} EDGES differ")
    if C["procs"] != U["procs"]:
        errs.append(f"procs differ: {C['procs']} vs {U['procs']}")
    if C["R"].shape != U["R"].shape:
        errs.append(f"R shape {C['R'].shape} vs {U['R'].shape}")
    if C["norm_sig"].shape != U["norm_sig"].shape:
        errs.append("norm shape differs")
    zc, zu = C["norm_sig"] == 0, U["norm_sig"] == 0
    if not np.array_equal(zc, zu):
        errs.append(f"zero pattern differs in {int((zc ^ zu).sum())} bins")
    for e in errs:
        print(f"  FAIL  {e}")
    if errs:
        raise SystemExit("cards are not comparable; T2 would be meaningless")
    print(
        f"  OK  axes {[f'{n}({len(e)-1})' for n, e in C['fit_axes']]}, "
        f"R {C['R'].shape}, procs {C['procs']}"
    )
    print(f"  lumi: corrected {C['lumi']}, uncorrected {U['lumi']}")
    shape = tuple(len(e) - 1 for _, e in C["fit_axes"])

    # ---- T1: row-sum audit, per card ---------------------------------------
    print("\n=== T1  row-sum audit (theory cancels; plumbing only) ===")
    for tag, D in (("corrected  ", C), ("uncorrected", U)):
        rowsum = D["R"].sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = rowsum / np.where(D["norm_sig"] == 0, np.nan, D["norm_sig"])
        print(
            f"\n  [{tag}] sum R_rowsum={rowsum.sum():.8g}  "
            f"sum norm_sig={D['norm_sig'].sum():.8g}  "
            f"global={rowsum.sum()/D['norm_sig'].sum():.8f}"
        )
        summarize(tag, ratio, D["norm_sig"], D["fit_axes"], shape)

    # ---- T2: granularity ----------------------------------------------------
    print("\n=== T2  granularity: matrix route vs event route (DECISIVE) ===")
    ng_c, ng_u = C["N_gen"], U["N_gen"]
    pos = ng_u > 0
    if not pos.all():
        print(f"  note: {int((~pos).sum())} gen bins have N_gen_unc <= 0; excluded")
    cbar = np.ones_like(ng_u)
    cbar[pos] = ng_c[pos] / ng_u[pos]
    print(
        f"  cbar = N_gen_corr/N_gen_unc over {int(pos.sum())} gen bins: "
        f"min {cbar[pos].min():.6f}  max {cbar[pos].max():.6f}  "
        f"mean {cbar[pos].mean():.6f}"
    )
    print("  (that IS the correction, coarse-grained onto the response gen grid)")

    Ru = U["R"]
    num = Ru @ cbar
    den = Ru.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        M = num / np.where(den == 0, np.nan, den)
        E = C["norm_sig"] / np.where(U["norm_sig"] == 0, np.nan, U["norm_sig"])
        MoverE = M / np.where(E == 0, np.nan, E)
    print(f"\n  M (matrix) min/max: {np.nanmin(M):.6f} / {np.nanmax(M):.6f}")
    print(f"  E (event)  min/max: {np.nanmin(E):.6f} / {np.nanmax(E):.6f}")
    ywm = summarize("M/E", MoverE, U["norm_sig"], C["fit_axes"], shape)

    # R invariance: R = R_raw/N_gen must be identical between the two cards,
    # since a gen-level reweighting multiplies numerator and denominator alike.
    # The whole construction rests on this, so check it rather than assume it.
    print("\n=== R invariance (R = R_raw/N_gen must be card-independent) ===")
    both = (C["N_gen"] > 0) & (U["N_gen"] > 0)
    Pc = C["R"][:, both] / C["N_gen"][both][None, :]
    Pu = U["R"][:, both] / U["N_gen"][both][None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.abs(Pc / np.where(Pu == 0, np.nan, Pu) - 1.0)
    ok = np.isfinite(rel)
    print(
        f"  |P_corr/P_unc - 1| over {int(ok.sum())} entries: "
        f"max {rel[ok].max():.3e}  median {np.median(rel[ok]):.3e}"
    )
    if Pu[ok].sum() > 0:
        print(
            f"  response-weighted mean = " f"{np.average(rel[ok], weights=Pu[ok]):.3e}"
        )
    print(
        "  (big values in near-empty response entries are MC noise, not a "
        "violation; the weighted mean is the meaningful number)"
    )

    names = [n for n, _ in C["fit_axes"]]
    # Which axis drives the residual? The axis whose M/E profile carries the
    # most structure is the one whose gen binning is too coarse.
    for ai, nm in enumerate(names):
        pm = profile(f"M (matrix) [{nm}]", num, den, ai, C["fit_axes"], shape)
        pe = profile(
            f"E (event)  [{nm}]",
            C["norm_sig"],
            U["norm_sig"],
            ai,
            C["fit_axes"],
            shape,
        )
        d = pm / pe - 1.0
        print(
            f"\n  ==> M/E - 1 profiled in {nm}: "
            f"max|.| = {np.nanmax(np.abs(d)):.3e}, "
            f"rms = {np.sqrt(np.nanmean(d**2)):.3e}"
        )

    print("\n=== verdict ===")
    print(f"  yield-weighted |M/E - 1| = {100*ywm:.4f}%")
    if ywm < 1e-3:
        print("  -> below the 0.1% decision threshold: the matrix route reproduces")
        print("     the event route on this gen grid. Refining it is unnecessary.")
    else:
        print("  -> ABOVE 0.1%: the response gen grid is too coarse for the")
        print("     correction it has to carry. Consider refining it toward the")
        print("     correction's own grid (qT 70 x absY 17), weighing the extra")
        print("     statistical noise per gen bin in R.")


if __name__ == "__main__":
    main()
