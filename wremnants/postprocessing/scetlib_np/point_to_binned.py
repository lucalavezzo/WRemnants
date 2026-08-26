#!/usr/bin/env python3
"""Convert a POINT-spectrum SCETlib pickle -> a binned ``{hist: hist.Hist}``
pickle (NAMED StrCategory vars axis) that ``make_theory_corr`` reads (via
``input_tools.read_scetlib_hist`` / ``read_matched_scetlib_dyturbo_hist``).

Does no integration of its own; reuses the validated NP-param-model machinery
(:mod:`btgrid_integrate`):

  * Q  : arctan-Q^2 Simpson over the Z mass window   (``q_integrate_weights``)
         -> integrated to ONE Q bin [qlo, qhi].
  * Y  : 3-point Simpson per experimental |Y| bin    (``rebin_weights``)
  * qT : 3-point Simpson per experimental qT bin     (``rebin_weights``)

The point run must be on the btgrid grid = union of the experimental bin EDGES
and CENTERS, so each bin holds exactly (edge, center, edge) = 3 source samples ->
3-point Simpson per bin (the <0.05% rebin the btgrid runcard targets). Bin edges
are recovered as every-other grid point (grid = interleaved edges+centers); a
midpoint check guards that ``[0::2]`` are the edges.

Two benign integrity items of a `calculation_piece = sing` point run:
  * qT = 0 is NaN (differential spectrum ill-defined there; physical limit 0).
    Zeroed so the first qT bin's Simpson is well-defined.
  * Negative ``sing``-only sigma (nonsingular-dominated region) are PHYSICAL and
    pass through unchanged; ``make_theory_corr``'s matching (DYTurbo - FO_sing)
    corrects them. Simpson is linear, so rebinning negatives is fine.

Y is kept SIGNED (the production resummed pkl is signed Y); make_theory_corr does
the |Y| fold, as for the real input.

Usage (inside the WRemnants apptainer with setup.sh sourced), from WREM_BASE:
  python3 wremnants/postprocessing/scetlib_np/point_to_binned.py \\
      <point_spectrum.pkl> -o <out_binned.pkl>
"""
import argparse
import pickle

import numpy as np

from wremnants.postprocessing.scetlib_np import btgrid_integrate as fz


def edges_from_grid(g, name, sub=2):
    """Recover experimental bin edges from an (edges + interior) source grid.

    ``sub`` = sub-intervals per experimental bin (``sub+1`` samples per bin,
    edges shared between neighbours): ``sub=2`` = default edges+1-centre layout
    (3-point Simpson per bin); ``sub=4`` = 2x-refined edges+3-interior layout
    (5-point composite Simpson); etc. Bin edges are then every ``sub``-th grid
    point. rebin_weights does composite Simpson over the samples in each bin, so
    finer ``sub`` lowers the rebin floor (~h^4).
    """
    g = np.asarray(sorted(set(g)), dtype=np.float64)
    if (g.size - 1) % sub != 0:
        raise SystemExit(
            f"{name}: {g.size} grid points incompatible with --subsample {sub} "
            f"(need (N_points - 1) divisible by {sub})"
        )
    edges = g[0::sub]
    # For even sub the bin centre is a grid point; verify edges really are edges
    # (uniform sub-sampling within each bin -> midpoint == centre sample).
    if sub % 2 == 0:
        centers = g[sub // 2::sub]
        mid = 0.5 * (edges[:-1] + edges[1:])
        if not np.allclose(mid, centers, atol=1e-6):
            raise SystemExit(
                f"{name}: grid[0::{sub}] are not bin edges (edge-midpoints != bin "
                f"centres). Wrong --subsample, or not an evenly-refined grid?"
            )
    return edges


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile", help="point-spectrum pkl: {spectra:{var:{(Q,Y,qT,lep):sigma}}}")
    ap.add_argument("-o", "--outfile", required=True)
    ap.add_argument("--qlo", type=float, default=60.0, help="mass-window low edge")
    ap.add_argument("--qhi", type=float, default=120.0, help="mass-window high edge")
    ap.add_argument("--no-selfcheck", action="store_true",
                    help="skip parsing the output back through read_scetlib_hist")
    ap.add_argument("--subsample", type=int, default=2,
                    help="sub-intervals per experimental Y/qT bin (sub+1 samples "
                         "per bin): 2 = edges+centre (3-pt Simpson, default); "
                         "4 = edges+3 interior (5-pt, 2x-refined); etc.")
    args = ap.parse_args()

    with open(args.infile, "rb") as f:
        d = pickle.load(f)
    sp = d["spectra"]
    nvars = len(sp)
    names = ({k: d["vars"][k].get("name", str(k)) for k in d["vars"]}
             if "vars" in d else {})

    # --- axes from the point grid ---
    keys0 = list(sp[0])
    Qg = np.array(sorted({k[0] for k in keys0}), dtype=np.float64)
    Yg = np.array(sorted({k[1] for k in keys0}), dtype=np.float64)
    Tg = np.array(sorted({k[2] for k in keys0}), dtype=np.float64)
    NQ, NY, NT = Qg.size, Yg.size, Tg.size
    print(f"point grid: {NQ} Q x {NY} Y x {NT} qT = {NQ*NY*NT} pts, {nvars} vars")

    Yedges = edges_from_grid(Yg, "Y", sub=args.subsample)
    Tedges = edges_from_grid(Tg, "qT", sub=args.subsample)
    print(f"target: 1 Q bin [{args.qlo},{args.qhi}] x {Yedges.size-1} Y x "
          f"{Tedges.size-1} qT")

    # --- validated param-model Simpson weights ---
    wQ = fz.q_integrate_weights(Qg, args.qlo, args.qhi)   # (NQ,)
    WY = fz.rebin_weights(Yg, Yedges, name="Y")           # (NY_t, NY)
    WT = fz.rebin_weights(Tg, Tedges, name="qT")          # (NT_t, NT)

    # --- dense fill: (nvars, NQ, NY, NT); zero non-finite (only the qT=0 edge) ---
    qpos = {q: i for i, q in enumerate(Qg)}
    ypos = {y: i for i, y in enumerate(Yg)}
    tpos = {t: i for i, t in enumerate(Tg)}
    Qi = np.fromiter((qpos[k[0]] for k in keys0), int, len(keys0))
    Yi = np.fromiter((ypos[k[1]] for k in keys0), int, len(keys0))
    Ti = np.fromiter((tpos[k[2]] for k in keys0), int, len(keys0))
    H = np.zeros((nvars, NQ, NY, NT), dtype=np.float64)
    for v in range(nvars):
        H[v, Qi, Yi, Ti] = np.fromiter((sp[v][k] for k in keys0), float, len(keys0))
    nbad = int((~np.isfinite(H)).sum())
    H = np.nan_to_num(H, nan=0.0, posinf=0.0, neginf=0.0)
    print(f"filled; zeroed {nbad} non-finite samples (expected: the qT=0 edge)")

    # --- Q-integrate, then Simpson-rebin Y and qT (all linear; param-model order) ---
    A = np.einsum("q,vqyt->vyt", wQ, H)        # (v, NY, NT)
    A = np.einsum("ry,vyt->vrt", WY, A)        # (v, NY_t, NT)
    A = np.einsum("st,vrt->vrs", WT, A)        # (v, NY_t, NT_t)
    hist = A[:, None, :, :]                    # (vars, Q=1, Y, qT)
    print(f"output hist (vars,Q,Y,qT) = {hist.shape}")

    # Emit a hist.Hist with a NAMED StrCategory vars axis (as production scetlib
    # hists have). read_scetlib_hist uses the hist as-is so names survive;
    # read_matched_scetlib_hist then maps each (non-scale) lambda var to the
    # central nonsingular via string ops -- an Integer vars axis would TypeError.
    # Named vars also let feedRabbitSigmaUL pick truth/templates by name. Match
    # the production scetlib hist axes EXACTLY: Q/Y/qT carry under+overflow
    # (flow=True, the hist default), vars is a StrCategory with overflow. flow=False
    # here left the real bins aligned but the flow bins absent, so
    # read_matched_scetlib_hist's addHists couldn't broadcast against the
    # (flow=True) nonsingular -> shape (1,82,70,..) vs (3,84,72,..).
    import hist as _hist
    var_labels = [names.get(v, str(v)) for v in range(nvars)]
    # WEIGHT storage (value + variance), NOT Double. make_theory_corr propagates
    # our storage type to the matched `_hist`: a Double scetlib -> Double `_hist`,
    # whose .variances() is None in newer boost_histogram, which the rabbit tensor
    # writer coerces to a single NaN and rejects ("1 NaN ... in variances"). The
    # standard scetlib_dyturbo `_hist` is Weight. SCETlib is analytic resummed with
    # NO MC statistics, so set the variance to exactly 0; the matched `_hist`'s
    # MC-stat variance then comes purely from the DYTurbo/MiNNLO pieces, identical
    # to the standard path. Also fixes the pseudodata, read from the same file.
    h = _hist.Hist(
        _hist.axis.Variable([float(args.qlo), float(args.qhi)], name="Q"),
        _hist.axis.Variable(Yedges, name="Y"),
        _hist.axis.Variable(Tedges, name="qT"),
        _hist.axis.StrCategory(var_labels, name="vars"),
        storage=_hist.storage.Weight(),
    )
    h.view()["value"][...] = np.moveaxis(hist, 0, -1)  # (Q,Y,qT,vars); variance stays 0
    # Carry the SCETlib config/meta from the point pkl: make_theory_corr stores
    # them as provenance and `get_scetlib_config` hard-requires "config".
    out = {"hist": h, "var_names": names,
           "config": d.get("config"), "meta_data": d.get("meta_data")}
    with open(args.outfile, "wb") as f:
        pickle.dump(out, f)
    print(f"wrote {args.outfile}  (vars: {var_labels[0]} ... {var_labels[-1]})")

    # sanity: per-var window-integrated xsec (central vs joint should differ)
    tot = hist.sum(axis=(1, 2, 3))
    print("per-var integrated sigma (var: name = total):")
    for v in range(nvars):
        print(f"  {v}: {names.get(v,'?'):<40} {tot[v]:.5g}")

    if not args.no_selfcheck:
        try:
            from wremnants.utilities.io_tools.input_tools import read_scetlib_hist
            h = read_scetlib_hist(args.outfile)
            print(f"self-check OK: read_scetlib_hist -> axes={list(h.axes.name)}, "
                  f"shape={h.shape}")
        except Exception as e:
            print(f"[warn] self-check via read_scetlib_hist failed: {e!r}")


if __name__ == "__main__":
    main()
