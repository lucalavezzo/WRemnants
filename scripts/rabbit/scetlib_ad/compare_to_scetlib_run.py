#!/usr/bin/env python3
"""Validate the autodiff cache against a native SCETlib production run.

This is the apples-to-apples test of our configuration and quadrature: a
``calculation_piece = sing`` production pkl is the RESUMMED cross section,
bin-integrated, and replaying only our cache's compressed bin rules gives the
same object. No matching, no DYTurbo, no MiNNLO, no correction file in between --
so any disagreement is ours (runcard, quadrature, Q integration, rule
compression), not a difference between two predictions.

The reference is summed from its own fine bins onto our grid, which is exact
because both sides are bin-integrated; the script refuses to run unless our bin
edges really are a subset of the reference's.

    source <scetlib-cms>/setup.sh
    python scripts/rabbit/scetlib_ad/compare_to_scetlib_run.py \
        --conf <cache>.conf --cache <cache>.npz --reference <...>_combined.pkl
"""

import argparse
import os
import pickle
import sys

import numpy as np

sys.path.insert(
    0,
    os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ),
)

from wremnants.postprocessing.scetlib_ad.xsec_backend import (  # noqa: E402
    ScetlibADXsec,
)

TOL = 1e-9


def _sum_map(fine_edges, coarse_edges, what):
    """Index of the coarse bin each fine bin falls in, -1 if outside.

    Requires every coarse edge to be a fine edge, so the sum is exact rather
    than an interpolation.
    """
    missing = [e for e in coarse_edges if not np.any(np.abs(fine_edges - e) < TOL)]
    if missing:
        raise SystemExit(
            f"our {what} edges {missing} are not edges of the reference grid, so "
            f"the reference cannot be summed onto our bins exactly. Rebuild the "
            f"cache on a subset of the reference grid."
        )
    idx = np.full(fine_edges.size - 1, -1, dtype=np.int64)
    for i in range(fine_edges.size - 1):
        lo, hi = fine_edges[i], fine_edges[i + 1]
        for j in range(coarse_edges.size - 1):
            if lo >= coarse_edges[j] - TOL and hi <= coarse_edges[j + 1] + TOL:
                idx[i] = j
                break
    return idx


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--conf", required=True)
    ap.add_argument("--cache", required=True)
    ap.add_argument("--reference", required=True, help="SCETlib *_combined.pkl (sing)")
    ap.add_argument("--var", default="pdf0", help="reference vars entry to compare")
    ap.add_argument("--threads", type=int, default=0)
    args = ap.parse_args()

    with open(args.reference, "rb") as fh:
        ref = pickle.load(fh)
    h = ref["hist"]
    cfg = ref.get("config", {})
    piece = cfg.get("Calculation_settings", {}).get("calculation_piece")
    if piece != "sing":
        print(
            f"WARNING: the reference is calculation_piece={piece!r}, not 'sing'. "
            f"This script compares the RESUMMED piece; a matched reference would "
            f"also contain the nonsingular and will not agree.",
            flush=True,
        )

    core = ScetlibADXsec(args.conf, args.cache, threads=args.threads)
    ours = core.resummed_only(core.anchor)

    # --- cross-check the configuration we are validating against
    print("configuration, reference vs our runcard:")
    ref_np = {
        k: v
        for k, v in cfg.get("Nonperturbative", {}).items()
        if k.startswith(("lambda", "delta_lambda", "b0_over_bmax", "np_model"))
    }
    our_np = (
        dict(core.conf["Nonperturbative"])
        if core.conf.has_section("Nonperturbative")
        else {}
    )
    bad = [
        (k, v, our_np[k])
        for k, v in ref_np.items()
        if k in our_np and not _same(v, our_np[k])
    ]
    only_ref = [
        k
        for k in ref_np
        if k not in our_np and _num(ref_np[k]) and _f(ref_np[k]) != 0.0
    ]
    for k in (
        "lambda",
        "transition_points",
        "mu0_min",
        "muf_min",
        "compensate_fo",
        "form_np_prescription",
        "muf_follows_muB",
        "disable_asymmetry",
        "run_order",
        "fixed_order",
    ):
        rv = cfg.get("Calculation_settings", {}).get(k)
        ov = core.conf["Calculation_settings"].get(k)
        flag = "" if _same(rv, ov) else "   <-- DIFFERS"
        print(f"   {k:<24} ref {str(rv):<22} ours {str(ov):<22}{flag}")
    if bad:
        print("   NP ANCHOR DIFFERS:", bad)
    else:
        print(f"   NP anchor: all {len(ref_np)} shared entries agree")
    if only_ref:
        print(f"   reference sets NP keys we do not, at non-zero values: {only_ref}")
    print(f"   reference has [TNPs]: {'TNPs' in cfg}")

    # --- sum the reference onto our bins
    Q = np.asarray(h.axes["Q"].edges)
    Y = np.asarray(h.axes["Y"].edges)
    T = np.asarray(h.axes["qT"].edges)
    b = core.bins
    ourQ = np.unique(b[:, 0:2], axis=0)
    if ourQ.shape[0] != 1:
        raise SystemExit("this script expects a single Q bin in the cache")
    ourY = np.unique(b[:, 2:4], axis=0)
    ourT = np.unique(b[:, 4:6], axis=0)
    ourY = ourY[np.argsort(ourY[:, 0])]
    ourT = ourT[np.argsort(ourT[:, 0])]
    Ye = np.concatenate([ourY[:, 0], ourY[-1:, 1]])
    Te = np.concatenate([ourT[:, 0], ourT[-1:, 1]])
    if abs(Q[0] - ourQ[0, 0]) > TOL or abs(Q[-1] - ourQ[0, 1]) > TOL:
        raise SystemExit(
            f"Q window differs: reference {Q[0]}..{Q[-1]} vs ours {ourQ[0]}"
        )

    vals = (
        h[{"vars": args.var}].values()
        if "vars" in [a.name for a in h.axes]
        else h.values()
    )
    vals = np.asarray(vals)
    # collapse the single Q axis
    vals = vals.sum(axis=0) if vals.ndim == 3 else vals  # (Y, qT)

    iy = _sum_map(Y, Ye, "Y")
    it = _sum_map(T, Te, "qT")
    # Both sides are bin-integrated, so summing the reference's fine bins into
    # ours is exact. Fine bins outside our range (iy/it < 0) are dropped.
    R = np.zeros((Ye.size - 1, Te.size - 1), dtype=np.float64)
    for i in np.nonzero(iy >= 0)[0]:
        for j in np.nonzero(it >= 0)[0]:
            R[iy[i], it[j]] += vals[i, j]

    # ours is flattened qT-major, then Y (bins_from_gen_axes order)
    O = ours.reshape(Te.size - 1, Ye.size - 1).T  # -> (Y, qT)

    print(
        f"\ncomparing {O.size} bins  (Y {Ye.size-1} x qT {Te.size-1}), var={args.var}"
    )
    print(
        f"   totals: ours {O.sum():.8g}   reference {R.sum():.8g}   "
        f"ours/ref = {O.sum()/R.sum():.8f}"
    )
    rel = O / np.where(R != 0, R, np.nan) - 1.0
    print(
        f"   per-bin ours/ref - 1:  max |.| = {np.nanmax(np.abs(rel)):.3e}   "
        f"median |.| = {np.nanmedian(np.abs(rel)):.3e}"
    )
    print(f"\n   worst bins:")
    order = np.argsort(-np.abs(np.nan_to_num(rel)).ravel())[:6]
    for k in order:
        i, j = np.unravel_index(k, rel.shape)
        print(
            f"      Y [{Ye[i]:6.2f},{Ye[i+1]:6.2f}]  qT [{Te[j]:5.1f},{Te[j+1]:5.1f}]  "
            f"ours {O[i,j]:12.6g}  ref {R[i,j]:12.6g}  rel {rel[i,j]:+.3e}"
        )

    print("\n   by qT (max over Y):")
    for j in range(Te.size - 1):
        print(
            f"      qT [{Te[j]:5.1f},{Te[j+1]:5.1f}]  max |rel| = "
            f"{np.nanmax(np.abs(rel[:, j])):.3e}"
        )

    # Free extra: our own +Y / -Y symmetry. The folded-|Y| production cache
    # assumes it exactly, so it is worth knowing.
    if np.allclose(Ye, -Ye[::-1]):
        asym = np.abs(O / O[::-1, :] - 1.0)
        print(
            f"\n   our own +Y vs -Y symmetry: max |sigma(+Y)/sigma(-Y) - 1| = "
            f"{np.nanmax(asym):.3e}  "
            f"(the folded-|Y| production cache assumes this is 0)"
        )


def _num(x):
    try:
        float(str(x).strip())
        return True
    except Exception:
        return False


def _f(x):
    return float(str(x).strip())


def _same(a, b):
    if a is None or b is None:
        return a == b
    if _num(a) and _num(b):
        return abs(_f(a) - _f(b)) < 1e-9
    return str(a).strip() == str(b).strip()


if __name__ == "__main__":
    main()
