#!/usr/bin/env python3
"""Build a SCETlib autodiff cache for the gen binning of a specific rabbit card.

The cache is only valid for the bins it was built for, so the binning should
come from the card rather than being kept in sync by hand. This reads the gen
axes out of a datacard -- the fit channel's own axes for a gen-level sigmaUL
card, or the response auxiliary's gen axes for a reco card -- writes the
matching SCETlib runcard next to the output, and runs the expensive build:

    compressed bin rules (resummed)  +  frozen fixed-order grid (nonsingular)

Cost scales with the number of gen bins. Measured on this SCETlib build:
~0.34 s/bin of rule building and ~2 s/bin of fixed-order warming, and ~0.84 MB
of cache per bin. A 200-bin gen-level card is therefore minutes and ~170 MB; the
5740-bin correction grid is hours and several GB.

    source <scetlib-cms>/setup.sh
    python scripts/rabbit/scetlib_ad/prepare_cache_for_card.py \
        --card <card>.hdf5 --base-conf <scetlib-cms>/examples/matched_ad/matched.conf \
        -o /path/to/cachedir
"""

import argparse
import configparser
import os
import sys
import time

import h5py
import numpy as np

sys.path.insert(
    0,
    os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ),
)

from wremnants.postprocessing.scetlib_ad.response import (  # noqa: E402
    DEFAULT_RESPONSE_GROUP,
)
from wremnants.postprocessing.scetlib_ad.xsec_backend import (  # noqa: E402
    _scetlib_src,
    bins_from_gen_axes,
    configure,
)


def gen_axes_from_card(path, gen_level):
    """[(qT name, edges), (|Y| name, edges)] for the card's gen grid."""
    from wums.ioutils import pickle_load_h5py

    with h5py.File(path, "r") as f:
        if gen_level:
            meta = pickle_load_h5py(f["meta"])
            channels = {
                n: i
                for n, i in meta["channel_info"].items()
                if not i.get("masked", False)
            }
            if len(channels) != 1:
                raise SystemExit(
                    f"expected a single non-masked channel, got {list(channels)}"
                )
            info = next(iter(channels.values()))
            axes = [
                (ax.name, np.asarray(ax.edges, dtype=np.float64)) for ax in info["axes"]
            ]
            if len(axes) != 2:
                raise SystemExit(
                    f"expected 2 gen axes (qT, |Y|), got {[n for n, _ in axes]}"
                )
            return axes
        from rabbit.auxiliary import read_auxiliary_from_h5

        aux = read_auxiliary_from_h5(f.get("auxiliary")) or {}
        if DEFAULT_RESPONSE_GROUP not in aux:
            raise SystemExit(
                f"the card has no {DEFAULT_RESPONSE_GROUP!r} auxiliary, so its "
                "gen binning is "
                "unknown. Rebuild it with setupRabbit --storeResponseMatrix, or "
                "pass --gen-level for a gen-level sigmaUL card."
            )
        b = aux[DEFAULT_RESPONSE_GROUP]
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in b["gen_axes"]]
        return [(n, np.asarray(b[f"edges__{n}"], dtype=np.float64)) for n in names]


def _upstream_prepare_cache():
    """The upstream example module, for its PDF-variation helpers.

    ``alphas_of`` / ``find_alphas_pair`` / ``pdf_set_size`` /
    ``ensure_beamfunc_grids`` are non-trivial (the last one fans out one
    single-threaded ~3.5 min process per PDF member and works around a shared
    .info race), and they live in the SCETlib checkout we already depend on.
    Loading them by path keeps one implementation, so an upstream fix arrives
    for free -- importing is safe because its ``main()`` is guarded.
    """
    import importlib.util

    path = os.path.join(_scetlib_src(), "examples", "matched_ad", "prepare_cache.py")
    if not os.path.exists(path):
        raise SystemExit(
            f"scetlib_ad: cannot find {path}, needed for the PDF/alphaS/muF "
            f"variation helpers. Pass --no-pdf to build a physics-only cache."
        )
    spec = importlib.util.spec_from_file_location("_scetlib_prepare_cache", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_variations(sing, nons, bins, p0, names, conf, args):
    """PDF eigenvector + alphaS-pair + muF variation members.

    Returns ``(n_eig, has_as, has_muf)`` for ``ScetlibCachedXsecTF``. These are
    what turn the model from "alpha_s at fixed PDF plus the NP lambdas" into a
    continuous parametrisation of the whole SCETlib theory uncertainty, so the
    corresponding card templates can be dropped instead of double-counted.
    """
    up = _upstream_prepare_cache()
    repo = _scetlib_src()
    pdf_set = conf["QCD"]["pdf_set"]
    n_eig = args.pdf_eig if args.pdf_eig >= 0 else (up.pdf_set_size(pdf_set) - 1) // 2
    nf = conf["QCD"].getint("nf", fallback=5)
    as_cen = float(p0[names.index("alphas")]) if "alphas" in names else 0.0
    as_step = 0.0
    pair = up.find_alphas_pair(pdf_set, args.as_pair, as_cen)
    # muF is not live in the kernel (the beam convolutions are frozen at their
    # own muF), so it rides on two extra members at kappa_F = 0.5 / 2.0,
    # interpolated in t = ln(kappa_F)/ln(muf_hi): exact at 0.5, 1, 2.
    muf_lo, muf_hi = (0.0, 0.0) if args.no_muf else (0.5, 2.0)
    if not (n_eig or pair or muf_hi):
        return 0, False, False

    members = list(range(1, 2 * n_eig + 1))
    sets = [pdf_set] * len(members)
    if n_eig:
        up.ensure_beamfunc_grids(repo, pdf_set, members, args.grid_jobs)
    if pair:
        # LAST and in (down, up) order -- both builders index the pair from the
        # end. Its effect is added to the EXISTING alphas slot, so one parameter
        # moves the calculation and the PDF together; without it, alphas is a
        # partial derivative at fixed PDF.
        down, up_set, as_step = pair
        sets += [down, up_set]
        members += [0, 0]
        up.ensure_beamfunc_grids(repo, down, [0], args.grid_jobs)
        up.ensure_beamfunc_grids(repo, up_set, [0], args.grid_jobs)
        print(
            f"  alphaS pair: {down} / {up_set}, central {as_cen:.4f} "
            f"+- {as_step:.4f}",
            flush=True,
        )
    if muf_hi:
        # After the alphaS pair, in (lo, hi) order. The set/member entries are
        # unused for these two -- the scale moves, not the PDF -- but one entry
        # per member is still required.
        sets += [pdf_set, pdf_set]
        members += [0, 0]
        print(f"  muF pair: kappa_F = {muf_lo} / {muf_hi}", flush=True)

    mem = np.array(members, dtype=np.int32)
    t0 = time.time()
    sing.build_pdf_variations(
        sets,
        mem,
        nf,
        p0,
        n_train_var=3,
        n_eig=n_eig,
        as_cen=as_cen,
        as_step=as_step,
        muf_lo=muf_lo,
        muf_hi=muf_hi,
    )
    print(
        f"  {n_eig} PDF eigenvector pairs for the resummed piece in "
        f"{(time.time() - t0) / 60:.1f} min",
        flush=True,
    )
    t0 = time.time()
    nons.build_fo_pdf_variations(
        sets,
        mem,
        nf,
        bins,
        np.asarray(nons.gradient_central()),
        n_eig=n_eig,
        as_cen=as_cen,
        as_step=as_step,
        muf_lo=muf_lo,
        muf_hi=muf_hi,
    )
    print(
        f"  ... and for the fixed-order piece in " f"{(time.time() - t0) / 60:.1f} min",
        flush=True,
    )
    return n_eig, as_step > 0.0, muf_hi > 0.0


def write_runcard(base_conf, out_path, gen_axes, Q_lo, Q_hi):
    """Base runcard + the card's grids, written where the cache can find it."""
    conf = configparser.ConfigParser(inline_comment_prefixes="#")
    conf.optionxform = str  # SCETlib option names are case-sensitive
    if not conf.read(base_conf):
        raise SystemExit(f"cannot read base runcard {base_conf!r}")
    for key, values in (
        ("Q", [Q_lo, Q_hi]),
        ("Y", list(gen_axes[1][1])),
        ("qT", list(gen_axes[0][1])),
    ):
        sec = f"Grid_{key}"
        if not conf.has_section(sec):
            conf.add_section(sec)
        conf[sec]["custom_grid"] = "true"
        conf[sec]["bins"] = "true"
        conf[sec]["values"] = "[" + ", ".join(f"{v:g}" for v in values) + "]"
    header = (
        "# Generated by scripts/rabbit/scetlib_ad/prepare_cache_for_card.py.\n"
        f"# Base runcard: {os.path.abspath(base_conf)}\n"
        "# The Grid_* sections are the card's gen binning; everything else is\n"
        "# inherited. Keep this file next to the cache -- the fit needs it to\n"
        "# rebuild the identical calculation the rules attach to.\n"
    )
    with open(out_path, "w") as f:
        f.write(header)
        conf.write(f)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--card", default=None, help="rabbit datacard hdf5")
    ap.add_argument(
        "--grid-json",
        default=None,
        help="explicit gen grid instead of a card, as JSON: "
        '\'{"Q": [60, 120], "Y": [...], "qT": [...]}\'. The Y edges may be '
        "signed (no |Y| folding assumed) -- useful for validating against a "
        "SCETlib reference run on its own signed grid, and required for W.",
    )
    ap.add_argument(
        "--base-conf",
        required=True,
        help="SCETlib runcard supplying the physics settings "
        "(orders, PDF, nonperturbative model); its Grid_* sections are replaced",
    )
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--outname", default="cache")
    ap.add_argument(
        "-g",
        "--gen-level",
        action="store_true",
        help="take the gen binning from the fit channel's own axes "
        "(a gen-level sigmaUL card) instead of the response "
        "auxiliary",
    )
    ap.add_argument("--Q-lo", type=float, default=60.0)
    ap.add_argument("--Q-hi", type=float, default=120.0)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument(
        "--n-train",
        type=int,
        default=9,
        help="training points for the rule compression. Accuracy "
        "tracks n_train/n_params, but the non-negative "
        "least-squares solve grows roughly like n_train^2, so "
        "raising it with the parameter count is the expensive "
        "knob (see doc/autodiff-design.md)",
    )
    ap.add_argument(
        "--pdf-eig",
        type=int,
        default=-1,
        help="number of PDF eigenvector pairs (default: all of the set). 0 "
        "still keeps the alphaS pair, which is a separate direction.",
    )
    ap.add_argument(
        "--as-pair",
        default="auto",
        help="'auto' finds <set>_as_0116/_as_0120, 'off' disables (alphas then "
        "moves the calculation but NOT the PDF), or 'down,up' explicitly.",
    )
    ap.add_argument(
        "--no-muf",
        action="store_true",
        help="skip the muF member pair (then resumScaleMuF does nothing and the "
        "card's resumFOScale* must be kept).",
    )
    ap.add_argument(
        "--no-pdf",
        action="store_true",
        help="physics-only cache: no PDF eigenvectors, no alphaS pair, no muF. "
        "alphaS is then a fixed-PDF derivative -- do not quote it.",
    )
    ap.add_argument(
        "--grid-jobs",
        type=int,
        default=0,
        help="parallel beamfunc-grid generation jobs (0 = one per core).",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="write the runcard and report the bin count and the "
        "projected cost, without building",
    )
    args = ap.parse_args()

    if (args.card is None) == (args.grid_json is None):
        raise SystemExit("give exactly one of --card and --grid-json")
    if args.grid_json:
        import json

        g = json.loads(args.grid_json)
        args.Q_lo, args.Q_hi = float(g["Q"][0]), float(g["Q"][-1])
        gen_axes = [
            ("qT", np.asarray(g["qT"], dtype=np.float64)),
            ("Y", np.asarray(g["Y"], dtype=np.float64)),
        ]
    else:
        gen_axes = gen_axes_from_card(args.card, args.gen_level)
    n_bins = int(np.prod([len(e) - 1 for _, e in gen_axes]))
    os.makedirs(args.outdir, exist_ok=True)
    runcard = os.path.join(args.outdir, args.outname + ".conf")
    write_runcard(args.base_conf, runcard, gen_axes, args.Q_lo, args.Q_hi)
    print(f"gen binning from {args.card or 'explicit --grid-json'}:")
    for name, edges in gen_axes:
        print(f"   {name:<12} {len(edges) - 1:4d} bins  [{edges[0]:g}, {edges[-1]:g}]")
    print(f"   Q            1 bin    [{args.Q_lo:g}, {args.Q_hi:g}]")
    print(f"   -> {n_bins} SCETlib bins; runcard written to {runcard}")
    # Per-bin costs measured on a 5740-bin build with all cores busy. A small
    # cache does not reach that: with 30 bins the same build ran ~10x slower per
    # bin, because there is not enough work to fill the pool.
    print(
        f"   projected (at full parallelism): ~{n_bins * 0.34 / 60:.0f} min of "
        f"rules, ~{n_bins * 2.0 / 60:.0f} min of fixed-order warming, "
        f"~{n_bins * 0.84:.0f} MB of cache"
    )
    if args.dry_run:
        return

    conf, sigma = configure(runcard, args.threads)
    bins = bins_from_gen_axes(gen_axes, args.Q_lo, args.Q_hi)
    p0 = np.asarray(sigma.gradient_central(), dtype=np.float64)
    names = list(sigma.gradient_param_names())
    print(f"\n{len(p0)} differentiable parameters:")
    for n, v in zip(names, p0):
        print(f"   {n:<24} {v:.6g}")

    sing, nons = sigma.sub_pieces()
    t0 = time.time()
    info = sing.build_bin_rules(
        bins,
        p0,
        n_train=args.n_train,
        n_hvp=1,
        seed=4242,
        n_jobs=args.threads or 0,
    )
    nodes = [d["nodes"] for d in info]
    resid = max(d["resid"] for d in info)
    # The upstream docstring quotes ~1e-15 for a fully converged solve, but the
    # active-set iteration cap deliberately trades residual for build time
    # (see the build_bin_rules comment in py/qT/qT.cpp): ~1e-7..1e-9 is normal.
    print(
        f"\nrules built in {(time.time() - t0) / 60:.1f} min "
        f"(median {int(np.median(nodes))} nodes/bin, worst training residual "
        f"{resid:.1e})",
        flush=True,
    )
    if resid > 1e-6:
        print(
            f"   WARNING: residual {resid:.1e} is large; the rules may not "
            f"reproduce the direct calculation. Cross-check before using.",
            flush=True,
        )

    # Populate the fixed-order grid over EVERY bin before saving. It fills lazily
    # per evaluated region, so warming one bin would leave the rest to be computed
    # on first use -- i.e. minutes of fixed-order work inside what is supposed to
    # be the cheap stage.
    t0 = time.time()
    sigma.sigma_binned_batch(bins, p0)
    print(f"fixed-order grid warmed in {(time.time() - t0) / 60:.1f} min", flush=True)

    if args.no_pdf:
        n_eig, has_as, has_muf = 0, False, False
        print(
            "\n--no-pdf: physics-only cache. alphaS will be a derivative at "
            "FIXED PDF, and the card's pdf*/pdfAlphaS/resumFOScale* templates "
            "must be kept.",
            flush=True,
        )
    else:
        n_eig, has_as, has_muf = build_variations(
            sing, nons, bins, p0, names, conf, args
        )

    from scetlib_tf import ScetlibCachedXsecTF

    out = os.path.join(args.outdir, args.outname)
    ScetlibCachedXsecTF(
        sing, nons, bins=bins, n_eig=n_eig, has_as=has_as, has_muf=has_muf
    ).save(out)
    path = out + ".npz"
    print(f"\nwrote {path} ({os.path.getsize(path) / 1e6:.1f} MB)")
    print(
        "\nNow check it:\n"
        f"   python scripts/rabbit/scetlib_ad/backend_check.py "
        f"--conf {runcard} --cache {path}"
    )


if __name__ == "__main__":
    main()
