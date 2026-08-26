#!/usr/bin/env python3
"""Build a SCETlib autodiff cache for the gen binning of a specific rabbit card.

The cache is only valid for the bins it was built for, so the binning should
come from the card rather than being kept in sync by hand. This reads the gen
axes out of a datacard -- the fit channel's own axes for a gen-level sigmaUL
card, or the response auxiliary's gen axes for a reco card -- writes the
matching SCETlib runcard next to the output, and runs the expensive build:

    compressed bin rules (resummed)  +  frozen fixed-order grid (nonsingular)

WHAT IS HERE AND WHAT IS NOT. The build itself is nothing but SCETlib calls, so
it lives in SCETlib, in ``examples/matched_ad/prepare_cache.py``: the variation
plan, the node set, the rules, the member loop, the cache file. Node sets, bin
rules, member variations and the on-disk format are all theirs, and keeping our
copy of them meant learning about a layout change by getting wrong answers.
What is ours, and what this file is, is the rabbit side: reading the gen axes
off a card, writing the runcard from them, cutting the card's bins to a
``--subset``, and driving the upstream steps in the order they have to happen.

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

# F401 on `configure`: it is used, but through _resolve_steps below rather than
# by name, so that the equivalence harness can stub it.
from wremnants.postprocessing.scetlib_ad.xsec_backend import (  # noqa: E402,F401
    _scetlib_src,
    bins_from_gen_axes,
    configure,
)

# The build steps, all of them SCETlib's, reached through the example module.
# They are exposed as attributes of THIS module (see ``__getattr__``) and
# ``main`` calls them through it, so a study that needs one step changed can
# rebind it -- ``mod.plan_variations = ...`` -- and get the upstream behaviour
# for everything else. That is how the muF knot-spacing scans work.
_UPSTREAM_STEPS = (
    "plan_variations",
    "check_member_range",
    "build_variations",
    "build_prologue",
    "write_cache",
    "fork_member_build",
    "alphas_of",
    "find_alphas_pair",
    "pdf_set_size",
    "ensure_beamfunc_grids",
    "load_bins",
)

_UPSTREAM = None

# What ``main`` resolves through _resolve_steps: the upstream steps plus
# ``configure``, which is ours but is stubbed by the equivalence harness.
_STEPS = _UPSTREAM_STEPS + ("configure",)


def _upstream_prepare_cache():
    """The upstream builder module, ``examples/matched_ad/prepare_cache.py``.

    Loaded by path rather than imported, because the SCETlib examples are not a
    package. Deliberately LAZY: ``--dry-run`` and the card readers below work
    without SCETlib on the path at all, and importing this pulls in the
    compiled extension.
    """
    global _UPSTREAM
    if _UPSTREAM is not None:
        return _UPSTREAM
    import importlib.util

    path = os.path.join(_scetlib_src(), "examples", "matched_ad", "prepare_cache.py")
    if not os.path.exists(path):
        raise SystemExit(
            f"scetlib_ad: cannot find {path}, which is the builder. Point "
            "SCETLIB_SRC at a checkout of the autodiff-sigmaul branch (source "
            "its setup.sh)."
        )
    spec = importlib.util.spec_from_file_location("_scetlib_prepare_cache", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    missing = [n for n in _UPSTREAM_STEPS if not hasattr(mod, n)]
    if missing:
        raise SystemExit(
            f"{path} is missing {missing}. That checkout predates the build "
            "steps moving into the example; update SCETlib."
        )
    _UPSTREAM = mod
    return mod


def __getattr__(name):
    """Expose the upstream build steps as attributes of this module (PEP 562).

    Only reached for names NOT in the module dict, so an assignment from a
    caller wins -- which is what makes the steps overridable.
    """
    if name in _UPSTREAM_STEPS:
        return getattr(_upstream_prepare_cache(), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def _resolve_steps():
    """The build steps this run will use, an override taking precedence.

    Resolved ONCE, at the top of ``main``, from this module's own globals and
    only then from the upstream module -- so a caller that did
    ``mod.plan_variations = ...`` after importing this file gets its own
    version, and everything it did not touch stays upstream's. Going through
    ``globals()`` rather than ``sys.modules[__name__]`` is deliberate: this
    file is usually loaded by path (``spec_from_file_location``), which does
    not register a module name.
    """
    import types

    g = globals()
    return types.SimpleNamespace(
        **{n: g.get(n) or getattr(_upstream_prepare_cache(), n) for n in _STEPS}
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


def select_subset(bins, gen_axes, spec):
    """The rows of ``bins`` named by a '<absY indices>/<ptV indices>' subset.

    ``gen_axes`` is [(qT_name, edges), (Y_name, edges)] -- a LIST of tuples, not
    a dict -- and ``bins_from_gen_axes`` flattens qT-MAJOR (i over qT outer, j
    over absY inner). Index accordingly, so a subset cache's bins are a literal
    sub-list of the full one and nothing downstream cares.
    """
    n_qt = np.asarray(gen_axes[0][1]).size - 1
    n_y = np.asarray(gen_axes[1][1]).size - 1
    ysel, qsel = (x.strip() for x in spec.split("/"))
    iy = list(range(n_y)) if ysel == "*" else [int(v) for v in ysel.split(",")]
    iq = list(range(n_qt)) if qsel == "*" else [int(v) for v in qsel.split(",")]
    for v, n, what in ((iy, n_y, "absY"), (iq, n_qt, "ptV")):
        bad = [i for i in v if not 0 <= i < n]
        if bad:
            raise SystemExit(f"--subset: {what} index {bad} out of range 0..{n - 1}")
    keep = [i * n_y + j for i in iq for j in iy]
    print(
        f"   --subset: {len(keep)} of {n_y * n_qt} bins "
        f"(absY {list(iy)}, ptV {list(iq)}) -- a TEST cache, not for a fit",
        flush=True,
    )
    return bins[keep]


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
    ap.add_argument(
        "--subset",
        default=None,
        help="build only a SUBSET of the card's bins, as "
        "'iy0,iy1,.../iqt0,iqt1,...' of indices into the card's absY and ptV "
        "axes (either side may be '*' for all). The point is to make the REAL "
        "test -- validate_variations against the production corrections, every "
        "variation -- cheap enough to iterate on: 12 bins instead of 210 is "
        "minutes instead of hours, and it is the same code path, so a subset "
        "cache validates exactly what a full one does over the bins it covers. "
        "\n"
        "The indices must be CONTIGUOUS in both axes: the gen fold requires the "
        "cache to tile a rectangle exactly, so scattered picks are refused "
        "downstream with 'gen bin(s) are not exactly tiled by the cache'. Choose "
        "by COST, not by count -- the parallel axis is the bins, so wall time is "
        "roughly the slowest bin, and the lowest ptV bin is far more expensive "
        "than all the others together.\n"
        "Do NOT use a subset cache for a fit.",
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
        "--members",
        default=None,
        help="build only the member range 'LO:HI' (half-open, into the "
        "canonical variation-member list [eig pairs..., alphaS pair, muF "
        "pair]) and write a PARTIAL cache plus a .shard.json sidecar. The "
        "point is to split the one stage whose axis is serial -- the member "
        "loop -- across PROCESSES, since members share global AD state behind "
        "a mutex and cannot be threaded. Use "
        "scripts/rabbit/scetlib_ad/build_cache_parallel.py, which picks the "
        "ranges, runs the workers and merges them; a partial cache on its own "
        "is NOT usable.",
    )
    ap.add_argument(
        "--fork-members",
        type=int,
        default=1,
        help="split the PDF-member loop over this many forked children (1 = "
        "the serial loop, the default). EXACT but SLOW, and almost never what "
        "you want: forking after the rules are built is the only way to give "
        "two processes the SAME rules (they share the parent's node cache by "
        "copy-on-write, where independent processes each adapt their own node "
        "set and cannot be merged at all), but a forked child loses the TBB "
        "worker pool -- measured 99% CPU per child against the parent's 1900% "
        "-- so each child is single-threaded and a real bin count is ~100x "
        "slower per member. To parallelise a real build, split BINS across "
        "processes instead (--subset, then build_cache_parallel.py "
        "--merge-bins) or simply raise --threads: the member stage is parallel "
        "over NODES, not bins, and scales past the bin count. This exists for "
        "--fork-selftest, which is what proves the member merge exact.",
    )
    ap.add_argument(
        "--fork-selftest",
        action="store_true",
        help="after the forked build, build EVERY member serially in the "
        "parent as well and write <outname>.serial.npz. Same process, same "
        "node cache, same rules, so the only difference is forked-and-merged "
        "against serial: comparing the two is the decisive test of the merge "
        "(compare_caches.py --bytes). Costs one extra member pass.",
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

    # Every build step below is reached through this namespace, not called
    # directly, so that a rebound step is the one that runs. See _resolve_steps.
    me = _resolve_steps()

    conf, sigma = me.configure(runcard, args.threads)
    bins = bins_from_gen_axes(gen_axes, args.Q_lo, args.Q_hi)
    if args.subset:
        bins = select_subset(bins, gen_axes, args.subset)
    p0 = np.asarray(sigma.gradient_central(), dtype=np.float64)
    names = list(sigma.gradient_param_names())

    sing, nons = sigma.sub_pieces()

    # The variation members are DECIDED here, before the rules, because the PDF
    # eigenvector coefficients are ordinary AD parameters: registering them
    # changes the parameter vector, and so the anchor, the training points and
    # the rule fingerprint. Doing it afterwards would leave the rules built for
    # a shorter vector than the members are interpolated in, which is the
    # "call set_pdf_eig_params before building" error the extension raises.
    plan = None if args.no_pdf else me.plan_variations(p0, names, conf, args)
    mlo, mhi = 0, None
    if args.members and args.fork_members > 1:
        raise SystemExit(
            "--members builds ONE shard and --fork-members splits the members "
            "over children; use one or the other"
        )
    if args.members:
        if plan is None:
            raise SystemExit("--members needs variation members; not with --no-pdf")
        mlo, mhi = (int(x) for x in args.members.split(":"))
        mhi = me.check_member_range(plan, mlo, mhi)
    if plan and plan["n_eig"]:
        # BOTH pieces: each kernel interpolates its own members from the
        # coefficients, so they are ordinary parameters on both sides and the
        # gradient columns line up by name.
        sing.set_pdf_eig_params(plan["n_eig"])
        nons.set_pdf_eig_params(plan["n_eig"])
        p0 = np.asarray(sing.gradient_central(), dtype=np.float64)
        names = list(sing.gradient_param_names())
    print(f"\n{len(p0)} differentiable parameters:")
    for n, v in zip(names, p0):
        print(f"   {n:<24} {v:.6g}")

    me.build_prologue(sigma, sing, bins, p0, args.n_train, args.threads or 0)

    out = os.path.join(args.outdir, args.outname)
    if plan is None:
        print(
            "\n--no-pdf: physics-only cache. alphaS will be a derivative at "
            "FIXED PDF, and the card's pdf*/pdfAlphaS/resumFOScale* templates "
            "must be kept.",
            flush=True,
        )
        path = me.write_cache(sing, nons, bins, plan, out)
    elif args.fork_members > 1:
        # The children inherit the state as it is RIGHT NOW -- fork before the
        # parent touches the calculation again, so a --fork-selftest pass
        # cannot influence what they build. build_fn/write_fn are passed so a
        # rebound step reaches the children too.
        path = me.fork_member_build(
            sing,
            nons,
            bins,
            p0,
            plan,
            args,
            out,
            build_fn=me.build_variations,
            write_fn=me.write_cache,
        )
        if args.fork_selftest:
            t0 = time.time()
            me.build_variations(sing, nons, bins, p0, plan, args)
            print(
                f"selftest: all {len(plan['members'])} members serially in "
                f"{(time.time() - t0) / 60:.1f} min",
                flush=True,
            )
            me.write_cache(sing, nons, bins, plan, out + ".serial")
    else:
        me.build_variations(sing, nons, bins, p0, plan, args, mlo, mhi)
        path = me.write_cache(
            sing,
            nons,
            bins,
            plan,
            out,
            *((mlo, mhi) if args.members else (None, None)),
            args=args,
        )
        if args.members:
            print(
                "Merge the shards with\n"
                "   python scripts/rabbit/scetlib_ad/build_cache_parallel.py "
                "--merge-only -o <outdir>",
                flush=True,
            )
            return
    print(
        "\nNow check it:\n"
        f"   python scripts/rabbit/scetlib_ad/backend_check.py "
        f"--conf {runcard} --cache {path}"
    )


if __name__ == "__main__":
    main()
