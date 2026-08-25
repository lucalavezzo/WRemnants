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


def plan_variations(p0, names, conf, args):
    """The variation-member list, decided before anything is built.

    Returns ``None`` when there is nothing to vary, else a dict with the
    canonical member list and the metadata both builders need. Pure: it reads
    the runcard and LHAPDF and builds nothing, which is what lets a member
    range be selected (``--members``) and a split build reproduce the order a
    single-process build would have produced.

    The ORDER is load-bearing and is the one both builders index:
    ``[eig0 up, eig0 down, ..., alphaS down, alphaS up, muF lo, muF hi]``. The
    alphaS pair is read from the END (``dn = nvar-1``, ``up = nvar`` in the
    fixed-order evaluation) and the muF pair sits after it, so anything that
    reorders or drops a pair silently reinterprets the members.

    ``pairs`` gives the atomic [lo, hi) member ranges with a label each: a
    member range may only be cut on those boundaries, since ``[up, down]``
    pairs are what the interpolation is built from and the muF pair must stay
    LAST in whatever call builds it (``build_pdf_variations`` decides the
    Vary.muf leg from the position within the call).
    """
    up = _upstream_prepare_cache()
    pdf_set = conf["QCD"]["pdf_set"]
    n_eig = args.pdf_eig if args.pdf_eig >= 0 else (up.pdf_set_size(pdf_set) - 1) // 2
    as_cen = float(p0[names.index("alphas")]) if "alphas" in names else 0.0
    pair = up.find_alphas_pair(pdf_set, args.as_pair, as_cen)
    # muF is not live in the kernel (the beam convolutions are frozen at their
    # own muF), so it rides on two extra members at kappa_F = 0.5 / 2.0,
    # interpolated in t = ln(kappa_F)/ln(muf_hi): exact at 0.5, 1, 2.
    muf_lo, muf_hi = (0.0, 0.0) if args.no_muf else (0.5, 2.0)
    if not (n_eig or pair or muf_hi):
        return None

    members = list(range(1, 2 * n_eig + 1))
    sets = [pdf_set] * len(members)
    pairs = [(2 * e, 2 * e + 2, f"eig{e}") for e in range(n_eig)]
    as_step = 0.0
    if pair:
        # LAST and in (down, up) order -- both builders index the pair from the
        # end. Its effect is added to the EXISTING alphas slot, so one parameter
        # moves the calculation and the PDF together; without it, alphas is a
        # partial derivative at fixed PDF.
        down, up_set, as_step = pair
        pairs.append((len(members), len(members) + 2, "alphaS"))
        sets += [down, up_set]
        members += [0, 0]
    if muf_hi:
        # After the alphaS pair, in (lo, hi) order. The set/member entries are
        # unused for these two -- the scale moves, not the PDF -- but one entry
        # per member is still required.
        pairs.append((len(members), len(members) + 2, "muF"))
        sets += [pdf_set, pdf_set]
        members += [0, 0]
    return dict(
        pdf_set=pdf_set,
        nf=conf["QCD"].getint("nf", fallback=5),
        n_eig=n_eig,
        as_cen=as_cen,
        as_step=as_step,
        as_pair=None if not pair else (pair[0], pair[1]),
        muf_lo=muf_lo,
        muf_hi=muf_hi,
        sets=sets,
        members=members,
        pairs=pairs,
    )


def check_member_range(plan, lo, hi):
    """Validate a ``--members`` slice against the pair boundaries; return hi.

    Called early (from ``main``) so a bad range costs a second rather than the
    26 minutes of node set and rules that come before the member loop, and
    again where it is used.
    """
    n_all = len(plan["members"])
    hi = n_all if hi is None else hi
    if (lo, hi) == (0, n_all):
        return hi
    starts = {a for a, _, _ in plan["pairs"]} | {n_all}
    ends = {b for _, b, _ in plan["pairs"]} | {0}
    if lo not in ends or hi not in starts or not 0 <= lo < hi <= n_all:
        raise SystemExit(
            f"--members {lo}:{hi} does not fall on the variation pair "
            f"boundaries {[(a, b, t) for a, b, t in plan['pairs']]} of the "
            f"{n_all}-member list. A pair cannot be split: the interpolation is "
            "built from [up, down] pairs and the muF pair must stay last in the "
            "call that builds it."
        )
    return hi


def build_variations(sing, nons, bins, p0, plan, args, lo=0, hi=None):
    """PDF eigenvector + alphaS-pair + muF variation members.

    Returns ``(n_eig, has_as, has_muf)`` for ``ScetlibCachedXsecTF``. These are
    what turn the model from "alpha_s at fixed PDF plus the NP lambdas" into a
    continuous parametrisation of the whole SCETlib theory uncertainty, so the
    corresponding card templates can be dropped instead of double-counted.

    ``[lo, hi)`` selects a contiguous slice of ``plan['members']``, which must
    fall on the pair boundaries ``plan['pairs']`` gives. The default is the
    whole list, and then this is exactly what it always was. A slice builds a
    PARTIAL cache -- valid only after ``build_cache_parallel.py`` has merged it
    with the shards holding the other members.

    The returned triple describes the slice that was BUILT (the flags the
    C++ builders were given), not the full plan; ``main`` records the full
    plan's numbers in a shard's sidecar so the merge can restore them.
    """
    up = _upstream_prepare_cache()
    repo = _scetlib_src()
    hi = check_member_range(plan, lo, hi)
    sets = plan["sets"][lo:hi]
    members = plan["members"][lo:hi]
    # What THIS slice carries. n_eig is the number of eigenvector pairs inside
    # it, which is what build_pdf_variations needs to split its own member list
    # (has_as/has_muf are then implied by the count, exactly as for a full
    # build); the FULL n_eig is what the merged cache header must carry.
    n_eig = sum(1 for a, b, t in plan["pairs"] if t.startswith("eig") and lo <= a < hi)
    has_as = any(t == "alphaS" and lo <= a < hi for a, b, t in plan["pairs"])
    has_muf = any(t == "muF" and lo <= a < hi for a, b, t in plan["pairs"])
    as_step = plan["as_step"] if has_as else 0.0
    muf_lo, muf_hi = (plan["muf_lo"], plan["muf_hi"]) if has_muf else (0.0, 0.0)
    nf = plan["nf"]

    if n_eig:
        eig = [
            plan["members"][i]
            for a, b, t in plan["pairs"]
            if t.startswith("eig") and lo <= a < hi
            for i in range(a, b)
        ]
        up.ensure_beamfunc_grids(repo, plan["pdf_set"], eig, args.grid_jobs)
    if has_as:
        down, up_set = plan["as_pair"]
        up.ensure_beamfunc_grids(repo, down, [0], args.grid_jobs)
        up.ensure_beamfunc_grids(repo, up_set, [0], args.grid_jobs)
        print(
            f"  alphaS pair: {down} / {up_set}, central {plan['as_cen']:.4f} "
            f"+- {as_step:.4f}",
            flush=True,
        )
    if has_muf:
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
        as_cen=plan["as_cen"],
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
        as_cen=plan["as_cen"],
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


def write_cache(sing, nons, bins, plan, out, lo=None, hi=None, args=None):
    """Write the one-file cache, plus a sidecar when it holds only a member slice.

    The header always describes the WHOLE cache -- ``n_eig``, ``has_as``,
    ``has_muf`` come from the plan, not from the slice -- because a shard's
    member data is merged into that layout and the merge takes the header from
    any shard. The sidecar is what makes the slice self-describing: the merge
    needs to know which members these are, and that they tile [0, n_members)
    exactly once, to reproduce the single-process member ORDER.
    """
    from scetlib_tf import ScetlibCachedXsecTF

    n_eig = 0 if plan is None else plan["n_eig"]
    has_as = plan is not None and plan["as_step"] > 0.0
    has_muf = plan is not None and plan["muf_hi"] > 0.0
    ScetlibCachedXsecTF(
        sing, nons, bins=bins, n_eig=n_eig, has_as=has_as, has_muf=has_muf
    ).save(out)
    path = out + ".npz"
    print(f"\nwrote {path} ({os.path.getsize(path) / 1e6:.1f} MB)", flush=True)
    if lo is None:
        return path
    import json

    n_all = len(plan["members"])
    hi = n_all if hi is None else hi
    side = out + ".shard.json"
    with open(side, "w") as f:
        json.dump(
            dict(
                lo=lo,
                hi=hi,
                n_members=n_all,
                n_eig=n_eig,
                has_as=has_as,
                has_muf=has_muf,
                pairs=[[a, b, t] for a, b, t in plan["pairs"]],
                sets=plan["sets"],
                pdf_members=plan["members"],
                as_cen=plan["as_cen"],
                as_step=plan["as_step"],
                subset=None if args is None else args.subset,
                n_train=None if args is None else args.n_train,
            ),
            f,
            indent=1,
        )
    print(
        f"wrote {side}\nPARTIAL cache: members [{lo}, {hi}) of {n_all}. It is "
        "NOT usable on its own.",
        flush=True,
    )
    return path


def fork_member_build(sing, nons, bins, p0, plan, args, out):
    """Run the member loop in forked children and merge their shards.

    WHY FORK, and not N independent builds. The member data is stored as
    DIFFERENCES against the nominal rule and the frozen fixed-order grid, so
    merging members from two processes is only valid if those are identical.
    They are NOT reproducible across processes: the bin loop is a
    ``tbb::parallel_for`` whose range splitting depends on the workers actually
    available, and the integrator objects it hands out keep internal buffers,
    so which bins share a thread changes the adaptive outcome. Measured on a
    10-bin subset, four independent builds of the SAME configuration gave
    matched bin sums of 28.8515 / 28.8517 / 28.8518 / 28.8518 pb and 357 / 359
    retained nodes per bin -- a different NUMBER of sites, which would make
    ``Var::w`` shorter than ``sites`` and the replay read past the end.

    Forking after ``build_bin_rules`` removes the question: every child sees
    the parent's node cache and rules in its own copy-on-write image of the
    same memory, so "the same rules" is structural rather than hoped for. It
    also means the expensive prologue (the outer node set, then the rules) is
    paid ONCE rather than once per process.

    The children write shards and exit; the parent merges. Nothing in the
    parent touches the calculation while they run.

    THE CATCH, measured: a forked child does not get the TBB worker pool back.
    Each child runs at 99% CPU against the parent's 1900%, i.e. single-threaded,
    and the member stage is parallel over the NODES of every bin, so on a real
    bin count that is a ~100x loss per member. Use this for --fork-selftest,
    which is what proves the merge exact, not to make a build faster: for that,
    split BINS across processes (--subset + build_cache_parallel.py
    --merge-bins) or raise --threads.
    """
    import traceback

    import build_cache_parallel as bcp

    ranges = bcp.split_pairs(plan["pairs"], args.fork_members)
    shard_dir = out + "_shards"
    os.makedirs(shard_dir, exist_ok=True)
    print(
        f"\nforking {len(ranges)} children for {len(plan['members'])} members:",
        flush=True,
    )
    paths, kids = [], {}
    t0 = time.time()
    for lo, hi in ranges:
        shard = os.path.join(shard_dir, f"m{lo:04d}-{hi:04d}")
        paths.append(shard + ".npz")
        print(f"   members [{lo:3d}, {hi:3d}) -> {shard}.log", flush=True)
        sys.stdout.flush()
        sys.stderr.flush()
        pid = os.fork()
        if pid:
            kids[pid] = (lo, hi)
            continue
        # ---- child. os._exit, never a normal return: the parent's atexit
        # handlers, TF teardown and buffered streams must not run twice.
        rc = 1
        try:
            fd = os.open(shard + ".log", os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
            os.dup2(fd, 1)
            os.dup2(fd, 2)
            print(f"child {os.getpid()}: members [{lo}, {hi})", flush=True)
            build_variations(sing, nons, bins, p0, plan, args, lo, hi)
            write_cache(sing, nons, bins, plan, shard, lo, hi, args)
            rc = 0
        except BaseException:  # noqa: BLE001 -- report and exit, never propagate
            traceback.print_exc()
        finally:
            try:
                sys.stdout.flush()
                sys.stderr.flush()
            finally:
                os._exit(rc)
    bad = []
    while kids:
        pid, status = os.wait()
        if pid not in kids:
            continue  # someone else's child (a beamfunc-grid job, say)
        lo, hi = kids.pop(pid)
        code = os.waitstatus_to_exitcode(status)
        print(
            f"   child for [{lo}, {hi}) finished with {code} "
            f"({(time.time() - t0) / 60:.1f} min in)",
            flush=True,
        )
        if code:
            bad.append((lo, hi))
    print(f"member loop in {(time.time() - t0) / 60:.1f} min", flush=True)
    if bad:
        raise SystemExit(f"member ranges {bad} failed; see {shard_dir}/m*.log")
    return bcp.merge_shards(paths, out)


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

    conf, sigma = configure(runcard, args.threads)
    bins = bins_from_gen_axes(gen_axes, args.Q_lo, args.Q_hi)
    if args.subset:
        # gen_axes is [(qT_name, edges), (Y_name, edges)] -- a LIST of tuples,
        # not a dict -- and bins_from_gen_axes flattens qT-MAJOR (i over qT
        # outer, j over absY inner). Index accordingly, so a subset cache's bins
        # are a literal sub-list of the full one and nothing downstream cares.
        n_qt = np.asarray(gen_axes[0][1]).size - 1
        n_y = np.asarray(gen_axes[1][1]).size - 1
        ysel, qsel = (x.strip() for x in args.subset.split("/"))
        iy = list(range(n_y)) if ysel == "*" else [int(v) for v in ysel.split(",")]
        iq = list(range(n_qt)) if qsel == "*" else [int(v) for v in qsel.split(",")]
        for v, n, what in ((iy, n_y, "absY"), (iq, n_qt, "ptV")):
            bad = [i for i in v if not 0 <= i < n]
            if bad:
                raise SystemExit(
                    f"--subset: {what} index {bad} out of range 0..{n - 1}"
                )
        keep = [i * n_y + j for i in iq for j in iy]
        bins = bins[keep]
        print(
            f"   --subset: {len(bins)} of {n_y * n_qt} bins "
            f"(absY {list(iy)}, ptV {list(iq)}) -- a TEST cache, not for a fit",
            flush=True,
        )
    p0 = np.asarray(sigma.gradient_central(), dtype=np.float64)
    names = list(sigma.gradient_param_names())

    sing, nons = sigma.sub_pieces()

    # The variation members are DECIDED here, before the rules, because the PDF
    # eigenvector coefficients are ordinary AD parameters: registering them
    # changes the parameter vector, and so the anchor, the training points and
    # the rule fingerprint. Doing it afterwards would leave the rules built for
    # a shorter vector than the members are interpolated in, which is the
    # "call set_pdf_eig_params before building" error the extension raises.
    plan = None if args.no_pdf else plan_variations(p0, names, conf, args)
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
        mhi = check_member_range(plan, mlo, mhi)
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

    # STEP ONE, and the ordering is not cosmetic (scetlib 95c9f7e): the outer
    # node set must be built deliberately, not as a side effect of whichever
    # call first reaches _ad_bin_grid. Paired, the fixed-order half OWNS that
    # set and adapts it on the matched integrand; leave it implicit and
    # build_bin_rules gets there first and freezes the grid before the precision
    # target that was supposed to configure it exists. Symptom when we skipped
    # it: the cache builds, then every evaluation dies with "the fixed-order
    # column set no longer matches the rule; rebuild the cache".
    t0 = time.time()
    m0 = sigma.prepare(bins, p0)
    print(
        f"outer node set + matched cross sections in "
        f"{(time.time() - t0) / 60:.1f} min (sum {float(np.sum(m0)):.6g} pb)",
        flush=True,
    )

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

    # The fixed-order grid is already populated over every bin -- step one built
    # it. This is the check that makes that claim falsifiable rather than
    # assumed: if the matched values moved, something below rebuilt a grid under
    # different conditions and the cache about to be written is not the one that
    # was measured.
    t0 = time.time()
    m1 = np.asarray(sigma.sigma_binned_batch(bins, p0)[0], dtype=np.float64)
    drift = float(np.max(np.abs(m1 - m0) / np.maximum(np.abs(m0), 1e-300)))
    print(
        f"fixed-order grid verified in {time.time() - t0:.1f} s "
        f"(max drift from step one {drift:.2e})",
        flush=True,
    )
    if drift > 1.0e-12:
        raise RuntimeError(
            f"matched cross sections moved by {drift:.3e} between the node-set "
            "build and the rule build; a grid was rebuilt under different "
            "conditions and the cache would not match what was measured."
        )

    out = os.path.join(args.outdir, args.outname)
    if plan is None:
        print(
            "\n--no-pdf: physics-only cache. alphaS will be a derivative at "
            "FIXED PDF, and the card's pdf*/pdfAlphaS/resumFOScale* templates "
            "must be kept.",
            flush=True,
        )
        path = write_cache(sing, nons, bins, plan, out)
    elif args.fork_members > 1:
        # The children inherit the state as it is RIGHT NOW -- fork before the
        # parent touches the calculation again, so a --fork-selftest pass
        # cannot influence what they build.
        path = fork_member_build(sing, nons, bins, p0, plan, args, out)
        if args.fork_selftest:
            t0 = time.time()
            build_variations(sing, nons, bins, p0, plan, args)
            print(
                f"selftest: all {len(plan['members'])} members serially in "
                f"{(time.time() - t0) / 60:.1f} min",
                flush=True,
            )
            write_cache(sing, nons, bins, plan, out + ".serial")
    else:
        build_variations(sing, nons, bins, p0, plan, args, mlo, mhi)
        path = write_cache(
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
