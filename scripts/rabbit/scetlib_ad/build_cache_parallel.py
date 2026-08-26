#!/usr/bin/env python3
"""Parallel SCETlib autodiff cache builds for a rabbit card: split the BINS.

WHICH AXIS, and it is not symmetric. The expensive stage of a full cache is the
PDF-member loop (measured on the 210-bin Z card at target_precision_rel 1e-3:
21.9 min for the outer node set, 4.4 min for the rules, and the fixed-order
member sweep dominating the rest -- 7-13 h for the 62 members that 29
eigenvector pairs need). The member loop is serial within a process, which
invites splitting it over processes. Do not. Three measured reasons, in order of
importance:

1. THE MEMBER LOOP IS NOT STARVED OF PARALLELISM. What a member costs is
   `set_pdf_keep_nodes`, which refills every frozen fixed-order node of every
   bin at the new PDF and is parallel over NODES ("Parallel over ALL nodes of
   ALL bins at once, which is the half that scales"), not over bins. On a
   10-bin subset, where a bin-parallel-only stage would saturate at 10 threads,
   the member stage still runs 3-4x faster at --threads 200 than at 32 (0.22
   against 0.7-1.0 min per member), and a live 210-bin build sits at 145 busy
   cores of the 200 it asked for. So the lever is threads (or nodes), not
   processes.

2. MEMBERS FROM INDEPENDENT PROCESSES CANNOT BE MERGED. The per-member data is
   stored as DIFFERENCES against the nominal rule and the frozen fixed-order
   grid, and the nominal side is not reproducible: four independent builds of
   the same configuration at the same --threads gave matched bin sums 28.8515 /
   28.8517 / 28.8518 / 28.8518 pb and 357 / 359 / 359 / 371 median nodes per
   bin. A different SITE COUNT is not a rounding difference -- `Var::w` is one
   weight per site, so the merged rule would read past the end of a shorter
   vector. `scetlib_cache.merge_shards` compares every non-member field of
   every rule byte for byte and refuses; that guard is the only thing standing
   between a cross-build member merge and a silently wrong cache, because the
   C++ loaders check settings, struct sizes and version and would accept one.

3. A FORKED CHILD LOSES THE TBB WORKER POOL. Forking after `build_bin_rules`
   does give the children the parent's node cache and rules by copy-on-write,
   which makes the member merge exact (`prepare_cache_for_card --fork-members`,
   and `--fork-selftest` proves it against the serial build in the same
   process). But each child then runs at 99% CPU against the parent's 1900%:
   correct, and a net loss for any real bin count.

WHAT DOES WORK: SPLIT BINS. A bin's rule is self-contained -- its own outer
grid, sites, node data, members and fixed-order deltas -- so nothing in it is a
difference against another bin, and processes that adapted their own node sets
can each contribute whole bins. Only the global header has to agree
(configuration fingerprint, rule options, anchor, parameter names, variation
metadata), and every one of those is a deterministic function of the runcard and
the flags rather than of the quadrature. `--bin-groups` builds them and
`--merge-bins` assembles them; `prepare_cache_for_card --subset` is the same
thing done by hand.

    # one node, bins split four ways (ptV 0 alone: it costs more than the rest)
    singularity exec ... ./incontainer.sh python3 \
        scripts/rabbit/scetlib_ad/build_cache_parallel.py \
        --card <card>.hdf5 --base-conf <base>.conf -o <cachedir> --threads 150 \
        --bin-groups '*/0;*/1,2,3,4;*/5,...,12;*/13,...,20'

    # or merge caches built separately (condor, several nodes)
    ... build_cache_parallel.py --base-conf <base>.conf -o <cachedir> \
        --merge-bins shard*/cache.npz

WHERE THE MERGES LIVE. Splicing the two blobs is an operation on SCETlib's own
serialisation format -- which per-rule fields must match bit for bit, that the
options struct carries uninitialised padding, that a shard's stored eigenvector
count is per-shard and not per-cache -- so it lives in SCETlib, in
`py/scetlib_cache.py`, next to the C++ that defines the layout. This module is
the scheduler: it decides the bin groups, runs the workers, budgets the threads
and calls those merges. The names are re-exported here so the tools that read a
cache keep one import.
"""

import argparse
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(
    0,
    os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ),
)

BUILDER = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "prepare_cache_for_card.py"
)

# SCETlib's cache-format module: the blob readers and the two merges. Pure
# python and numpy -- no TensorFlow, no compiled extension -- but it is only on
# the path once <scetlib-cms>/setup.sh has been sourced, so it is imported
# lazily and re-exported, which keeps `--bin-groups` (which needs no merge until
# its workers finish) able to start with a clear error rather than a traceback.
_MERGE_API = (
    "parse_rule_blob",
    "merge_rule_blobs",
    "parse_fo_blob",
    "merge_fo_blobs",
    "copy_runcard",
    "merge_shards",
    "merge_bin_caches",
    "split_pairs",
    # The blob primitives too, not out of tidiness: the validation tools that
    # prove the merges right (split_merge_selftest.py, compare_caches.py) have
    # to WRITE a blob to split one, so they need the same struct formats and
    # emitters the merges use. Two readers of one format is how a merge starts
    # passing against a bug in itself.
    "_U64",
    "_I32",
    "_F64",
    "_RULE_MAGIC",
    "_FO_MAGIC",
    "_FO_NODE",
    "_opts_fields",
    "_parse_fo_grid",
    "_emit_fo_grid",
)


def _cache_module():
    try:
        import scetlib_cache
    except ImportError as e:
        raise SystemExit(
            "scetlib_ad: cannot import scetlib_cache, which holds the cache "
            "format and the merges. Run\n"
            "    source <path-to>/scetlib-cms/setup.sh\n"
            "inside the container first (it puts <src>/py on PYTHONPATH).\n"
            f"Original error: {e}"
        ) from e
    return scetlib_cache


def __getattr__(name):
    """Re-export SCETlib's merge API from this module (PEP 562, so it is lazy)."""
    if name in _MERGE_API:
        return getattr(_cache_module(), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--base-conf", required=True)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--outname", default="cache")
    ap.add_argument(
        "--bin-groups",
        default=None,
        help="build these bin subsets as separate processes and "
        "merge them: a ';'-separated list of --subset strings, "
        "e.g. '*/0;*/1,2,3;*/4,...'. Each process builds ALL "
        "the members for its bins, which is what makes the "
        "merge safe. Pick the groups by COST, not count -- the "
        "lowest ptV bin costs more than all the others "
        "together, so it wants a group of its own.",
    )
    ap.add_argument("--threads", type=int, default=120)
    ap.add_argument("--pdf-eig", type=int, default=-1)
    ap.add_argument("--as-pair", default="auto")
    ap.add_argument("--no-muf", action="store_true")
    ap.add_argument(
        "--stagger",
        type=float,
        default=20.0,
        help="seconds between worker launches; they all read the "
        "same PDF and beamfunc .info files at startup",
    )
    ap.add_argument(
        "--merge-bins",
        nargs="+",
        metavar="CACHE.npz",
        help="merge caches built over DISJOINT BIN SUBSETS "
        "(--subset), each carrying the complete member list, "
        "into -o/--outname. Use it for shards built by hand or "
        "on separate nodes; --bin-groups does the same thing "
        "and runs them for you.",
    )
    ap.add_argument(
        "--merge-only",
        action="store_true",
        help="merge the shards already under <outdir>/shards "
        "without building anything",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="print the bin groups and the commands, build nothing",
    )
    args, rest = ap.parse_known_args()

    os.makedirs(args.outdir, exist_ok=True)
    shard_dir = os.path.join(args.outdir, "shards")
    out = os.path.join(args.outdir, args.outname)

    if args.merge_bins:
        _cache_module().merge_bin_caches(args.merge_bins, out)
        return
    if args.merge_only:
        # Both layouts: <outdir>/shards/*/<outname>.npz from --bin-groups, and
        # <outdir>/<outname>_shards/m*.npz from --fork-members (whose parent
        # normally merges them itself -- this is for picking up the pieces when
        # it died after the children succeeded). Bin groups if they carry no
        # shard sidecar, member shards if they do.
        paths = []
        if os.path.isdir(shard_dir):
            paths += [
                os.path.join(shard_dir, d, args.outname + ".npz")
                for d in os.listdir(shard_dir)
            ]
        fork_dir = out + "_shards"
        if os.path.isdir(fork_dir):
            paths += [
                os.path.join(fork_dir, f)
                for f in os.listdir(fork_dir)
                if f.endswith(".npz")
            ]
        paths = sorted(p for p in paths if os.path.exists(p))
        if not paths:
            raise SystemExit(f"no {args.outname}.npz under {shard_dir}")
        member_shards = [
            p for p in paths if os.path.exists(p[: -len(".npz")] + ".shard.json")
        ]
        if member_shards and len(member_shards) != len(paths):
            raise SystemExit(
                "some shards are member slices (.shard.json) and some are not; "
                "merge the member slices of each bin group first, then "
                "--merge-bins the results"
            )
        sc = _cache_module()
        (sc.merge_shards if member_shards else sc.merge_bin_caches)(paths, out)
        return

    if not args.bin_groups:
        raise SystemExit(
            "nothing to do. There are two ways to parallelise a cache build, "
            "and they are not symmetric:\n"
            "  * BINS across processes: --bin-groups '<subset>;<subset>;...' "
            "here, or build them yourself with prepare_cache_for_card "
            "--subset and merge with --merge-bins. Safe: a bin's rule is "
            "self-contained.\n"
            "  * MEMBERS inside one process: prepare_cache_for_card "
            "--fork-members N. Exact, but a forked child loses the TBB worker "
            "pool (measured 99% CPU against the parent's 1900%), so it is a "
            "net loss for any real bin count.\n"
            "Members across INDEPENDENT processes is not on the list because "
            "it cannot be merged: two processes adapt different node sets and "
            "even retain a different NUMBER of sites, and the member data is "
            "stored as differences against those. See the module docstring."
        )

    groups = [g.strip() for g in args.bin_groups.split(";") if g.strip()]
    print(f"{len(groups)} bin groups, each with the complete member list:")
    for g in groups:
        print(f"   --subset {g!r}")
    procs, logs = [], []
    t0 = time.time()
    for i, sub in enumerate(groups):
        d = os.path.join(shard_dir, f"b{i:02d}")
        os.makedirs(d, exist_ok=True)
        cmd = (
            [
                sys.executable,
                "-u",
                BUILDER,
                "--base-conf",
                args.base_conf,
                "-o",
                d,
                "--outname",
                args.outname,
                "--threads",
                str(args.threads),
                "--pdf-eig",
                str(args.pdf_eig),
                "--as-pair",
                args.as_pair,
                "--subset",
                sub,
            ]
            + (["--no-muf"] if args.no_muf else [])
            + rest
        )
        log = os.path.join(d, "build.log")
        print(f"\n[{i}] {' '.join(cmd)}\n    -> {log}", flush=True)
        if args.dry_run:
            continue
        f = open(log, "w")
        logs.append(f)
        procs.append((subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT), d))
        if i + 1 < len(groups):
            time.sleep(args.stagger)
    if args.dry_run:
        return

    bad = []
    for p, d in procs:
        if p.wait() != 0:
            bad.append(d)
    for f in logs:
        f.close()
    print(f"\nworkers done in {(time.time() - t0) / 60:.1f} min", flush=True)
    if bad:
        raise SystemExit(
            "these bin groups failed; see their build.log:\n  " + "\n  ".join(bad)
        )

    paths = [os.path.join(d, args.outname + ".npz") for _, d in procs]
    merged = _cache_module().merge_bin_caches(paths, out)
    print(
        f"\nNow check it:\n"
        f"   python scripts/rabbit/scetlib_ad/backend_check.py "
        f"--conf {out}.conf --cache {merged}"
    )


if __name__ == "__main__":
    main()
