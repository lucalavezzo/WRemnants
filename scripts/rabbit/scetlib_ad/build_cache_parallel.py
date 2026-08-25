#!/usr/bin/env python3
"""Parallel SCETlib autodiff cache builds, and the two merges they need.

WHICH AXIS. The expensive stage of a full cache is the PDF-member loop
(measured on the 210-bin Z card at target_precision_rel 1e-3: 21.9 min for the
outer node set, 4.4 min for the rules, 1.6 min for four members of the resummed
piece and 54.8 min for four of the fixed-order one, so ~14 h for the 62 members
that 29 eigenvector pairs need). The member loop is serial -- members share
global AD state behind a mutex -- which invites splitting it over processes.
Do not. Three measured reasons, in order of importance:

1. THE MEMBER LOOP IS NOT STARVED OF PARALLELISM. What a member costs is
   `set_pdf_keep_nodes`, which refills every frozen fixed-order node of every
   bin at the new PDF and is parallel over NODES ("Parallel over ALL nodes of
   ALL bins at once, which is the half that scales"), not over bins. On a
   10-bin subset, where a bin-parallel-only stage would saturate at 10 threads,
   the member stage still runs 3-4x faster at --threads 200 than at 32 (0.22
   against 0.7-1.0 min per member), and a live 210-bin build sits at 145 busy
   cores of the 200 it asked for. So the lever is threads (or nodes), not
   processes; splitting members can only recover the ~28% the tail wastes.

2. MEMBERS FROM INDEPENDENT PROCESSES CANNOT BE MERGED. The per-member data is
   stored as DIFFERENCES against the nominal rule and the frozen fixed-order
   grid -- `_fo_var_d[m]` is literally (member sweep - central sweep), and the
   resummed `Var::w`/`Var::c_val` are anchored on sites and a c_val the nominal
   rule owns. That is only mergeable if the nominal side is identical, and it
   is not reproducible: `_parallel_run` is a tbb::parallel_for whose range
   splitting depends on the workers actually available, and the integrators it
   hands out keep internal buffers, so which bins share a thread changes the
   adaptive outcome. Four independent builds of the same configuration at the
   same --threads gave matched bin sums 28.8515 / 28.8517 / 28.8518 / 28.8518 pb
   and 357 / 359 / 359 / 371 median nodes per bin; bin by bin, the structure
   (n_grid, n_sites, n_fo_w) differed in 9 of 10 bins. A different SITE COUNT is
   not a rounding difference: `Var::w` is one weight per site, so the merged
   rule would read past the end of a shorter vector. `merge_shards` therefore
   compares every non-member field of every rule byte for byte and refuses.

3. A FORKED CHILD LOSES THE TBB WORKER POOL. Forking after `build_bin_rules`
   does give the children the parent's node cache and rules by copy-on-write,
   which makes the member merge exact (that is what
   `prepare_cache_for_card --fork-members` does, and `--fork-selftest` proves
   it against the serial build in the same process). But each child then runs
   at 99% CPU against the parent's 1900%: single-threaded. Correct, and a net
   loss for any real bin count.

WHAT DOES WORK: SPLIT BINS. A bin's rule is self-contained -- its own outer
grid, sites, node data, members and fixed-order deltas -- so nothing in it is a
difference against another bin, and processes that adapted their own node sets
can each contribute whole bins. Only the global header has to agree
(configuration fingerprint, rule options, anchor, parameter names, variation
metadata), and every one of those is a deterministic function of the runcard
and the flags rather than of the quadrature. `--bin-groups` builds them and
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

MEMBER ORDER, when a member merge is used at all. The order is load-bearing:
the fixed-order evaluation indexes the alphaS pair from the END
(`i_as_dn = size - n_muf - 2`) and the muF pair after it, while the resummed
kernel indexes eigenvectors forward from 0. It is preserved by construction:
the canonical list comes from `prepare_cache_for_card.plan_variations`, the
same function the workers use; shards are contiguous slices of it cut only on
pair boundaries; and they are concatenated in `lo` order with a check that the
slices tile [0, n_members) exactly once.
"""

import argparse
import json
import os
import struct
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(
    0,
    os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ),
)

# NB prepare_cache_for_card imports THIS module (for the merge, in its
# --fork-members mode), so nothing here may import it at module level.

BUILDER = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "prepare_cache_for_card.py"
)

# The blobs are written with ostream::write of native PODs, so this reader is
# valid on the machine that wrote them (little-endian x86 here). Nothing in the
# cache is portable across architectures anyway -- the C++ loader refuses a
# different sizeof(GlobalData).
_U64 = struct.Struct("<Q")
_I32 = struct.Struct("<i")
_F64 = struct.Struct("<d")
_RULE_MAGIC = b"SCTRULE"
_FO_MAGIC = b"SCETFOG"
# _Fo_node, field by field: Q, Y, qT, weight, C[3][3], mu_ref, err,
# int32 rule, uint64 n_orders, uint8 muR_resolved.
_FO_NODE = 4 * 8 + 9 * 8 + 8 + 8 + 4 + 8 + 1


def _opts_fields(raw):
    """The Bin_rule_opts fields, at their natural x86-64 offsets.

    Only the ones that decide what the rules ARE: n_jobs is concurrency and
    append is a build mode, and neither changes the answer.
    """
    return dict(
        n_train=_I32.unpack_from(raw, 0)[0],
        n_hvp=_I32.unpack_from(raw, 4)[0],
        with_hess=raw[8],
        scale=_F64.unpack_from(raw, 16)[0],
        seed=_I32.unpack_from(raw, 24)[0],
        orthogonal=raw[28],
        zero_absolute=raw[29],
        tol=_F64.unpack_from(raw, 32)[0],
        max_iter=_U64.unpack_from(raw, 40)[0],
    )


class _R:
    """Cursor over a blob, so every field consumed is accounted for."""

    def __init__(self, buf, what):
        self.b = memoryview(buf)
        self.p = 0
        self.what = what

    def raw(self, n):
        if self.p + n > len(self.b):
            raise ValueError(f"{self.what}: truncated at {self.p} (+{n})")
        out = self.b[self.p : self.p + n]
        self.p += n
        return out

    def u64(self):
        return _U64.unpack(self.raw(8))[0]

    def i32(self):
        return _I32.unpack(self.raw(4))[0]

    def f64(self):
        return _F64.unpack(self.raw(8))[0]

    def vec(self, elem):
        """A rule_put_vec / fo vector: count then count*elem bytes."""
        n = self.u64()
        self.raw(n * elem)
        return n

    def at_end(self):
        return self.p == len(self.b)


# --------------------------------------------------------------------------- #
# the compressed bin rules                                                    #
# --------------------------------------------------------------------------- #
def _parse_rules(buf, opts_size):
    """Byte offsets of everything the merge has to compare or splice.

    ``opts_size`` is sizeof(Bin_rule_opts), which the format does NOT record --
    it is written as one raw POD. It is determined by trying the plausible
    values and keeping the one that parses the whole blob exactly (see
    :func:`parse_rule_blob`), which is a stronger check than reading the struct
    definition would be.
    """
    r = _R(buf, "rule blob")
    magic = bytes(r.raw(8))
    if magic[:7] != _RULE_MAGIC:
        raise ValueError(f"not a compressed-rule blob: {magic!r}")
    version = magic[7:8]
    sizes = [
        r.raw(4).tobytes() for _ in range(4)
    ]  # Site, GlobalData, HardData, NodeData
    sz_site, sz_g, sz_h, sz_nd = (_I32.unpack(x)[0] for x in sizes)
    fn = r.u64()
    fingerprint = bytes(r.raw(fn))
    n_anchor = r.u64()
    anchor = np.frombuffer(r.raw(8 * n_anchor).tobytes(), dtype="<f8")
    opts_raw = bytes(r.raw(opts_size))
    # Bin_rule_opts is written as ONE raw POD, padding included, and the
    # padding is not initialised -- two processes that set identical options
    # still differ in those bytes. Compare the FIELDS.
    opts = _opts_fields(opts_raw) if opts_size == 56 else None
    meta_at = r.p
    meta = dict(
        n_eig=r.i32(),
        as_index=r.i32(),
        as_cen=r.f64(),
        as_step=r.f64(),
        muf_index=r.i32(),
        muf_lnstep=r.f64(),
    )
    head = bytes(r.b[:meta_at])
    n_bins = r.u64()
    rules = []
    for _ in range(n_bins):
        start = r.p
        key = bytes(r.raw(48))
        n_grid = r.vec(32)  # grid: array<double, 4>
        n_sites = r.vec(sz_site)  # sites
        r.vec(sz_g)  # g
        r.vec(sz_h)  # h_incl
        r.vec(sz_h)  # h_asym
        r.vec(1)  # has_asym
        r.vec(sz_nd)  # nd
        c_val = r.f64()
        r.vec(8)  # c_grad
        r.vec(8)  # c_hess
        n_fo_w = r.vec(8)  # fo_w
        r.u64()  # n_sites_full
        resid = r.f64()
        r.raw(4)  # n_iter
        r.raw(1)  # capped
        nv_at = r.p
        nv = r.u64()
        var = []
        for _ in range(nv):
            v0 = r.p
            r.vec(8)  # w
            r.vec(sz_nd)  # nd
            r.f64()  # c_val
            r.vec(8)  # c_grad
            r.f64()  # g_muf_ratio
            r.f64()  # g_v_muf
            r.raw(4)  # is_muf
            var.append((v0, r.p))
        rules.append(
            dict(
                key=key,
                shared=(start, nv_at),
                n_var=nv,
                var=var,
                end=r.p,
                n_grid=n_grid,
                n_sites=n_sites,
                n_fo_w=n_fo_w,
                c_val=c_val,
                resid=resid,
            )
        )
    if not r.at_end():
        raise ValueError(f"rule blob: {len(r.b) - r.p} trailing bytes")
    return dict(
        version=version,
        head=head,
        sizes=(sz_site, sz_g, sz_h, sz_nd),
        fingerprint=fingerprint,
        anchor=anchor,
        opts=opts,
        opts_raw=opts_raw,
        meta=meta,
        rules=rules,
        buf=buf,
    )


def parse_rule_blob(buf):
    """As :func:`_parse_rules`, determining sizeof(Bin_rule_opts) by parsing.

    A wrong size shifts every subsequent field, so the vector counts come out
    as garbage and the parse dies or ends in the wrong place. Exactly one
    candidate can survive a full parse of a multi-MB blob.
    """
    ok, err = [], None
    for cand in range(24, 129, 4):
        try:
            ok.append((cand, _parse_rules(buf, cand)))
        except Exception as e:  # noqa: BLE001 -- candidate rejected, keep the last
            err = e
    if len(ok) != 1:
        raise ValueError(
            "cannot pin down sizeof(Bin_rule_opts) from the rule blob "
            f"({[c for c, _ in ok]} all parse); last error: {err}"
        )
    return ok[0][1]


def merge_rule_blobs(parsed, meta):
    """One rule blob with every shard's members, in shard order.

    Everything but the members is taken from the first shard VERBATIM, after
    checking the others match it byte for byte -- so a merged blob of one shard
    is that shard, and a merged blob of a full member set is exactly what a
    single-process build would have written.
    """
    a = parsed[0]
    if a["opts"] is None:
        raise SystemExit(
            f"sizeof(Bin_rule_opts) is {len(a['opts_raw'])}, not the 56 this "
            "reader knows how to interpret; the field offsets in _opts_fields "
            "need updating before a merge can be trusted"
        )
    for b in parsed[1:]:
        for field in ("version", "sizes", "fingerprint", "opts"):
            if a[field] != b[field]:
                raise SystemExit(f"shards disagree on the rule {field}")
        if not np.array_equal(a["anchor"], b["anchor"]):
            raise SystemExit("shards disagree on the rule anchor")
        if len(a["rules"]) != len(b["rules"]):
            raise SystemExit("shards have different bin counts")
    out = [a["head"]]  # magic..opts, byte for byte from the first shard
    out.append(
        _I32.pack(meta["n_eig"])
        + _I32.pack(meta["as_index"])
        + _F64.pack(meta["as_cen"])
        + _F64.pack(meta["as_step"])
        + _I32.pack(meta["muf_index"])
        + _F64.pack(meta["muf_lnstep"])
    )
    out.append(_U64.pack(len(a["rules"])))
    for ib in range(len(a["rules"])):
        ra = a["rules"][ib]
        s0, s1 = ra["shared"]
        shared = bytes(a["buf"][s0:s1])
        n_var = 0
        chunks = []
        for sh in parsed:
            rb = sh["rules"][ib]
            if rb["key"] != ra["key"]:
                raise SystemExit(
                    f"bin {ib}: shards built different bins ({rb['key']!r})"
                )
            t0, t1 = rb["shared"]
            if bytes(sh["buf"][t0:t1]) != shared:
                raise SystemExit(
                    f"bin {ib}: the shards' NOMINAL rule differs byte for byte, so "
                    "their members are differences against different "
                    "calculations and must not be merged. The rules are only "
                    "reproducible from the same bins, anchor, seed, n_train and "
                    "parameter set."
                )
            n_var += rb["n_var"]
            for v0, v1 in rb["var"]:
                chunks.append(bytes(sh["buf"][v0:v1]))
        out.append(shared)
        out.append(_U64.pack(n_var))
        out.extend(chunks)
    return b"".join(out)


# --------------------------------------------------------------------------- #
# the frozen fixed-order grid                                                 #
# --------------------------------------------------------------------------- #
def _parse_fo_grid(r):
    """One ``put_bins`` block -> {bin key bytes: node block bytes}.

    A dict, not a byte range: ``_Fo_cache::bins`` is an unordered_map filled by
    the parallel bin loop, so the ORDER a process writes its bins in is not
    reproducible even when every byte of content is. Comparing per key is the
    check that means something; comparing the raw block would fail on thread
    scheduling alone.
    """
    n = r.u64()
    out = {}
    for _ in range(n):
        key = bytes(r.raw(48))
        p0 = r.p
        nn = r.u64()
        r.raw(nn * _FO_NODE)
        if key in out:
            raise ValueError("duplicate bin key in a fixed-order grid")
        out[key] = bytes(r.b[p0 : r.p])
    return out


def _emit_fo_grid(grid):
    return b"".join([_U64.pack(len(grid))] + [k + v for k, v in grid.items()])


def parse_fo_blob(buf):
    r = _R(buf, "fo blob")
    magic = bytes(r.raw(8))
    if magic[:7] != _FO_MAGIC:
        raise ValueError(f"not a frozen fixed-order blob: {magic!r}")
    fn = r.u64()
    fingerprint = bytes(r.raw(fn))
    grid = _parse_fo_grid(r)
    n_mem = r.u64()
    meta, var_bins, deltas, muf = None, None, [], []
    if n_mem:
        meta = dict(
            n_eig=r.i32(),
            as_cen=r.f64(),
            as_step=r.f64(),
            as_anchor=r.f64(),
            muf_lnstep=r.f64(),
            muf_index=r.i32(),
        )
        nb = r.u64()
        var_bins = np.frombuffer(r.raw(8 * nb).tobytes(), dtype="<f8")
        for _ in range(n_mem):
            nd = r.u64()
            deltas.append(bytes(r.raw(8 * nd)))
        for _ in range(2):
            have = r.raw(1)[0]
            muf.append(_parse_fo_grid(r) if have else None)
    if not r.at_end():
        raise ValueError(f"fo blob: {len(r.b) - r.p} trailing bytes")
    return dict(
        version=magic[7:8],
        fingerprint=fingerprint,
        grid=grid,
        meta=meta,
        var_bins=var_bins,
        deltas=deltas,
        muf=muf,
    )


def merge_fo_blobs(parsed, meta):
    a = parsed[0]
    for b in parsed[1:]:
        if (a["version"], a["fingerprint"]) != (b["version"], b["fingerprint"]):
            raise SystemExit("shards disagree on the fixed-order fingerprint")
        if set(a["grid"]) != set(b["grid"]):
            raise SystemExit("shards froze different fixed-order bins")
        bad = [k for k in a["grid"] if a["grid"][k] != b["grid"][k]]
        if bad:
            raise SystemExit(
                f"{len(bad)} of {len(a['grid'])} frozen fixed-order bins differ "
                "byte for byte between shards. Every member is stored as "
                "(member sweep - central sweep) on this grid, so the deltas do "
                "not belong to a common central and must not be merged."
            )
        if not np.array_equal(a["var_bins"], b["var_bins"]):
            raise SystemExit("shards disagree on the fixed-order variation bins")
    muf = [None, None]
    for sh in parsed:
        for i, g in enumerate(sh["muf"]):
            if g is not None:
                if muf[i] is not None:
                    raise SystemExit(f"two shards carry the muF grid {i}")
                muf[i] = g
    out = [
        _FO_MAGIC + a["version"],
        _U64.pack(len(a["fingerprint"])),
        a["fingerprint"],
        _emit_fo_grid(a["grid"]),
    ]
    deltas = [d for sh in parsed for d in sh["deltas"]]
    out.append(_U64.pack(len(deltas)))
    if deltas:
        out.append(
            _I32.pack(meta["n_eig"])
            + _F64.pack(meta["as_cen"])
            + _F64.pack(meta["as_step"])
            + _F64.pack(meta["as_anchor"])
            + _F64.pack(meta["muf_lnstep"])
            + _I32.pack(meta["muf_index"])
        )
        out.append(_U64.pack(a["var_bins"].size))
        out.append(a["var_bins"].astype("<f8").tobytes())
        for d in deltas:
            out.append(_U64.pack(len(d) // 8))
            out.append(d)
        for g in muf:
            out.append(b"\x01" + _emit_fo_grid(g) if g is not None else b"\x00")
    return b"".join(out)


# --------------------------------------------------------------------------- #
# the merge                                                                   #
# --------------------------------------------------------------------------- #
def copy_runcard(paths, out, verbose=True):
    """Put the shards' runcard next to the merged cache, if they have one.

    A cache without the runcard it was built from is unusable -- the fit has to
    rebuild the identical calculation for the rules to attach to -- and the
    shards all wrote the same one (``write_runcard`` uses the CARD's full grid,
    not the subset). Copied rather than assumed identical: a difference here
    means the shards are not the same calculation.
    """
    confs = [p[: -len(".npz")] + ".conf" for p in paths]
    if not all(os.path.exists(c) for c in confs):
        return None
    first = open(confs[0]).read()
    for c in confs[1:]:
        if open(c).read() != first:
            raise SystemExit(f"{c} differs from {confs[0]}: not one calculation")
    dest = (out[: -len(".npz")] if out.endswith(".npz") else out) + ".conf"
    with open(dest, "w") as f:
        f.write(first)
    if verbose:
        print(f"wrote {dest}")
    return dest


def merge_shards(shard_paths, out_path, verbose=True):
    """Merge partial caches into one, in canonical member order."""
    shards = []
    for path in shard_paths:
        if not path.endswith(".npz"):
            raise SystemExit(f"{path}: expected a .npz partial cache")
        side = path[: -len(".npz")] + ".shard.json"
        if not os.path.exists(side):
            raise SystemExit(
                f"{path}: no {os.path.basename(side)} sidecar; a "
                "partial cache must come from --members"
            )
        with open(side) as f:
            info = json.load(f)
        shards.append((info, path))
    shards.sort(key=lambda x: x[0]["lo"])

    # The slices must tile the canonical member list exactly once -- this is
    # what guarantees the merged order IS the single-process order.
    n_members = shards[0][0]["n_members"]
    want = 0
    for info, path in shards:
        if info["lo"] != want:
            raise SystemExit(
                f"the shards do not tile the member list: expected a shard "
                f"starting at {want}, got {info['lo']} ({path})"
            )
        want = info["hi"]
        for k in (
            "n_members",
            "n_eig",
            "has_as",
            "has_muf",
            "pairs",
            "sets",
            "pdf_members",
        ):
            if info[k] != shards[0][0][k]:
                raise SystemExit(f"{path}: shard plans disagree on {k!r}")
    if want != n_members:
        raise SystemExit(
            f"the shards cover members [0, {want}) of {n_members}; a cache "
            "missing members would silently mis-index the alphaS pair, which "
            "is counted from the END"
        )

    npz = [np.load(p, allow_pickle=False) for _, p in shards]
    d0 = npz[0]
    for d, (_, p) in zip(npz[1:], shards[1:]):
        for k in ("format", "bins", "anchor", "names", "n_eig", "has_as", "has_muf"):
            if not np.array_equal(d0[k], d[k]):
                raise SystemExit(f"{p}: shards disagree on {k!r}")

    rules = [parse_rule_blob(d["rules"].tobytes()) for d in npz]
    for (info, path), pr in zip(shards, rules):
        n_var = {r["n_var"] for r in pr["rules"]}
        if n_var != {info["hi"] - info["lo"]}:
            raise SystemExit(
                f"{path}: sidecar says {info['hi'] - info['lo']} members but "
                f"the rules carry {sorted(n_var)}"
            )
    fos = [parse_fo_blob(d["fo"].tobytes()) for d in npz]
    for (info, path), pf in zip(shards, fos):
        if len(pf["deltas"]) != info["hi"] - info["lo"]:
            raise SystemExit(
                f"{path}: sidecar says {info['hi'] - info['lo']} members but "
                f"the fixed-order cache carries {len(pf['deltas'])}"
            )

    # The header describes the WHOLE cache: the full eigenvector count, and the
    # alphaS / muF metadata from whichever shard built those pairs (the others
    # wrote the "absent" values, -1 / 0).
    n_eig_total = int(shards[0][0]["n_eig"])
    rmeta = dict(
        n_eig=n_eig_total,
        as_index=max(p["meta"]["as_index"] for p in rules),
        as_cen=rules[0]["meta"]["as_cen"],
        as_step=max(p["meta"]["as_step"] for p in rules),
        muf_index=max(p["meta"]["muf_index"] for p in rules),
        muf_lnstep=max(p["meta"]["muf_lnstep"] for p in rules),
    )
    fmeta = dict(
        n_eig=n_eig_total,
        as_cen=fos[0]["meta"]["as_cen"],
        as_step=max(p["meta"]["as_step"] for p in fos),
        as_anchor=fos[0]["meta"]["as_anchor"],
        muf_lnstep=max(p["meta"]["muf_lnstep"] for p in fos),
        muf_index=max(p["meta"]["muf_index"] for p in fos),
    )
    for p in fos:
        if p["meta"]["as_anchor"] != fmeta["as_anchor"]:
            raise SystemExit("shards disagree on the fixed-order alphaS anchor")
    # as_cen is written by every shard whether or not it built the pair, so it
    # is a consistency check rather than something to pick a winner for.
    for tag, got in (
        ("rule", [p["meta"]["as_cen"] for p in rules]),
        ("fixed-order", [p["meta"]["as_cen"] for p in fos]),
    ):
        if len(set(got)) != 1:
            raise SystemExit(f"shards disagree on the {tag} alphaS centre: {got}")
    if bool(shards[0][0]["has_as"]) != (rmeta["as_step"] > 0.0):
        raise SystemExit("the plan wants an alphaS pair but no shard built one")
    if bool(shards[0][0]["has_muf"]) != (rmeta["muf_lnstep"] > 0.0):
        raise SystemExit("the plan wants a muF pair but no shard built one")

    merged_rules = merge_rule_blobs(rules, rmeta)
    merged_fo = merge_fo_blobs(fos, fmeta)
    out = out_path[: -len(".npz")] if out_path.endswith(".npz") else out_path
    np.savez_compressed(
        out,
        format=d0["format"],
        n_eig=d0["n_eig"],
        has_as=d0["has_as"],
        has_muf=d0["has_muf"],
        rules=np.frombuffer(merged_rules, dtype=np.uint8),
        fo=np.frombuffer(merged_fo, dtype=np.uint8),
        bins=d0["bins"],
        anchor=d0["anchor"],
        names=d0["names"],
    )
    for d in npz:
        d.close()
    copy_runcard([p for _, p in shards], out, verbose)
    if verbose:
        print(
            f"merged {len(shards)} shards -> {out}.npz "
            f"({os.path.getsize(out + '.npz') / 1e6:.1f} MB): {n_members} members "
            f"({n_eig_total} eigenvector pairs, alphaS pair "
            f"{'yes' if rmeta['as_step'] > 0 else 'no'}, muF pair "
            f"{'yes' if rmeta['muf_lnstep'] > 0 else 'no'})"
        )
    return out + ".npz"


# --------------------------------------------------------------------------- #
# the other axis: merging caches built over DIFFERENT BINS                     #
# --------------------------------------------------------------------------- #
def merge_bin_caches(paths, out_path, verbose=True):
    """Assemble one cache from several built over disjoint bin subsets.

    This merge is much weaker than the member one, and that is the point: a
    bin's rule is SELF-CONTAINED -- its own outer grid, sites, node data, PDF
    members and fixed-order deltas -- so bins from processes that adapted their
    own node sets are each internally consistent. Nothing is a difference
    against anything in another bin. Two processes therefore CAN contribute
    different bins to one cache, where they cannot contribute different members
    of the same bin.

    What still has to agree is only what is global: the configuration
    fingerprint, the rule options, the anchor, the parameter names and the
    variation metadata (which member is the alphaS pair, the muF step). All of
    those are deterministic functions of the runcard and the flags, not of the
    quadrature.

    Bins come out in the canonical order a full build would have written them,
    qT-major with |Y| inner, so the merged cache is indistinguishable in layout
    from a single-process one. (``GenFold`` matches by value, so the order is
    not load-bearing, but a surprise here would be hard to spot later.)
    """
    npz = [np.load(p, allow_pickle=False) for p in paths]
    d0 = npz[0]
    for d, p in zip(npz[1:], paths[1:]):
        for k in ("format", "n_eig", "has_as", "has_muf", "anchor", "names"):
            if not np.array_equal(d0[k], d[k]):
                raise SystemExit(f"{p}: bin-shards disagree on {k!r}")
    rules = [parse_rule_blob(d["rules"].tobytes()) for d in npz]
    fos = [parse_fo_blob(d["fo"].tobytes()) for d in npz]
    a, fa = rules[0], fos[0]
    if a["opts"] is None:
        raise SystemExit("unexpected sizeof(Bin_rule_opts); see _opts_fields")
    for b, fb, p in zip(rules[1:], fos[1:], paths[1:]):
        for field in ("version", "sizes", "fingerprint", "opts", "meta"):
            if a[field] != b[field]:
                raise SystemExit(
                    f"{p}: bin-shards disagree on the rule {field}\n"
                    f"  {a[field]}\n  {b[field]}"
                )
        if not np.array_equal(a["anchor"], b["anchor"]):
            raise SystemExit(f"{p}: bin-shards disagree on the anchor")
        if (fa["version"], fa["fingerprint"], fa["meta"]) != (
            fb["version"],
            fb["fingerprint"],
            fb["meta"],
        ):
            raise SystemExit(f"{p}: bin-shards disagree on the fixed-order header")
        if len(fa["deltas"]) != len(fb["deltas"]):
            raise SystemExit(
                f"{p}: {len(fb['deltas'])} members against {len(fa['deltas'])}; "
                "every bin-shard must carry the COMPLETE member list"
            )
    n_mem = len(fa["deltas"])

    # bin -> (which shard, which record). The rule key is the same six doubles
    # as a row of `bins`, so one byte key indexes all three structures.
    owner, order = {}, []
    for i, (d, pr, pf) in enumerate(zip(npz, rules, fos)):
        bins = np.ascontiguousarray(d["bins"], dtype="<f8")
        if not np.array_equal(bins, np.asarray(pf["var_bins"]).reshape(-1, 6)):
            raise SystemExit(
                f"{paths[i]}: the fixed-order variation bins are not the "
                "cache's bins; the shard was not built in one piece"
            )
        by_key = {r["key"]: r for r in pr["rules"]}
        for j in range(bins.shape[0]):
            key = bins[j].tobytes()
            if key not in by_key:
                raise SystemExit(f"{paths[i]}: bin {j} has no rule")
            if key in owner:
                raise SystemExit(
                    f"bin {bins[j]} is in more than one shard; the subsets must "
                    "be disjoint or the merged cache would double-count it"
                )
            owner[key] = (i, by_key[key], j)
            order.append((bins[j][4], bins[j][2], key))
    order.sort()  # qT-major, |Y| inner: the build order
    keys = [k for _, _, k in order]
    # ONE order for all three structures, and it is not cosmetic. The
    # fixed-order replay refuses bins that are not element-for-element
    # `_fo_var_bins` (fo_binned_pdf_batch), and the matched rule replay indexes
    # the fixed-order member deltas by the RULE's position
    # (_fo_columns_mapped(which[nb], ...) -> _fo_stage_members(ibin)). So
    # `bins`, the rule records and `_fo_var_bins`/`_fo_var_d` must all be
    # emitted in the same order, which is why everything below iterates `keys`.

    out = [a["head"]]  # magic..opts, byte for byte from the first shard
    m = a["meta"]
    out.append(
        _I32.pack(m["n_eig"])
        + _I32.pack(m["as_index"])
        + _F64.pack(m["as_cen"])
        + _F64.pack(m["as_step"])
        + _I32.pack(m["muf_index"])
        + _F64.pack(m["muf_lnstep"])
    )
    out.append(_U64.pack(len(keys)))
    for k in keys:
        i, rec, _ = owner[k]
        out.append(bytes(rules[i]["buf"][rec["shared"][0] : rec["end"]]))
    merged_rules = b"".join(out)

    grid, muf = {}, [None, None]
    for pf in fos:
        for k, v in pf["grid"].items():
            grid[k] = v
        for t in range(2):
            if pf["muf"][t] is not None:
                muf[t] = {**(muf[t] or {}), **pf["muf"][t]}
    if set(grid) != set(keys):
        raise SystemExit("the frozen fixed-order grid does not cover the bins")
    grid = {k: grid[k] for k in keys}
    for t in range(2):
        if muf[t] is not None:
            if set(muf[t]) != set(keys):
                raise SystemExit(f"the muF grid {t} does not cover the bins")
            muf[t] = {k: muf[t][k] for k in keys}
    fout = [
        _FO_MAGIC + fa["version"],
        _U64.pack(len(fa["fingerprint"])),
        fa["fingerprint"],
        _emit_fo_grid(grid),
        _U64.pack(n_mem),
    ]
    if n_mem:
        fm = fa["meta"]
        fout.append(
            _I32.pack(fm["n_eig"])
            + _F64.pack(fm["as_cen"])
            + _F64.pack(fm["as_step"])
            + _F64.pack(fm["as_anchor"])
            + _F64.pack(fm["muf_lnstep"])
            + _I32.pack(fm["muf_index"])
        )
        fout.append(_U64.pack(6 * len(keys)))
        fout.append(b"".join(keys))
        stride = 3 * 8  # kFoVarStride doubles per bin
        for mi in range(n_mem):
            blocks = []
            for k in keys:
                i, _, j = owner[k]
                d = fos[i]["deltas"][mi]
                blocks.append(d[j * stride : (j + 1) * stride])
            joined = b"".join(blocks)
            fout.append(_U64.pack(len(joined) // 8))
            fout.append(joined)
        for g in muf:
            fout.append(b"\x01" + _emit_fo_grid(g) if g is not None else b"\x00")
    merged_fo = b"".join(fout)

    bins = np.concatenate(
        [np.frombuffer(k, dtype="<f8").reshape(1, 6) for k in keys], axis=0
    )
    outp = out_path[: -len(".npz")] if out_path.endswith(".npz") else out_path
    np.savez_compressed(
        outp,
        format=d0["format"],
        n_eig=d0["n_eig"],
        has_as=d0["has_as"],
        has_muf=d0["has_muf"],
        rules=np.frombuffer(merged_rules, dtype=np.uint8),
        fo=np.frombuffer(merged_fo, dtype=np.uint8),
        bins=bins,
        anchor=d0["anchor"],
        names=d0["names"],
    )
    for d in npz:
        d.close()
    copy_runcard(paths, outp, verbose)
    if verbose:
        print(
            f"merged {len(paths)} bin-shards -> {outp}.npz "
            f"({os.path.getsize(outp + '.npz') / 1e6:.1f} MB): {len(keys)} bins "
            f"x {n_mem} members"
        )
    return outp + ".npz"


# --------------------------------------------------------------------------- #
# the split build                                                             #
# --------------------------------------------------------------------------- #
def split_pairs(pairs, n_procs):
    """Contiguous member ranges, cut only on pair boundaries.

    Balanced by PAIR COUNT, which is the right cost measure here: every member
    costs the same (one sweep over every bin), unlike the bins, which do not.
    """
    n = len(pairs)
    k = max(1, min(int(n_procs), n))
    out = []
    for i in range(k):
        lo = pairs[(i * n) // k][0]
        hi = pairs[((i + 1) * n) // k - 1][1]
        out.append((lo, hi))
    return out


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
        merge_bin_caches(args.merge_bins, out)
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
        (merge_shards if member_shards else merge_bin_caches)(paths, out)
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
    merged = merge_bin_caches(paths, out)
    print(
        f"\nNow check it:\n"
        f"   python scripts/rabbit/scetlib_ad/backend_check.py "
        f"--conf {out}.conf --cache {merged}"
    )


if __name__ == "__main__":
    main()
