"""Generate the CVH-refit efficiency scale-factor map for muon_efficiencies_cvh.hpp.

Reads the 'cvhEfficiency' histogram produced by
'scripts/histmakers/mz_dilepton.py --cvhEfficiencyHists' (axes
pt, eta, phi, charge, passCVH, in uncorrected muon kinematics) and writes the
measured SF = eps_data / eps_MC as a map in (eta, phi'), where

    phi' = phi - q * C / pt

is the muon azimuth at the module rather than at the vertex (C is the track
bending, see muon_efficiencies_cvh.hpp). Working in phi' makes the map charge-
and pt-independent: the two charges' holes, which sit at opposite +-C/pt
offsets in vertex phi, line up, so a single map applies to everything without
smearing the sharp module edges.

The affected cells (a couple of small rectangles) are written into the
GENERATED block of the header in place; everywhere else the SF is exactly 1.

Usage:
    python scripts/analysisTools/make_cvh_efficiency_sf.py -i <cvhEfficiency.hdf5>
"""

import argparse
import os

import h5py
import numpy as np

from wums import ioutils, logging

logger = logging.child_logger(__name__)

HEADER = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "wremnants/production/include/muon_efficiencies_cvh.hpp",
)

# track bending to the module radius, C = 0.3*B*r/2 [rad*GeV]; must match the
# value hard-coded in muon_efficiencies_cvh.hpp (cvh_bending_C)
BENDING_C = 0.214

# search boxes (module frame) in which to keep the correction; everything
# outside stays at SF = 1. (eta_lo, eta_hi, phiprime_lo, phiprime_hi)
SEARCH_BOXES = {
    "hotspot1 (TIB-L2, detId 369141860)": (-0.10, 0.85, 0.60, 1.20),
    "hotspot2 (uncorrected)": (-1.95, -1.40, 4.90, 5.55),
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input", required=True, help="cvhEfficiency hdf5")
    p.add_argument(
        "--nSigma",
        type=float,
        default=3.0,
        help="Keep a cell's SF only if it is below 1 by at least this many sigma",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the generated block but do not modify the header",
    )
    p.add_argument("-v", "--verbose", type=int, default=3)
    return p.parse_args()


def load(fname):
    """Return the (data, total MC) 'cvhEfficiency' histograms."""
    hdata, hmc = None, None
    with h5py.File(fname, "r") as f:
        for proc in f.keys():
            if proc == "meta_info":
                continue
            res = ioutils.pickle_load_h5py(f[proc])
            h = res["output"]["cvhEfficiency"].get()
            if res["dataset"].get("is_data", False):
                hdata = h if hdata is None else hdata + h
            else:
                hmc = h if hmc is None else hmc + h
    return hdata, hmc


def to_module_phi(h):
    """Collapse (pt, eta, phi, charge, pass) -> (eta, phi', pass).

    Each (pt, charge) slice is shifted in phi by -q*C/pt with a linear sub-bin
    redistribution (the shift is well below one phi bin over 25-65 GeV), then
    summed. This is the histogram-level equivalent of filling phi' per muon.
    """
    v = h.values()  # pt, eta, phi, charge, pass
    pt_c = h.axes["pt"].centers
    phi_w = h.axes["phi"].edges[1] - h.axes["phi"].edges[0]
    nphi = len(h.axes["phi"])
    charges = {0: -1.0, 1: +1.0}  # charge-axis index -> muon charge

    out = np.zeros((len(h.axes["eta"]), nphi, 2))
    for iq, q in charges.items():
        for ipt, pt in enumerate(pt_c):
            f = -q * BENDING_C / pt / phi_w  # shift in units of phi bins
            assert abs(f) < 1.0, f"sub-bin shift assumption violated: |f|={abs(f)}"
            step = 1 if f > 0 else -1
            src = v[ipt, :, :, iq, :]  # eta, phi, pass
            out += (1.0 - abs(f)) * src
            out += abs(f) * np.roll(src, step, axis=1)
    return out


def efficiency(v):
    """(eff, err) from an array with the pass axis last (index 1 = pass)."""
    n = v.sum(axis=-1)
    eff = np.divide(v[..., 1], n, out=np.ones_like(n), where=n > 0)
    err = np.divide(
        np.sqrt(np.maximum(eff * (1 - eff), 0) * n),
        n,
        out=np.zeros_like(n),
        where=n > 0,
    )
    return eff, err


def scale_factor(vd, vm):
    """SF = eff_data/eff_MC and its (data-dominated) uncertainty."""
    ed, dd = efficiency(vd)
    em, _ = efficiency(vm)
    em = np.where(em > 0, em, 1.0)
    return ed / em, dd / em


def build_regions(hd, hm, nsigma):
    """Return a list of (name, i0, i1, j0, j1, sf[nEta,nPhi]) rectangles."""
    md = to_module_phi(hd)
    mm = to_module_phi(hm)
    sf, err = scale_factor(md, mm)
    # holes only: never upweight, and drop cells not significantly below 1
    keep = (sf < 1.0 - nsigma * err) & (err > 0)
    sf = np.where(keep, np.clip(sf, 0.0, 1.0), 1.0)

    eta_e = hd.axes["eta"].edges
    phi_e = hd.axes["phi"].edges
    eta_c = hd.axes["eta"].centers
    phi_c = hd.axes["phi"].centers

    regions = []
    for name, (e0, e1, p0, p1) in SEARCH_BOXES.items():
        box = (
            (eta_c[:, None] > e0)
            & (eta_c[:, None] < e1)
            & (phi_c[None, :] > p0)
            & (phi_c[None, :] < p1)
        )
        sig = keep & box
        if not sig.any():
            logger.warning(f"no significant cells in {name}, skipping")
            continue
        ii, jj = np.where(sig)
        i0, i1, j0, j1 = ii.min(), ii.max(), jj.min(), jj.max()
        block = sf[i0 : i1 + 1, j0 : j1 + 1].copy()
        regions.append((name, i0, i1, j0, j1, block, eta_e, phi_e))
        logger.info(
            f"{name}: eta [{eta_e[i0]:.2f},{eta_e[i1+1]:.2f}] "
            f"phi' [{phi_e[j0]:.3f},{phi_e[j1+1]:.3f}], "
            f"{block.shape[0]}x{block.shape[1]} cells, min SF = {block.min():.3f}"
        )
    return regions


def render(regions):
    """Render the C++ GENERATED block."""
    lines = [
        "// BEGIN GENERATED -- do not edit by hand, see make_cvh_efficiency_sf.py",
        "",
    ]
    for k, (name, i0, i1, j0, j1, block, eta_e, phi_e) in enumerate(regions):
        neta, nphi = block.shape
        eta_lo, eta_w = eta_e[i0], eta_e[1] - eta_e[0]
        phi_lo, phi_w = phi_e[j0], phi_e[1] - phi_e[0]
        lines.append(f"// {name}")
        lines.append(
            f"//   eta in [{eta_e[i0]:.2f}, {eta_e[i1+1]:.2f}), "
            f"phi' in [{phi_e[j0]:.3f}, {phi_e[j1+1]:.3f})"
        )
        lines.append(f"inline const double cvh_sf_region{k}[{neta} * {nphi}] = {{")
        for ie in range(neta):
            row = ", ".join(f"{block[ie, ip]:.3f}" for ip in range(nphi))
            eta_center = eta_lo + (ie + 0.5) * eta_w
            lines.append(f"    {row}, // eta ~ {eta_center:+.3f}")
        # (trailing comma on the last row is legal in a C++ aggregate initializer)
        lines.append("};")
        lines.append("")

    lines.append("inline const CvhSFRegion cvh_sf_regions[] = {")
    for k, (name, i0, i1, j0, j1, block, eta_e, phi_e) in enumerate(regions):
        neta, nphi = block.shape
        eta_lo, eta_w = eta_e[i0], eta_e[1] - eta_e[0]
        phi_lo, phi_w = phi_e[j0], phi_e[1] - phi_e[0]
        lines.append(
            f"    {{{eta_lo:.4f}, {eta_w:.4f}, {neta}, "
            f"{phi_lo:.4f}, {phi_w:.4f}, {nphi}, cvh_sf_region{k}}},"
        )
    lines.append("};")
    lines.append("")
    lines.append("// END GENERATED")
    return "\n".join(lines)


def patch_header(block):
    with open(HEADER) as f:
        text = f.read()
    b0 = text.index("// BEGIN GENERATED")
    e0 = text.index("// END GENERATED")
    e0 = text.index("\n", e0) + 1
    # keep the surrounding clang-format guards, which sit just outside the markers
    new = text[:b0] + block + "\n" + text[e0:]
    with open(HEADER, "w") as f:
        f.write(new)
    logger.info(f"Patched {HEADER}")


def closure(hd, hm, regions):
    """Per charge, compare the map (evaluated per cell in module phi) against
    the directly measured vertex-phi SF.

    In the analysis only passing muons are kept and each is downweighted by the
    map SF, so the map must reproduce the truth SF = eps_data/eps_MC cell by
    cell. Here that truth is measured in vertex phi, separately per charge, and
    the map is looked up at phi' = phi - q*C/pt. The residual, averaged over
    (phi, pt) weighted by the muons being corrected, should be ~1 in every eta
    bin -- unlike a charge-averaged vertex-phi correction, which is off because
    the hole sits at opposite phi offsets for the two charges.
    """

    def sf_lookup(eta, phiprime):
        for _, i0, i1, j0, j1, block, eta_e, phi_e in regions:
            ie = int(np.floor((eta - eta_e[i0]) / (eta_e[1] - eta_e[0])))
            ip = int(np.floor((phiprime - phi_e[j0]) / (phi_e[1] - phi_e[0])))
            if 0 <= ie < block.shape[0] and 0 <= ip < block.shape[1]:
                return block[ie, ip]
        return 1.0

    eta_c = hd.axes["eta"].centers
    phi_c = hd.axes["phi"].centers
    pt_c = hd.axes["pt"].centers
    vd, vm = hd.values(), hm.values()
    twopi = 2 * np.pi
    eta_boxes = [(e0, e1) for (e0, e1, _, _) in SEARCH_BOXES.values()]

    logger.info(
        "Closure: effective correction, measured (needed) -> residual after map"
    )
    for iq, q in ((0, -1), (1, +1)):
        # per-cell truth SF = eff_data/eff_MC (an efficiency ratio, so non-CVH
        # data/MC differences cancel), and the map SF at phi'(q,pt)
        ed, _ = efficiency(vd[:, :, :, iq, :])  # pt, eta, phi
        em, _ = efficiency(vm[:, :, :, iq, :])
        em = np.where(em > 0, em, 1.0)
        sf_meas = ed / em
        w = vm[:, :, :, iq, 1]  # passing MC, the muons being corrected

        for je, eta in enumerate(eta_c):
            if not any(e0 < eta < e1 for e0, e1 in eta_boxes):
                continue
            wsum = w[:, je, :].sum()
            if wsum < 500:
                continue
            need = (w[:, je, :] * sf_meas[:, je, :]).sum() / wsum
            applied = np.zeros(1)
            sf_map = np.array(
                [
                    [
                        sf_lookup(eta, (phi - q * BENDING_C / pt) % twopi)
                        for phi in phi_c
                    ]
                    for pt in pt_c
                ]
            )
            appl = (w[:, je, :] * sf_map).sum() / wsum
            resid = need / appl
            if abs(need - 1) > 1e-3 or abs(resid - 1) > 1e-3:
                logger.info(
                    f"  q={q:+d} eta~{eta:+.2f}: needed {need:.4f}, "
                    f"applied {appl:.4f}, residual {resid:.4f}"
                )


def main():
    args = parse_args()
    logging.setup_logger(__file__, args.verbose)
    hd, hm = load(args.input)
    regions = build_regions(hd, hm, args.nSigma)
    block = render(regions)
    if args.dry_run:
        print(block)
    else:
        patch_header(block)
    closure(hd, hm, regions)


if __name__ == "__main__":
    main()
