"""Generate the CVH-refit efficiency scale-factor map for muon_efficiencies_cvh.hpp.

Reads the 'cvhEfficiency' histogram produced by
'scripts/histmakers/mz_dilepton.py --cvhEfficiencyHists' (axes
pt, eta, phi, charge, passCVH, in uncorrected muon kinematics) and writes the
measured SF = eps_data / eps_MC as a map in (eta, phi'), where

    phi' = phi - q * C / pt

is the muon azimuth at the module rather than at the vertex (C = 0.3*B*r/2 is
the track bending, see muon_efficiencies_cvh.hpp). Working in phi' makes the map
charge- and pt-independent: the two charges' holes, which sit at opposite +-C/pt
offsets in vertex phi, line up, so a single map applies to everything without
smearing the sharp module edges.

C belongs to the module, not to the analysis -- it is set by the radius the
track has reached when it crosses it -- so every hotspot carries its own, and
--fitBending measures each one from the charge splitting of its own hole.

The affected cells (a couple of small rectangles) are written into the
GENERATED block of the header in place; everywhere else the SF is exactly 1.
The map is what the histmakers apply by default ('--cvhBadModules sf'). The
bounding rectangles are written out as well, for the cross-check that drops the
affected phase space geometrically instead ('--cvhBadModules veto').

Usage:
    python scripts/analysisTools/w_mass_13TeV/make_cvh_efficiency_sf.py -i <cvhEfficiency.hdf5>
"""

import argparse
import os
from collections import namedtuple

import h5py
import numpy as np

from wums import ioutils, logging

logger = logging.child_logger(__name__)

# this file sits in scripts/analysisTools/w_mass_13TeV/, so the repo root is
# four levels up
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), *[os.pardir] * 3))
HEADER = os.path.join(
    REPO_ROOT, "wremnants/production/include/muon_efficiencies_cvh.hpp"
)

# C = 0.3*B*r/2 [rad*GeV] with B = 3.8 T, so C = 0.57 * r[m]
BENDING_PER_METRE = 0.3 * 3.8 / 2.0

# Search boxes (module frame) in which to keep the correction; everything
# outside stays at SF = 1. Each entry is the box (eta_lo, eta_hi, phiprime_lo,
# phiprime_hi) plus 'bendC', the track bending to THAT module's radius, which is
# what maps vertex phi to module phi'. The two hotspots need not sit at the same
# radius, so they get their own value.
#
# Both are taken from the tracker geometry, C = 0.3*B*rho/2 at the module's own
# rho (runtree of the alignment study, ZMass/calibration_studies/
# module_level_corrections). That beats fitting them: --fitBending returns
# 0.214 +- 0.014 for hotspot1 and 0.119 +- 0.081 for hotspot2, both consistent
# with the geometry but the second one carrying no information.
#
# Caveat for hotspot2: it is a disk module, so the crossing radius varies across
# its radial extent (~32-42 cm, which is what gives the hole its eta width) and
# with it C, by about +-13%. bendC is the value at the module centre; the
# residual eta dependence is far below the phi binning.
SEARCH_BOXES = {
    "hotspot1 (TIB-L2, detId 369141860)": {
        "box": (-0.10, 0.85, 0.60, 1.20),
        "bendC": 0.199,  # rho = 34.96 cm (z = +14.8, eta +0.41, phi 0.898)
    },
    "hotspot2 (TID-2, detId 402666798)": {
        "box": (-1.95, -1.40, 4.90, 5.55),
        "bendC": 0.211,  # rho = 36.97 cm (z = -95.4, eta -1.68, phi 5.230)
    },
}

Region = namedtuple("Region", "name i0 i1 j0 j1 sf eta_e phi_e bendC")


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
        "--fitBending",
        action="store_true",
        help="""Use the bending constant fitted from each hole's charge splitting
        instead of the geometric one configured in SEARCH_BOXES. The fit is
        reported either way, as a cross-check of the module identification and as
        the only handle on the radius if a new hotspot shows up unidentified.""",
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


def to_module_phi(h, bendC):
    """Collapse (pt, eta, phi, charge, pass) -> (eta, phi', pass).

    Each (pt, charge) slice is shifted in phi by -q*C/pt with a linear sub-bin
    redistribution (the shift is well below one phi bin over 25-65 GeV), then
    summed. This is the histogram-level equivalent of filling phi' per muon.

    bendC is the bending to the radius of the module being looked at, so a map
    built here is only valid for the region it was built for.
    """
    v = h.values()  # pt, eta, phi, charge, pass
    pt_c = h.axes["pt"].centers
    phi_w = h.axes["phi"].edges[1] - h.axes["phi"].edges[0]
    nphi = len(h.axes["phi"])
    charges = {0: -1.0, 1: +1.0}  # charge-axis index -> muon charge

    out = np.zeros((len(h.axes["eta"]), nphi, 2))
    for iq, q in charges.items():
        for ipt, pt in enumerate(pt_c):
            f = -q * bendC / pt / phi_w  # shift in units of phi bins
            assert abs(f) < 1.0, f"sub-bin shift assumption violated: |f|={abs(f)}"
            step = 1 if f > 0 else -1
            src = v[ipt, :, :, iq, :]  # eta, phi, pass
            out += (1.0 - abs(f)) * src
            out += abs(f) * np.roll(src, step, axis=1)
    return out


def fit_bending_C(h, box, min_excess=100.0):
    """Measure C = 0.3*B*r/2 for one hotspot from the charge splitting of its hole.

    The dead module sits at a fixed azimuth on the detector; a muon of charge q
    and momentum pt reaches it at vertex phi = phi_module + q*C/pt. So the mean
    vertex phi of the muons whose refit failed, taken in a window around the
    hole, moves linearly in q/pt with slope C, and a straight-line fit over the
    (charge, pt) cells returns it.

    The flat component of refit failures is subtracted first, using the failure
    rate measured at the same eta outside the hole. It does not move with q/pt,
    so leaving it in would dilute the slope by the signal purity -- which
    matters for a shallow hole even if it is negligible for a deep one.

    Takes the DATA histogram (the holes are a data-only alignment effect).
    Returns (C, sigma_C, chi2/ndf, npoints), or None if there is too little
    excess to fit.
    """
    e0, e1, p0, p1 = box
    v = h.values()  # pt, eta, phi, charge, pass
    eta_c, phi_c, pt_c = (h.axes[a].centers for a in ("eta", "phi", "pt"))
    dphi = phi_c[1] - phi_c[0]

    ineta = (eta_c > e0) & (eta_c < e1)
    # window around the hole, one bin wider so the shifted tails stay inside
    inphi = (phi_c > p0 - dphi) & (phi_c < p1 + dphi)

    sel = v[:, ineta]  # pt, eta(sel), phi, charge, pass
    n_fail = sel[..., 0].sum(axis=1)  # pt, phi, charge
    n_tot = sel.sum(axis=(1, 4))  # pt, phi, charge

    xs, ys, es = [], [], []
    for iq, q in ((0, -1.0), (1, +1.0)):
        for ipt, pt in enumerate(pt_c):
            f, n = n_fail[ipt, :, iq], n_tot[ipt, :, iq]
            out = ~inphi
            rate = f[out].sum() / max(n[out].sum(), 1.0)
            excess = np.clip(f[inphi] - rate * n[inphi], 0.0, None)
            tot = excess.sum()
            if tot < min_excess:
                continue
            w = excess / tot
            mean = (w * phi_c[inphi]).sum()
            var = (w * (phi_c[inphi] - mean) ** 2).sum()
            # error on the mean of the excess; fall back to a flat bin if the
            # excess sits in a single phi bin
            err = np.sqrt(var / tot) if var > 0 else dphi / np.sqrt(12.0 * tot)
            xs.append(q / pt)
            ys.append(mean)
            es.append(err)

    x, y, e = (np.asarray(a) for a in (xs, ys, es))
    if len(x) < 4 or np.ptp(x) == 0 or np.any(e <= 0):
        return None

    w = 1.0 / e**2
    Sw, Sx, Sy = w.sum(), (w * x).sum(), (w * y).sum()
    Sxx, Sxy = (w * x * x).sum(), (w * x * y).sum()
    den = Sw * Sxx - Sx * Sx
    if den <= 0:
        return None
    C = (Sw * Sxy - Sx * Sy) / den
    intercept = (Sxx * Sy - Sx * Sxy) / den
    sigma = np.sqrt(Sw / den)
    chi2 = (w * (y - intercept - C * x) ** 2).sum()
    return C, sigma, chi2 / (len(x) - 2), len(x)


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


def bending_constants(hd, use_fit):
    """Bending constant per region: the configured one, and the fitted one.

    The fit is always run and reported -- it is the only handle on the module
    radius when the detId is not known -- but it is only used when asked for,
    so that regenerating the map does not silently change the frame it is
    measured in.
    """
    out = {}
    for name, cfg in SEARCH_BOXES.items():
        bendC = cfg["bendC"]
        fit = fit_bending_C(hd, cfg["box"])
        if fit is None:
            logger.warning(
                f"{name}: not enough failing muons to fit the bending, "
                f"keeping the configured C = {bendC:.3f}"
            )
        else:
            C, sigma, chi2ndf, npts = fit
            logger.info(
                f"{name}: fitted C = {C:.3f} +- {sigma:.3f} rad*GeV "
                f"(r = {100*C/BENDING_PER_METRE:.0f} +- "
                f"{100*sigma/BENDING_PER_METRE:.0f} cm), "
                f"chi2/ndf = {chi2ndf:.2f}, {npts} points; configured {bendC:.3f}"
            )
            if use_fit:
                if sigma > abs(C):
                    logger.warning(
                        f"{name}: fitted bending is compatible with zero, "
                        f"keeping the configured C = {bendC:.3f}"
                    )
                else:
                    bendC = C
        out[name] = bendC
    return out


def build_regions(hd, hm, nsigma, bendCs):
    """Return the list of Region rectangles, each in its own module frame."""
    eta_e = hd.axes["eta"].edges
    phi_e = hd.axes["phi"].edges
    eta_c = hd.axes["eta"].centers
    phi_c = hd.axes["phi"].centers

    regions = []
    for name, cfg in SEARCH_BOXES.items():
        e0, e1, p0, p1 = cfg["box"]
        bendC = bendCs[name]
        sf, err = scale_factor(to_module_phi(hd, bendC), to_module_phi(hm, bendC))
        # holes only: never upweight, and drop cells not significantly below 1
        keep = (sf < 1.0 - nsigma * err) & (err > 0)
        sf = np.where(keep, np.clip(sf, 0.0, 1.0), 1.0)

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
        regions.append(Region(name, i0, i1, j0, j1, block, eta_e, phi_e, bendC))
        logger.info(
            f"{name}: eta [{eta_e[i0]:.2f},{eta_e[i1+1]:.2f}] "
            f"phi' [{phi_e[j0]:.3f},{phi_e[j1+1]:.3f}] at C = {bendC:.3f}, "
            f"{block.shape[0]}x{block.shape[1]} cells, min SF = {block.min():.3f}"
        )
    return regions


def render(regions):
    """Render the C++ GENERATED block."""
    lines = [
        "// BEGIN GENERATED -- do not edit by hand, see make_cvh_efficiency_sf.py",
        "",
    ]
    for k, r in enumerate(regions):
        name, i0, i1, j0, j1, block, eta_e, phi_e = (
            r.name,
            r.i0,
            r.i1,
            r.j0,
            r.j1,
            r.sf,
            r.eta_e,
            r.phi_e,
        )
        neta, nphi = block.shape
        eta_lo, eta_w = eta_e[i0], eta_e[1] - eta_e[0]
        phi_lo, phi_w = phi_e[j0], phi_e[1] - phi_e[0]
        lines.append(f"// {name}")
        lines.append(
            f"//   eta in [{eta_e[i0]:.2f}, {eta_e[i1+1]:.2f}), "
            f"phi' in [{phi_e[j0]:.3f}, {phi_e[j1+1]:.3f}), "
            f"C = {r.bendC:.3f} (r = {100*r.bendC/BENDING_PER_METRE:.0f} cm)"
        )
        lines.append(f"inline const double cvh_sf_region{k}[{neta} * {nphi}] = {{")
        for ie in range(neta):
            row = ", ".join(f"{block[ie, ip]:.3f}" for ip in range(nphi))
            eta_center = eta_lo + (ie + 0.5) * eta_w
            lines.append(f"    {row}, // eta ~ {eta_center:+.3f}")
        # (trailing comma on the last row is legal in a C++ aggregate initializer)
        lines.append("};")
        lines.append("")

    lines.append(
        "// last column of each row is bendC = 0.3*B*r/2, the bending to that module"
    )
    lines.append("inline const CvhSFRegion cvh_sf_regions[] = {")
    for k, r in enumerate(regions):
        neta, nphi = r.sf.shape
        eta_lo, eta_w = r.eta_e[r.i0], r.eta_e[1] - r.eta_e[0]
        phi_lo, phi_w = r.phi_e[r.j0], r.phi_e[1] - r.phi_e[0]
        lines.append(
            f"    {{{eta_lo:.4f}, {eta_w:.4f}, {neta}, "
            f"{phi_lo:.4f}, {phi_w:.4f}, {nphi}, {r.bendC:.3f}, cvh_sf_region{k}}},"
        )
    lines.append("};")
    lines.append("")

    lines.append(
        "// Bounding rectangles of the same regions, used by the geometric veto."
    )
    lines.append("inline const CvhModuleRegion cvh_bad_modules[] = {")
    for r in regions:
        lines.append(
            f"    {{{r.eta_e[r.i0]:.4f}, {r.eta_e[r.i1+1]:.4f}, "
            f"{r.phi_e[r.j0]:.4f}, {r.phi_e[r.j1+1]:.4f}, {r.bendC:.3f}}}, // {r.name}"
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


def veto_summary(hd, hm, regions):
    """Cost and closure of the geometric veto.

    The alternative to weighting MC by the SF is to drop, in data and MC alike,
    every event with a muon inside one of the rectangles. Report what that costs
    in acceptance and what data/MC efficiency mismatch is left outside them.

    Each rectangle lives in its own module frame; the frames differ by less than
    a tenth of a phi bin, so the plane-wide sums below are taken in the first
    region's frame, while the in-box sums use each region's own.
    """
    if not regions:
        return
    ref = regions[0].bendC
    md, mm = to_module_phi(hd, ref), to_module_phi(hm, ref)
    sf, _ = scale_factor(md, mm)
    nd_pass, nm_pass = md[..., 1], mm[..., 1]

    inbox = np.zeros(sf.shape, dtype=bool)
    nd_box = nm_box = 0.0
    for r in regions:
        cell = (slice(r.i0, r.i1 + 1), slice(r.j0, r.j1 + 1))
        own_d = md if r.bendC == ref else to_module_phi(hd, r.bendC)
        own_m = mm if r.bendC == ref else to_module_phi(hm, r.bendC)
        nd_box += own_d[..., 1][cell].sum()
        nm_box += own_m[..., 1][cell].sum()
        inbox[cell] = True

    w = nm_pass * ~inbox
    # MC-yield-weighted mean of (1 - SF): the efficiency mismatch that the veto
    # leaves behind, to be compared with the same number over the full plane
    resid = (w * (1.0 - sf)).sum() / w.sum()
    total = (nm_pass * (1.0 - sf)).sum() / nm_pass.sum()
    logger.info(
        f"Veto: removes {100*nd_box/nd_pass.sum():.2f}% of the data "
        f"muons and {100*nm_box/nm_pass.sum():.2f}% of the MC ones"
    )
    logger.info(
        f"Veto: mean data/MC efficiency mismatch <1-SF> = {total:.2e} over the full "
        f"plane, {resid:.2e} outside the vetoed rectangles"
    )


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

    twopi = 2 * np.pi

    def sf_lookup(eta, phi, q, pt):
        """Mirror of wrem::cvh_muon_sf, including the per-region unbending."""
        for r in regions:
            ie = int(np.floor((eta - r.eta_e[r.i0]) / (r.eta_e[1] - r.eta_e[0])))
            if not 0 <= ie < r.sf.shape[0]:
                continue
            phiprime = (phi - q * r.bendC / pt) % twopi
            ip = int(np.floor((phiprime - r.phi_e[r.j0]) / (r.phi_e[1] - r.phi_e[0])))
            if 0 <= ip < r.sf.shape[1]:
                return r.sf[ie, ip]
        return 1.0

    eta_c = hd.axes["eta"].centers
    phi_c = hd.axes["phi"].centers
    pt_c = hd.axes["pt"].centers
    vd, vm = hd.values(), hm.values()
    eta_boxes = [cfg["box"][:2] for cfg in SEARCH_BOXES.values()]

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
            sf_map = np.array(
                [[sf_lookup(eta, phi, q, pt) for phi in phi_c] for pt in pt_c]
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
    bendCs = bending_constants(hd, args.fitBending)
    regions = build_regions(hd, hm, args.nSigma, bendCs)
    block = render(regions)
    if args.dry_run:
        print(block)
    else:
        patch_header(block)
    veto_summary(hd, hm, regions)
    closure(hd, hm, regions)


if __name__ == "__main__":
    main()
