"""Plot the single-muon CVH refit efficiency in data and MC.

Input is the 'cvhEfficiency' histogram produced by
'scripts/histmakers/mz_dilepton.py --cvhEfficiencyHists' (see there for the
required options), with axes (pt, eta, phi, charge, passCVH) in uncorrected
muon kinematics.

The measured scale factor SF = eff(data)/eff(MC) per (eta,phi) cell is what
'wremnants/production/include/muon_efficiencies_cvh.hpp' has to reproduce; the
hard-coded version of it is overlaid for comparison.
"""

import argparse
import datetime
import os

import h5py
import make_cvh_efficiency_sf as cvhgen  # sibling module, same directory
import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np

from wums import ioutils, logging, output_tools, plot_tools

hep.style.use(hep.style.ROOT)
logger = logging.child_logger(__name__)

# the two localized CVH refit efficiency holes seen in data, as (eta, phi) boxes
HOTSPOTS = {
    "hotspot1": dict(
        eta=(-0.10, 0.85),
        phi=(0.65, 1.15),
        label=r"TIB-L2 glued module (detId 369141860)",
    ),
    "hotspot2": dict(
        eta=(-1.95, -1.40),
        phi=(4.95, 5.50),
        label=r"TID-2 out-of-plane module (detId 402666798)",
    ),
}

# the original 1D-in-eta approximation (flat over phi 0.80-1.00), now superseded
# by the 2D module-phi map in muon_efficiencies_cvh.hpp; kept here only as a
# reference curve on the vertex-phi diagnostic
SF_CODE_ETA_EDGES = np.arange(0.05, 0.701, 0.05)
SF_CODE_VALUES = np.array(
    [
        0.982,
        0.933,
        0.910,
        0.823,
        0.718,
        0.645,
        0.585,
        0.606,
        0.653,
        0.777,
        0.847,
        0.955,
        0.980,
    ]
)
SF_CODE_PHI_GATE = (0.80, 1.00)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="hdf5 output of mz_dilepton.py --cvhEfficiencyHists",
    )
    parser.add_argument(
        "--hist",
        type=str,
        default="cvhEfficiency",
        help="Name of the histogram to read",
    )
    parser.add_argument("--outpath", type=str, default=None, help="Output base dir")
    parser.add_argument("--postfix", type=str, default="", help="Suffix for filenames")
    parser.add_argument("--title", type=str, default="CMS")
    parser.add_argument("--subtitle", type=str, default="Work in progress")
    parser.add_argument("--lumi", type=float, default=16.8, help="Luminosity in /fb")
    parser.add_argument(
        "--nSigma",
        type=float,
        default=3.0,
        help="Significance threshold for the applied-map plot (matches the generator)",
    )
    parser.add_argument("--titlePos", type=int, default=0)
    parser.add_argument("-v", "--verbose", type=int, default=3)
    return parser.parse_args()


def default_outpath(script_path):
    stem = os.path.splitext(os.path.basename(script_path))[0]
    repo_root = os.path.dirname(os.path.abspath(script_path))
    while repo_root != "/" and not os.path.isdir(os.path.join(repo_root, ".git")):
        repo_root = os.path.dirname(repo_root)
    repo_name = os.path.basename(repo_root) if repo_root != "/" else "plots"
    today = datetime.date.today().strftime("%y%m%d")
    return os.path.expanduser(f"~/public_html/{repo_name}/{today}_{stem}/")


def load(fname, hname):
    """Return the (data, total MC prediction) 'cvhEfficiency' histograms."""
    hdata, hmc = None, None
    with h5py.File(fname, "r") as f:
        for proc in f.keys():
            if proc == "meta_info":
                continue
            res = ioutils.pickle_load_h5py(f[proc])
            h = res["output"][hname].get()
            # aggregated groups drop the 'is_data' flag, they are all MC
            if res["dataset"].get("is_data", False):
                hdata = h if hdata is None else hdata + h
            else:
                hmc = h if hmc is None else hmc + h
    logger.info(f"Loaded {hname} with axes {[a.name for a in hdata.axes]}")
    return hdata, hmc


def efficiency(v):
    """(eff, err) from an array with the passCVH axis last."""
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
    """SF = eff_data/eff_MC with the (dominant) data uncertainty."""
    ed, dd = efficiency(vd)
    em, _ = efficiency(vm)
    em = np.where(em > 0, em, 1.0)
    return ed / em, dd / em


def sf_hardcoded(eta, phi):
    """The correction currently applied to MC, evaluated on a grid."""
    eta, phi = np.broadcast_arrays(np.asarray(eta), np.asarray(phi))
    sf = np.ones(eta.shape)
    ibin = np.digitize(eta, SF_CODE_ETA_EDGES) - 1
    inside = (
        (phi >= SF_CODE_PHI_GATE[0])
        & (phi < SF_CODE_PHI_GATE[1])
        & (ibin >= 0)
        & (ibin < len(SF_CODE_VALUES))
    )
    sf[inside] = SF_CODE_VALUES[ibin[inside]]
    return sf


def decorate(ax, args, fontsize=17):
    """CMS label, with an explicit font size so it fits next to the lumi text."""
    hep.cms.label(
        ax=ax,
        label=args.subtitle,
        data=True,
        lumi=args.lumi,
        loc=args.titlePos,
        fontsize=fontsize,
        exp=args.title,
    )


def plot_map(outdir, args, hd, hm, name):
    """2D map of the measured SF over the full (eta, phi) plane."""
    sf, _ = scale_factor(hd.values().sum(axis=(0, 3)), hm.values().sum(axis=(0, 3)))
    eta_e, phi_e = hd.axes["eta"].edges, hd.axes["phi"].edges

    fig, ax = plt.subplots(figsize=(10, 6))
    mesh = ax.pcolormesh(
        eta_e, phi_e, sf.T, cmap="viridis", vmin=0.35, vmax=1.0, rasterized=True
    )
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$\epsilon_\mathrm{data}\,/\,\epsilon_\mathrm{MC}$", labelpad=10)
    for hs in HOTSPOTS.values():
        ax.add_patch(
            mpl.patches.Rectangle(
                (hs["eta"][0], hs["phi"][0]),
                hs["eta"][1] - hs["eta"][0],
                hs["phi"][1] - hs["phi"][0],
                fill=False,
                edgecolor="red",
                lw=1.5,
            )
        )
    ax.set_xlabel(r"muon $\eta$")
    ax.set_ylabel(r"muon $\phi$")
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
    decorate(ax, args)
    plot_tools.save_pdf_and_png(outdir, name, fig=fig)
    plt.close(fig)


def plot_applied_map(outdir, args, hd, hm, name):
    """2D map of the applied SF in (eta, phi'), i.e. the actual correction that
    goes into muon_efficiencies_cvh.hpp. Built with the generator's own logic
    (module-phi folding + significance clamp) so it cannot drift from the header.
    """
    # each region is folded with the bending to its own module radius, as the
    # header does; pass the configured constants rather than refitting them here
    bendCs = {name: cfg["bendC"] for name, cfg in cvhgen.SEARCH_BOXES.items()}
    regions = cvhgen.build_regions(hd, hm, args.nSigma, bendCs)
    eta_e, phi_e = hd.axes["eta"].edges, hd.axes["phi"].edges
    sf = np.ones((len(hd.axes["eta"]), len(hd.axes["phi"])))
    for r in regions:
        sf[r.i0 : r.i1 + 1, r.j0 : r.j1 + 1] = r.sf

    fig, ax = plt.subplots(figsize=(10, 6))
    mesh = ax.pcolormesh(
        eta_e,
        phi_e,
        np.ma.masked_equal(sf.T, 1.0),  # leave SF=1 cells blank
        cmap="viridis",
        vmin=0.35,
        vmax=1.0,
        rasterized=True,
    )
    mesh.cmap.set_bad("#f0f0f0")
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(
        r"applied SF  $\epsilon_\mathrm{data}\,/\,\epsilon_\mathrm{MC}$", labelpad=10
    )
    ax.set_xlabel(r"muon $\eta$")
    ax.set_ylabel(r"muon $\phi'$ (module frame, $\phi - qC/p_\mathrm{T}$)")
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
    decorate(ax, args)
    plot_tools.save_pdf_and_png(outdir, name, fig=fig)
    plt.close(fig)

    # zoomed panels, one per region, with the SF value printed in each cell
    for k, r in enumerate(regions):
        rname, i0, j0, block, ee, pe = r.name, r.i0, r.j0, r.sf, r.eta_e, r.phi_e
        i1, j1 = r.i1, r.j1
        fig, ax = plt.subplots(figsize=(1.6 + 1.1 * block.shape[0], 6))
        mesh = ax.pcolormesh(
            ee[i0 : i1 + 2],
            pe[j0 : j1 + 2],
            np.ma.masked_equal(block.T, 1.0),
            cmap="viridis",
            vmin=0.35,
            vmax=1.0,
        )
        mesh.cmap.set_bad("#f0f0f0")
        for ie in range(block.shape[0]):
            for ip in range(block.shape[1]):
                if block[ie, ip] == 1.0:
                    continue
                ax.text(
                    ee[i0 + ie] + 0.5 * (ee[1] - ee[0]),
                    pe[j0 + ip] + 0.5 * (pe[1] - pe[0]),
                    f"{block[ie, ip]:.2f}",
                    ha="center",
                    va="center",
                    fontsize=9,
                    color="white" if block[ie, ip] < 0.75 else "black",
                )
        fig.colorbar(mesh, ax=ax, pad=0.02, label="applied SF")
        ax.set_xlabel(r"muon $\eta$")
        ax.set_ylabel(r"muon $\phi'$ (module frame)")
        ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
        ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
        decorate(ax, args)
        fig.suptitle(rname, fontsize=14)
        znm = "_".join(filter(None, [name, f"region{k}", args.postfix]))
        plot_tools.save_pdf_and_png(outdir, znm, fig=fig)
        output_tools.write_logfile(
            outdir, znm, args=args, wd=os.path.dirname(os.path.abspath(__file__))
        )
        plt.close(fig)


def plot_hotspot(outdir, args, hd, hm, key, name):
    """Zoomed SF map and the SF vs eta in each phi slice of one hotspot."""
    hs = HOTSPOTS[key]
    eta_c, phi_c = hd.axes["eta"].centers, hd.axes["phi"].centers
    eta_e, phi_e = hd.axes["eta"].edges, hd.axes["phi"].edges
    ie = np.where((eta_c > hs["eta"][0]) & (eta_c < hs["eta"][1]))[0]
    ip = np.where((phi_c > hs["phi"][0]) & (phi_c < hs["phi"][1]))[0]

    vd = hd.values().sum(axis=(0, 3))
    vm = hm.values().sum(axis=(0, 3))
    sf, err = scale_factor(vd, vm)

    fig, axs = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw=dict(wspace=0.45))

    ax = axs[0]
    mesh = ax.pcolormesh(
        eta_e[ie[0] : ie[-1] + 2],
        phi_e[ip[0] : ip[-1] + 2],
        sf[np.ix_(ie, ip)].T,
        cmap="viridis",
        vmin=0.35,
        vmax=1.0,
    )
    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(r"$\epsilon_\mathrm{data}\,/\,\epsilon_\mathrm{MC}$", labelpad=10)
    ax.set_xlabel(r"muon $\eta$")
    ax.set_ylabel(r"muon $\phi$")
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())

    ax = axs[1]
    for k, j in enumerate(ip):
        if sf[ie, j].min() > 0.995:
            continue
        ax.errorbar(
            eta_c[ie],
            sf[ie, j],
            yerr=err[ie, j],
            marker="o",
            ms=3,
            lw=1.2,
            label=rf"$\phi \in [{phi_e[j]:.2f}, {phi_e[j+1]:.2f}]$",
        )
    if key == "hotspot1":
        # the hard-coded correction is flat in phi inside its gate
        ax.plot(
            eta_c[ie],
            sf_hardcoded(eta_c[ie], np.full(len(ie), 0.9)),
            color="k",
            ls="--",
            lw=2,
            label="original 1D approx. (superseded)",
        )
    ax.axhline(1.0, color="gray", lw=1, ls=":")
    ax.set_xlabel(r"muon $\eta$")
    ax.set_ylabel(r"$\epsilon_\mathrm{data}\,/\,\epsilon_\mathrm{MC}$")
    ax.set_ylim(0.2, 1.1)
    ax.legend(fontsize=13, ncol=1, loc="lower left")
    decorate(ax, args)
    fig.suptitle(hs["label"], fontsize=16)
    plot_tools.save_pdf_and_png(outdir, name, fig=fig)
    plt.close(fig)


def plot_stability(outdir, args, hd, hm, name):
    """SF integrated over each hotspot vs pt, separately per charge."""
    eta_c, phi_c = hd.axes["eta"].centers, hd.axes["phi"].centers
    pt_c, pt_e = hd.axes["pt"].centers, hd.axes["pt"].edges
    E, P = np.meshgrid(eta_c, phi_c, indexing="ij")

    fig, ax = plt.subplots(figsize=(8, 6))
    for key, hs in HOTSPOTS.items():
        msk = (
            (E > hs["eta"][0])
            & (E < hs["eta"][1])
            & (P > hs["phi"][0])
            & (P < hs["phi"][1])
        )
        for q, (ql, ls) in enumerate((("$\\mu^-$", "-"), ("$\\mu^+$", "--"))):
            vd = hd.values()[:, :, :, q, :][:, msk].sum(axis=1)
            vm = hm.values()[:, :, :, q, :][:, msk].sum(axis=1)
            sf, err = scale_factor(vd, vm)
            ax.errorbar(
                pt_c, sf, yerr=err, marker="o", ms=4, ls=ls, label=f"{key}, {ql}"
            )
    ax.set_xlim(pt_e[0], pt_e[-1])
    ax.set_xlabel(r"muon $p_\mathrm{T}$ (GeV)")
    ax.set_ylabel(r"$\epsilon_\mathrm{data}\,/\,\epsilon_\mathrm{MC}$ in hotspot")
    ax.legend(fontsize=14, ncol=2)
    decorate(ax, args)
    plot_tools.save_pdf_and_png(outdir, name, fig=fig)
    plt.close(fig)


def hole_position(hd, key, ipt=None, q=None):
    """Failure-weighted mean phi of a hotspot, i.e. where the hole sits."""
    hs = HOTSPOTS[key]
    eta_c, phi_c = hd.axes["eta"].centers, hd.axes["phi"].centers
    E, P = np.meshgrid(eta_c, phi_c, indexing="ij")
    msk = (
        (E > hs["eta"][0])
        & (E < hs["eta"][1])
        & (P > hs["phi"][0])
        & (P < hs["phi"][1])
    )
    v = hd.values()[..., 0]  # failing muons
    v = v[ipt] if ipt is not None else v.sum(axis=0)
    v = v[..., q] if q is not None else v.sum(axis=-1)
    w = v[msk]
    n = w.sum()
    mean = (w * P[msk]).sum() / n
    rms = np.sqrt((w * (P[msk] - mean) ** 2).sum() / n)
    return mean, rms / np.sqrt(n), n


def plot_position(outdir, args, hd, key, name):
    """Position of the hole in phi vs pt and charge.

    The muon phi in the histogram is the one at the vertex, while the module
    sits at a fixed phi in the detector, so the hole is seen displaced by the
    track bending, +-C/pt for the two charges. This is what makes the hotspot
    scale factor look charge dependent when parametrized in the vertex phi.
    """
    pt_c = hd.axes["pt"].centers

    fig, ax = plt.subplots(figsize=(8, 6))
    x, y, ey, qs = [], [], [], []
    for q, (ql, col) in enumerate(((r"$\mu^-$", "C0"), (r"$\mu^+$", "C1"))):
        m = np.array([hole_position(hd, key, ipt=i, q=q) for i in range(len(pt_c))])
        ok = m[:, 2] > 50
        ax.errorbar(
            pt_c[ok],
            m[ok, 0],
            yerr=m[ok, 1],
            marker="o",
            ms=5,
            ls="none",
            color=col,
            label=ql,
        )
        x.append(pt_c[ok])
        y.append(m[ok, 0])
        ey.append(m[ok, 1])
        qs.append(np.full(ok.sum(), -1 if q == 0 else 1))
    x, y, ey, qs = (np.concatenate(a) for a in (x, y, ey, qs))

    # phi_vertex = phi_module + q*C/pt, C = 0.3*B*r/2 for a track reaching radius r
    def chi2(p):
        return (((y - (p[0] + qs * p[1] / x)) / ey) ** 2).sum()

    from scipy import optimize

    fit = optimize.minimize(chi2, [y.mean(), 0.2])
    phi0, C = fit.x
    radius = 2 * C / (0.3 * 3.8)
    xx = np.linspace(x.min(), x.max(), 100)
    for sign, col in ((-1, "C0"), (1, "C1")):
        ax.plot(xx, phi0 + sign * C / xx, color=col, lw=1.5)
    ax.axhline(phi0, color="k", ls=":", lw=1)
    ax.text(
        0.04,
        0.06,
        rf"$\phi = \phi_0 + q\,C/p_\mathrm{{T}}$, $C = {C:.3f}$ rad GeV"
        + "\n"
        + rf"$\rightarrow r = 2C/(0.3B) = {100*radius:.0f}$ cm, $\chi^2/\mathrm{{ndf}} = {fit.fun:.1f}/{len(x)-2}$",
        transform=ax.transAxes,
        fontsize=14,
    )
    ax.set_xlabel(r"muon $p_\mathrm{T}$ (GeV)")
    ax.set_ylabel(r"$\langle \phi \rangle$ of failed refits")
    ax.legend(fontsize=15, loc="center right")
    decorate(ax, args)
    fig.suptitle(HOTSPOTS[key]["label"], fontsize=15)
    plot_tools.save_pdf_and_png(outdir, name, fig=fig)
    plt.close(fig)
    logger.info(
        f"{key}: hole at phi0 = {phi0:.4f}, bending C = {C:.4f} rad GeV "
        f"(r = {100*radius:.1f} cm), chi2/ndf = {fit.fun:.1f}/{len(x)-2}"
    )


def print_summary(hd, hm):
    """Integrated efficiencies inside and outside the hotspots."""
    eta_c, phi_c = hd.axes["eta"].centers, hd.axes["phi"].centers
    E, P = np.meshgrid(eta_c, phi_c, indexing="ij")
    vd, vm = hd.values().sum(axis=(0, 3)), hm.values().sum(axis=(0, 3))
    rest = np.ones(E.shape, dtype=bool)
    regions = {}
    for key, hs in HOTSPOTS.items():
        msk = (
            (E > hs["eta"][0])
            & (E < hs["eta"][1])
            & (P > hs["phi"][0])
            & (P < hs["phi"][1])
        )
        regions[key] = msk
        rest &= ~msk
    regions["outside"] = rest
    regions["inclusive"] = np.ones(E.shape, dtype=bool)
    for key, msk in regions.items():
        d, m = vd[msk].sum(axis=0), vm[msk].sum(axis=0)
        sf, err = scale_factor(d, m)
        logger.info(
            f"{key:10s}: data fail {d[0]:8.0f}/{d.sum():10.0f} (eff {efficiency(d)[0]:.6f}), "
            f"MC fail {m[0]:7.1f}/{m.sum():10.0f} (eff {efficiency(m)[0]:.6f}), "
            f"SF = {sf:.5f} +/- {err:.5f}"
        )


def main():
    args = parse_args()
    logging.setup_logger(__file__, args.verbose)
    outdir = output_tools.make_plot_dir(args.outpath or default_outpath(__file__))
    logger.info(f"Writing plots to {outdir}")

    hd, hm = load(args.input, args.hist)
    print_summary(hd, hm)

    for name, fn in (
        ("cvhSF_map", lambda o, n: plot_map(o, args, hd, hm, n)),
        ("cvhSF_appliedMap", lambda o, n: plot_applied_map(o, args, hd, hm, n)),
        ("cvhSF_hotspot1", lambda o, n: plot_hotspot(o, args, hd, hm, "hotspot1", n)),
        ("cvhSF_hotspot2", lambda o, n: plot_hotspot(o, args, hd, hm, "hotspot2", n)),
        ("cvhSF_stability", lambda o, n: plot_stability(o, args, hd, hm, n)),
        ("cvhSF_position", lambda o, n: plot_position(o, args, hd, "hotspot1", n)),
    ):
        name = "_".join(filter(None, [name, args.postfix]))
        fn(outdir, name)
        output_tools.write_logfile(
            outdir, name, args=args, wd=os.path.dirname(os.path.abspath(__file__))
        )


if __name__ == "__main__":
    main()
