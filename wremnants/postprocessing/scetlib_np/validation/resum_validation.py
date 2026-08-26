"""Validate the SCETlib *resummed* σ_gen(λ_central) from two sources.

This isolates ONE comparison: the resummed (TMD) piece only — no fixed-order
matching, no response matrix, no MiNNLO. Two predictions of the same object,
the SCETlib resummed cross section at the datacard's λ_central, projected to
the model's gen grid (ptVGen, absYVGen), Q∈[Q_lo,Q_hi], |Y|-folded:

  (A) ParamModel  : our bt-grid integration, ``model.sigma_gen_central``
                    (Hankel transform of the cached I_pert/C_ν with λ_central
                    NP applied, then Q-integrated and rebinned).
  (B) correction  : the SCETlib spectrum-mode resummed histogram that the
                    theory correction is built from (the ``…_combined.pkl``
                    resummed input to make_theory_corr), read with
                    ``input_tools.read_scetlib_hist``.

Output: a ptVGen-projection ratio plot (B as reference, A as the red curve),
like scetlib_np.validation.histmaker_validation.

NOTE the resummed inputs live under /work — run inside a container that binds
it (the WRemnantsHelpers default singularity command binds /work; the WRemnants
ci/run_with_singularity.sh does NOT, so add /work to APPTAINER_BIND there).

Run, e.g.:
    source setup.sh
    python3 -m wremnants.postprocessing.scetlib_np.validation.resum_validation
"""

import argparse
import sys
import time

import numpy as np

from wremnants.postprocessing.scetlib_np.sigma_gen import _default_btgrid_dir

from rabbit.inputdata import FitInputData
from wums import boostHistHelpers as hh

from wremnants.utilities.io_tools import input_tools
from wremnants.postprocessing.scetlib_np.param_model import SCETlibNPParamModel
from wremnants.postprocessing.scetlib_np.validation_plots import (
    _merge_matrix,
    plot_ptll_ratio,
    summarize,
    tf_to_hist,
)

# ---- Defaults (FranksVals CT18Z N3p0LL setup) ----
# The resummed SCETlib spectrum input to the FranksVals theory correction:
RESUM_PKL = (
    "/work/submit/areimers/wmass/TheoryCorrections/SCETlib/"
    "com13_ct18z_newnps_n3+0ll_lattice_lambda4bugfix_franksvals_fine/"
    "inclusive_Z_COM13_CT18Z_N3+0LL_lattice_lambda4bugfix_franksvals_fine_combined.pkl"
)
SIGNAL_PROC = "Zmumu"
CHARGE = 0  # Z
Q_LO, Q_HI = 60.0, 120.0


def _central_var(h):
    """Select the central variation of a SCETlib hist (vars axis), if present."""
    if "vars" in h.axes.name:
        names = list(h.axes["vars"])
        idx = 0
        for cand in ("central", "pdf0", "nominal"):
            if cand in names:
                idx = names.index(cand)
                break
        h = h[{"vars": idx}]
    return h


def resum_from_correction(resum_pkl, gen_axes_meta, q_lo, q_hi, charge, tol=1e-6):
    """Read the SCETlib spectrum resum pkl and project onto the model gen grid.

    Returns an ndarray (NptVGen, NabsYVGen) matching ``gen_axes_meta`` order,
    integrated over Q∈[q_lo,q_hi] and |Y|-folded to the model's absYVGen range.
    """
    h = input_tools.read_scetlib_hist(resum_pkl, charge=charge)
    h = _central_var(h)
    h = hh.makeAbsHist(h, "Y")  # signed Y -> |Y|

    # integrate the Q window
    Qe = np.asarray(h.axes["Q"].edges, dtype=np.float64)
    qi = np.where(np.isclose(Qe, q_lo, atol=tol))[0]
    qj = np.where(np.isclose(Qe, q_hi, atol=tol))[0]
    # slice(...,sum) sums ONLY the in-range Q bins; a plain slice + {"Q":sum}
    # re-adds the sliced-out low-mass [10,60] bin via underflow (a ~4× leak).
    if qi.size and qj.size:
        h = h[{"Q": slice(int(qi[0]), int(qj[0]), sum)}]
    elif "Q" in h.axes.name:
        h = h[{"Q": sum}]
    if "charge" in h.axes.name:
        h = h[{"charge": sum}]

    v = h.project("qT", "absY").values(flow=False)  # (qT, absY), model-axis order
    ptv_edges = np.asarray(gen_axes_meta[0][1], dtype=np.float64)
    absy_edges = np.asarray(gen_axes_meta[1][1], dtype=np.float64)
    # _merge_matrix raises if a model edge isn't present in the pkl binning.
    Wp = _merge_matrix(h.axes["qT"].edges, ptv_edges, name="qT->ptVGen", tol=tol)
    Wa = _merge_matrix(h.axes["absY"].edges, absy_edges, name="absY->absYVGen", tol=tol)
    return Wp @ v @ Wa.T  # (ptVGen, absYVGen)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--datacard", required=True, help="fit-input hdf5 (λ_central + binning)")
    p.add_argument("--btgrid", default=_default_btgrid_dir(), help="SCETlib bT-grid directory")
    p.add_argument("--resum-pkl", default=RESUM_PKL, help="SCETlib spectrum resummed pkl (correction input)")
    p.add_argument("--signal-proc", default=SIGNAL_PROC)
    p.add_argument("--charge", type=int, default=CHARGE, help="0 for Z")
    p.add_argument("--q-lo", type=float, default=Q_LO)
    p.add_argument("--q-hi", type=float, default=Q_HI)
    p.add_argument(
        "--plot-out",
        default="/home/submit/lavezzo/alphaS/WRemnantsHelpers/260601_debug/resum_comp.png",
        help="output path for the ptVGen ratio plot ('' to skip)",
    )
    args = p.parse_args(argv)

    print("=" * 70)
    print("SCETlib resummed σ_gen(λ_central): ParamModel bt-grid  vs  correction pkl")
    print("=" * 70)
    print(f"  datacard  : {args.datacard}")
    print(f"  btgrid    : {args.btgrid}")
    print(f"  resum pkl : {args.resum_pkl}")

    print("\nConstructing SCETlibNPParamModel (runs the bt-grid integral at λ_central) …")
    t0 = time.time()
    indata = FitInputData(args.datacard)
    model = SCETlibNPParamModel(
        indata,
        btgrid_dir=args.btgrid,
        signal_proc=args.signal_proc,
    )
    print(f"  constructed in {time.time()-t0:.1f}s; gen grid {model.gen_shape}")

    # (A) our bt-grid resummed σ_gen at λ_central (no response matrix)
    paramodel_resum = model.sigma_gen_central  # tf, (NptVGen, NabsYVGen)

    # (B) SCETlib spectrum resummed from the correction input pkl
    print(f"\nReading correction resum pkl and projecting to the gen grid …")
    corr_resum = resum_from_correction(
        args.resum_pkl, model._gen_axes_meta, args.q_lo, args.q_hi, args.charge
    )

    h_param = tf_to_hist(paramodel_resum, model._gen_axes_meta)
    h_corr = tf_to_hist(corr_resum, model._gen_axes_meta)

    print("\n[resum] ParamModel (bt-grid) vs correction pkl, on the model gen grid:")
    summarize(h_param.values(flow=False), h_corr.values(flow=False), model._gen_axes_meta)

    if args.plot_out:
        scale = float(h_corr.values().sum() / h_param.values().sum())
        plot_ptll_ratio(
            h_param,
            h_corr,
            axis="ptVGen",
            out_path=args.plot_out,
            scale=scale,
            ref_label="SCETlib resum (correction pkl)",
            model_label=r"ParamModel resum (bt-grid)",
            rlabel="bt-grid / correction",
            title="SCETlib resummed σ_gen — bt-grid vs correction pkl",
            rrange=(0.5, 1.5),
        )

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
