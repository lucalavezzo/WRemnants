#!/usr/bin/env python3
"""Dump per-region smoothing polynomial coefficients for SmoothExtendedABCD.

Refit the linearized smoothing polynomial from a nominal fake histogram and
save the per-region polynomial coefficients in the SmoothExtendedABCD
power-series basis.  The output can be loaded as initial parameter values in a
rabbit fit via ``params:PATH`` in the ``--paramModel SmoothExtendedABCDIsoMT``
CLI token.

This is usually not needed: setupRabbit.py computes the same coefficients and
stores them as auxiliary data in the fit input file (see --noSmoothingParams),
where the param model picks them up by itself.  Use this script to derive the
starting values from a different input file than the one that is fit, or to
inspect them standalone.

Typical use::

    python scripts/rabbit/regen_smoothing_params.py \\
        -i mw_with_mu_eta_pt_scetlib_dyturbo.hdf5 \\
        -o /path/to/params.hdf5
"""

import os

import h5py
import numpy as np

from wremnants.postprocessing.datagroups.datagroups import Datagroups
from wremnants.postprocessing.histselections import (
    FakeSelectorSimpleABCD,
    compute_extended_abcd_initial_params,
)
from wremnants.postprocessing.regression import Regressor
from wremnants.utilities import common, parsing
from wums import ioutils, logging, output_tools

logger = logging.child_logger(__name__)


def make_parser():
    parser = parsing.base_parser()
    parser.description = __doc__
    parser.add_argument(
        "-i",
        "--inputFile",
        required=True,
        type=str,
        help="Input HDF5 histogram file (output of a histmaker).",
    )
    parser.add_argument(
        "-o",
        "--outpath",
        required=True,
        type=str,
        help="Output path for the params HDF5 file (.hdf5 appended if missing).",
    )
    parser.add_argument(
        "--inputBaseName",
        default="nominal",
        type=str,
        help="Name of the nominal histogram inside the input file.",
    )
    parser.add_argument(
        "--fakerateAxes",
        nargs="+",
        default=["eta", "pt", "charge"],
        help="Axes for the fakerate binning.",
    )
    parser.add_argument(
        "--fakeEstimation",
        type=str,
        default="extended1D",
        choices=["simple", "extrapolate", "extended1D", "extended2D"],
        help="Fake estimation mode (must match what will be used in setupRabbit).",
    )
    parser.add_argument(
        "--fakeSmoothingMode",
        type=str,
        default="full",
        choices=FakeSelectorSimpleABCD.smoothing_modes,
        help="Smoothing mode for fake estimate.",
    )
    parser.add_argument(
        "--fakeSmoothingOrder",
        type=int,
        default=3,
        help="Polynomial order for the spectrum smoothing.",
    )
    parser.add_argument(
        "--fakeSmoothingPolynomial",
        type=str,
        default="chebyshev",
        choices=Regressor.polynomials,
        help="Polynomial type for the spectrum smoothing.",
    )
    parser.add_argument(
        "--lumiScale",
        type=float,
        default=1.0,
        help="Rescale equivalent luminosity by this value (must match the value used in setupRabbit).",
    )
    parser.add_argument(
        "--flowToExplicitBins",
        type=str,
        nargs="*",
        default=["mt", "relIso", "iso"],
        help="""
        Axes for which under-/overflow bins are converted into explicit bins with an infinite outer edge,
        must match the value used in setupRabbit. Needed for histograms from older histmakers, where the
        open ended ABCD regions were stored in the overflow bin, instead of an explicit last bin.
        Axes that are not in the histogram, or that have no flow bins, are ignored. Pass no argument to disable.
        """,
    )
    parser.add_argument(
        "--excludeProcGroups",
        type=str,
        nargs="*",
        default=["QCD"],
        help="Process groups to exclude when building Datagroups.",
    )
    parser.add_argument(
        "--filterProcGroups",
        type=str,
        nargs="*",
        default=None,
        help="If set, keep only these process groups when building Datagroups.",
    )
    return parser


def dump_smoothing_params(
    outpath, fakeselector, datagroups, inputBaseName, meta_data_dict=None, postfix=""
):
    """
    Dump the per-region Chebyshev polynomial coefficients of the nominal fake
    histogram smoothing fit to a standalone file, in the layout expected by
    SmoothExtendedABCD as ``initial_params``.

    The coefficients are computed by
    ``histselections.compute_extended_abcd_initial_params``, which is also used by
    setupRabbit.py to store them directly in the fit input file. See there for
    the details of the layout and the region ordering.
    """
    datasets = compute_extended_abcd_initial_params(
        fakeselector, datagroups, inputBaseName
    )

    if outpath and not os.path.isdir(outpath):
        os.makedirs(outpath)

    outfile = f"{outpath}/params"
    if postfix:
        outfile += f"_{postfix}"

    outfile += ".hdf5"

    with h5py.File(outfile, mode="w") as f:
        f.create_dataset("params", data=datasets["params"])
        f.create_dataset("order", data=datasets["order"])
        f.create_dataset(
            "smoothing_axis_name",
            data=np.array(
                datasets["smoothing_axis_name"][0], dtype=h5py.string_dtype()
            ),
        )
        f.create_dataset("n_outer", data=datasets["n_outer"])
        f.create_dataset("outer_shape", data=datasets["outer_shape"])
        if meta_data_dict is not None:
            ioutils.pickle_dump_h5py("meta", meta_data_dict, f)

    logger.info(f"Saved smoothing initial params to {outfile}")


def main():
    parser = make_parser()
    args = parser.parse_args()
    logging.setup_logger(__file__, args.verbose, args.noColorLogger)

    dg = Datagroups(
        args.inputFile,
        excludeGroups=args.excludeProcGroups if args.excludeProcGroups else None,
        filterGroups=args.filterProcGroups if args.filterProcGroups else None,
    )

    dg.lumiScale = args.lumiScale
    dg.fakerate_axes = args.fakerateAxes
    # older histmakers stored the open ended ABCD regions in the overflow bins,
    # turn them into explicit bins with an infinite outer edge, as done in setupRabbit
    dg.flowToExplicitBinsAxes = args.flowToExplicitBins or []
    dg.set_histselectors(
        dg.getNames(),
        args.inputBaseName,
        mode=args.fakeEstimation,
        smoothing_mode=args.fakeSmoothingMode,
        smoothingOrderSpectrum=args.fakeSmoothingOrder,
        smoothingPolynomialSpectrum=args.fakeSmoothingPolynomial,
        mcCorr=None,
        integrate_x=True,
        forceGlobalScaleFakes=False,
        abcdExplicitAxisEdges={},
        fakeTransferAxis="",
        fakeTransferCorrFileName=None,
        histAxesRemovedBeforeFakes=[],
    )

    fakeselector = dg.groups[dg.fakeName].histselector
    logger.info(f"fakeselector type: {type(fakeselector).__name__}")

    meta_data_dict = {
        "meta_info": output_tools.make_meta_info_dict(
            args=args,
            wd=common.base_dir,
        ),
        "meta_info_input": dg.getMetaInfo(),
    }

    dump_smoothing_params(
        args.outpath,
        fakeselector,
        dg,
        args.inputBaseName,
        meta_data_dict=meta_data_dict,
    )


if __name__ == "__main__":
    main()
