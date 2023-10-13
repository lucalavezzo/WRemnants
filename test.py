#!/usr/bin/env python3

import argparse
import os
import sys

import lz4.frame
import pickle
import hist
from utilities import boostHistHelpers as hh, logging
#import wremnants
import hdf5plugin
import h5py
import narf
from narf import ioutils
import ROOT
from utilities import input_tools as it

logger = logging.child_logger(__name__)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", type=str, nargs=1, help="Input file")
    args = parser.parse_args()

    h5file = h5py.File(args.inputfile[0], "r")
    results = narf.ioutils.pickle_load_h5py(h5file["results"])

    h = it.read_and_scale(args.inputfile[0], "ZmumuPostVFP", "nominal_ptll")