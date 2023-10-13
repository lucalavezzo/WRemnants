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

logger = logging.child_logger(__name__)

if __name__ == "__main__":

     parser = argparse.ArgumentParser()
     parser.add_argument("inputfile", type=str, nargs=1, help="Input file")
     parser.add_argument("-w", "--what", type=str, default="all",
                         choices=["all", "procs", "hists", "axes"],
                         help="What to print")
     parser.add_argument("-p", "--process", type=str, default=None,
                         help="Process to print")
     parser.add_argument("-n", "--histo", type=str, default=None,
                         help="Select specific histogram, and printout will have more information")
     args = parser.parse_args()

     h5file = h5py.File(args.inputfile[0], "r")
     results = narf.ioutils.pickle_load_h5py(h5file["results"])

     if args.what == "all":
          print(results)
     else:
          if args.what == "procs":
               print(results.keys())
          elif args.what in ["hists", "axes"]:
               print("="*30)
               print("GOING TO PRINT HISTOGRAMS")
               print("="*30)
               #print(results.keys())
               for p in results.keys():
                    if p == "meta_info":
                         continue
                    if args.process and p != args.process:
                         continue
                    print("-"*30)
                    print(f"Process {p}")
                    space = " "*5
                    for k in results[p]["output"].keys():
                         if args.histo and k != args.histo:
                              continue
                         print(f"{space}{k}")
                         if args.what == "axes":
                              print(f"{space}  Axes = {results[p]['output'][k].axes.name}")
                              if k == args.histo:
                                   for n in results[p]['output'][k].axes:
                                        print(f"{space}{space} {n}")