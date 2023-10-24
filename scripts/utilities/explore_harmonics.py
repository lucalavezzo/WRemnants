import sys
import argparse
import os
import lz4.frame
import pickle
import hist
from utilities import boostHistHelpers as hh, logging
import hdf5plugin
import h5py
import narf
from narf import ioutils
#from utilities import input_tools as it
import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np
hep.style.use("CMS")
plt.style.use(hep.style.CMS)
plt.rcParams['figure.dpi'] = 50

def integrate_pT(h, lower_bound=None, upper_bound=None, pt_axis='ptVgen'):
    if lower_bound is None: lower_bound = h.axes[-1].edges[0]
    if upper_bound is None: upper_bound = h.axes[-1].edges[-1]
    h = h[{pt_axis: slice(lower_bound*1.0j, upper_bound*1.0j)}]
    h = h[{pt_axis: sum}]
    return h

def flip_phi_axis(h, options):
    assert h.axes[1].name == 'phiStarll', 'phiStarll must be second axis if not integrating over phi, found ' + str(h.axes[1].name)
    if h.axes[1].centers[0] >= 0: return h # don't need to flip if already in [0, 2pi]
    h_new = hist.Hist(h.axes[0], hist.axis.Regular(len(h.axes[1].centers), 0, 2*np.pi, circular = True, name = "phiStarll"), storage = hist.storage.Double())
    original_values = h.values()
    new_values = np.copy(original_values)
    N = new_values.shape[1]
    new_values[:, 0:int(N/2)] = original_values[:, int(N/2):int(N)]
    new_values[:, int(N/2):int(N)] = original_values[:, 0:int(N/2)]        
    h_new.values()[...] = new_values
    return h_new

def make_Pi_plot(h_pi, hel, options):
    fig = plt.figure()
    ax = fig.subplots()
    if options.integrateTheta or options.integratePhi: _ = hep.histplot(h_pi, ax=ax)
    else: _ = hep.hist2dplot(h_pi, ax=ax)
    hep.cms.label(llabel='Preliminary',data=False, lumi=LUMI, ax=ax)
    helicityLabel = "P_{}".format(hel) if hel != -1 else "UL"
    dataType = 'Templated' if not options.gen else 'Gen'
    var_labels = ''
    if not options.integrateTheta: var_labels += r'\cos \theta_{CS}'
    if not options.integrateTheta and not options.integratePhi: var_labels += ', '
    if not options.integratePhi: var_labels += r'\phi_{CS}'
    fig.suptitle("{dataType} ${helicityLabel}$ ( ${var_labels}$ ), $p_T^Z = {lb}-{ub}$ GeV".format(dataType=dataType, helicityLabel = helicityLabel, var_labels = var_labels, lb = options.lower_bound, ub = options.upper_bound), fontsize=24)
    if options.integrateTheta: ax.set_xlabel(r"$\phi_{CS}$")
    elif options.integratePhi: ax.set_xlabel(r"$\cos\theta_{CS}$")
    else:
        ax.set_xlabel(r"$\cos\theta_{CS}$")
        ax.set_ylabel(r"$\phi_{CS}$")
    plotLabel = dataType + "_" + helicityLabel
    if options.integrateTheta: plotLabel += "_phi"
    elif options.integratePhi: plotLabel += "_cosTheta"
    fig.savefig(OUTPUT_LABEL.format(plotLabel))
    fig.savefig(OUTPUT_LABEL.format(plotLabel).replace("pdf", "png"))
    print("Saving to", OUTPUT_LABEL.format(plotLabel))

def main():
    # argument parser for input file and output directory
    parser = argparse.ArgumentParser(description='Process input and output files.')
    parser.add_argument('-f', '--file', type=str, required=True, help='Input file path')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output directory path')
    parser.add_argument('-lb', '--lower_bound', type=float, default=None, required=False, help="Integrating lower bound window of pT")
    parser.add_argument('-ub', '--upper_bound', type=float, default=None, required=False, help="Integrating upper bound window of pT")
    parser.add_argument('-it', '--integrateTheta', action='store_true', default=False, required=False, help='Integrate Theta')
    parser.add_argument('-ip', '--integratePhi', action='store_true', default=False, required=False, help='Integrate Phi')
    parser.add_argument('-g', '--gen', action='store_true', default=False, required=False, help='Plot gen level (default expected reco)')
    options = parser.parse_args()

    if options.integrateTheta and options.integratePhi:
        raise ValueError("Cannot integrate both theta and phi!")

    global LUMI 
    LUMI = '16.8'
    global INPUT_FILE 
    INPUT_FILE = options.file
    global OUTPUT_LABEL
    OUTPUT_LABEL = options.output + '/{}.pdf'

    if not os.path.isdir(options.output):
        os.mkdir(options.output)
        print("Made output directory", options.output)

    # load input file
    h5file = h5py.File(INPUT_FILE, "r")
    results = narf.ioutils.pickle_load_h5py(h5file["results"])

    # grab the plots
    if options.gen:
        h = results['ZmumuPostVFP']['output']['helicity_moments_scale'].get()
        pt_axis = 'ptVgen'
    else:
        h = results['ZmumuPostVFP']['output']['nominal'].get()
        pt_axis = 'ptVgenSig'

    # plot the P_i's
    for hel in range(-1, len(h.axes['helicitySig'].centers)-1):

        # grab the ith harmonic
        h_pi = h[{'helicitySig':hel*1.0j}]

        # integrate harmonic across other axes
        h_pi = h_pi.project("cosThetaStarll", "phiStarll", pt_axis)

        # integrate over specific pt window
        h_pi = integrate_pT(h_pi, options.lower_bound, options.upper_bound, pt_axis)

        # flip the phi axis from [-pi, pi] to [0, 2pi]
        # (to match ATLAS plots)
        h_pi = flip_phi_axis(h_pi, options)

        # if requested, integrate over a particular kinematic angle
        if options.integrateTheta: h_pi = h_pi[{'cosThetaStarll': sum}]
        elif options.integratePhi: h_pi = h_pi[{'phiStarll': sum}]

        # now plot each Pi
        make_Pi_plot(h_pi, hel, options)

if __name__ == '__main__':
    main()

