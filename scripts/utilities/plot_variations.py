import uproot
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse

def unroll_th2d_to_th1d(hist2D, hist1D):
    if hist2D is None or hist1D is None:
        print("Input histograms are null!")
        return

    x_bins = hist2D.shape[0]
    y_bins = hist2D.shape[1]

    for i in range(0, x_bins):
        for j in range(0, y_bins):
            bin_1D = (i - 1) * y_bins + j
            content = hist2D[i,j]
            if type(content) == float:
                hist1D[bin_1D] = content
            else:
                hist1D[bin_1D] = content.value

    return hist1D

def plot_histograms(sys_name, nominal, up, down):
    
    fig, ax = plt.subplots(figsize=(18, 6))

    ax.plot(nominal, label="Nominal", color="black")
    ax.plot(up, label=f"{sys_name} Up", color="red")
    ax.plot(down, label=f"{sys_name} Down", color="blue")

    ax.legend()
    ax.set_title("Histograms")
    ax.set_xlabel("Bin")
    ax.set_ylabel("Content")

    fig.savefig(f"{out_dir}/{sys_name}_histograms.png")
    plt.close()

def plot_ratios(sys_name, up, down, nominal):
    
    fig, ax = plt.subplots(figsize=(8, 6))

    ratio_up = np.where(nominal > 1e-4, up / nominal, 1)
    ratio_down = np.where(nominal > 1e-4, down / nominal, 1)

    ax.plot(ratio_up, label=f"Ratio {sys_name} Up", color="red")
    ax.plot(ratio_down, label=f"Ratio {sys_name} Down", color="blue")

    ax.axhline(1, color="black", linestyle="--", linewidth=1)

    ax.legend()
    ax.set_title("Ratios")
    ax.set_xlabel("Bin")
    ax.set_ylabel("Ratio")

    fig.savefig(f"{out_dir}/{sys_name}_ratios.png")
    plt.close()

def plot_histograms_and_ratios(sys_name, nominal, up, down):
    
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    axs[0].plot(nominal, label="Nominal", color="black")
    axs[0].plot(up, label=f"{sys_name} Up", color="red")
    axs[0].plot(down, label=f"{sys_name} Down", color="blue")
    axs[0].legend()
    axs[0].set_title("Histograms")
    axs[0].set_xlabel("Bin")
    axs[0].set_ylabel("Content")

    ratio_up = np.where(nominal > 1e-4, up / nominal, 1)
    ratio_down = np.where(nominal > 1e-4, down / nominal, 1)
    axs[1].plot(ratio_up, label="Ratiof {sys_name} Up", color="red")
    axs[1].plot(ratio_down, label="Ratiof {sys_name} Down", color="blue")
    axs[1].axhline(1, color="black", linestyle="--", linewidth=1)
    axs[1].legend()
    axs[1].set_title("Ratios")
    axs[1].set_xlabel("Bin")
    axs[1].set_ylabel("Ratio")

    fig.tight_layout()
    fig.savefig(f"{out_dir}/{sys_name}_histograms_and_ratios.png")
    plt.close()

def proces_systematic(file, process, nominal, sys_name, systematic_up, systematic_down):

    nominal_hist = file[process][nominal].to_hist()
    up_hist = file[process][systematic_up].to_hist()
    down_hist = file[process][systematic_down].to_hist()

    nominal1D = np.zeros(nominal_hist.size)
    up1D = np.zeros(up_hist.size)
    down1D = np.zeros(down_hist.size)

    nominal1D = unroll_th2d_to_th1d(nominal_hist, nominal1D)
    up1D = unroll_th2d_to_th1d(up_hist, up1D)
    down1D = unroll_th2d_to_th1d(down_hist, down1D)

    plot_histograms(sys_name, nominal1D, up1D, down1D)
    plot_ratios(sys_name, up1D, down1D, nominal1D)
    plot_histograms_and_ratios(sys_name, nominal1D, up1D, down1D)

def parse_hists_for_systematic(systematics, base_name, sys_name):
    all_hists = file[process].keys()
    for hist in all_hists:
        if sys_name not in hist: continue
        if base_name not in hist: continue
        if 'up' not in hist.lower(): continue
        if hist.replace('Up', 'Down') not in all_hists: continue
        up = hist
        down = hist.replace('Up', 'Down')
        up = up.replace(";1", '')
        down = down.replace(";1", '')

        sys_label = up.replace("up", '')

        systematics[sys_label] = (up, down)
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-p', '--process', required=False, default='Zmumu', type=str, help='the process name')
    parser.add_argument('-f', '--file', type=str, help='the input file name')
    parser.add_argument('-o', '--out_dir', type=str, help='the output directory')

    args = parser.parse_args()

    process = args.process
    infile = args.file
    out_dir = args.out_dir

    file = uproot.open(infile)

    # either add the up and down variation for each systematic by hand here...
    systematics = {
        'pdfAlphaS' : ('nominal_Zmumu_pdfAlphaSDown_inclusive', 'nominal_Zmumu_pdfAlphaSUp_inclusive'),
    }  

    # ... or automatically add all systmatics matching a certain pattern here
    parse_hists_for_systematic(systematics, 'nominal_Zmumu', 'resumTNP')
    parse_hists_for_systematic(systematics, 'nominal_Zmumu', 'Resolution_correction_smearing_variation')
    parse_hists_for_systematic(systematics, 'nominal_Zmumu', 'Z_nonClosure_parametrized_A')
    parse_hists_for_systematic(systematics, 'nominal_Zmumu', 'CMS_prefire_stat')
    parse_hists_for_systematic(systematics, 'nominal_Zmumu', 'massShift')

    print("Will output the following systematics:")
    print(list(systematics.keys()))

    # process each systematic
    for sys_name, (up, down) in systematics.items():
        proces_systematic(file, process, 'nominal_Zmumu_inclusive', sys_name, up, down)