__author__ = 'jgwall'

import argparse
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

debug = False


def main():
    args = parse_args()

    # Load data
    data=None
    if debug:
        # data = np.genfromtxt(args.infile, dtype=str, max_rows=100)
        data = pd.read_csv(args.infile, sep='\t', index_col=0, nrows=10000)
    else:
        data = pd.read_csv(args.infile, sep='\t', index_col=0)
    print("Loaded",len(data),"sites of data")

    # Set up graphic
    nrow=len(args.subpops) + 1
    ncol=2
    fig = plt.figure(figsize=(ncol * 4, nrow * 4))
    grid = gridspec.GridSpec(nrows=nrow, ncols=ncol, hspace=0.5, wspace=0.5)
    myrow=0
    for subpop in [None] + args.subpops:
        add_plots(data, subpop, fig, grid, args.min, args.max, args.flip, myrow, args.outfile)
        myrow+=1
    fig.savefig(args.graphfile, dpi=100)





def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="TASSEL numerical genotype transform / reference probability output (0-1 scaled)")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-g", "--graphfile")
    parser.add_argument("--min", default=0, type=float, help="Minimum allowable frequency")
    parser.add_argument("--max", default=1, type=float, help="Maximum allowable frequency")
    parser.add_argument("-f", "--flip", default=False, action="store_true", help="Whether to flip allele frequencies (= subtract them all from 1)")
    parser.add_argument("--subpops", nargs="*", help="Subpopulations to split the results into for graphing; each will be searched among sample names using regex")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def add_plots(data, pop, fig, grid, min, max, flip, row, outfile):
    groupname = pop if pop is not None else "Everything"
    print("Adding plot for group",groupname)

    # Convert to numpy arrays for faster processing
    counts = np.array(data)
    samples = np.array(data.columns)
    sites = np.array(data.index)
    if np.any(counts > 1):
        print("\tError. Some genotypes are >1")
        sys.exit(1)

    # Slice out subgroup if provided
    if pop is not None:
        tokeep = np.array([pop in s for s in samples])
        counts = counts[:,tokeep]
    print("\tFound",counts.shape,"samples in this group")

    # Determine which sites pass muster
    freqs = np.nanmean(counts, axis=1)
    if flip: freqs = 1 - freqs
    goodsites = (freqs >= min) & (freqs <= max)
    print("\tFound", np.sum(goodsites),"out of",len(goodsites),"sites between frequencies",min,"and",max)

    if (pop is None) and (outfile is not None):
        bad_sites = sites[~goodsites]
        print("\tWriting",len(bad_sites),"bad sites to",outfile)
        OUT = open(outfile, "wt")
        OUT.writelines([s + "\n" for s in bad_sites])
        OUT.close()

    # Set up axes
    ax_before = fig.add_subplot(grid[row, 0], title="Distribution pre-filter ("+groupname+")", xlabel="Allele Frequency", ylabel="Count")
    ax_after = fig.add_subplot(grid[row, 1], title="Distribution post-filter ("+groupname+")", xlabel="Allele Frequency", ylabel="Count")

    # Graph
    bins = np.linspace(0, np.nanmax(freqs), num=50)
    ax_before.hist(freqs[goodsites], bins=bins, color="blue")
    ax_before.hist(freqs[~goodsites], bins=bins, color="red")
    ax_after.hist(freqs[goodsites], bins=50, color="blue")

    # Pretty formatting
    for ax in [ax_before, ax_after]:
        for t in ax.get_xticklabels():
            t.set_rotation("vertical")



if __name__ == '__main__': main()