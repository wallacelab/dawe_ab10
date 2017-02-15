__author__ = 'jgwall'

import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#from scipy.stats import ttest_ind as ttest
from scipy.stats import chisquare #fisher_exact, 
from scipy.stats.distributions import chi2

debug = False


def main():
    args = parse_args()

    # Load data
    print("Loading input data")
    data=pd.read_csv(args.infile, sep='\t', nrows = 1000 if debug else None)
    print("\t Loaded data on",len(data),"sites")
    data=data.loc[data['chrom']!=0,:]
    print("\t\t",len(data), "sites after removing chromosome 0")
    chromlengths = load_chromlengths(args.chromlengths)
    print("\t", len(chromlengths), "chromosomes")

    # Calculate allele frequencies
    data['ab10_freq'] = data['ab10_B'] / (data['ab10_A'] + data['ab10_B'])
    data['n10_freq'] = data['n10_B'] / (data['n10_A'] + data['n10_B'])
    # Get plot position of raw data
    data['plot_x'] = [ chromlengths[str(c)]['cum'] + p  for c,p in zip(data['chrom'], data['pos']) ]
    
    # Remove any with too low allele frequencies
    print("Removing any sites with minor alleles with freqs below",args.min_freq)
    data = data.loc[(data['ab10_freq'] >= args.min_freq) & (data['n10_freq'] >= args.min_freq),:]	# Lower tail
    data = data.loc[(data['ab10_freq'] <= 1- args.min_freq) & (data['n10_freq'] <= 1- args.min_freq),:] # upper tail
    print("\t",len(data), "sites remain after removing these")

    # plot raw data if asked
    if(args.plot_raw_freqs): graph_raw_freqs(data, chromlengths, args.outprefix)

    ## calculate and plot windows
    data['plot_x'] = [chromlengths[str(c)]['cum'] + p for c, p in zip(data['chrom'], data['pos'])] # Get plotting coorindates
    results = calc_differences(data, chromlengths, args.window, args.step)
    graph_differences(data, results, chromlengths, args.window, args.step, args.outprefix, args.quantile)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Summarized file of allele counts")
    parser.add_argument("--chromlengths", help="Chromosome lenths input file")
    parser.add_argument("--outprefix", help="Prefix for all the output files")
    parser.add_argument("-m", "--min_freq", type=float, default=0, help="Minimum allele frequency to be included")
    parser.add_argument("-r","--plot-raw-freqs", default=False, action="store_true", help="Whether to plot the raw values across the chromosomes")
    parser.add_argument("-q", "--quantile", type=float, default=80, help="What quantile of spread to plot for frequencies (so, a value of 50 gets the 25th and 75th percentile)")

    # Window-based arguments
    parser.add_argument("--window", default=100, type=int, help="Number of sites to include in the sliding window")
    parser.add_argument("--step", default=20, type=int, help="Number of sites to step each sliding window")

    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def load_chromlengths(infile):
    lengths=dict()
    for line in open(infile, "r"):
        chrom, length, cum = line.strip().split('\t')[:3]   # ':2' in case there are additional columns for some reas
        lengths[chrom] = {"length":int(length), "cum":int(cum)}
    return lengths

def graph_raw_freqs(data, chromlengths, outprefix):
    outfile=outprefix + ".raw_freqs.png"
    print("Graphing raw fequency data to", outfile)
    
     # Plot results
    chroms = sorted(set(data['chrom']))
    nrow = math.ceil(len(chroms)/2)
    fig = plt.figure(figsize=(10, 3 * nrow))
    grid = gridspec.GridSpec(nrows=nrow, ncols=2, hspace=0.5, wspace=0.5)
    chrom_max = max([chromlengths[str(c)]['length'] for c in chroms])

    # Plot 1 chrom at a time
    row, col = 0, 0
    for mychrom in chroms:
      subdata = data.loc[data['chrom'] == mychrom,:]
      # Make axes
      ax_freqs = fig.add_subplot(grid[row,col], title="Chromosome " + str(mychrom), xlabel="Position (megabases)", ylabel="Frequency of W23")
         
      # Plot
      adata = pd.DataFrame({"pos":subdata['pos'], 'freq':subdata['ab10_freq'], 'color':"red"})
      ndata = pd.DataFrame({"pos":subdata['pos'], 'freq':subdata['n10_freq'], 'color':"blue"})
      plotdata = pd.concat([adata, ndata])
      plotdata = plotdata.sort('pos')
      #print(plotdata.head())
      
      scatter = ax_freqs.scatter(x=plotdata['pos']/1000000, y=plotdata['freq'], c=plotdata['color'], alpha=0.05, linewidths=0)
      scatter.set_rasterized(True)
      
      # Make pretty
      xlim = [0, chrom_max]
      #ax_freqs.set_xlim(xlim)
      ax_freqs.legend(fontsize="xx-small", bbox_to_anchor=(1.045, 1.0))
         
      row+=1
      if row >= nrow:
        row = 0
        col +=1
  
    fig.savefig(outfile, dpi=150)
    

def calc_differences(data, chromlengths, window, step):
    # print(data.head())
    chroms = np.array(data['chrom'])
    poses = np.array(data['pos'])
    plot_pos = np.array(data['plot_x'])

    # Go through windows
    start = 0
    results=list()
    while start < len(data):
            stop=start + window
            if stop > len(data): stop = len(data)
            
            # If break into another chromosome, skip forward to next one
            if chroms[stop-1] != chroms[start]:
                # print("Retrocessing stop because chrom",chroms[stop-1],"at",stop,"does not equal chrom",chroms[start],"at",start)
                while chroms[stop-1] != chroms[start]: stop-=1
                start=stop
                continue

            # print(start, stop, poses[start])
            target = data.iloc[start:stop,:]
            n_freqs = target['n10_B'] / (target['n10_B'] + target['n10_A']) # Frequency of allele B in N10
            a_freqs = target['ab10_B'] / (target['ab10_B'] + target['ab10_A'])# Frequency of allele B in Ab10
            observations = [np.sum(target['ab10_A']), np.sum(target['ab10_B'])]
            expected = [np.sum(target['n10_A']), np.sum(target['n10_B'])]
            chisq, pval = chisquare(f_obs=observations, f_exp=expected)
            # print(observations, expected, pval)

            mychrom = chroms[start]
            label= str(mychrom) + ":" + str(poses[start]) + "-" + str(poses[stop-1])
            plot_position = (plot_pos[start] + plot_pos[stop-1]) / 2
            tmp = {'label':label,'chrom':mychrom,'plot_x':plot_position,"pval":pval, "a_freqs":a_freqs,"n_freqs":n_freqs}
            tmp['ab10_A'], tmp['ab10_B'] = observations
            tmp['n10_A'], tmp['n10_B'] = expected
            #results.append([label, mychrom, plot_position, pval, a_freqs, n_freqs])
            results.append(tmp)

            start+=step
            #if debug and start > 100: break
    print("Extracted data on",len(results),"windows")
    return(results)


def graph_differences(data, results, chromlengths, window, step, outprefix, quantile):

    # Unpack data
    poses = np.array(data['plot_x'])
    labels = [r['label'] for r in results]
    chroms = [r['chrom'] for r in results]
    plot_positions = np.array([r['plot_x'] for r in results])
    pvals = [r['pval'] for r in results]
    a_freqs = [r['a_freqs'] for r in results]
    n_freqs = [r['n_freqs'] for r in results]
    data['chrom_string'] = [str(c) for c in data['chrom']]
    
    # Get quantiles
    a_quants = get_quantiles(a_freqs, quantile)
    n_quants = get_quantiles(n_freqs, quantile)
    
    # Output text
    output = pd.DataFrame({'region':labels})
    output['pval'] =pvals
    for column in "ab10_A", "ab10_B", "n10_A", "n10_B":
      output[column] = [r[column] for r in results]
    output.to_csv(outprefix + ".txt", sep='\t')
    
    # Correct pvalues that went to zero to be 10x less than other minimum
    pvals = np.array(pvals)
    old_pvals = pvals.copy()
    min_p = np.nanmin(pvals[pvals > 0])
    pvals[pvals == 0] = min_p * 0.1
    print("\tAltering pvalues of 0 resulted in",np.nansum(pvals != old_pvals),"changes set to 10x less than nonzero values")

    # Plot results
    print("Outputting graphic")
    fig = plt.figure(figsize=(12,3))
    #grid = gridspec.GridSpec(nrows=2, ncols=1, hspace=0.5, wspace=0.5)

    # Allele frequency quantiles
    ax_quants =   fig.add_axes([0.05, 0.4, 0.9, 0.5], ylabel="Frequency of W23")
    ax_meanpval = fig.add_axes([0.05, 0.1, 0.9, 0.3], ylabel="-log10 p-value")
    
    #colors
    ab10color, n10color = "red", "dodgerblue"
    
    # PLot 1 chrom at a time
    for mychrom in sorted(chromlengths.keys()):
      start = chromlengths[mychrom]['cum']
      stop = chromlengths[mychrom]['length'] + start
      toplot = (plot_positions >= start) & (plot_positions <= stop)
      #print("\tPlotting",np.sum(toplot),"points on chrom",mychrom,"between",start,"and",stop)
      
      # Ghostly original values behind
      subdata = data.loc[data['chrom_string'] == mychrom,:]      
      # Plot
      adata = pd.DataFrame({"plot_x":subdata['plot_x'], 'freq':subdata['ab10_freq'], 'color':ab10color})
      ndata = pd.DataFrame({"plot_x":subdata['plot_x'], 'freq':subdata['n10_freq'], 'color':n10color})
      plotdata = pd.concat([adata, ndata])
      plotdata = plotdata.sort('plot_x')
      #print(plotdata.head())
      
      scatter = ax_quants.scatter(x=plotdata['plot_x'], y=plotdata['freq'], c=plotdata['color'], alpha=0.025, linewidths=0, label='')
      scatter.set_rasterized(True)

      # Line of mean values
      ax_quants.plot(plot_positions[toplot], a_quants['mean'][toplot], color='black', linewidth=3, alpha=0.5)   # outline
      ax_quants.plot(plot_positions[toplot], n_quants['mean'][toplot], color='black', linewidth=3, alpha=0.5)   # outline
      ax_quants.plot(plot_positions[toplot], a_quants['mean'][toplot], color=ab10color, label="Ab10" if mychrom=='1' else '')
      ax_quants.plot(plot_positions[toplot], n_quants['mean'][toplot], color=n10color, label="N10" if mychrom=='1' else '')


      # # Pvalues from means
      ax_meanpval.plot(plot_positions[toplot], -np.log10(pvals)[toplot], color="dimgray")
    
    # Set desired padding
    xmin, xmax = np.min(poses), np.max(poses)
    xpad = (xmax-xmin) / 100
    ax_quants.set_xlim(left=xmin-xpad, right=xmax+xpad)
    ax_meanpval.set_xlim(left=xmin-xpad, right=xmax+xpad)
    ax_quants.set_ylim([0,1])
    
    # Add legend and chromosome divisions
    ax_quants.legend(fontsize="xx-small", bbox_to_anchor=(1.055, 1.0))
    add_chromdivs(ax_quants, chromlengths)
    add_chromdivs(ax_meanpval, chromlengths)

    # Alter axis lines
    ax_quants.spines['right'].set_visible(False)
    ax_quants.spines['top'].set_visible(False)
    ax_quants.spines['bottom'].set_visible(False)
    ax_meanpval.spines['left'].set_visible(False)
    ax_meanpval.spines['top'].set_visible(False) 
    ax_meanpval.spines['bottom'].set_visible(False) 
    
    # Alter tick marks
    ax_quants.tick_params(axis='y', which='both', right='off')
    ax_quants.tick_params(axis='x', which='both', top='off', bottom='off')
    ax_meanpval.tick_params(axis='x', which='both', top='off', bottom='off')

    # Plot chromosome names
    cnames, cpos, chromdiv = list(), list(), set()
    for mychrom in sorted(chromlengths.keys()):
      if mychrom in {'UNKNOWN', 'Pt', 'Mt', '0'}:continue # skip inconsequential ones
      cnames.append("chr" + mychrom)
      cpos.append(chromlengths[mychrom]['cum'] + chromlengths[mychrom]['length']/2)
    ax_meanpval.set_xticks(cpos)
    ax_meanpval.set_xticklabels(cnames, fontsize="small")
    ax_quants.get_xaxis().set_visible(False)
    
    # Move mean pval axis labels to right hand side
    ax_meanpval.yaxis.set_label_position("right")
    ax_meanpval.yaxis.tick_right()
    
    # Alter text sizes
    title="Transmission distortion: window size " + str(window) + " & step size " + str(step) + ", Chi-square test of equivalence"
    ax_quants.set_title(title, fontsize="medium")
    for l in ax_quants.get_yticklabels() + ax_meanpval.get_yticklabels():
       l.set_size("xx-small")
    
    ax_meanpval.patch.set_alpha(0)

    fig.savefig(outprefix + ".png", dpi=600)
    fig.savefig(outprefix + ".svg", dpi=600)

def get_quantiles(freqs, quantile):
    quantiles=dict()
    spread=quantile/2
    quantiles['q1'] = np.array([np.percentile(f, 50+spread) for f in freqs])
    quantiles['q2'] = np.array([np.percentile(f, 50-spread) for f in freqs])
    quantiles['mean'] = np.array([np.mean(f) for f in freqs])
    # print(quantiles)
    return quantiles


def add_chromdivs(ax, chromlengths):
    for chr in chromlengths.keys():
        ax.axvline(chromlengths[chr]['cum'], color="lightgray", linestyle="dashed")

def plotchroms(ax, chromlengths, xvals, yvals, chroms, color, label):
    xvals, yvals = np.array(xvals), np.array(yvals)
    for mychrom in np.unique(chroms):
        toplot = chroms == mychrom
        if mychrom == chroms[0]:    # Hack to have only 1 legend value
            ax.plot(xvals[toplot], yvals[toplot], color=color, label=label)
        else:
            ax.plot(xvals[toplot], yvals[toplot], color=color)
    add_chromdivs(ax, chromlengths)

if __name__ == '__main__': main()
