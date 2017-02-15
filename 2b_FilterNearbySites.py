__author__ = 'jgwall'

import argparse
import numpy as np
import re

debug = False

sitekey = re.compile("S(.+)_(.+)")  # Regular expression to extract chrom and positions

def main():
    args = parse_args()
    print("Filtering SNPs in", args.infile,"based on windows of",args.winsize,"bp")
    IN = open(args.infile, "r")
    OUT = open(args.outfile, "w")

    # Header
    OUT.write(IN.readline())

    # first line
    line = IN.readline()
    chrom, pos = get_position(line)
    oldsnp = make_snp(chrom, pos, line)

    # Go through sites
    n, kept, removed = 0, 0, 0
    for line in IN:
        n+=1
        if debug and n > 1000: break
        chrom, pos = get_position(line)
        mysnp = make_snp(chrom, pos, line)

        # If on a new chromosome or if position is far enough away, output old SNP and continue
        if (mysnp['chrom'] != oldsnp['chrom']) or (mysnp['pos'] - oldsnp['pos'] > args.winsize):
            # print("Keeping ", oldsnp['chrom'], ":", oldsnp['pos'])
            OUT.write(oldsnp['data'])
            oldsnp = mysnp
            kept+=1
            continue

        # If SNPs are within window size, keep the one with less missing data. (If a tie, keep the first one)
        if sum_missing(mysnp['data']) < sum_missing(oldsnp['data']):
            oldsnp = mysnp
            # print("Skipping",oldsnp['chrom'],":",oldsnp['pos'])
        else:
            # print("Skipping", mysnp['chrom'], ":", mysnp['pos'])
            pass    # If previous SNP has less missing data, no need to update. Kept statement for clarity

        removed +=1

    OUT.write(oldsnp['data'])   # final cleanup
    print("Processed",n,"sites")
    print("\tRemoved",removed,"sites and kept",kept,"sites")

    IN.close()
    OUT.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("-w", "--winsize", type=int, default=100, help="Window size; all SNPs within this window size will be reduced to 1")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()

def get_position(line):
    site = line.split('\t')[0]
    sitedata = sitekey.search(site)
    chrom, pos = sitedata.group(1), sitedata.group(2)
    # print(site, "becomes",chrom,":",pos)
    return chrom, int(pos)

def make_snp(chrom, pos, data):
    snp=dict()
    snp['chrom']=chrom
    snp['pos']=pos
    snp['data']=data
    return snp

def sum_missing(line):
    data=np.array(line.strip().split('\t'))
    missing_count = np.sum(data == '-')
    # print(data[0],"has", missing_count,"missing sites")
    return missing_count

if __name__ == '__main__': main()