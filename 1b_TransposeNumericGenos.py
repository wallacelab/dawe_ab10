__author__ = 'jgwall'

import argparse
import numpy as np

debug = False


def main():
    args = parse_args()

    print("Transposing numeric genotypes from",args.infile)
    if debug:
        data = np.genfromtxt(args.infile, skip_header=1, dtype=str, max_rows=5)
    else:
        data = np.genfromtxt(args.infile, skip_header=1, dtype=str)
    print("Data has shape",data.shape,"; transposing")
    print("Writing transposed file to",args.outfile)
    data = data.T

    OUT = open(args.outfile, "w")
    n=0
    for i in range(len(data)):
        OUT.write("\t".join(data[i,:]) + "\n")
        n+=1
    OUT.close()
    print("\tWrote",n,"lines of data")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile")
    parser.add_argument("-o", "--outfile")
    parser.add_argument("--debug", default=False, action="store_true")
    args = parser.parse_args()

    # Handle debug flag
    global debug
    debug = args.debug

    return parser.parse_args()


if __name__ == '__main__': main()