#! /usr/bin/Rscript

# Plot MDS from a TASSEL distance file

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile")
parser$add_argument("-c", "--chrom", type="integer")
parser$add_argument("--start", type="integer")
parser$add_argument("--end", type="integer")
args=parser$parse_args()
# setwd('/home/jgwall/Projects/DaweRecombination/publication_graphic/')
#args=parser$parse_args(c('-i','4f_genes.gff3','-c','4','--start','157507372','--end','9999999999'))

data=read.delim(args$infile, header=F, skip=1)
names(data)[c(1,4,5)] = c("chrom","start","end")
subdata=subset(data, data$chrom==args$chrom & data$start >= args$start & data$end <= args$end)
cat(nrow(subdata),"genes are located in",args$chrom,":",args$start,"-",args$end,'\n')
