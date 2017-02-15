#! /usr/bin/Rscript

library(argparse)
parser=ArgumentParser()
parser$add_argument("-i", "--infile")
parser$add_argument("-o", "--outfile")
args=parser$parse_args()
#setwd('/home/jgwall/Projects/DaweRecombination/publication_graphic/')
#args=parser$parse_args(c('-i','4b_raw_all.ab_transposed.txt','-o','99_tmp'))

# Load data
cat("Summarizing Ab10/N10 data in",args$infile,"\n")
data=read.delim(args$infile, row.names=1, colClasses="character")

# Split by subpop
is.ab10 = grepl(names(data), pattern="W23Acr.B73_Ab10")
is.n10 = grepl(names(data), pattern="W23Acr.B73_N10")
ab10 = subset(data, select=is.ab10)
n10 = subset(data, select=is.n10)
if(ncol(ab10) + ncol(n10) != ncol(data)){
  cat("WARNING! Not all taxa split out:",ncol(data),"originally but only",ncol(ab10),"in ab10 and",ncol(n10),"in n10")
}

# Get allele counts
count_alleles=function(x){
  x=as.matrix(x)
  AA = x=="AA"
  AB = x=="AB"
  BB = x=="BB"
  
  countA = rowSums(AA) *2 + rowSums(AB)
  countB = rowSums(BB) *2 + rowSums(AB)
  return(data.frame(countA, countB))
}
ab10_stats=count_alleles(ab10)
n10_stats=count_alleles(n10)

# Make output data
chroms = sub(rownames(data), pattern="S(.+)_.+", repl="\\1")
poses = sub(rownames(data), pattern="S.+_(.+)", repl="\\1")
output=data.frame(site=rownames(data), chrom=chroms, pos=poses, ab10_A=ab10_stats$countA,  ab10_B=ab10_stats$countB, n10_A=n10_stats$countA, n10_B=n10_stats$countB)

# Write out
write.table(output, file=args$outfile, sep='\t', quote=F, row.names=F, col.names=T)
