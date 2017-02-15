#Calculate and show # of rare alleles by taxon


cat("Testing for bad taxa via too many rare alleles\n")
args=commandArgs(TRUE)	
infile=args[1]
outfile=args[2]
namesfile=args[3]	# File with full taxa names so can be output correctly
parents=args[4:length(args)]	# Any parental genotypes at the end
#setwd("/home/jgwall/Documents/Papers/DAWE_ab10recombination/Figures/File SX - Distortion Scripts/")
#infile='1_Analysis/1c.rare_alleles.numeric.txt';  outfile='1_Analysis/1c.bad_taxa.rare_alleles' ; namesfile='1_Analysis/1c.full_taxa_names.txt'; parents=c('B73:250052761","W23-Acr:250052635')


nrows=1000	#Overestimate of how many rows need

#Read in data
taxanames = scan(namesfile, what=character())
header=scan(infile, what=character(), nlines=1, skip=1)
data=read.delim(infile, skip=1, row.names=1, header=T, colClasses = c('character', rep('numeric',length(header)-1)))
data=1-as.matrix(data)	# TASSEL claims its ReferenceProbability is outputting the probability of a minor allele, but looking at the file it's actually the major allele
data=ceiling(data)	#Round hets up to 1 to just count "instances of a rare allele"

#Do significance test
rarecount = rowSums(data, na.rm=T)	#Count of # rare allele calls in each line
totalsnps = rowSums(!is.na(data))	#Count of total called SNPs in each line
rarefreq=rarecount/totalsnps
sig.test = pnorm(q=rarefreq, mean=mean(rarefreq, na.rm=T), sd=sd(rarefreq, na.rm=T), lower.tail=F)

#TODO: Switch this so take anything more than 2 SD away from the mean?


fdr = p.adjust(sig.test, method="fdr")	#False discovery rate 
col = ifelse(fdr < 0.1, yes="blue", no="black")
col = ifelse(fdr < 0.01, yes="red", no=col)


png(file=paste(outfile, ".png", sep=""), width=2000, height=1000)
	par(mfrow=c(1,2), cex=1.2)
	plot(x=1:length(rarefreq), y=rarefreq, pch=20, col=col, main="Rare alleles by line", xlab="line", ylab="Rare allele counts")
	problems = fdr < 0.1
	if(sum(problems, na.rm=T)>0){
		text(x=(1:length(rarefreq))[problems], y=rarefreq[problems], labels=rownames(data)[problems], pos=4, col=col[col!="black"], cex=0.7)	
		output=rownames(data)[problems]
	}else{
		text(x=1, y=min(rarefreq, na.rm=T), labels="No problem taxa identified", pos=4, col="darkgreen", cex=2)	
		output=""
	}
	plot(x=1:length(sort(rarefreq)), y=sort(rarefreq), pch=20, col=col[order(rarefreq)], main="Rare alleles by line (ordered)", xlab="line", ylab="Rare allele freq")
dev.off()

# fix taxa names
cat("Putting taxa names back to full\n")
tofind=paste("^",output,":",sep="")
for(i in 1:length(output)){
  matched = grepl(taxanames, pattern=tofind[i])
  if(sum(matched)==1){
	output[i] = taxanames[matched]
  }else if(sum(matched)>1){
	cat("WARNING! More than one taxon matched for",output[i],":",taxanames[matched],"\n")
  }
  # Don't filter out parental lines
  if(output[i] %in% parents){
	output[i]=""
  }
}
  
write.table(file=paste(outfile,".txt",sep=""), x=output, row.names=F, col.names=F, sep="\t", quote=F)