#! /usr/bin/env bash

# Analyze Kelly Dawe's GBS data on Ab10 recombination to see if Ab10 dramatically increases recombination rates in lines

# TODO List tASSEL version
TASSEL5="perl $HOME/Software/TASSEL/tassel-5-standalone/run_pipeline.pl -Xms1g -Xmx4g" # TASSEL version 5.5.26

# Source files
genos=1_genos.hmp.txt.gz

# Analysis variables
parentB73="B73:250052761"
parentW23="W23-Acr:250052635"	#Names of the two parents

# Working directory
workdir="1_Analysis"
if [ ! -e $workdir ]; then mkdir $workdir; fi


#######
# Fitlering raw data to reliable sites and taxa
#######

#   Filter sites based on excess heterozygostiy (= paralogs)
site_het_cutoff=0.05 # Cutoff for calling that a site has too many heterozygous calls (which generally result from paralogs)
$TASSEL5 -h $genos -genotypeSummary site -export $workdir/1a_genos.sitesummary.txt 
Rscript 1a_FindProbableParalogousSites.r $workdir/1a_genos.sitesummary.txt  $workdir/1a_genos.bad_sites_too_heterozygous.txt $site_het_cutoff
$TASSEL5 -h $genos -filterAlign -excludeSiteNamesInFile $workdir/1a_genos.bad_sites_too_heterozygous.txt -export $workdir/1a_genos.filter_paralogs.hmp.txt.gz 

#   # Filter based on expected allele frequencies
$TASSEL5 -h $workdir/1a_genos.filter_paralogs.hmp.txt.gz -NumericalGenotypePlugin -endPlugin -export $workdir/1a_genos.filter_paralogs.numeric.txt -exportType ReferenceProbablity
python3 1b_TransposeNumericGenos.py -i $workdir/1a_genos.filter_paralogs.numeric.txt -o $workdir/1a_genos.filter_paralogs.numeric.transposed.txt
python3 1b_FindSitesWithBadAlleleFrequencies.py -i $workdir/1a_genos.filter_paralogs.numeric.transposed.txt -o $workdir/1b_genos.bad_sites_allele_freqs.txt \
  --min 0.05 --max 1.0 -g $workdir/1b_genos.bad_sites_allele_freqs.png --subpops Ab10 N10 --flip #--debug 
$TASSEL5 -h $workdir/1a_genos.filter_paralogs.hmp.txt.gz -excludeSiteNamesInFile $workdir/1b_genos.bad_sites_allele_freqs.txt -export $workdir/1b_genos.filter_sites.hmp.txt.gz

#   Check for outcrosses via number of rare alleles 
$TASSEL5 -h $genos -filterAlign -filterAlignMinFreq 0.001 -filterAlignMaxFreq 0.05 -export $workdir/1c.rare_alleles.hmp.txt.gz -exportType Hapmap
$TASSEL5 -h $workdir/1c.rare_alleles.hmp.txt.gz -NumericalGenotypePlugin -endPlugin -export $workdir/1c.rare_alleles.numeric.txt -exportType ReferenceProbablity
zcat $workdir/1c.rare_alleles.hmp.txt.gz | head -n1 | cut -f12- | tr '\t' '\n' >  $workdir/1c.full_taxa_names.txt
Rscript 1c_TestForOutcrossedTaxa_RareAlleles.r $workdir/1c.rare_alleles.numeric.txt $workdir/1c.bad_taxa.rare_alleles $workdir/1c.full_taxa_names.txt $parentB73 $parentW23

#   #Test for excess heterozygosity and missingness
missing_cutoff=0.8
het_cutoff=0.5
$TASSEL5 -fork1 -h $workdir/1b_genos.filter_sites.hmp.txt.gz -genotypeSummary taxa -export $workdir/1d_genos.filter_sites.taxasummary.txt -runfork1
Rscript 1d_TestForOutcrossedTaxa_HetsAndMissing.r $workdir/1d_genos.filter_sites.taxasummary.txt $workdir/1d.bad_taxa.hets_and_missing.txt $het_cutoff $missing_cutoff 

#   # Filter out bad taxa so far and remove minor SNP states
$TASSEL5 -h  $workdir/1b_genos.filter_sites.hmp.txt.gz -excludeTaxaInFile $workdir/1c.bad_taxa.rare_alleles.txt -excludeTaxaInFile $workdir/1d.bad_taxa.hets_and_missing.txt  \
  -filterAlign -filterAlignRemMinor -export $workdir/1e_genos.initial_filter.hmp.txt.gz
  
#############
# Transmission Distortion Graphic
#############


# Convert genos to ABH, where B73=A and W23=B
echo $parentB73 > $workdir/2_parentA.txt
echo $parentW23 > $workdir/2_parentB.txt
$TASSEL5 -h $workdir/1e_genos.initial_filter.hmp.txt.gz -GenosToABHPlugin -parentA $workdir/2_parentA.txt -parentB $workdir/2_parentB.txt -o $workdir/2a_genos.abh.txt
   
# Transpose so is easier to deal with; also reformat ABH to AA-AB-BB
Rscript -e "x=read.csv('$workdir/2a_genos.abh.txt', colClasses='character', header=F)"\
  -e "x=t(x); x[1,1:2]=c('site','chrom');" \
  -e "x[x=='A'] = 'AA'; x[x=='B'] = 'BB'; x[x=='H'] = 'AB'; x[is.na(x)] = '-';" \
  -e "write.table(x[,-2], file='$workdir/2b_genos.ab_transposed.txt', sep='\t', quote=F, row.names=F, col.names=F)"    # Delete column 2 (chromosome) b/c is part of site name
   
# Remove all but 1 site within 75 bp so as to not double-count the same GBS tag
python3 2b_FilterNearbySites.py -i $workdir/2b_genos.ab_transposed.txt -o $workdir/2b_genos.ab_transposed.window_filtered.txt -w 75 #--debug
   
# Summarize allele counts
Rscript 2c_SummarizeAlleleCounts.r -i $workdir/2b_genos.ab_transposed.window_filtered.txt -o $workdir/2c_genos.allele_summary.txt
   
# Make graphics
chromlengths=0_Zmays_agpv2_chrom_lengths.txt
winsize=50
step=50
# Main paper graphic, exlcuding sites with freqs before 0.05
python3 2d_GraphTransmissionDistortion.py -i $workdir/2c_genos.allele_summary.txt --window $winsize --step $step \
  --outprefix $workdir/2d_transmission_distortion --chromlengths $chromlengths  --min_freq 0.05  #--debug

# Supplemental graphic, including all sites
python3 2d_GraphTransmissionDistortion.py -i $workdir/2c_genos.allele_summary.txt --window $winsize --step $step \
  --outprefix $workdir/2e_transmission_distortion_supplemental --chromlengths $chromlengths --plot-raw-freqs  #--debug


#######
# Find # of genes that were affected; late addition, so some manual stuff here
#######

# Redo plot with a smaller step size to get more refined regions that are significantly affected by distortion
chromlengths=0_Zmays_agpv2_chrom_lengths.txt
winsize=50
step=10
python3 2d_GraphTransmissionDistortion.py -i $workdir/2c_genos.allele_summary.txt --window $winsize --step $step \
  --outprefix $workdir/2f_transmission_distortion_smaller_step --chromlengths $chromlengths  --min_freq 0.05  #--debug


# Taking the arbitrary division of pval 10e-20 as the "distorted" cutoff, these are the windows for the 3 knobs in this population
# Chrom 4: 4:157630250-end
# chrom 6: 6:129472045-162727151
# chrom 10: 10:95119023-end # Note: the very last window technically falls outside the cutoff, but that's because it's sliding off the chromosome and so has only a few sites in it. The previous few windows have the same end point and much higher significance

gff=0_ZmB73_5b_FGS.genes_only.gff
Rscript 2f_TabulateGenes.r -i  $gff -c 4 --start 157630250 --end 9999999999
Rscript 2f_TabulateGenes.r -i  $gff -c 6 --start 129472045 --end 162727151
Rscript 2f_TabulateGenes.r -i  $gff -c 10 --start 95119023 --end 9999999999
