#!/usr/bin/env bash

# FILE=05_samples_filtered_stats2
# FILE=07_samples_filtered_stats2
# FILE=18_samples_filtered_stats2
# FILE=38_samples_filtered_stats2

# Generate the input file in plink format
plink --vcf ~/Documents/Github/Crypto/Genome_Analysis/products/vcfFiltered/38_samples_filtered_stats2.vcf \
--make-bed --out $FILE --allow-extra-chr --const-fid 0

# ADMIXTURE does not accept chromosome names that are not human chromosomes. We will thus just exchange the first column by 0
awk '{$1="0";print $0}' $FILE.bim > $FILE.bim.tmp
mv $FILE.bim.tmp $FILE.bim
