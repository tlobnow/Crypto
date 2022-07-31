#!/usr/bin/env bash

# set filters
MAF=0.1
MISS=0.85
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=135.01


VCF_IN=~/Documents/Github/Crypto/Genome_Analysis/products/vcf/07_samples.vcf.gz
VCF_OUT=~/Documents/Github/Crypto/Genome_Analysis/products/vcfFiltered/07_samples_filtered_stats2.vcf.gz

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT

