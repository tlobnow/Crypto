#!/usr/bin/env bash

VCF_IN=~/Documents/Github/Crypto/Genome_Analysis/vcf/HZ_UGA_TYGZ1_PAR.vcf.gz
VCF_OUT=~/Documents/Github/Crypto/Genome_Analysis/vcf/HZ_UGA_TYGZ1_PAR_filtered.vcf.gz

# set filters
MAF=0.1
MISS=0.9
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=50

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT
