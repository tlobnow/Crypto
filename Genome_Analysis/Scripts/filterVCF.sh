#!/usr/bin/env bash

VCF_IN=~/Desktop/vcf/06_samples/HZ_UGA_TYGZ1_PAR.vcf
VCF_OUT=~/Desktop/vcfFiltered/06_samples_filtered.vcf.gz

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
