#!/usr/bin/env bash

# set filters
MAF=0.05
MISS=0.75
QUAL=30
MIN_DEPTH=10
MAX_DEPTH=150


VCF_IN=~/Documents/Github/Crypto/Genome_Analysis/products/vcf/15_samples.vcf.gz
VCF_OUT=~/Documents/Github/Crypto/Genome_Analysis/products/vcfFiltered/15_samples_filtered.vcf.gz

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT

VCF_IN=~/Documents/Github/Crypto/Genome_Analysis/products/vcf/12_samples.vcf.gz
VCF_OUT=~/Documents/Github/Crypto/Genome_Analysis/products/vcfFiltered/12_samples_filtered.vcf.gz

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT

VCF_IN=~/Documents/Github/Crypto/Genome_Analysis/products/vcf/07_samples.vcf.gz
VCF_OUT=~/Documents/Github/Crypto/Genome_Analysis/products/vcfFiltered/07_samples_filtered.vcf.gz

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT

VCF_IN=~/Documents/Github/Crypto/Genome_Analysis/products/vcf/06_samples.vcf.gz
VCF_OUT=~/Documents/Github/Crypto/Genome_Analysis/products/vcfFiltered/06_samples_filtered.vcf.gz

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT

VCF_IN=~/Documents/Github/Crypto/Genome_Analysis/products/vcf/05_samples.vcf.gz
VCF_OUT=~/Documents/Github/Crypto/Genome_Analysis/products/vcfFiltered/05_samples_filtered.vcf.gz

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT

VCF_IN=~/Documents/Github/Crypto/Genome_Analysis/products/vcf/04_samples.vcf.gz
VCF_OUT=~/Documents/Github/Crypto/Genome_Analysis/products/vcfFiltered/04_samples_filtered.vcf.gz

# perform the filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL \
--min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT
