#!/usr/bin/env bash

### Generate statistics on minimum depth, genotype quality, minor allele frequency, and missing data 
### to set filters for the vcf in the next script

VCF=~/Documents/Github/Crypto/Genome_Analysis/products/vcf/38_samples.vcf.gz
OUT=~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/38_samples

# Calculate allele frequency
vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2

# Calculate mean depth of coverage per individual
vcftools --gzvcf $VCF --depth --out $OUT

# Calculate mean depth of coverage for  each site
vcftools --gzvcf $VCF --site-mean-depth --out $OUT

# Calculate the site quality score for each site
vcftools --gzvcf $VCF --site-quality --out $OUT

# Calculate the proportion of missing data per individual
vcftools --gzvcf $VCF --missing-indv --out $OUT

# Calculate the propoortion of missing data per site
vcftools --gzvcf $VCF --missing-site --out $OUT

# Calculate heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf $VCF --het --out $OUT
