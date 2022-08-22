#!/usr/bin/env bash

#### Calculate allele frequency
vcftools --gzvcf ~/Desktop/tyz_final.vcf.gz --freq2 --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/tyz_final

#### Calculate mean depth of coverage per individual
vcftools --gzvcf ~/Desktop/tyz_final.vcf.gz --depth --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/tyz_final

#### Calculate mean depth of coverage for  each site
vcftools --gzvcf ~/Desktop/tyz_final.vcf.gz --site-mean-depth --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/tyz_final

#### Calculate the site quality score for each site
vcftools --gzvcf ~/Desktop/tyz_final.vcf.gz --site-quality --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/tyz_final

#### Calculate the proportion of missing data per individual
vcftools --gzvcf ~/Desktop/tyz_final.vcf.gz --missing-indv --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/tyz_final

#### Calculate the propoortion of missing data per site
vcftools --gzvcf ~/Desktop/tyz_final.vcf.gz --missing-site --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/tyz_final

#### Calculate heterozygosity and inbreeding coefficient per individual
vcftools --gzvcf ~/Desktop/tyz_final.vcf.gz --het --out      ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/tyz_final
