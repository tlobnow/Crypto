#!/usr/bin/env bash

FILE=tyz_final

vcftools --gzvcf /Users/finnlo/Documents/Github/Crypto/Genome_Analysis/products/vcf/$FILE.vcf.gz --max-missing 1 --recode --stdout | gzip > $FILE.noN.vcf.gz

