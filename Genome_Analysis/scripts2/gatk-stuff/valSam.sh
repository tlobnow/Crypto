#!/usr/bin/env bash

gatk ValidateVariants \
-R ~/Documents/Github/Crypto/Genome_Analysis/resources/CryptoDB-57_CtyzzeriUGA55_Genome.fasta \
-V ~/Documents/Github/Crypto/Genome_Analysis/products/vcfAnnotated/38_samples_filtered_stats2_annotated.vcf.gz
