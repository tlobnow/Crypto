#!/usr/bin/env bash

~/Documents/Github/gatk/build/gatk-4.2.6.1-15-ga3141f7-SNAPSHOT/gatk FastaAlternateReferenceMaker \
-R ~/Documents/Github/Crypto/Genome_Analysis/resources/CryptoDB-57_CtyzzeriUGA55_Genome.fasta \
-V ~/Documents/Github/Crypto/Genome_Analysis/products/vcfFiltered/05_samples_filtered_stats2.vcf.gz \
-O 05_samples.fasta
