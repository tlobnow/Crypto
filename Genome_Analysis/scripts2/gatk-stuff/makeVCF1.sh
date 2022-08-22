#!/usr/bin/env bash

gatk HaplotypeCaller \

--sample-ploidy 1 \
--input ~/Desktop/USA_T_GA_fixed.bam \
--reference ~/Desktop/CryptoDB-57_CtyzzeriUGA55_Genome.fasta \
--output ~/Desktop/USA_T_GA.vcf
