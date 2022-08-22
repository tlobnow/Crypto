#!/usr/bin/env bash

# 05_samples
~/Documents/Github/gatk/build/gatk-4.2.6.1-15-ga3141f7-SNAPSHOT/gatk VariantsToTable \
-V /Users/finnlo/Documents/Github/Crypto/Genome_Analysis/products/vcfAnnotated/05_samples_filtered_stats2_annotated.vcf \
-F CHROM -F POS -F TYPE -GF AD \
-O 05_variants.txt


# 07_samples
~/Documents/Github/gatk/build/gatk-4.2.6.1-15-ga3141f7-SNAPSHOT/gatk VariantsToTable \
-V /Users/finnlo/Documents/Github/Crypto/Genome_Analysis/products/vcfAnnotated/07_samples_filtered_stats2_annotated.vcf \
-F CHROM -F POS -F TYPE -GF AD \
-O 07_variants.txt

# 18_samples
~/Documents/Github/gatk/build/gatk-4.2.6.1-15-ga3141f7-SNAPSHOT/gatk VariantsToTable \
-V /Users/finnlo/Documents/Github/Crypto/Genome_Analysis/products/vcfAnnotated/18_samples_filtered_stats2_annotated.vcf \
-F CHROM -F POS -F TYPE -GF AD \
-O 18_variants.txt

# 38_samples
~/Documents/Github/gatk/build/gatk-4.2.6.1-15-ga3141f7-SNAPSHOT/gatk VariantsToTable \
-V /Users/finnlo/Documents/Github/Crypto/Genome_Analysis/products/vcfAnnotated/38_samples_filtered_stats2_annotated.vcf \
-F CHROM -F POS -F TYPE -GF AD \
-O 38_variants.txt
