#!/usr/bin/env bash

~/Documents/Github/gatk/build/gatk-4.2.6.1-15-ga3141f7-SNAPSHOT/gatk HaplotypeCaller \
--sample-ploidy 1 --native-pair-hmm-threads 8 \
-I ~/Desktop/USA_T_GA_sort.rmd.bam \
-R ~/Desktop/CryptoDB-57_CtyzzeriUGA55_Genome.fasta \
-O ~/Documents/Github/Crypto/Genome_Analysis/products/gatkVCF/USA_T_GA.vcf \

### --annotate-with-num-discovered-alleles false 
### --heterozygosity 0.001 
### --indel-heterozygosity 1.25E-4 --heterozygosity-stdev 0.01 \
### --standard-min-confidence-threshold-for-calling 30.0 
--max-alternate-alleles 6 --max-genotype-count 1024 --num-reference-samples-if-no-call 0 \
--contamination-fraction-to-filter 0.0 --output-mode EMIT_VARIANTS_ONLY --all-site-pls false --gvcf-gq-bands 1 --gvcf-gq-bands 2 --gvcf-gq-bands 3 \
--gvcf-gq-bands 4 --gvcf-gq-bands 5 --gvcf-gq-bands 6 --gvcf-gq-bands 7 --gvcf-gq-bands 8 --gvcf-gq-bands 9 --gvcf-gq-bands 10 --gvcf-gq-bands 11 \
--gvcf-gq-bands 12 --gvcf-gq-bands 13 --gvcf-gq-bands 14 --gvcf-gq-bands 15 --gvcf-gq-bands 16 --gvcf-gq-bands 17 --gvcf-gq-bands 18 --gvcf-gq-bands 19 \
--gvcf-gq-bands 20 --gvcf-gq-bands 21 --gvcf-gq-bands 22 --gvcf-gq-bands 23 --gvcf-gq-bands 24 --gvcf-gq-bands 25 --gvcf-gq-bands 26 --gvcf-gq-bands 27 \
--gvcf-gq-bands 28 --gvcf-gq-bands 29 --gvcf-gq-bands 30 --gvcf-gq-bands 31 --gvcf-gq-bands 32 --gvcf-gq-bands 33 --gvcf-gq-bands 34 --gvcf-gq-bands 35 \
--gvcf-gq-bands 36 --gvcf-gq-bands 37 --gvcf-gq-bands 38 --gvcf-gq-bands 39 --gvcf-gq-bands 40 --gvcf-gq-bands 41 --gvcf-gq-bands 42 --gvcf-gq-bands 43 \
--gvcf-gq-bands 44 --gvcf-gq-bands 45 --gvcf-gq-bands 46 --gvcf-gq-bands 47 --gvcf-gq-bands 48 --gvcf-gq-bands 49 --gvcf-gq-bands 50 --gvcf-gq-bands 51 \
--gvcf-gq-bands 52 --gvcf-gq-bands 53 --gvcf-gq-bands 54 --gvcf-gq-bands 55 --gvcf-gq-bands 56 --gvcf-gq-bands 57 --gvcf-gq-bands 58 --gvcf-gq-bands 59 \
--gvcf-gq-bands 60 --gvcf-gq-bands 70 --gvcf-gq-bands 80 --gvcf-gq-bands 90 --gvcf-gq-bands 99 --floor-blocks false --indel-size-to-eliminate-in-ref-model 10 \
--disable-optimizations false --just-determine-active-regions false --dont-genotype false --do-not-run-physical-phasing false --do-not-correct-overlapping-quality false \
--use-filtered-reads-for-annotations false --adaptive-pruning false --do-not-recover-dangling-branches false --recover-dangling-heads false --kmer-size 10 \
--kmer-size 25 --dont-increase-kmer-sizes-for-cycles false --allow-non-unique-kmers-in-ref false --num-pruning-samples 1 --min-dangling-branch-length 4 \
--recover-all-dangling-branches false --max-num-haplotypes-in-population 128 --min-pruning 2 --adaptive-pruning-initial-error-rate 0.001  \
--pruning-lod-threshold 2.302585092994046 --max-unpruned-variants 100 --linked-de-bruijn-graph false --disable-artificial-haplotype-recovery false  \
--debug-assembly false --debug-graph-transformations false --capture-assembly-failure-bam false --error-correction-log-odds -Infinity --error-correct-reads false \
--kmer-length-for-read-error-correction 25 --min-observations-for-kmer-to-be-solid 20 --base-quality-score-threshold 18 --pair-hmm-gap-continuation-penalty 10  \
--pair-hmm-implementation FASTEST_AVAILABLE --pcr-indel-model CONSERVATIVE --phred-scaled-global-read-mismapping-rate 45 --native-pair-hmm-use-double-precision false \
--bam-writer-type CALLED_HAPLOTYPES --dont-use-soft-clipped-bases false --min-base-quality-score 10 \
--smith-waterman JAVA --emit-ref-confidence NONE --max-mnp-distance 0 --force-call-filtered-alleles false --allele-informative-reads-overlap-margin 2 \
--min-assembly-region-size 50 --max-assembly-region-size 300 --active-probability-threshold 0.002 --max-prob-propagation-distance 50 --force-active false \
--assembly-region-padding 100 --padding-around-indels 75 --padding-around-snps 20 --padding-around-strs 75 \
--max-reads-per-alignment-start 50 --interval-set-rule UNION --interval-padding 0 --interval-exclusion-padding 0 --interval-merging-rule ALL \
--read-validation-stringency SILENT --seconds-between-progress-updates 10.0 \
--disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true \
--create-output-variant-md5 false --lenient false --add-output-sam-program-record true \
--add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 \
--disable-bam-index-caching false --sites-only-vcf-output false --help false \
--version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 \
--gcs-project-for-requester-pays  --disable-tool-default-read-filters false --minimum-mapping-quality 20 --disable-tool-default-annotations false \
--enable-all-annotations false --allow-old-rms-mapping-quality-annotation-data false
