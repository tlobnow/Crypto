#!/usr/bin/env bash

INDS=($(for i in /SAN/Ctyzzeri/all/align/*.bam; do echo $(basename -a -s .bam $i); done))

for IND in ${INDS[@]};
do

	gatk AddOrReplaceReadGroups \
	-I /SAN/Ctyzzeri/all/results/align/${IND}.bam \
	-O /SAN/Ctyzzeri/all/results/fixed/${IND}_fixed.bam \
	-RGID 4 \
	-RGLB lib1 \
	-RGPL ILLUMINA \
	-RGPU unit1 \
	-RGSM 20

done
