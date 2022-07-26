#!/usr/bin/env bash

~/plink-1.9/plink \
--bfile /SAN/Ctyzzeri/genomics_bootcamp/resources/ADAPTmap_genotypeTOP_20160222_full \
--cow --nonfounders --allow-no-sex --recode \
--out ADAPTmap_TOP \
--geno 0.1 \
--mind 0.1 \
--maf 0.05
