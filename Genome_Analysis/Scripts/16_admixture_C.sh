#!/usr/bin/env bash

FILE=05_samples_filtered_stats2

Rscript plotADMIXTURE.r \
-p ~/Documents/Github/Crypto/Genome_Analysis/products/admixture/$FILE \
-i ~/Documents/Github/Crypto/Genome_Analysis/products/admixture/$FILE.list \
-k 5 -l TYZ1,TYZ2,TYZ3,TYZ4,TYZ5


FILE=07_samples_filtered_stats2

Rscript plotADMIXTURE.r \
-p ~/Documents/Github/Crypto/Genome_Analysis/products/admixture/$FILE \
-i ~/Documents/Github/Crypto/Genome_Analysis/products/admixture/$FILE.list \
-k 5 -l TYZ1,TYZ2,TYZ3,HOM,TYZ4,TYZ5,PAR


FILE=18_samples_filtered_stats2

Rscript plotADMIXTURE.r \
-p ~/Documents/Github/Crypto/Genome_Analysis/products/admixture/$FILE \
-i ~/Documents/Github/Crypto/Genome_Analysis/products/admixture/$FILE.list \
-k 5 -l CHN_C_GU1,CHN_C_GU2,CHN_C_GU3,CHN_C_SH1,CHN_C_SH2,CHN_M_GU,EUR_C_CZ,EUR_H_CZ,EUR_M_866,EUR_M_900,EUR_M_942,USA_C_WA1,USA_C_WA2,USA_C_WI,USA_H_AL1,USA_H_AL2,USA_H_OR,USA_M_GA
