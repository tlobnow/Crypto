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
-k 5 -l CHN_C_GU1,CHN_C_GU2,CHN_C_GU3,CHN_C_SH1,CHN_C_SH2,CHN_M_GU,EUR_C_CZ,EUR_H_CZ,\
EUR_M_866,EUR_M_900,EUR_M_942,USA_C_WA1,USA_C_WA2,USA_C_WI,USA_H_AL1,USA_H_AL2,USA_H_OR,USA_M_GA


FILE=38_samples_filtered_stats2

Rscript plotADMIXTURE.r \
-p ~/Documents/Github/Crypto/Genome_Analysis/products/admixture/$FILE \
-i ~/Documents/Github/Crypto/Genome_Analysis/products/admixture/$FILE.list \
-k 5 -l EUR_T_866,EUR_T_942,EUR_T_900,USA_T_GA,CHN_T_GU,EUR_P_CZ2,USA_P_OR,USA_P_AL1,USA_P_AL2,EUR_P_CZ1,\
CHN_P_GU1,CHN_P_GU2,CHN_P_GU3,CHN_P_SH1,CHN_P_SH2,USA_P_WI,USA_P_WA1,USA_P_WA2,\
COL_H_MD,USA_H_MO,USA_H_ID,UGA_H_KA1,UGA_H_KA2,EUR_H_WL1,EUR_H_WL2,EUR_H_WL3,EUR_H_WL4,EUR_H_WL5,\
AFR_H_GH1,AFR_H_GH2,AFR_H_GH3,AFR_H_TZ1,AFR_H_TZ2,AFR_H_TZ3,AFR_H_TZ4,MDG_H_1,MDG_H_2,MDG_H_3
