### Generate statistics on minimum depth, genotype quality, minor allele frequency, and missing data
### to set filters for the vcf in the next script

INDS=($(for i in ~/Documents/Github/Crypto/Genome_Analysis/products/vcf/*.vcf.gz; do echo $(basename -a -s *.vcf.gz $i); done))

for IND in ${INDS[@]};
do

        # Calculate allele frequency
        vcftools --gzvcf ~/Documents/Github/Crypto/Genome_Analysis/products/vcf/${IND} --freq2 --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/${IND}

        # Calculate mean depth of coverage per individual
        vcftools --gzvcf ~/Documents/Github/Crypto/Genome_Analysis/products/vcf/${IND} --depth --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/${IND}

        # Calculate mean depth of coverage for  each site
        vcftools --gzvcf ~/Documents/Github/Crypto/Genome_Analysis/products/vcf/${IND} --site-mean-depth --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/${IND}

        # Calculate the site quality score for each site
        vcftools --gzvcf ~/Documents/Github/Crypto/Genome_Analysis/products/vcf/${IND} --site-quality --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/${IND}

        # Calculate the proportion of missing data per individual
        vcftools --gzvcf ~/Documents/Github/Crypto/Genome_Analysis/products/vcf/${IND} --missing-indv --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/${IND}

        # Calculate the propoortion of missing data per site
        vcftools --gzvcf ~/Documents/Github/Crypto/Genome_Analysis/products/vcf/${IND} --missing-site --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/${IND}

        # Calculate heterozygosity and inbreeding coefficient per individual
        vcftools --gzvcf ~/Documents/Github/Crypto/Genome_Analysis/products/vcf/${IND} --het --out ~/Documents/Github/Crypto/Genome_Analysis/products/vcftools/${IND}

done
