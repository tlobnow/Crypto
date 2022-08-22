#!/usr/bin/env bash

### RUN THIS SCRIPT IN THE ACTUAL COMMAND LINE FOR EACH OF THE FILES.. IT DOES NOT WORK YET..

# FILE=/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/products/admixture/05_samples_filtered_stats2
# FILE=/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/products/admixture/07_samples_filtered_stats2
# FILE=/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/products/admixture/18_samples_filtered_stats2
FILE=/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/products/admixture/38_samples_filtered_stats2


grep "CV" *out | awk '{print $3,$4}' | cut -c 4,7-20 > $FILE.cv.error
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > $FILE.cv.error
