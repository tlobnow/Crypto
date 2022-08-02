#!/usr/bin/env bash

FILE=~/Documents/Github/Crypto/Genome_Analysis/products/admixture/05_samples_filtered_stats2

admixture --cv $FILE.bed 2 > log2.out
admixture --cv $FILE.bed 3 > log3.out
admixture --cv $FILE.bed 4 > log4.out
admixture --cv $FILE.bed 5 > log5.out

FILE=~/Documents/Github/Crypto/Genome_Analysis/products/admixture/07_samples_filtered_stats2

admixture --cv $FILE.bed 2 > log2.out
admixture --cv $FILE.bed 3 > log3.out
admixture --cv $FILE.bed 4 > log4.out
admixture --cv $FILE.bed 5 > log5.out

FILE=~/Documents/Github/Crypto/Genome_Analysis/products/admixture/18_samples_filtered_stats2

admixture --cv $FILE.bed 2 > log2.out
admixture --cv $FILE.bed 3 > log3.out
admixture --cv $FILE.bed 4 > log4.out
admixture --cv $FILE.bed 5 > log5.out
