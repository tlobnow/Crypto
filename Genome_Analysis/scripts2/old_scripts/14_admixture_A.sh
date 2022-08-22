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

FILE=~/Documents/Github/Crypto/Genome_Analysis/products/admixture/38_samples_filtered_stats2

admixture --cv $FILE.bed 2 > log2.out
admixture --cv $FILE.bed 3 > log3.out
admixture --cv $FILE.bed 4 > log4.out
admixture --cv $FILE.bed 5 > log5.out
admixture --cv $FILE.bed 6 > log6.out
admixture --cv $FILE.bed 7 > log7.out
admixture --cv $FILE.bed 8 > log8.out
admixture --cv $FILE.bed 9 > log9.out
admixture --cv $FILE.bed 10 > log10.out
admixture --cv $FILE.bed 11 > log11.out
admixture --cv $FILE.bed 12 > log12.out

