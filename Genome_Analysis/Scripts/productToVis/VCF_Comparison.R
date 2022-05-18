####### Look at the differences between samples ################################
# load libraries
library(tidyverse)
library(data.table)
library(VariantAnnotation)

# set the working directory
setwd("/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/Products/")

# read the csv files
AA_0866 <- read.csv("Ctyz_AA_0866.csv")
AA_0900 <- read.csv("Ctyz_AA_0900.csv")
AA_0942 <- read.csv("Ctyz_AA_0942.csv")
IXa <- read.csv("Ctyz_IXa.csv")

