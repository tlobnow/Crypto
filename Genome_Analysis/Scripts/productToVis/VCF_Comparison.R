####### Look at the differences between samples ################################
# load libraries
library(tidyverse)
library(data.table)

# read the csv files
AA_0866 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0866.csv")
AA_0900 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0900.csv")
AA_0942 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0942.csv")
IXa     <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_IXa.csv")

# calculate variants per GeneID
pull_AA_0866 <- AA_0866 %>% filter(GeneID != ".") %>% group_by(GeneID) %>% dplyr::count(GeneID)
AA_0866 <- left_join(AA_0866, pull_AA_0866, by = 'GeneID')
rm(pull_AA_0866)

# calculate variants per GeneID
pull_AA_0900 <- AA_0900 %>% filter(GeneID != ".") %>% group_by(GeneID) %>% dplyr::count(GeneID)
AA_0900 <- left_join(AA_0900, pull_AA_0900, by = 'GeneID')
rm(pull_AA_0900)

# calculate variants per GeneID
pull_AA_0942 <- AA_0942 %>% filter(GeneID != ".") %>% group_by(GeneID) %>% dplyr::count(GeneID)
AA_0942 <- left_join(AA_0942, pull_AA_0942, by = 'GeneID')
rm(pull_AA_0942)

# calculate variants per GeneID
pull_IXa <- IXa %>% filter(GeneID != ".") %>% group_by(GeneID) %>% dplyr::count(GeneID)
IXa <- left_join(IXa, pull_IXa, by = 'GeneID')
rm(pull_IXa)


################################################################################
#### make subsets of the different samples that exclude intergenic areas
# Rename the REF & ALT Columns to match the respective sample origins
sub866 <- AA_0866 %>% filter(GeneID != '.') %>% dplyr::select(GeneID, n, REF, ALT)
sub900 <- AA_0900 %>% filter(GeneID != '.') %>% dplyr::select(GeneID, n, REF, ALT)
sub942 <- AA_0942 %>% filter(GeneID != '.') %>% dplyr::select(GeneID, n, REF, ALT)
subIXa <-     IXa %>% filter(GeneID != '.') %>% dplyr::select(GeneID, n, REF, ALT)

setnames(sub866, old = c("REF", "ALT"), new = c("REF.866", "ALT.866"), skip_absent = T)
setnames(sub900, old = c("REF", "ALT"), new = c("REF.900", "ALT.900"), skip_absent = T)
setnames(sub942, old = c("REF", "ALT"), new = c("REF.942", "ALT.942"), skip_absent = T)
setnames(subIXa, old = c("REF", "ALT"), new = c("REF.IXa", "ALT.IXa"), skip_absent = T)

#### Make a csv file that contains the variants between different samples
divVariants <- full_join(sub866, sub900)
divVariants <- full_join(divGenes, sub942)
divVariants <- full_join(divGenes, subIXa)

################################################################################
#### Make a csv file that contains the most divergent genes (n=?)

genesub866 <- AA_0866 %>% filter(GeneID != '.') %>% dplyr::select(GeneID, n) %>% distinct(GeneID, n)
genesub900 <- AA_0900 %>% filter(GeneID != '.') %>% dplyr::select(GeneID, n) %>% distinct(GeneID, n)
genesub942 <- AA_0942 %>% filter(GeneID != '.') %>% dplyr::select(GeneID, n) %>% distinct(GeneID, n)
genesubIXa <-     IXa %>% filter(GeneID != '.') %>% dplyr::select(GeneID, n) %>% distinct(GeneID, n)
















