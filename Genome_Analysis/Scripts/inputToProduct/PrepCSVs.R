# load libraries
library(tidyverse)
library(data.table)
library(VariantAnnotation)

# set the working directory
setwd(dir = "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/InputFiles/")

## Input = annotated VCF file
#CollapsedVCF <- readVcf("/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/InputFiles/IXa_annotated.vcf")
#CollapsedVCF <- readVcf("/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/InputFiles/866_annotated.vcf")
#CollapsedVCF <- readVcf("/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/InputFiles/900_annotated.vcf")
#CollapsedVCF <- readVcf("/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/InputFiles/942_annotated.vcf")

# Extract interesting data
AminoAcidChange <- as.data.frame(CollapsedVCF@info@listData[["AminoAcidChange"]])
Locus           <- as.data.frame(CollapsedVCF@info@listData[["LocusTag"]]@unlistData)
Product         <- as.data.frame(CollapsedVCF@info@listData[["Product"]]@unlistData)
ProteinID       <- as.data.frame(CollapsedVCF@info@listData[["ProteinID"]]@unlistData)
VariantType     <- as.data.frame(CollapsedVCF@info@listData[["VariantType"]]@unlistData)
FeatureType     <- as.data.frame(CollapsedVCF@info@listData[["FeatureType"]]@unlistData)
NAME            <- as.data.frame(CollapsedVCF@rowRanges@ranges@NAMES)
START           <- as.data.frame(CollapsedVCF@rowRanges@ranges@start)
WIDTH           <- as.data.frame(CollapsedVCF@rowRanges@ranges@width)
QUAL            <- as.data.frame(CollapsedVCF@fixed@listData[["QUAL"]])

# separate the NAMES column
NAME <- setDT(NAME)[, paste0("CollapsedVCF@rowRanges@ranges@NAMES", 1:2) := tstrsplit(CollapsedVCF@rowRanges@ranges@NAMES, ":")]
setnames(NAME, old = c("CollapsedVCF@rowRanges@ranges@NAMES1", "CollapsedVCF@rowRanges@ranges@NAMES2"), new = c("Chr", "Pos.REFALT"), skip_absent = T)
NAME <- setDT(NAME)[, paste0("Pos.REFALT", 1:2) := tstrsplit(Pos.REFALT, "_")]
setnames(NAME, old = c("Pos.REFALT1", "Pos.REFALT2"), new = c("Pos", "REF.ALT"), skip_absent = T)
NAME <- NAME %>% dplyr::select(Chr, Pos, REF.ALT)

AminoAcidChange <- AminoAcidChange %>% dplyr::select(value)

# combine the vectors
CollapsedVCF <- cbind(Locus, Product, ProteinID, VariantType, FeatureType, NAME, AminoAcidChange, START, WIDTH, QUAL)

# change column names
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["LocusTag"]]@unlistData'] <- 'GeneID'
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["ProteinID"]]@unlistData'] <- 'ProteinID'
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["VariantType"]]@unlistData'] <- 'VariantType'
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["FeatureType"]]@unlistData'] <- 'FeatureType'
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["Product"]]@unlistData'] <- 'Product'
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'value'] <- 'AminoAcidChange'
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@rowRanges@ranges@start'] <- 'START'
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@rowRanges@ranges@width'] <- 'WIDTH'
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@fixed@listData[["QUAL"]]'] <- 'QUAL'

CollapsedVCF$Chr <- gsub(pattern = "PYHZ01000001.1", replacement = "Ctyz_1", x = CollapsedVCF$Chr)
CollapsedVCF$Chr <- gsub(pattern = "PYHZ01000002.1", replacement = "Ctyz_2", x = CollapsedVCF$Chr)
CollapsedVCF$Chr <- gsub(pattern = "PYHZ01000003.1", replacement = "Ctyz_3", x = CollapsedVCF$Chr)
CollapsedVCF$Chr <- gsub(pattern = "PYHZ01000004.1", replacement = "Ctyz_4", x = CollapsedVCF$Chr)
CollapsedVCF$Chr <- gsub(pattern = "PYHZ01000005.1", replacement = "Ctyz_5", x = CollapsedVCF$Chr)
CollapsedVCF$Chr <- gsub(pattern = "PYHZ01000006.1", replacement = "Ctyz_6", x = CollapsedVCF$Chr)
CollapsedVCF$Chr <- gsub(pattern = "PYHZ01000007.1", replacement = "Ctyz_7", x = CollapsedVCF$Chr)
CollapsedVCF$Chr <- gsub(pattern = "PYHZ01000008.1", replacement = "Ctyz_8", x = CollapsedVCF$Chr)
CollapsedVCF$Chr <- gsub(pattern = "PYHZ01000009.1", replacement = "Ctyz_00_1", x = CollapsedVCF$Chr)
CollapsedVCF$Chr <- gsub(pattern = "PYHZ01000010.1", replacement = "Ctyz_00_2", x = CollapsedVCF$Chr)
CollapsedVCF$Chr <- gsub(pattern = "PYHZ01000011.1", replacement = "Ctyz_00_3", x = CollapsedVCF$Chr)


## Output = CSV file
#write.csv(CollapsedVCF, "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/Products/Ctyz_IXa.csv")
#write.csv(CollapsedVCF, "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/Products/Ctyz_AA_0866.csv")
#write.csv(CollapsedVCF, "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/Products/Ctyz_AA_0900.csv")
#write.csv(CollapsedVCF, "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/Products/Ctyz_AA_0942.csv")







################################################################################################################################################################
################################################################################################################################################################

# read the csv files
AA_0866 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0866.csv")
AA_0900 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0900.csv")
AA_0942 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0942.csv")
IXa     <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_IXa.csv")

# calculate variants per GeneID
pull_AA_0866 <- AA_0866 %>% filter(GeneID != ".") %>% group_by(GeneID) %>% dplyr::count(GeneID)
AA_0866 <- left_join(AA_0866, pull_AA_0866, by = 'GeneID')

# calculate variants per GeneID
pull_AA_0900 <- AA_0900 %>% filter(GeneID != ".") %>% group_by(GeneID) %>% dplyr::count(GeneID)
AA_0900 <- left_join(AA_0900, pull_AA_0900, by = 'GeneID')

# calculate variants per GeneID
pull_AA_0942 <- AA_0942 %>% filter(GeneID != ".") %>% group_by(GeneID) %>% dplyr::count(GeneID)
AA_0942 <- left_join(AA_0942, pull_AA_0942, by = 'GeneID')

# calculate variants per GeneID
pull_IXa <- IXa %>% filter(GeneID != ".") %>% group_by(GeneID) %>% dplyr::count(GeneID)
IXa <- left_join(IXa, pull_IXa, by = 'GeneID')

rm(pull_IXa, pull_AA_0866, pull_AA_0900, pull_AA_0942)


################################################################################
################################################################################
#### Make a csv file that contains the most divergent genes (n=?)

genesub866 <- AA_0866 %>% filter(GeneID != '.') %>% dplyr::select(GeneID, Chr, n) %>% distinct(GeneID, Chr, n)
genesub900 <- AA_0900 %>% filter(GeneID != '.') %>% dplyr::select(GeneID, Chr, n) %>% distinct(GeneID, Chr, n)
genesub942 <- AA_0942 %>% filter(GeneID != '.') %>% dplyr::select(GeneID, Chr, n) %>% distinct(GeneID, Chr, n)
genesubIXa <-     IXa %>% filter(GeneID != '.') %>% dplyr::select(GeneID, Chr, n) %>% distinct(GeneID, Chr, n)

setnames(genesub866, old = c("n"), new = c("n.866"), skip_absent = T)
setnames(genesub942, old = c("n"), new = c("n.942"), skip_absent = T)
setnames(genesub900, old = c("n"), new = c("n.900"), skip_absent = T)
setnames(genesubIXa, old = c("n"), new = c("n.IXa"), skip_absent = T)

divGenes <- full_join(genesub866, genesub942)
divGenes <- full_join(divGenes, genesub900)
divGenes <- full_join(divGenes, genesubIXa)

# add mean Diversity column (without NAs)
meanDiv <- divGenes %>% dplyr::select(n.866, n.942, n.900, n.IXa) %>% rowMeans(na.rm=TRUE)
divGenes <- cbind(divGenes, meanDiv)

# add Standard Deviation for difference between samples (without NAs)
StDev <- divGenes %>% dplyr::select(n.866, n.942, n.900, n.IXa) 
StDev <- apply(StDev, 1, sd, na.rm = TRUE)
divGenes <- cbind(divGenes, StDev)

# round meanDiv and sd
divGenes$meanDiv <- round(divGenes$meanDiv, digits = 2)
divGenes$StDev <- round(divGenes$StDev, digits = 2)

# write csv file
write.csv(divGenes, "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/Products/DivGenes.csv")

