# load libraries
library(tidyverse)
library(data.table)
library(VariantAnnotation)

# set the working directory
setwd(dir = "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/InputFiles/")

## Input = annotated VCF file
#CollapsedVCF <- readVcf("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/InputFiles/IXa_annotated.vcf")
#CollapsedVCF <- readVcf("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/InputFiles/866_annotated.vcf")
#CollapsedVCF <- readVcf("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/InputFiles/900_annotated.vcf")
#CollapsedVCF <- readVcf("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/InputFiles/942_annotated.vcf")


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
NAME <- setDT(NAME)[, paste0("REF.ALT", 1:2) := tstrsplit(REF.ALT, "/")]
setnames(NAME, old = c("REF.ALT1", "REF.ALT2"), new = c("REF", "ALT"), skip_absent = T)
NAME <- NAME %>% dplyr::select(Chr, REF, ALT)


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


################################################################################

ct_AA_0866 <- AA_0866 %>% group_by(GeneID, Chr, FeatureType) %>% dplyr::count(VariantType)
ct_AA_0942 <- AA_0942 %>% group_by(GeneID, Chr, FeatureType) %>% dplyr::count(VariantType)
ct_AA_0900 <- AA_0900 %>% group_by(GeneID, Chr, FeatureType) %>% dplyr::count(VariantType)
ct_IXa     <-     IXa %>% group_by(GeneID, Chr, FeatureType) %>% dplyr::count(VariantType)

setnames(ct_AA_0866, old = c("VariantType", "n"), new = c("VariantType", "Variants.866"), skip_absent = T)
setnames(ct_AA_0942, old = c("VariantType", "n"), new = c("VariantType", "Variants.942"), skip_absent = T)
setnames(ct_AA_0900, old = c("VariantType", "n"), new = c("VariantType", "Variants.900"), skip_absent = T)
setnames(ct_IXa,     old = c("VariantType", "n"), new = c("VariantType", "Variants.IXa"), skip_absent = T)

divGenes2 <- full_join(ct_AA_0866, ct_AA_0942)
divGenes2 <- full_join(divGenes2, ct_AA_0900)
divGenes2 <- full_join(divGenes2, ct_IXa)


divGenes3 <- full_join(divGenes, divGenes2)


### SUBSET ALL SNPs
sub_SNP <- divGenes3 %>% filter(VariantType == "SNP")
colnames(sub_SNP)[colnames(sub_SNP)%in%"Variants.866"] <- "SNPs.866"
colnames(sub_SNP)[colnames(sub_SNP)%in%"Variants.942"] <- "SNPs.942"
colnames(sub_SNP)[colnames(sub_SNP)%in%"Variants.900"] <- "SNPs.900"
colnames(sub_SNP)[colnames(sub_SNP)%in%"Variants.IXa"] <- "SNPs.IXa"


### SUBSET AMBIGUOUS SNPs
sub_AmbigSNP <- divGenes3 %>% filter(VariantType == "Ambiguous_SNP")
colnames(sub_AmbigSNP)[colnames(sub_AmbigSNP)%in%"Variants.866"] <- "AmSNPs.866"
colnames(sub_AmbigSNP)[colnames(sub_AmbigSNP)%in%"Variants.942"] <- "AmSNPs.942"
colnames(sub_AmbigSNP)[colnames(sub_AmbigSNP)%in%"Variants.900"] <- "AmSNPs.900"
colnames(sub_AmbigSNP)[colnames(sub_AmbigSNP)%in%"Variants.IXa"] <- "AmSNPs.IXa"

### SUBSET ALL INSERTIONS
sub_INS <- divGenes3 %>% filter(VariantType == "Insertion")
colnames(sub_INS)[colnames(sub_INS)%in%"Variants.866"] <- "INS.866"
colnames(sub_INS)[colnames(sub_INS)%in%"Variants.942"] <- "INS.942"
colnames(sub_INS)[colnames(sub_INS)%in%"Variants.900"] <- "INS.900"
colnames(sub_INS)[colnames(sub_INS)%in%"Variants.IXa"] <- "INS.IXa"


### SUBSET ALL DELETIONS
sub_DEL <- divGenes3 %>% filter(VariantType == "Deletion")
colnames(sub_DEL)[colnames(sub_DEL)%in%"Variants.866"] <- "DEL.866"
colnames(sub_DEL)[colnames(sub_DEL)%in%"Variants.942"] <- "DEL.942"
colnames(sub_DEL)[colnames(sub_DEL)%in%"Variants.900"] <- "DEL.900"
colnames(sub_DEL)[colnames(sub_DEL)%in%"Variants.IXa"] <- "DEL.IXa"

### JOIN THE SUBSETS INTO THE FINISHED TABLE
divGenes4 <- full_join(sub_SNP, sub_AmbigSNP)
divGenes4 <- full_join(divGenes4, sub_INS)
divGenes4 <- full_join(divGenes4, sub_DEL)
divGenes4 <- divGenes4 %>% select(-VariantType)


divGenes4 <- divGenes4 %>% 
  filter(GeneID != ".") %>%
  arrange(GeneID) %>% 
  group_by(GeneID) %>% 
  fill(c(everything()), .direction = "downup") %>% 
  ungroup() %>% 
  distinct(GeneID, .keep_all = T)

colnames(divGenes4)[colnames(divGenes4)%in%"n.866"] <- "ChangesSum.866"
colnames(divGenes4)[colnames(divGenes4)%in%"n.942"] <- "ChangesSum.942"
colnames(divGenes4)[colnames(divGenes4)%in%"n.900"] <- "ChangesSum.900"
colnames(divGenes4)[colnames(divGenes4)%in%"n.IXa"] <- "ChangesSum.IXa"

# add Product Column
Product <- AA_0866 %>% filter(GeneID != ".") %>% select(GeneID, Product) %>% arrange(GeneID) %>% group_by(GeneID) %>% fill(c(everything()), .direction = "downup") %>% ungroup() %>% distinct(GeneID, .keep_all = T)
divGenes4 <- full_join(divGenes4, Product)


# add mean Diversity column (without NAs)
meanDiversity <- divGenes4 %>% dplyr::select(ChangesSum.866, ChangesSum.942, ChangesSum.900, ChangesSum.IXa) %>% rowMeans(na.rm=TRUE)
divGenes4 <- cbind(divGenes4, meanDiversity)

# add Standard Deviation for difference between samples (without NAs)
StDev <- divGenes4 %>% dplyr::select(ChangesSum.866, ChangesSum.942, ChangesSum.900, ChangesSum.IXa) 
StDev <- apply(StDev, 1, sd, na.rm = TRUE)
divGenes4 <- cbind(divGenes4, StDev)

# round meanDiv and sd
divGenes4$meanDiversity <- round(divGenes4$meanDiversity, digits = 2)
divGenes4$StDev <- round(divGenes4$StDev, digits = 2)

# arrange cols in specific order
DivGenes <- divGenes4 %>% select(GeneID, Chr, FeatureType, Product,
                                 meanDiversity, StDev, 3:6, 8:23)


################################################################################
################################################################################
# write csv file
write.csv(DivGenes, "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/Products/divGenes.csv")

DivGenes <- read.csv("/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/Products/divGenes.csv")








