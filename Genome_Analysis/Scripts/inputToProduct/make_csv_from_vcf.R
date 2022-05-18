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
CollapsedVCF <- readVcf("/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/InputFiles/942_annotated.vcf")

# Extract interesting data
#RefCodon        <- as.data.frame(CollapsedVCF@info@listData[["RefCodon"]]@unlistData)
#AltCodon        <- as.data.frame(CollapsedVCF@info@listData[["AltCodon"]]@unlistData)
#RefAminoAcid    <- as.data.frame(CollapsedVCF@info@listData[["RefAminoAcid"]]@unlistData)
#AltAminoAcid    <- as.data.frame(CollapsedVCF@info@listData[["AltAminoAcid"]]@unlistData)
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
CollapsedVCF <- cbind(Locus, Product, ProteinID, VariantType, FeatureType, NAME, 
                      #RefCodon, AltCodon, 
                      AminoAcidChange, START, WIDTH, QUAL)

# change column names
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["LocusTag"]]@unlistData'] <- 'GeneID'
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["ProteinID"]]@unlistData'] <- 'ProteinID'
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["VariantType"]]@unlistData'] <- 'VariantType'
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["FeatureType"]]@unlistData'] <- 'FeatureType'
colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["Product"]]@unlistData'] <- 'Product'
#colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["RefCodon"]]@unlistData'] <- 'RefCodon'
#colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["AltCodon"]]@unlistData'] <- 'AltCodon'
#colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["RefAminoAcid"]]@unlistData'] <- 'RefAminoAcid'
#colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["AltAminoAcid"]]@unlistData'] <- 'AltAminoAcid'
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
write.csv(CollapsedVCF, "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/Products/Ctyz_AA_0942.csv")



