# load libraries
library(tidyverse)
library(data.table)
library(VariantAnnotation)

makeVcfToCsv <- function(x){
  CollapsedVCF <- readVcf(x)
  # Extract interesting data
  Locus       <- as.data.frame(CollapsedVCF@info@listData[["LocusTag"]]@unlistData)
  ProteinID   <- as.data.frame(CollapsedVCF@info@listData[["ProteinID"]]@unlistData)
  VariantType <- as.data.frame(CollapsedVCF@info@listData[["VariantType"]]@unlistData)
  FeatureType <- as.data.frame(CollapsedVCF@info@listData[["FeatureType"]]@unlistData)
  Product     <- as.data.frame(CollapsedVCF@info@listData[["Product"]]@unlistData)
  QUAL        <- as.data.frame(CollapsedVCF@fixed@listData[["QUAL"]])
  NAME        <- as.data.frame(CollapsedVCF@rowRanges@ranges@NAMES)
  
  # separate the NAMES column
  NAME <- setDT(NAME)[, paste0("CollapsedVCF@rowRanges@ranges@NAMES", 1:2) := tstrsplit(CollapsedVCF@rowRanges@ranges@NAMES, ":")]
  setnames(NAME, old = c("CollapsedVCF@rowRanges@ranges@NAMES1", "CollapsedVCF@rowRanges@ranges@NAMES2"), new = c("chr", "Pos.REFALT"), skip_absent = T)
  NAME <- setDT(NAME)[, paste0("Pos.REFALT", 1:2) := tstrsplit(Pos.REFALT, "_")]
  setnames(NAME, old = c("Pos.REFALT1", "Pos.REFALT2"), new = c("Pos", "REFALT"), skip_absent = T)
  NAME <- setDT(NAME)[, paste0("REFALT", 1:2) := tstrsplit(REFALT, "/")]
  setnames(NAME, old = c("REFALT1", "REFALT2"), new = c("REF", "ALT"), skip_absent = T)
  NAME <- NAME %>% dplyr::select(chr, Pos, REF, ALT)
  
  # combine the vectors
  CollapsedVCF <- cbind(Locus, ProteinID, VariantType, FeatureType, QUAL, NAME)
  
  # change column names
  colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["LocusTag"]]@unlistData'] <- 'GeneID'
  colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["ProteinID"]]@unlistData'] <- 'ProteinID'
  colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["VariantType"]]@unlistData'] <- 'VariantType'
  colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["FeatureType"]]@unlistData'] <- 'FeatureType'
  colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@info@listData[["Product"]]@unlistData'] <- 'Product'
  colnames(CollapsedVCF)[colnames(CollapsedVCF) %in% 'CollapsedVCF@fixed@listData[["QUAL"]]'] <- 'QUAL'
}


# IXa
makeVcfToCsv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/InputFiles/IXa_annotated.vcf")
CollapsedVCF$chr[CollapsedVCF$chr %in% c("PYHZ01000001.1", "PYHZ01000002.1", "PYHZ01000003.1", "PYHZ01000004.1", "PYHZ01000005.1", "PYHZ01000006.1", "PYHZ01000007.1", "PYHZ01000008.1", "PYHZ01000009.1", "PYHZ01000010.1", "PYHZ01000011.1")] <- c("Ctyz_1", "Ctyz_2", "Ctyz_3", "Ctyz_4", "Ctyz_5", "Ctyz_6", "Ctyz_7", "Ctyz_8", "Ctyz_00_1", "Ctyz_00_2", "Ctyz_00_3")
write.csv(CollapsedVCF, "~/Genome_Analysis/Products/Ctyz_IXa.csv")

# AA_0866
makeVcfToCsv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/InputFiles/866_annotated.vcf")
write.csv(CollapsedVCF, "~/Genome_Analysis/Products/Ctyz_AA_0866.csv")

# AA_0900
makeVcfToCsv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/InputFiles/900_annotated.vcf")
write.csv(CollapsedVCF, "~/Genome_Analysis/Products/Ctyz_AA_0900.csv")

# AA_0942
makeVcfToCsv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/InputFiles/942_annotated.vcf")
write.csv(CollapsedVCF, "~/Genome_Analysis/Products/Ctyz_AA_0942.csv")
