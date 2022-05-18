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
