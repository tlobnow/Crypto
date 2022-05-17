# load libraries
library(tidyverse)
library(data.table)
library(VariantAnnotation)

# set the working directory
setwd(dir = "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/InputFiles/")

# load the annotated vcf file
IXa <- readVcf("IXa_annotated.vcf")

# Extract interesting data
Locus       <- as.data.frame(IXa@info@listData[["LocusTag"]]@unlistData)
ProteinID   <- as.data.frame(IXa@info@listData[["ProteinID"]]@unlistData)
VariantType <- as.data.frame(IXa@info@listData[["VariantType"]]@unlistData)
QUAL        <- as.data.frame(IXa@fixed@listData[["QUAL"]])
NAME        <- as.data.frame(IXa@rowRanges@ranges@NAMES)

# separate the NAMES column
NAME <- setDT(NAME)[, paste0("IXa@rowRanges@ranges@NAMES", 1:2) := tstrsplit(IXa@rowRanges@ranges@NAMES, ":")]
setnames(NAME, old = c("IXa@rowRanges@ranges@NAMES1", "IXa@rowRanges@ranges@NAMES2"), new = c("chr", "Pos.REFALT"), skip_absent = T)
NAME <- setDT(NAME)[, paste0("Pos.REFALT", 1:2) := tstrsplit(Pos.REFALT, "_")]
setnames(NAME, old = c("Pos.REFALT1", "Pos.REFALT2"), new = c("Pos", "REFALT"), skip_absent = T)
NAME <- setDT(NAME)[, paste0("REFALT", 1:2) := tstrsplit(REFALT, "/")]
setnames(NAME, old = c("REFALT1", "REFALT2"), new = c("REF", "ALT"), skip_absent = T)
NAME <- NAME %>% dplyr::select(chr, Pos, REF, ALT)

# combine the vectors
Ctyz_IXa <- cbind(Locus, ProteinID, VariantType, QUAL, NAME)

# change column names
colnames(Ctyz_IXa)[colnames(Ctyz_IXa) %in% 'IXa@info@listData[["LocusTag"]]@unlistData'] <- 'LocusTag'
colnames(Ctyz_IXa)[colnames(Ctyz_IXa) %in% 'IXa@info@listData[["ProteinID"]]@unlistData'] <- 'ProteinID'
colnames(Ctyz_IXa)[colnames(Ctyz_IXa) %in% 'IXa@info@listData[["VariantType"]]@unlistData'] <- 'VariantType'
colnames(Ctyz_IXa)[colnames(Ctyz_IXa) %in% 'IXa@fixed@listData[["QUAL"]]'] <- 'QUAL'

Ctyz_IXa$chr[Ctyz_IXa$chr %in% c("PYHZ01000001.1", "PYHZ01000002.1", "PYHZ01000003.1",
                                 "PYHZ01000004.1", "PYHZ01000005.1", "PYHZ01000006.1",
                                 "PYHZ01000007.1", "PYHZ01000008.1", "PYHZ01000009.1",
                                 "PYHZ01000010.1", "PYHZ01000011.1")] <- c("Ctyz_1", "Ctyz_2", "Ctyz_3",
                                                                           "Ctyz_4", "Ctyz_5", "Ctyz_6", 
                                                                           "Ctyz_7", "Ctyz_8", "Ctyz_00_1",
                                                                           "Ctyz_00_2", "Ctyz_00_3")
# write csv
write.csv(Ctyz_IXa, "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/Products/Ctyz_IXa.csv")
