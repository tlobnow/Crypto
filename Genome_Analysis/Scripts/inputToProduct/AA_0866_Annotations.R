# load libraries
library(tidyverse)
library(data.table)
library(VariantAnnotation)

# set the working directory
setwd(dir = "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/InputFiles/")

# load the annotated vcf file
AA_0866 <- readVcf("AA_0866_annotated.vcf")

# Extract interesting data
Locus       <- as.data.frame(AA_0866@info@listData[["LocusTag"]]@unlistData)
ProteinID   <- as.data.frame(AA_0866@info@listData[["ProteinID"]]@unlistData)
VariantType <- as.data.frame(AA_0866@info@listData[["VariantType"]]@unlistData)
QUAL        <- as.data.frame(AA_0866@fixed@listData[["QUAL"]])
NAME        <- as.data.frame(AA_0866@rowRanges@ranges@NAMES)

# separate the NAMES column
NAME <- setDT(NAME)[, paste0("AA_0866@rowRanges@ranges@NAMES", 1:2) := tstrsplit(AA_0866@rowRanges@ranges@NAMES, ":")]
setnames(NAME, old = c("AA_0866@rowRanges@ranges@NAMES1", "AA_0866@rowRanges@ranges@NAMES2"), new = c("chr", "Pos.REFALT"), skip_absent = T)
NAME <- setDT(NAME)[, paste0("Pos.REFALT", 1:2) := tstrsplit(Pos.REFALT, "_")]
setnames(NAME, old = c("Pos.REFALT1", "Pos.REFALT2"), new = c("Pos", "REFALT"), skip_absent = T)
NAME <- setDT(NAME)[, paste0("REFALT", 1:2) := tstrsplit(REFALT, "/")]
setnames(NAME, old = c("REFALT1", "REFALT2"), new = c("REF", "ALT"), skip_absent = T)
NAME <- NAME %>% dplyr::select(chr, Pos, REF, ALT)

# combine the vectors
Ctyz_AA_0866 <- cbind(Locus, ProteinID, VariantType, QUAL, NAME)

# change column names
colnames(Ctyz_AA_0866)[colnames(Ctyz_AA_0866) %in% 'AA_0866@info@listData[["LocusTag"]]@unlistData'] <- 'LocusTag'
colnames(Ctyz_AA_0866)[colnames(Ctyz_AA_0866) %in% 'AA_0866@info@listData[["ProteinID"]]@unlistData'] <- 'ProteinID'
colnames(Ctyz_AA_0866)[colnames(Ctyz_AA_0866) %in% 'AA_0866@info@listData[["VariantType"]]@unlistData'] <- 'VariantType'
colnames(Ctyz_AA_0866)[colnames(Ctyz_AA_0866) %in% 'AA_0866@fixed@listData[["QUAL"]]'] <- 'QUAL'
                                                                        
# write csv
write.csv(Ctyz_AA_0866, "/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/Products/Ctyz_AA_0866.csv")
