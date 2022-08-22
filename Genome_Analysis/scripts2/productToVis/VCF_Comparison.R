####### Look at the differences between samples ################################
# load libraries
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(cowplot)
library(adegenet)

# read the csv files
divGenes <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/divGenes.csv")

#AA_0866 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0866.csv")
#AA_0900 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0900.csv")
#AA_0942 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0942.csv")
#IXa     <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_IXa.csv")


################################################################################
################################################################################

#### calculate the sum of each variant type per sample per Chromosome
ChrSums <- divGenes %>% group_by(Chr) %>% summarise(SNPsum866   = sum(SNPs.866, na.rm = T),
                                                    SNPsum942   = sum(SNPs.942, na.rm = T),
                                                    SNPsum900   = sum(SNPs.900, na.rm = T),
                                                    SNPsumIXa   = sum(SNPs.IXa, na.rm = T),
                                                    AmSNPsum866 = sum(AmSNPs.866, na.rm = T),
                                                    AmSNPsum942 = sum(AmSNPs.942, na.rm = T),
                                                    AmSNPsum900 = sum(AmSNPs.900, na.rm = T),
                                                    AmSNPsumIXa = sum(AmSNPs.IXa, na.rm = T),
                                                    INSsum866   = sum(INS.866, na.rm = T),
                                                    INSsum942   = sum(INS.942, na.rm = T),
                                                    INSsum900   = sum(INS.900, na.rm = T),
                                                    INSsumIXa   = sum(INS.IXa, na.rm = T),
                                                    DELsum866   = sum(DEL.866, na.rm = T),
                                                    DELsum942   = sum(DEL.942, na.rm = T),
                                                    DELsum900   = sum(DEL.900, na.rm = T),
                                                    DELsumIXa   = sum(DEL.IXa, na.rm = T))

################################################################################
################################################################################

#### 
# pivot longer each Variant and the corresponding sum
SNPs   <- ChrSums %>% select(1:5)      %>% pivot_longer(cols = c(SNPsum866,SNPsum942,SNPsum900,SNPsumIXa),         names_to = "SNPs",   values_to = "SNP")
AmSNPs <- ChrSums %>% select(1, 6:9)   %>% pivot_longer(cols = c(AmSNPsum866,AmSNPsum942,AmSNPsum900,AmSNPsumIXa), names_to = "AmSNPs", values_to = "AmSNP")
INSs   <- ChrSums %>% select(1, 10:13) %>% pivot_longer(cols = c(INSsum866,INSsum942,INSsum900,INSsumIXa),         names_to = "INSs",    values_to = "INS")
DELs   <- ChrSums %>% select(1, 14:17) %>% pivot_longer(cols = c(DELsum866,DELsum942,DELsum900,DELsumIXa),         names_to = "DELs",    values_to = "DEL")

# bind the 
pivChrSums <- cbind(SNPs, AmSNPs, INSs, DELs)
pivChrSums <- pivChrSums[!duplicated(as.list(pivChrSums))]

pivChrSums$SNPs <- gsub(pattern = "SNPsum", replacement = "", x = pivChrSums$SNPs)
colnames(pivChrSums)[colnames(pivChrSums)%in%"SNPs"] <- "Sample"

# deselect duplicated sample columns
pivChrSums <- pivChrSums %>% select(-c(4,6,8))

# pivot longer the Variants
pivChrSums <- pivChrSums %>% pivot_longer(cols = c(3:6), names_to = "Variant", values_to = "n")

# add Host Origin column
pivChrSums <- pivChrSums %>% mutate(Origin = case_when(Sample == '866' ~ 'W',
                                                       Sample == '942' ~ 'W',
                                                       Sample == '900' ~ 'E',
                                                       Sample == 'IXa' ~ 'E'))

################################################################################
################################################################################

# pure stacked per Chromosome, Sample fill, Variant wrap
a <- pivChrSums %>% 
  ggplot(aes(x = Chr, y = n,  fill = Sample)) +
  geom_col() +
  facet_wrap(~Variant, nrow = 1) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Variant Calls per Chromosome versus UGA55 RefSeq")
a # VCperChrSample.jpg


# per Chromosome, Variant fill, Sample wrap
b <- pivChrSums %>% 
  ggplot(aes(x = Chr, y = n,  fill = Variant)) +
  geom_col() +
  facet_wrap(~Sample, nrow = 1) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_brewer(palette = "Dark2", direction = -1) +
  ggtitle("Variant Calls per Chromosome versus UGA55 RefSeq")
b # VCperSampleChr.jpg

# per sample, Variant wrap
c <- pivChrSums %>% 
  ggplot(aes(x = Sample, y = n,  fill = Sample)) +
  geom_col() +
  facet_wrap(~Variant, nrow = 1) +
  scale_fill_brewer(palette = "Dark2", direction = -1) +
  ggtitle("Variant Calls per Chromosome versus UGA55 RefSeq")
c # VCperSampleVar.jpg

# all Variants per sample
d <- pivChrSums %>% 
  ggplot(aes(x = Sample, y = n,  fill = Variant)) +
  geom_col() +
  scale_fill_brewer(palette = "Dark2", direction = -1) +
  ggtitle("Variant Calls per Chromosome versus UGA55 RefSeq")
d # VCperSample.jpg

################################################################################
################################################################################
save_plot(plot = a, filename = '/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/figures/VCperChrSample.jpg', base_height = 5, base_width = 7.5) 
save_plot(plot = b, filename = '/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/figures/VCperSampleChr.jpg', base_height = 5, base_width = 7.5) 
save_plot(plot = c, filename = '/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/figures/VCperSampleVar.jpg', base_height = 5, base_width = 7.5) 
save_plot(plot = d, filename = '/Users/finnlo/Documents/Github/Crypto/Genome_Analysis/figures/VCperSample.jpg'   , base_height = 5, base_width = 7.5) 
################################################################################
################################################################################




