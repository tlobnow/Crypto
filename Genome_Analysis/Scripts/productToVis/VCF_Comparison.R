####### Look at the differences between samples ################################
# load libraries
library(tidyverse)
library(data.table)
library(RColorBrewer)

# read the csv files
AA_0866 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0866.csv")
AA_0900 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0900.csv")
AA_0942 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0942.csv")
IXa     <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_IXa.csv")
divGenes <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/divGenes.csv")


################################################################################
################################################################################

divGenes2 %>% 
  ggplot(aes(x = Chr)) +
  geom_col(aes(y = Variants.942, fill = '942')) +
  geom_col(aes(y = Variants.IXa, fill = 'IXa')) +
  geom_col(aes(y = Variants.866, fill = '866')) +
  geom_col(aes(y = Variants.900, fill = '900')) +
  facet_wrap(~VariantType, nrow = 1) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Variant Calls per Chromosome versus UGA55 RefSeq") +
  ylab("n")

################################################################################
################################################################################












