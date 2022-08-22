### VCF Distance Mapping #######################################################

# Load libraries
library(tidyverse)
library(data.table)
library(adegenet)

# Load data
AA_0866 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0866.csv") %>% dplyr::select(-X)
AA_0900 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0900.csv") %>% dplyr::select(-X)
AA_0942 <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_AA_0942.csv") %>% dplyr::select(-X)
IXa     <- read.csv("https://raw.githubusercontent.com/tlobnow/Crypto/main/Genome_Analysis/Products/Ctyz_IXa.csv")     %>% dplyr::select(-X)

# Join data
df1 <- left_join(AA_0866, AA_0900, suffix = c(".866", ".900"), by = 'START')
df2 <- left_join(AA_0942, IXa,     suffix = c(".942", ".IXa"), by = 'START')
df <- full_join(df1, df2, by = 'START')

GeneID          <- df %>% dplyr::select('START', 'GeneID.866',  'GeneID.900',  'GeneID.942',  'GeneID.IXa')
Product         <- df %>% dplyr::select('START', 'Product.866', 'Product.900', 'Product.942', 'Product.IXa')
WIDTH           <- df %>% dplyr::select('START', 'WIDTH.866', 'WIDTH.900', 'WIDTH.942', 'WIDTH.IXa')
AminoAcidChange <- df %>% dplyr::select('START', 'AminoAcidChange.866', 'AminoAcidChange.900', 'AminoAcidChange.942', 'AminoAcidChange.IXa')
REF.ALT         <- df %>% dplyr::select('START', 'REF.ALT.866', 'REF.ALT.900', 'REF.ALT.942', 'REF.ALT.IXa')




                               #'Product.866', 'Product.900', 'Product.942', 'Product.IXa',
                               #'START', 
                               #'WIDTH.866', 'WIDTH.900', 'WIDTH.942', 'WIDTH.IXa',
                               #'REF.ALT.866', 'REF.ALT.900', 'REF.ALT.942', 'REF.ALT.IXa',
                               #'QUAL.866', 'QUAL.900', 'QUAL.942', 'QUAL.IXa',
                               #'AminoAcidChange.866', 'AminoAcidChange.900', 'AminoAcidChange.942', 'AminoAcidChange.IXa')




