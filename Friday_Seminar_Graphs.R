library(tidyverse)  #install.packages('tidyverse')
library(data.table) #install.packages('data.table') # for the setnames() function
library(tibble)     #install.packages('tibble')     # to use tibbles
library(readr)      #install.packages('readr')      # to read tsv files
library(xlsx)       #install.packages('readxl')     # to read xlsx files
library(knitr)      #install.packages('knitr')      # for making tables in markdown
library(visdat)     #install.packages('visdat')     # for looking at missing data
library(scales)     #install.packages('scales')     # for the percent() function
library(tidyr)
library(RColorBrewer) #display.brewer.all()
library(leaflet)
library(htmltools)

SOTA <- read.csv('https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv')
PCR_Res <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/PCR_Results.csv")
Amp_Res <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Amplification_Success.csv")
Crypto_Detection <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/Crypto_Detection.csv") 
Crypto_Detection_tested <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/Crypto_Detection.csv")  %>% filter(ILWE_Crypto_Ct >= 0, !is.na(Latitude), !is.na(Longitude))
Crypto_Positive <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/Crypto_Detection.csv")  %>% filter(ILWE_Crypto_Ct > 0, !is.na(Latitude), !is.na(Longitude))
Crypto_Detection_21 <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/Crypto_Detection.csv")  %>% filter(Year == 2021, ILWE_Crypto_Ct > 0, !is.na(Latitude), !is.na(Longitude))
HMHZ <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/HMHZ_Samples_Locations.csv", na.strings=c(""," ","NA")) %>% filter(!is.na(Longitude)) %>% select(-X)
qual_div <- read.csv('https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Quality_Diversity_Assessment.csv')
pairs_designed <- read.csv('https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Primer_Pairs_designed.csv', na.strings = c('', ' ', '-'))



#SOTA <- read.csv('https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv')
SOTA <- SOTA[SOTA$Mouse_ID %like% "AA_", ] # we want Brandenburg HMHZ samples only

# we want to compare the catching rate overall to the mice that were Crypto-positive
All_Samples <- SOTA %>% group_by(Year) %>% count()
Pos_Samples <- SOTA %>% filter(ILWE_Crypto_Ct > 0, Year >= 2016) %>% group_by(Year) %>% count()

# let's combine both and calculate the prevalence
Samples_Yr <- full_join(Pos_Samples, All_Samples, by = "Year") %>% 
  mutate(Prevalence = (n.x / n.y)) %>% 
  filter(!is.na(Prevalence))

# let's visualize that with a bar plot (geom_col is like geom_bar, but can take both x AND y)

fig1 <- Samples_Yr %>%
  ggplot(aes(x = Year, label = Prevalence)) +
  geom_col(aes(y = n.y, fill = "blue")) +
  geom_col(aes(y = n.x, fill = "red")) +
  geom_text(aes(label = percent(Prevalence),
                y = (n.x / n.y)), nudge_y = 9) +
  labs(y = "Samples [n]") +
  theme(legend.position = "none") +
  ggtitle("Cryptosporidium spp. prevalence in the HMHZ since 2016")                                       
#fig1

#PCR_Res <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/PCR_Results.csv")

PCR_Res <- PCR_Res %>% group_by(Locus..Gene.Family) %>% mutate(n_Primers = sum(Successful.PCR.Amplification))
PCR_Res_lim <- PCR_Res %>% distinct(Locus..Gene.Family, n_Primers)
PCR_Res_lim <- PCR_Res_lim %>% mutate(sequenced = case_when(Locus..Gene.Family %in% c('GP60', 'CP56', 'GST', 'MEDLE', 'SKSR') ~ 1))

#PCR_Res_lim %>%
#  ggplot(aes(x = Locus..Gene.Family)) +
#  geom_col(aes(y = n_Primers, fill = 's')) +
#  geom_col(aes(y = sequenced, fill = 'a')) +
#  coord_polar()+
#  theme(legend.position = "none") +
#  theme(axis.text.x = element_text(angle = 0))



#Amp_Res <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Amplification_Success.csv")
Amp_Res$Primer.Pair[ Amp_Res$Locus..Gene.Family %in% "HSP90"] <- 'HSP90'
Amp_Res$Successful.PCR.Amplification[ Amp_Res$Locus..Gene.Family %in% "HSP90"] <- F

Amp_Res <- Amp_Res %>% mutate(Sanger_num = case_when(Sanger.Sequenced == T ~ 1,
                                                     Sanger.Sequenced == F ~ 0),
                              Amp_num = case_when(Successful.PCR.Amplification == T ~1,
                                                  Successful.PCR.Amplification == F ~0),
                              Primer.pairs.ordered = case_when(Locus..Gene.Family == 'HSP90' ~ 0,
                                                               Locus..Gene.Family != 'HSP90' ~ 1))

fig2 <- Amp_Res %>%
  ggplot(aes(x = Target.Gene.ID..Reference.Sequence.UGA55.)) +
  geom_col(aes(y = Primer.pairs.ordered)) +
  geom_col(aes(y = Amp_num), col = 'coral2', size = 1) +
  geom_col(aes(y = Sanger_num), fill = 'cyan4') +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "bottom") +
  ylab('Primer Pair') +
  xlab('Target Gene ID (RefSeq UGA55)') +
  #ggtitle('Sequencing Success of ordered primer pairs') +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 6, colour = 'black'))+
  theme(panel.grid = element_blank())+
  facet_wrap(~Locus..Gene.Family, scales = 'free_x', nrow = 1, as.table = T)
#fig2


#qual_div <- read.csv('https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Quality_Diversity_Assessment.csv')
colnames(qual_div)[colnames(qual_div) %in% 'Pairwise.Identity.with.C.parvum.for.product.region...PI..nt.'] <- 'Pairwise_Identity'
colnames(qual_div)[colnames(qual_div) %in% 'Diversity.High..80..PI.Med...80..PI.Low...95..PI'] <- 'Diversity'

qual_div$Pairwise_Identity <- qual_div$Pairwise_Identity %>% str_replace_all("%", "")
qual_div$Pairwise_Identity <- as.numeric(qual_div$Pairwise_Identity)

qual_div <- qual_div %>% mutate(Diversity = case_when(Pairwise_Identity < 85 ~ 'PI <80% (High)',
                                                      Pairwise_Identity >= 85 & Pairwise_Identity < 95 ~ 'PI <95% (Medium)',
                                                      Pairwise_Identity > 95 ~ 'PI >95% (Low)'),
                                PI_dbl = 100 - Pairwise_Identity)
fig3 <- qual_div %>%
  ggplot(aes(x = Primer.Pair)) +
  geom_col(aes(y = PI_dbl, fill = Diversity)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_brewer(palette = 'Dark2', direction = 1) +
  ylab('Diversity [%]') +
  xlab('Primer Pair') +
  #ggtitle('Diversity across ordered primer pairs (Ct vs Cp)') +
  #facet_wrap(~Gene...Gene.Family, scales = 'free_x', nrow = 1, as.table = T)
  facet_wrap(~reorder(Gene...Gene.Family, PI_dbl), scales = 'free_x', nrow = 1) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5))
#fig3

#pairs_designed <- read.csv('https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/Primer_Pairs_designed.csv', na.strings = c('', ' ', '-'))
colnames(pairs_designed)[colnames(pairs_designed) %in% 'Nucleotide.Pairwise.Identity....'] <- 'Pairwise_Identity'
pairs_designed$Passes <- rowSums(pairs_designed[51:59] == T)
pairs_designed$Fails <- rowSums(pairs_designed[51:59] == F)
pairs_designed <- pairs_designed %>% mutate(Diversity = case_when(Pairwise_Identity < 80 ~ 'PI <80% (High)',
                                                                  Pairwise_Identity >= 80 & Pairwise_Identity < 95 ~ 'PI <95% (Medium)',
                                                                  Pairwise_Identity > 95 ~ 'PI >95% (Low)'),
                                            PI_dbl = 100 - Pairwise_Identity,
                                            # Quality = Primer Pair Quality in grades (1 == excellent, 4 == poor)
                                            Quality = case_when(Passes == 4 ~ 4,
                                                                Passes == 5 ~ 3,
                                                                Passes == 6 ~ 2,
                                                                Passes == 7 ~ 1),
                                            # Hard Fails = Primer Quality Assessment reveals fails in Self-Annealing or Hairpin-Formation
                                            # leads to Quality repercussions (one grade worse than for 'light' fails)
                                            Hard_Fails = case_when(Self.Annealing.Pass == F ~ T,
                                                                   Hairpin.formation..Pass == F ~ T,
                                                                   Self.Annealing.Pass & Hairpin.formation..Pass == F ~ T,
                                                                   Self.Annealing.Pass == T ~ F,
                                                                   Hairpin.formation..Pass == T ~ F,
                                                                   Self.Annealing.Pass & Hairpin.formation..Pass == T ~ F),
                                            # adjusted Quality measurement
                                            Quality_adj = case_when(Hard_Fails == T & Quality == 1 ~ Quality + 1,
                                                                    Hard_Fails == T & Quality == 2 ~ Quality + 1,
                                                                    Hard_Fails == T & Quality == 3 ~ Quality,
                                                                    Hard_Fails == T & Quality == 4 ~ Quality,
                                                                    Hard_Fails == F ~ Quality))
pairs_designed <- pairs_designed %>% group_by(Primer.Pair) %>% mutate(Quality_mean = mean(Quality_adj),
                                                                      #Quality_chr = case_when(Quality_mean == 3.5 ~ 'terrible',
                                                                      #                        Quality_mean == 3.0 ~ 'poor',
                                                                      #                        Quality_mean == 2.5 | Quality_mean == 2.0 ~ 'medium',
                                                                      #                        Quality_mean == 1.5 ~ 'good',
                                                                      #                        Quality_mean == 1.0 ~ 'excellent'),
                                                                      Quality_chr = case_when(Quality_mean == 3.5 ~ 'poor',
                                                                                              Quality_mean == 3.0 ~ 'poor',
                                                                                              Quality_mean == 2.5 | Quality_mean == 2.0 ~ 'medium',
                                                                                              Quality_mean == 1.5 ~ 'medium',
                                                                                              Quality_mean == 1.0 ~ 'excellent'))

pairs_designed <- full_join(pairs_designed, Amp_Res)
pairs_designed <- pairs_designed %>% mutate(ordered = case_when(Primer.pairs.ordered == 1  ~ T,
                                                                Primer.pairs.ordered == 0  ~ F,
                                                                is.na(Primer.pairs.ordered) ~ F))

fig4 <- pairs_designed %>%
  group_by(Primer.Pair) %>%
  ggplot(aes(x = Primer.Pair)) +
  geom_hline(yintercept = 25, col = 'grey') +
  geom_col(aes(y = PI_dbl/2, fill = Quality_chr)) +
  #geom_label(aes(y = ordered)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_brewer(palette = 'Dark2', direction = 1) +
  ylab('Diversity [%]') +
  xlab('Primer Pair') +
  #ggtitle('Diversity across ordered primer pairs (Ct vs Cp)') +
  facet_wrap(~reorder(Gene...Gene.Family, PI_dbl), scales = 'free_x', nrow = 2) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        plot.title = element_text(hjust = 0.5))
#fig4


################################################################################
################################################################################
################################################################################
# Load Palette 
r <- c(0, 64, 128, 179, 217, 255)
g <- c(0, 12, 25, 25, 12,  0)
b <- c(255, 249, 243, 191,  95,   0)

myPal <- function (n, name = c("myPal.colors")) 
{
  myPal.colors = rgb(r,g,b,maxColorValue = 255)
  name = match.arg(name)
  orig = eval(parse(text = name))
  rgb = t(col2rgb(orig))
  temp = matrix(NA, ncol = 3, nrow = n)
  x = seq(0, 1, , length(orig))
  xg = seq(0, 1, , n)
  for (k in 1:3) {
    hold = spline(x, rgb[, k], n = n)$y
    hold[hold < 0] = 0
    hold[hold > 255] = 255
    temp[, k] = round(hold)
  }
  palette = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
  palette
}

# Load Palette 
s <- c(255, 147, 111)
t <- c(217, 196, 168)
u <- c(102, 125, 220)

YlGnBl <- function (n, name = c("YlGnBl.colors")) 
{
  YlGnBl.colors = rgb(s,t,u,maxColorValue = 255)
  name = match.arg(name)
  orig = eval(parse(text = name))
  rgb = t(col2rgb(orig))
  temp = matrix(NA, ncol = 3, nrow = n)
  x = seq(0, 1, , length(orig))
  xg = seq(0, 1, , n)
  for (k in 1:3) {
    hold = spline(x, rgb[, k], n = n)$y
    hold[hold < 0] = 0
    hold[hold > 255] = 255
    temp[, k] = round(hold)
  }
  palette = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
  palette
}




## HYBRID INDEX MAP ------------------------------------------------------------

library(RColorBrewer) #display.brewer.all()
library(leaflet)
library(htmltools)

SOTA <- read.csv('https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv') %>% filter(!is.na(HI))

High_Infection_Samples <- SOTA %>% 
  select(Mouse_ID, ILWE_Crypto_Ct, Oocyst_Predict_Crypto, Year, Latitude, Longitude, HI, Sex) %>% 
  filter(ILWE_Crypto_Ct > 0, Year >= 2016) %>% 
  arrange(ILWE_Crypto_Ct) %>% 
  head(30)

SOTA <- SOTA %>% 
  mutate(Eim_Species = ifelse(eimeriaSpecies == "E_falciformis", "E_falciformis",
                              ifelse(eimeriaSpecies == "E_ferrisi", "E_ferrisi",
                                     ifelse(eimeriaSpecies == "Eimeria_alorani", "Eimeria_alorani",
                                            ifelse(eimeriaSpecies == "Eimeria_apionodes", "Eimeria_apionodes",
                                                   ifelse(eimeriaSpecies == "Eimeria_falciformis", "Eimeria_falciformis",
                                                          ifelse(eimeriaSpecies == "Eimeria_sp_Apodemus", "Eimeria_sp_Apodemus",
                                                                 ifelse(eimeriaSpecies == "Eimeria_vermiformis", "Eimeria_vermiformis",
                                                                        ifelse(eimeriaSpecies == "Negative", "Negative",
                                                                               ifelse(NA))))))))),
         Eimeria_Positive = case_when(Eim_Species != "Negative" | Ct.Eimeria  > 0 | OPG > 0 ~ T,
                                      Eim_Species == "Negative" | Ct.Eimeria == 0 | OPG == 0 ~ F))


# Cutting the data into slices helps to make nice legends!
SOTA$HI <- as.numeric(SOTA$HI)
SOTA$HI_Level <-  cut(SOTA$HI, c(0, 0.001, 0.250, 0.500, 0.750, 0.999, 1), include.lowest = T , 
                      labels = c('HI = 0.00', 'HI < 0.25', 'HI < 0.50', 'HI < 0.75', 'HI < 1.00', 'HI = 1.00'))
HMHZ$HI_Level <-  cut(HMHZ$HI, c(0, 0.001, 0.250, 0.500, 0.750, 0.999, 1), include.lowest = T , 
                          labels = c('HI = 0.00', 'HI < 0.25', 'HI < 0.50', 'HI < 0.75', 'HI < 1.00', 'HI = 1.00'))

# We also want to look at samples that sent into sequencing (was interesting for my Bachelor thesis)
Sequenced <- Crypto_Detection %>%
  mutate(seq = Mouse_ID %in% c("AA_0144", "AA_0325", "AA_0689", "AA_0209", "AA_0282", "AA_0793", "AA_0667", "AA_0805", "AA_0900",
                               "AA_0523", "AA_0534", "AA_0537", "AA_0545", "AA_0546", "AA_0553", "AA_0554", "AA_0555", "AA_0557",
                               "AA_0559", "AA_0571", "AA_0578", "AA_0580", "AA_0585", "AA_0589",  "AA_0601", "AA_0660",  "AA_0666",
                               "AA_0667", "AA_0669", "AA_0679")) %>%
  filter(seq == T) # only retain samples that fulfill this condition

Illumina <- SOTA %>%
  mutate(illumina = Mouse_ID %in% c("AA_0900", "AA_0866", "AA_0942")) %>%
  filter(illumina == T)

HMHZ <- HMHZ %>% mutate(seq_by = case_when(Mouse_ID %in% c("AA_0144", "AA_0325", "AA_0689", "AA_0209", "AA_0282", "AA_0793", "AA_0805", "AA_0900", "AA_0523", "AA_0667") ~ 'Finn (2021)',
                                           Mouse_ID %in% c("AA_0534", "AA_0537", "AA_0545", "AA_0546", "AA_0553", "AA_0554", "AA_0555", "AA_0557", "AA_0559", "AA_0571",
                                                           "AA_0578", "AA_0580", "AA_0585", "AA_0589", "AA_0601", "AA_0660", "AA_0666", "AA_0669", "AA_0679") ~ 'Yasmin (2019)',
                                           !is.na(Location) ~ 'Kvac (2013)'))

HMHZ <- HMHZ %>% 
  select(Mouse_ID, GP60_Subtype, Actin_Subtype, COWP_Subtype, TRAP_C1_Subtype, 
         MSC6_Subtype, GST_Subtype, SKSR_Subtype, MEDLE_Subtype, CP56_Subtype, Location, seq_by, Longitude, Latitude, HI, Country)
HMHZ$summed <- rowSums(!is.na(HMHZ[2:10]))

Finn    <- HMHZ %>% filter(seq_by == "Finn (2021)")
Yasmin  <- HMHZ %>% filter(seq_by == "Yasmin (2019)")
Kvac    <- HMHZ %>% filter(seq_by == "Kvac (2013)")

data_col_HI        = colorFactor(myPal(6), SOTA$HI) # for our samples, uncut
data_col_HI_HMHZ   = colorFactor(myPal(6), HMHZ$HI) # for our samples, uncut
data_col_HI_Level  = colorFactor(myPal(6), HMHZ$HI_Level) # for our legend palette, cut
data_col_HI_Level  = colorFactor(myPal(6), SOTA$HI_Level) # for our legend palette, cut
data_col_seq       = colorFactor(YlGnBl(3), HMHZ$seq_by, reverse = F)


map <- HMHZ %>% leaflet() %>% addProviderTiles("CartoDB") %>% setView(lat = 52.520007, lng =13.404954, zoom = 6)
map1 <- map %>%
  addPolylines(lat = c(55.0000, 53.6000, 53.51885, 52.8875  , 52.6053, 51.8978, 45.0000), 
               lng = c(10.0000, 11.4563, 12.4464,13.8119 , 13.8756, 13.8103, 10.0000), 
               color = "purple", 
               weight = 55, 
               opacity = 0.1) %>%
  #addLabelOnlyMarkers(data = Finn, 
  #                    label = ~as.character(Mouse_ID),
  #                    labelOptions = labelOptions(noHide = T, direction = 'top', textOnly = F, opacity = 0.5),
  #                    group = 'Labels_F') %>%
  #addLabelOnlyMarkers(data = Yasmin, 
  #                    label = ~as.character(Mouse_ID),
  #                    labelOptions = labelOptions(noHide = T, direction = 'bottom', textOnly = F, opacity = 0.5),
  #                    group = 'Labels_Y') %>%
  #addLabelOnlyMarkers(data = Kvac, 
  #                    label = ~as.character(Mouse_ID),
  #                    labelOptions = labelOptions(noHide = T, direction = 'bottom', textOnly = F, opacity = 0.5),
  #                    group = 'Labels_K') %>%
  #addCircleMarkers(data = HMHZ, 
  #                 color = "black",
  #                 label = ~htmlEscape(Mouse_ID),
  #                 popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
  #                                "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
  #                                "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
  #                                sep=" "),
  #                 opacity = 5,
  #                 radius = 4,
  #                 group = "Underlay") %>%
  addCircleMarkers(data = SOTA,
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct Mean:<b>",      as.character(round(ILWE_Crypto_Ct, digits = 2)),"<br>",
                                  "<b>Oocyst Prediction:<b>", as.character(Oocyst_Predict_Crypto), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 1,
                   radius = 3,
                   group = "Samples (total)") %>%
  addCircleMarkers(data = Crypto_Detection_tested,
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct Mean:<b>",      as.character(round(ILWE_Crypto_Ct, digits = 2)),"<br>",
                                  "<b>Oocyst Prediction:<b>", as.character(Oocyst_Predict_Crypto), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 1,
                   radius = 3,
                   group = "Samples (Crypto-tested)") %>%
  addCircleMarkers(data = Crypto_Positive,
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct Mean:<b>",      as.character(round(ILWE_Crypto_Ct, digits = 2)),"<br>",
                                  "<b>Oocyst Prediction:<b>", as.character(Oocyst_Predict_Crypto), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 1,
                   radius = 3,
                   group = "Samples (Crypto-positive)") %>%
  addCircleMarkers(data = Sequenced, 
                   col = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct Mean:<b>",      as.character(round(ILWE_Crypto_Ct, digits = 2)),"<br>",
                                  "<b>Oocyst Prediction:<b>", as.character(Oocyst_Predict_Crypto), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 5,
                   radius = 3,
                   group = "Sequenced") %>%
  addCircleMarkers(data = High_Infection_Samples, 
                   color = ~data_col_HI(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  "<b>Year:<b>",    as.character(Year),"<br>",
                                  "<b>Ct Mean:<b>",      as.character(round(ILWE_Crypto_Ct, digits = 2)),"<br>",
                                  "<b>Oocyst Prediction:<b>", as.character(Oocyst_Predict_Crypto), "<br>",
                                  "<b>Sex:<b>", Sex, "<br>",
                                  sep=" "),
                   opacity = 3,
                   radius = 3,
                   group = "High_Infection_Samples") %>%
  addCircleMarkers(data = HMHZ, 
                   color = ~data_col_HI_HMHZ(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  sep=" "),
                   opacity = 3,
                   radius = 3,
                   group = "HMHZ") %>%
  addCircleMarkers(data = HMHZ, 
                   color = ~data_col_seq(seq_by),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  sep=" "),
                   opacity = 5,
                   radius = 3,
                   group = "sequenced") %>%
  addCircleMarkers(data = Illumina,
                   color = 'yellowgreen',
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>HI:<b>",      as.character(round(HI, digits = 2)), "<br>",
                                  sep=" "),
                   opacity = 5,
                   radius = 3,
                   group = 'Illumina') %>%
  addLegend("bottomleft", 
            pal = data_col_HI_Level, 
            title = "Hybrid Index",
            values = SOTA$HI_Level, 
            group = c('HI = 0.00', 'HI < 0.25', 'HI < 0.50', 'HI < 0.75', 'HI < 1.00', 'HI = 1.00'),
            opacity = 1) %>%
  addLegend("bottomright",
            pal = data_col_seq,
            title = 'Data Origin',
            values = HMHZ$seq_by,
            group = c('BR-HZ (2021)', 'BR-HZ (2019', 'CZ/BAV-HZ (2013'),
            opacity = 1) %>%
  addLayersControl(baseGroups = c("Samples (total)", "Samples (Crypto-tested)", "Samples (Crypto-positive)", "High_Infection_Samples", "Sequenced", 'HMHZ', 'sequenced'), #, "Eim_Infected"), 
                   overlayGroups = c(#"Underlay", 'Labels_F', 'Labels_Y', 'Labels_K', 
                                     'Illumina'),
                   options = layersControlOptions(collapsed = T))
#map1

HMHZ <- read.csv("https://raw.githubusercontent.com/tlobnow/Cryptosporidium-BSc/Main-Branch/Analysis/HMHZ_Samples_Locations.csv", na.strings=c(""," ","NA")) %>% filter(!is.na(Longitude)) %>% select(-X)

Yellow <- HMHZ %>% filter(GP60_Ssp        == "IXb.1"|
                            MSC6_Subtype    == "MS1"|
                            GST_Subtype     == "G1"|
                            CP56_Subtype    == "CP1" |
                            SKSR_Subtype    == "S1")

Light_green <- HMHZ %>% filter(GP60_Ssp    == "IXa"|
                                 Actin_Subtype   == "A1"|
                                 COWP_Subtype    == "C1"|
                                 TRAP_C1_Subtype == "T1")

Orange <- HMHZ %>% filter(GP60_Ssp    == "IXb.2"|
                            Actin_Subtype   == "A2"|
                            COWP_Subtype    == "C2"|
                            TRAP_C1_Subtype == "T2"|
                            MSC6_Subtype    == "MS2"|
                            GST_Subtype     == "G2"|
                            CP56_Subtype    == "CP2"|
                            SKSR_Subtype    == "S2")


Dark_Orange <- HMHZ %>% filter(GP60_Ssp == "IXb.3")


map <- Clades %>%
  leaflet() %>%
  addProviderTiles("CartoDB") %>%
  setView(lat = 51.2 , lng =13.679507, zoom = 6)

map2 <- map %>%
  addPolylines(lat = c(55.0000, 53.6000, 53.51885, 52.8875, 52.6053, 51.8978, 50.08506775 , 48.13659738768638, 45.0000), lng = c(10.0000, 11.4563, 12.4464,13.8119 , 13.8756, 13.8103, 12.5988909, 11.582597486663213, 11.50000), color = "purple", weight = 75, opacity = 0.1)  %>%
  addCircleMarkers(data = HMHZ, 
                   col = ~data_col_HI_HMHZ(HI),
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  sep=" "),
                   opacity = 0.1,
                   radius = 10,
                   group = "Colored by HI") %>%
  addCircleMarkers(data = Yellow,
                   col = 'gold',
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>GP60 Subtype:<b>",      as.character(GP60_Subtype), "<br>",
                                  "<b>MSC6-7 Subtype:<b>",      as.character(MSC6_Subtype), "<br>",
                                  "<b>GST Subtype:<b>",      as.character(GST_Subtype), "<br>",
                                  "<b>CP56 Subtype:<b>",      as.character(CP56_Subtype), "<br>",
                                  "<b>SKSR Subtype:<b>",      as.character(SKSR_Subtype), "<br>",
                                  "<b>MEDLE Subtype:<b>",      as.character(MEDLE_Subtype), "<br>",
                                  "<b>Actin Subtype:<b>",      as.character(Actin_Subtype), "<br>",
                                  "<b>TRAP-C1 Subtype:<b>",    as.character(TRAP_C1_Subtype), "<br>",
                                  "<b>COWP Subtype:<b>",      as.character(COWP_Subtype), "<br>",
                                  sep=" "),
                   opacity = 1,
                   radius = 3,
                   group = '1 .. East BR') %>%
  addCircleMarkers(data = Light_green,
                   col = 'yellowgreen',
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>GP60 Subtype:<b>",      as.character(GP60_Subtype), "<br>",
                                  "<b>MSC6-7 Subtype:<b>",      as.character(MSC6_Subtype), "<br>",
                                  "<b>GST Subtype:<b>",      as.character(GST_Subtype), "<br>",
                                  "<b>CP56 Subtype:<b>",      as.character(CP56_Subtype), "<br>",
                                  "<b>SKSR Subtype:<b>",      as.character(SKSR_Subtype), "<br>",
                                  "<b>MEDLE Subtype:<b>",      as.character(MEDLE_Subtype), "<br>",
                                  "<b>Actin Subtype:<b>",      as.character(Actin_Subtype), "<br>",
                                  "<b>TRAP-C1 Subtype:<b>",    as.character(TRAP_C1_Subtype), "<br>",
                                  "<b>COWP Subtype:<b>",      as.character(COWP_Subtype), "<br>",
                                  sep=" "),
                   opacity = 1,
                   radius = 3,
                   group = '1 .. East BR/CZ') %>%
  addCircleMarkers(data = Orange,
                   col = 'orange',
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>GP60 Subtype:<b>",      as.character(GP60_Subtype), "<br>",
                                  "<b>MSC6-7 Subtype:<b>",      as.character(MSC6_Subtype), "<br>",
                                  "<b>GST Subtype:<b>",      as.character(GST_Subtype), "<br>",
                                  "<b>CP56 Subtype:<b>",      as.character(CP56_Subtype), "<br>",
                                  "<b>SKSR Subtype:<b>",      as.character(SKSR_Subtype), "<br>",
                                  "<b>MEDLE Subtype:<b>",      as.character(MEDLE_Subtype), "<br>",
                                  "<b>Actin Subtype:<b>",      as.character(Actin_Subtype), "<br>",
                                  "<b>TRAP-C1 Subtype:<b>",    as.character(TRAP_C1_Subtype), "<br>",
                                  "<b>COWP Subtype:<b>",      as.character(COWP_Subtype), "<br>",
                                  sep=" "),
                   opacity = 1,
                   radius = 3,
                   group = '2 .. West BR/CZ') %>%
  addCircleMarkers(data = Dark_Orange,
                   col = 'brown',
                   label = ~htmlEscape(Mouse_ID),
                   popup = ~paste("<b>Mouse_ID:<b>",as.character(Mouse_ID), "<br>",
                                  "<b>HI:<b>",      as.character(HI), "<br>",
                                  "<b>Location:<b>", as.character(Latitude), "<b>,<b>", as.character(Longitude), "<br>",
                                  "<b>GP60 Subtype:<b>",      as.character(GP60_Subtype), "<br>",
                                  "<b>MSC6-7 Subtype:<b>",      as.character(MSC6_Subtype), "<br>",
                                  "<b>GST Subtype:<b>",      as.character(GST_Subtype), "<br>",
                                  "<b>CP56 Subtype:<b>",      as.character(CP56_Subtype), "<br>",
                                  "<b>SKSR Subtype:<b>",      as.character(SKSR_Subtype), "<br>",
                                  "<b>MEDLE Subtype:<b>",      as.character(MEDLE_Subtype), "<br>",
                                  "<b>Actin Subtype:<b>",      as.character(Actin_Subtype), "<br>",
                                  "<b>TRAP-C1 Subtype:<b>",    as.character(TRAP_C1_Subtype), "<br>",
                                  "<b>COWP Subtype:<b>",      as.character(COWP_Subtype), "<br>",
                                  sep=" "),
                   opacity = 1,
                   radius = 3,
                   group = '3 .. West BR') %>%
  addLegend("bottomleft", 
            pal = data_col_HI_Level, 
            title = "Hybrid Index",
            values = SOTA$HI_Level, 
            group = c('HI = 0.00', 'HI < 0.25', 'HI < 0.50', 'HI < 0.75', 'HI < 1.00', 'HI = 1.00'),
            opacity = 1) %>%
  addLegend("bottomright", 
            colors = c('yellow', 'yellowgreen', 'darkorange', 'brown'),
            #colors =c("#FFC125",  "#FFC125", "#8A4117", "#7D0552", "#571B7E"),
            labels = c('1 .. East BR',
                       '1 .. East BR/CZ',
                       '2 .. West BR/CZ', 
                       '3 .. West BR'),
            title = 'Clades',
            opacity = 1) %>%
  addLayersControl(baseGroups = "Colored by HI",
                   overlayGroups = c('1 .. East BR',
                                     '1 .. East BR/CZ',
                                     '2 .. West BR/CZ',
                                     '3 .. West BR'), 
                   options = layersControlOptions(collapsed = T))
#map2
