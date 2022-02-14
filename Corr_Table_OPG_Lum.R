library(tidyverse)
library(visdat)
library(data.table)
library(stringr)


ABI_Best_thSC     <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_input/Cryptosporidium/Crypto_Standard_Curve.csv")
ABI_Best_thSC     <-  filter(ABI_Best_thSC, Ct_mean > 0)
linear_model0     <- lm(log2(Amount_Oocysts) ~ Ct_mean, data = ABI_Best_thSC)
Oocyst_Predict    <- 2^predict(linear_model0, newdata = Crypto_qPCR)
Crypto_qPCR <- data.frame(Crypto_qPCR, Oocyst_Predict)
Crypto_qPCR <- Crypto_qPCR %>% mutate(Oocyst_Predict = replace(Oocyst_Predict, Oocyst_Predict == "4292821751815.77", "0"))
Crypto_qPCR$Oocyst_Predict <- as.integer(Crypto_qPCR$Oocyst_Predict)