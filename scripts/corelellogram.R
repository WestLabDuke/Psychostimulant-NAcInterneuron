# Author: Melyssa Minto
# Corelellogram for IP and UF ATAC samples


# loading libraries -------------------------------------------------------


library(corrplot)
library(tidyverse)
library(dichromat)


# reading in data ---------------------------------------------------------


VarianceStabilizedCounts <- read_csv("../data/All_ATAC_Counts.csv")


# wrangling data ----------------------------------------------------------


# extracting matrix
mat = VarianceStabilizedCounts[, 3:ncol(VarianceStabilizedCounts)]

# extracting labels 
x = colnames(mat)
PV_SST = 
  case_when( 
  grepl("PV", x) ~ "PV",
  grepl("SST", x) ~"SST")

IP_UF = 
  case_when( 
  grepl("IP", x) ~ "IP",
  grepl("UF", x) ~"UF")

l = paste(PV_SST, IP_UF, sep="_")

# perform correlation
colnames(mat) = l
cor.mat = cor(mat, method = "spearman")


# plotting data -----------------------------------------------------------

# plot
corrplot(cor.mat,
         type = "upper",
         title = "\n\n\n",
         is.corr = FALSE,
         order = "hclust",
         col.lim = c(0.3,1), 
         col = colorRampPalette(c("white", "deepskyblue", "blue4"))(100),
         addCoef.col = TRUE,
         method = "shade",
         tl.col = "black",
         diag = TRUE)

png("../output/correlograms/corr_atac.png", bg = "transparent", res=300,units = "in", height = 10, width =10)
corrplot(cor.mat,
         type = "upper",
         title = "\n\n\n",
         is.corr = FALSE,
         order = "hclust",
         col.lim = c(0.3,1), 
         col = colorRampPalette(c("white", "deepskyblue", "blue4"))(100),
         addCoef.col = TRUE,
         method = "shade",
         tl.col = "black",
         diag = TRUE)

dev.off()

