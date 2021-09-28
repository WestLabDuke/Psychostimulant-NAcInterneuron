# Author: Melyssa Minto
# This script will create the RRHO plot for the Differential Signal between PV Pull down MECP2 KI mice and PV Pull down in response to AMPH


# Loading Libraries -------------------------------------------------------

library(tidyverse)
library(readxl)
library(RRHO2)

# Reading in Data ---------------------------------------------------------


MECP2WTvKI <- read_excel("/media/west-lab-share/David_Gallegos/PVseqPaperALL/Psychostimulant-NAcInterneuron/FinalTables/Table S10-NAc PV+ RNA in MeCP2WTvsKI.xlsx")
Delayed_AMPH_nucRNAseq <- read_csv("../output/PV3HrSalvAMPH_RNA/diffExpressionMatrix.csv")
Rapid_AMPH_nucRNAseq <- read_csv("../output/PV35minSalvAMPH_RNA/diffExpressionMatrix.csv")



# Preparing data for RRHO -------------------------------------------------
# > using signed -log10 FDRs for plot
mecp2_diffsignal = MECP2WTvKI %>% 
  mutate(sgnlogpadj = -log10(q_value)* sign(log2fold_change)) %>% 
  dplyr::rename( SYMBOL = `Gene Symbol`) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

delayed_diffsignal = Delayed_AMPH_nucRNAseq %>% 
  mutate(sgnlogpadj = -log10(amphIPVsalIP_padj)*sign(amphIPVsalIP_lfc) ) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

rapid_diffsignal = Rapid_AMPH_nucRNAseq %>% 
  mutate(sgnlogpadj = -log10(amphIPVsalIP_padj)*sign(amphIPVsalIP_lfc) ) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

# > filtering and sorting by common genes
genelist = intersect(mecp2_diffsignal$SYMBOL, delayed_diffsignal$SYMBOL)
l1 = mecp2_diffsignal %>% 
  dplyr::filter(SYMBOL %in% genelist) %>% 
  arrange(SYMBOL)
l2 = delayed_diffsignal %>% 
  dplyr::filter(SYMBOL %in% genelist) %>% 
  arrange(SYMBOL)

# Performing RRHO  --------------------------------------------------------
RRHO(as.data.frame(l1),as.data.frame(l2), BY=TRUE, alternative='two.sided', log10.ind = TRUE, plot = TRUE, outputdir = "../output/RRHO", labels = c("MecP2  WT v KI (rank)", "PV (rank)"))


