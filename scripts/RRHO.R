# Author: Melyssa Minto
# This script will create the RRHO plot for the Differential Signal between experiments in the paper


# Loading Libraries -------------------------------------------------------

library(tidyverse)
library(RRHO)
library(pheatmap)
library(dichromat)
library(readxl)
library(ggpubr)
library(gridExtra)
library(xlsx)
# Defining functions ------------------------------------------------------

customRRHO = function(diff_signal1, diff_signal2, dir, labels, conditions, lims=NULL){
  
  # > setting up ourput directories
  rmcmd = paste0("rm -r ../output/RRHO/", dir, "/")
  makecmd = paste0("mkdir ../output/RRHO/", dir, "/")
  system(rmcmd)
  system(makecmd)
  
  
  # > filtering and sorting by common genes
  genelist = intersect(diff_signal1$SYMBOL, diff_signal2$SYMBOL)
  cat('Common Genes: ')
  cat(length(genelist))
  l1 = diff_signal1 %>% 
    dplyr::filter(SYMBOL %in% genelist) %>% 
    arrange(SYMBOL)
  l2 = diff_signal2 %>% 
    dplyr::filter(SYMBOL %in% genelist) %>% 
    arrange(SYMBOL)
  

  
  # > performing RRHO
  object = RRHO(as.data.frame(l1),
                as.data.frame(l2),
                BY=TRUE,
                stepsize = 100,
                alternative='two.sided',
                log10.ind = TRUE,
                plot = TRUE,
                outputdir = paste0("../output/RRHO/", dir, "/"),
                labels = labels)

  # extracting and rotating matrix
  mat = object$hypermat
  mat = t(mat)
  #mat = mat[nrow(mat):1,]
  
  # finding quadrants
  mat.signs = object$hypermat.signs
  
  
  mycol <- c("navy", "blue", "cyan", "green", "yellow", "red", "red4")

 p_map = mat %>% 
    as.data.frame() %>%
    rownames_to_column("diff_sig1") %>%
    pivot_longer(-c(diff_sig1), names_to = "diff_sig2", values_to = "hypergeometric_pvals") %>% 
    mutate(diff_sig1= fct_relevel(diff_sig1,rownames(as.data.frame(mat)))) %>%
    mutate(diff_sig2= fct_relevel(diff_sig2,colnames(as.data.frame(mat))))%>%
    ggplot(aes(x=diff_sig2, y=diff_sig1, fill=hypergeometric_pvals)) + 
    geom_raster() +
   scale_fill_gradientn(colours = mycol, limits = lims, n.breaks=7)+
   theme(axis.text=element_blank(),
         axis.ticks=element_blank(),
         legend.key.height= unit(2.7, 'cm'),
         legend.key.width= unit(1, 'cm')) +
   labs(x = "", y = "", fill = "-log10(p-value)") +
   ggtitle("Rank Rank Hypergeometric Overlap Map")+
   theme(plot.title = element_text(hjust = 0.5, face="bold"))
 


  
  png(paste0("../output/RRHO/", dir, "/RRHO.png"), res = 300, units = 'in', height = 6.62,width =7.02 )
  print( p_map  )
 
  grid.lines(x = unit(c(0.05,.8), "npc"), y = unit(c(0.03,0.03), "npc"), arrow = arrow(length = unit(.4, "cm"), type = "closed", ends = "both"))
  grid.text(conditions[2], x = unit(0.7, "npc"), y = unit(0.015, "npc"))
  grid.text(conditions[1], x = unit(0.15, "npc"), y = unit(0.015, "npc"))
  
  grid.lines(x = unit(c(0.03,.03), "npc"), y = unit(c(0.06,.94), "npc"), arrow = arrow(length = unit(.4, "cm"), type = "closed", ends = "both"), name = "Help")
  grid.text(conditions[4], x = unit(0.015, "npc"), y = unit(0.85, "npc"),rot = 90)
  grid.text(conditions[3], x = unit(0.015, "npc"), y = unit(0.17, "npc"),rot = 90)
  
  dev.off()
  
}
# Reading in Data ---------------------------------------------------------

mecp2WTvKI <- read_delim("../data/mecp2WTvKI.gene_exp_diff", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  dplyr::filter(sample_1 %in% "WT_IP" & sample_2 %in% "KO_IP") # samples were originally mislabeled as KO instead of KI
mecp2IPvInput <- read_delim("../data/mecp2WTvKI.gene_exp_diff", delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>% 
  dplyr::filter(sample_1 %in% "WT_IP" & sample_2 %in% "WT_INPUT") 

PV_Delayed_AMPH_nucRNAseq <- read_csv("../output/PV3HrSalvAMPH_RNA/diffExpressionMatrix.csv")
PV_Rapid_AMPH_nucRNAseq <- read_csv("../output/PV35minSalvAMPH_RNA/diffExpressionMatrix.csv")
SST_Delayed_AMPH_nucRNAseq <- read_csv("../output/SST3HrSalVAMPH_RNA/diffExpressionMatrix.csv")
SST_Rapid_AMPH_nucRNAseq <- read_csv("../output/SST3HrSalVAMPH_RNA/diffExpressionMatrix.csv")
UF_3hr <- read_csv("../output/CombinedUF3hr_RNA/diffExpressionMatrix.csv")
UF_35min <- read_csv("../output/CombinedUF35min_RNA/diffExpressionMatrix.csv")
PVIPvUF <- read_csv("../output/PVIPvUF_RNA/diffExpressionMatrix.csv")
SSTIPvUF <- read_excel("../FinalTables/Table S1-NAc PV and SST nucRNAseq.xlsx", sheet = "ALL SST IPvUF")
PV_24hr <-read_csv("../output/PV24HrPostRepeated_RNA/diffExpressionMatrix.csv")

# Preparing data for RRHO -------------------------------------------------
# > using signed -log10 FDRs for plot
mecp2WTvKI =  mecp2WTvKI[!duplicated(mecp2WTvKI$gene), ]
mecp2_diffsignal = mecp2WTvKI %>% 
  mutate(sgnlogpadj = log10(q_value)* sign(as.numeric(`log2(fold_change)`))) %>% 
  dplyr::rename( SYMBOL = gene) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

mecp2IPvInput =  mecp2IPvInput[!duplicated(mecp2IPvInput$gene), ]
mecp2IP_diffsignal = mecp2IPvInput %>% 
  mutate(sgnlogpadj = log10(q_value)* sign(as.numeric(`log2(fold_change)`))) %>% 
  dplyr::rename( SYMBOL = gene) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

mecp2IPvInput_opp = mecp2IPvInput %>% 
  mutate(sgnlogpadj = -log10(q_value)* sign(as.numeric(`log2(fold_change)`))) %>% 
  dplyr::rename( SYMBOL = gene) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

mecp2_diffsignal_opp = mecp2WTvKI %>% 
  mutate(sgnlogpadj = -log10(q_value)* sign(as.numeric(`log2(fold_change)`))) %>% 
  dplyr::rename( SYMBOL = gene) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

pv_delayed_diffsignal = PV_Delayed_AMPH_nucRNAseq %>% 
  mutate(sgnlogpadj = log10(amphIPVsalIP_padj)*sign(amphIPVsalIP_lfc) ) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

pv_delayed_diffsignal_opp = PV_Delayed_AMPH_nucRNAseq %>% 
  mutate(sgnlogpadj = -log10(amphIPVsalIP_padj)*sign(amphIPVsalIP_lfc) ) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

pv_rapid_diffsignal_opp = PV_Rapid_AMPH_nucRNAseq %>% 
  mutate(sgnlogpadj = -log10(amphIPVsalIP_padj)*sign(amphIPVsalIP_lfc) ) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

pv_rapid_diffsignal = PV_Rapid_AMPH_nucRNAseq %>% 
  mutate(sgnlogpadj = log10(amphIPVsalIP_padj)*sign(amphIPVsalIP_lfc) ) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

sst_delayed_diffsignal = SST_Delayed_AMPH_nucRNAseq %>% 
  mutate(sgnlogpadj = log10(amphIPVsalIP_padj)*sign(amphIPVsalIP_lfc) ) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

sst_rapid_diffsignal = SST_Rapid_AMPH_nucRNAseq %>% 
  mutate(sgnlogpadj = log10(amphIPVsalIP_padj)*sign(amphIPVsalIP_lfc) ) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

uf_3hr_diffsignal = UF_3hr %>% 
  mutate(sgnlogpadj = log10(amphUFVsalUF_padj)*sign(amphUFVsalUF_lfc) ) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

uf_35min_diffsignal = UF_35min %>% 
  mutate(sgnlogpadj = log10(amphUFVsalUF_padj)*sign(amphUFVsalUF_lfc) ) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()


pv_ipvuf_diffsignal = PVIPvUF %>% 
  mutate(sgnlogpadj = log10(salIPvUF_padj)*sign(salIPvUF_lfc) ) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

sst_ipvuf_diffsignal = SSTIPvUF %>% 
  mutate(sgnlogpadj = log10(`SAL IPvUF FDR`)*sign(`SAL IPvUF Log2 FC`) ) %>% 
  dplyr::rename(SYMBOL = `GENE SYMBOL`) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()

pv_24hr_diffsignal = PV_24hr %>% 
  mutate(sgnlogpadj = log10(amphIPVsalIP_padj)*sign(amphIPVsalIP_lfc)) %>% 
  dplyr::select(SYMBOL, sgnlogpadj) %>% 
  drop_na()


# Run RRHO --------------------------------------------------------------

customRRHO(mecp2_diffsignal, pv_delayed_diffsignal, "mecp2_PVDelayed_3hr", c("MecP2  KI v WT (rank)", "PV Delayed Response Amph v Sal (rank)"), conditions = c("MeCP2 KI", "MeCP2 WT", "PV 3hr Amph", "PV 3hr Sal"))
customRRHO(mecp2_diffsignal, pv_delayed_diffsignal_opp, "mecp2_PVDelayed_3hr_opp", c("MecP2  KI v WT (rank)", "PV Delayed Response Sal v Amph (rank)"), conditions = c("MeCP2 KI", "MeCP2 WT", "PV 3hr Sal", "PV 3hr Amph"))

customRRHO(mecp2_diffsignal, pv_rapid_diffsignal, "mecp2_PVRapid_35min", c("MecP2  KI v WT (rank)", "PV Rapid Response Amph v Sal (rank)"), conditions = c("MeCP2 KI", "MeCP2 WT", "PV 35min Amph", "PV 35min Sal"))
customRRHO(mecp2_diffsignal, pv_rapid_diffsignal_opp, "mecp2_PVRapid_35min_opp", c("MecP2  KI v WT (rank)", "PV Rapid Response Sal v Amph (rank)"), conditions = c("MeCP2 KI", "MeCP2 WT", "PV 35min Sal", "PV 35min Amph"))

customRRHO(pv_rapid_diffsignal, pv_delayed_diffsignal, "PV_Rapid_Delayed", c("Rapid Amph v Sal (rank)", "Delayed Amph v Sal (rank)"), conditions = c("PV 35min Amph", "PV 35min Sal", "PV 3hr Amph", "PV 3hr Sal"))

customRRHO(sst_delayed_diffsignal, pv_delayed_diffsignal_opp, "Delayed_PVvSST_opp",c("SST Amph v Sal (rank)", "PV Sal v Amph (rank)"), lims = c(-5, 230), conditions = c("SST 3hr Amph", "SST 3hr Sal", "PV 3hr Sal", "PV 3hr Amph"))
customRRHO(sst_delayed_diffsignal, pv_delayed_diffsignal, "Delayed_PVvSST",c("SST Amph v Sal (rank)", "PV Amph v Sal (rank)"), lims = c(-5, 230), conditions = c("SST 3hr Amph", "SST 3hr Sal", "PV 3hr Amph", "PV 3hr Sal"))

customRRHO(sst_rapid_diffsignal, pv_rapid_diffsignal, "Rapid_PVvSST",c("SST Amph v Sal (rank)", "PV Amph v Sal (rank)"), lims = c(-5, 230), conditions = c("SST 35min Amph", "SST 35min Sal", "PV 35min Amph", "PV 35min Sal"))

customRRHO(mecp2_diffsignal_opp, pv_ipvuf_diffsignal, "mecp2KIvWT_PVIPvUF", c("MecP2  WT v KI (rank)", "PV Sal IP v UF (rank)"), conditions = c("MeCP2 WT", "MeCP2 KI", "PV IP", "UF"))

customRRHO(mecp2IPvInput_opp, pv_ipvuf_diffsignal, "mecp2IPvInput_PVIPvUF", c("MecP2  IP v Input (rank)", "PV Sal IP v UF (rank)"), lim = c(-1,93), conditions = c("MeCP2 IP", "Input", "PV IP", "UF"))

customRRHO(sst_ipvuf_diffsignal, pv_ipvuf_diffsignal, "SST_PVIPvUF", c("SST Sal IP v UP (rank)", "PV Sal IP v UF (rank)"), conditions = c("SST IP", "UF", "PV IP", "UF"))

customRRHO( uf_3hr_diffsignal, uf_35min_diffsignal, "UF_3hrv35min",c("UF Delayed Amph v Sal (rank)", "UF Rapid Amph v Sal (rank)"), conditions = c("UF 3hr Amph", "UF 3hr Sal", "UF 35min Amph", "UF 35min Sal"))

customRRHO( pv_delayed_diffsignal, pv_24hr_diffsignal, "PV_3hrv24hr", c("3hr Amph v Sal (rank)", "24hr Amph v Sal (rank)"), conditions = c("PV 3hr Amph", "PV 3hr Sal", "PV 24hr Amph", "PV 24hr Sal"))

## writing final table of overlapping genes
wrangle_overlap <- function(dir){
  up_file = list.files(paste0("../output/RRHO/", dir, "/"), pattern = "Upregulated", full.names = T)
  down_file = list.files(paste0("../output/RRHO/", dir, "/"), pattern = "Downregulated", full.names = T)
  
  as.data.frame(bind_rows( read_csv(down_file, col_names = FALSE, col_types = cols()) %>% mutate(RRHO = "Downregulated"), 
             read_csv(up_file, col_names = FALSE, col_types = cols()) %>% mutate(RRHO = "Upregulated")) %>% 
    dplyr::rename(SYMBOL = X1))
  
}


write.xlsx(wrangle_overlap("Rapid_PVvSST"), "../output/RRHO/RRHO_Overlaps.xlsx", sheetName="Rapid PV v SST",  append=TRUE, row.names = FALSE)
write.xlsx(wrangle_overlap("Delayed_PVvSST"), "../output/RRHO/RRHO_Overlaps.xlsx", sheetName="Delayed PV v SST",  append=TRUE, row.names = FALSE)
write.xlsx(wrangle_overlap("mecp2IPvInput_PVIPvUF"), "../output/RRHO/RRHO_Overlaps.xlsx", sheetName="PV Ribotag v INTACT",  append=TRUE, row.names = FALSE)
write.xlsx(wrangle_overlap("SST_PVIPvUF"), "../output/RRHO/RRHO_Overlaps.xlsx", sheetName="IP SST v PV",  append=TRUE, row.names = FALSE)
write.xlsx(wrangle_overlap("PV_3hrv24hr"), "../output/RRHO/RRHO_Overlaps.xlsx", sheetName="PV 3hr v 24hr",  append=TRUE, row.names = FALSE)

