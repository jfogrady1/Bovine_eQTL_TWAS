library(data.table)
library(tidyverse)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library("ggsci")
args = commandArgs(trailingOnly = T)

# read in the data
#args[1] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_CH_TWAS_results.txt"
#args[2] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_LM_TWAS_results.txt"
#args[3] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_HOFR_TWAS_results.txt"
#args[4] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_ALL_TWAS_results.txt"
#args[5] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_CH_TWAS_results.txt"
#args[6] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_LM_TWAS_results.txt"
#args[7] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_HOFR_TWAS_results.txt"
#args[8] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_ALL_TWAS_results.txt"
#args[9] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_CH_TWAS_results.txt"
##args[10] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_LM_TWAS_results.txt"
#args[11] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_HOFR_TWAS_results.txt"
#args[12] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_ALL_TWAS_results.txt"
#args[13] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bovine_annotation_MF2.csv"
#args[14] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_Supplementary.txt"
#args[15] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_Supplementary.txt"
#args[16] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_Supplementary.txt"
#args[17] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/Plotting/INFECTED_volcano_TWAS.pdf"
#args[18] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/Plotting/CONTROL_volcano_TWAS.pdf"

ALL_CH <- fread(args[1]) %>% mutate(Group = "AA")
ALL_LM <- fread(args[2])  %>% mutate(Group = "AA")
ALL_HF <- fread(args[3])  %>% mutate(Group = "AA")
ALL_Multi <- fread(args[4])  %>% mutate(Group = "AA")
CON_CH <- fread(args[5])  %>% mutate(Group = "CT")
CON_LM <- fread(args[6])  %>% mutate(Group = "CT")
CON_HF <- fread(args[7])  %>% mutate(Group = "CT")
CON_Multi <- fread(args[8])  %>% mutate(Group = "CT")
REA_CH <- fread(args[9])  %>% mutate(Group = "RT")
REA_LM <- fread(args[10])  %>% mutate(Group = "RT")
REA_HF <- fread(args[11])  %>% mutate(Group = "RT")
REA_Multi <- fread(args[12])  %>% mutate(Group = "RT")

ALL_LM$GWAS <- "LM"
CON_LM$GWAS <- "LM"
REA_LM$GWAS <- "LM"

ALL_HF$GWAS <- "HF"
CON_HF$GWAS <- "HF"
REA_HF$GWAS <- "HF"

ALL_CH$GWAS <- "CH"
CON_CH$GWAS <- "CH"
REA_CH$GWAS <- "CH"

ALL_Multi$GWAS <- "Multi"
CON_Multi$GWAS <- "Multi"
REA_Multi$GWAS <- "Multi"


Infected_ALL <- rbind(REA_CH, REA_LM, REA_HF, REA_Multi)
Control_ALL <- rbind(CON_CH, CON_LM, CON_HF, CON_Multi)
ALL_ALL <- rbind(ALL_CH, ALL_LM, ALL_HF, ALL_Multi)


annotation <- fread(args[13]) %>% select(3,2)
annotation <- annotation[!(duplicated(annotation$Ensembl)),]
dim(annotation)
Infected_ALL <- left_join(Infected_ALL, annotation, by = c("Gene" = "Ensembl"))
Infected_ALL$Symbol <- if_else(is.na(Infected_ALL$Symbol), Infected_ALL$Gene, Infected_ALL$Symbol)
colnames(Infected_ALL)
Infected_ALL <- Infected_ALL %>% select(1,9,2,3,4,5,6,7,8)
head(Infected_ALL)

Control_ALL <- left_join(Control_ALL, annotation, by = c("Gene" = "Ensembl"))
Control_ALL$Symbol <- if_else(is.na(Control_ALL$Symbol), Control_ALL$Gene, Control_ALL$Symbol)
colnames(Control_ALL)
Control_ALL <- Control_ALL %>% select(1,9,2,3,4,5,6,7,8)
head(Control_ALL)

ALL_ALL <- left_join(ALL_ALL, annotation, by = c("Gene" = "Ensembl"))
ALL_ALL$Symbol <- if_else(is.na(ALL_ALL$Symbol), ALL_ALL$Gene, ALL_ALL$Symbol)
colnames(ALL_ALL)
ALL_ALL <- ALL_ALL %>% select(1,9,2,3,4,5,6,7,8)
head(ALL_ALL)


write.table(Infected_ALL, args[14], sep = "\t", quote = F, col.names = T, row.names = F)
write.table(Control_ALL, args[15], sep = "\t", quote = F, col.names = T, row.names = F)
write.table(ALL_ALL, args[16], sep = "\t", quote = F, col.names = T, row.names = F)


head(ALL_ALL) # keep for  plotting the circos

# Volcano plots of the other two

# cut offs
REA_cut_off <- 0.00003243 # number of expression models
CON_cut_off <- 0.00002048 # number of expression models

Infected_ALL$Z <- if_else(Infected_ALL$Z == "Cannot compute LD covariance matrix" | Infected_ALL$Z == "No Overlapping SNPs.", 0, as.numeric(Infected_ALL$Z))
Infected_ALL$P <- if_else(Infected_ALL$P == "Cannot compute LD covariance matrix" | Infected_ALL$P == "No Overlapping SNPs.", 1, as.numeric(Infected_ALL$P))
Infected_ALL$permute.P <- if_else(Infected_ALL$permute.P == "Cannot compute LD covariance matrix" | Infected_ALL$permute.P == "No Overlapping SNPs.", 1, as.numeric(Infected_ALL$permute.P))
# remember to convert to numeric
Infected_ALL$signif <- if_else(as.numeric(Infected_ALL$P) < REA_cut_off, "TRUE", "FALSE")
Infected_ALL$permute_signif <- if_else(Infected_ALL$signif == "TRUE" & Infected_ALL$permute.P < 0.05, "TRUE", "FALSE")
Infected_ALL <- Infected_ALL %>% mutate(Shape = case_when(
     signif == "TRUE" & permute_signif == "TRUE" ~ "Padj < 0.05 & Ppermute < 0.05",
    signif == "TRUE" & permute_signif == "FALSE" ~ "Padj < 0.05",
    signif == "FALSE" ~ "Not significnat"))

Infected_ALL <- Infected_ALL %>% mutate(Alpha = case_when(
     signif == "TRUE" & permute_signif == "TRUE" ~ 1,
    signif == "TRUE" & permute_signif == "FALSE" ~ 0.2,
    signif == "FALSE" ~ 0.0001))

Infected_ALL <- Infected_ALL %>% mutate(Colour = case_when(
     GWAS == "HF" ~ "#053061",
    GWAS == "CH" ~ "grey50", 
    GWAS == "LM" ~ "#A85307",
    GWAS == "Multi" ~ "black"))


Infected_ALL$label <- if_else(Infected_ALL$signif == "TRUE" & Infected_ALL$permute_signif == "TRUE", Infected_ALL$Symbol, NA)

ggplot(data = Infected_ALL, aes(x = as.numeric(Z), y = -log10(as.numeric(P)), shape = Shape, colour = GWAS, alpha = Alpha, label = label)) +
geom_point(size = 3) + scale_color_aaas(breaks = c("CH", "HF", "LM", "Multi")) +
theme_bw() + labs(x = "TWAS Z-score",
                  y = "-log10(P)",
                  shape = "Significance",
                  colour = "GWAS set") +
                  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12)) +
        geom_hline(yintercept = -log10(REA_cut_off), linetype = 2) +
         geom_text_repel(nudge_y = 0.5, size = 3, fontface = "bold")
ggsave(args[17], width = 12, height = 12, dpi = 600)



                  
Control_ALL$Z <- if_else(Control_ALL$Z == "Cannot compute LD covariance matrix" | Control_ALL$Z == "No Overlapping SNPs.", 0, as.numeric(Control_ALL$Z))
Control_ALL$P <- if_else(Control_ALL$P == "Cannot compute LD covariance matrix" | Control_ALL$P == "No Overlapping SNPs.", 1, as.numeric(Control_ALL$P))
Control_ALL$permute.P <- if_else(Control_ALL$permute.P == "Cannot compute LD covariance matrix" | Control_ALL$permute.P == "No Overlapping SNPs.", 1, as.numeric(Control_ALL$permute.P))
# remember to convert to numeric
Control_ALL$signif <- if_else(as.numeric(Control_ALL$P) < REA_cut_off, "TRUE", "FALSE")
Control_ALL$permute_signif <- if_else(Control_ALL$signif == "TRUE" & Control_ALL$permute.P < 0.05, "TRUE", "FALSE")
Control_ALL <- Control_ALL %>% mutate(Shape = case_when(
     signif == "TRUE" & permute_signif == "TRUE" ~ "Padj < 0.05 & Ppermute < 0.05",
    signif == "TRUE" & permute_signif == "FALSE" ~ "Padj < 0.05",
    signif == "FALSE" ~ "Not significnat"))

Control_ALL <- Control_ALL %>% mutate(Alpha = case_when(
     signif == "TRUE" & permute_signif == "TRUE" ~ 1,
    signif == "TRUE" & permute_signif == "FALSE" ~ 0.2,
    signif == "FALSE" ~ 0.0001))

Control_ALL <- Control_ALL %>% mutate(Colour = case_when(
     GWAS == "HF" ~ "#053061",
    GWAS == "CH" ~ "grey50", 
    GWAS == "LM" ~ "#A85307",
    GWAS == "Multi" ~ "black"))


Control_ALL$label <- if_else(Control_ALL$signif == "TRUE" & Control_ALL$permute_signif == "TRUE", Control_ALL$Symbol, NA)

ggplot(data = Control_ALL, aes(x = as.numeric(Z), y = -log10(as.numeric(P)), shape = Shape, colour = GWAS, alpha = Alpha, label = label)) +
geom_point(size = 3) + scale_color_aaas(breaks = c("CH", "HF", "LM", "Multi")) +
theme_bw() + labs(x = "TWAS Z-score",
                  y = "-log10(P)",
                  shape = "Significance",
                  colour = "GWAS set") +
                  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12)) +
        geom_hline(yintercept = -log10(CON_cut_off), linetype = 2) +
         geom_text_repel(nudge_y = 0.5, size = 3, fontface = "bold")
ggsave(args[18], width = 12, height = 12, dpi = 600)

