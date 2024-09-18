library(UpSetR)
library(tidyverse)
library(ggridges)
library(ggplot2)


args = commandArgs(trailing = T)

########################################################################################
########################################################################################
## 1. Upset plot of cis-eGenes
########################################################################################
########################################################################################

# Now onto upset plot of eGenes in reference groups
#args[1-3]
ALL_perm <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt", header = T) %>% 
filter(is_eGene == T) 
CONTROL_perm <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_fdr0.05.txt", header = T) %>% 
filter(is_eGene == T) 
INFECTED_perm <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_fdr0.05.txt", header = T) %>% 
filter(is_eGene == T) 

ALL_perm$group <- "ALL"
CONTROL_perm$group <- "CONTROL"
INFECTED_perm$group <- "REACTOR"

all_merged <- rbind(ALL_perm, CONTROL_perm, INFECTED_perm)


ggplot(data  = all_merged, aes(x=start_distance, fill = group)) + geom_histogram(bins = 40) +
scale_fill_manual(labels = c("All", "Control", "Reactor"), values = c("#542788","#2166ac", "#b2182b")) + 
  theme_bw() +
  labs(y = "Counts", x = "Cis-eQTL Distance to Transcriptional Start Site (TSS)") +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"))




dim(ALL_eGene)

ALL <- read.table(gzfile("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_independent_qtl.txt.gz"), sep = "\t", header = T)
CONTROL <- read.table(gzfile("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_independent_qtl.txt.gz"), sep = "\t", header = T)
INFECTED <- read.table(gzfile("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_independent_qtl.txt.gz"), sep = "\t", header = T)


ALL_perm <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt", header = T) %>% 
filter(is_eGene == T) %>% select(1)
CONTROL_perm <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_fdr0.05.txt", header = T) %>% 
filter(is_eGene == T) %>% select(1)
INFECTED_perm <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_fdr0.05.txt", header = T) %>% 
filter(is_eGene == T) %>% select(1)

ALL <- inner_join(ALL, ALL_perm) %>% filter(rank >=2)
CONTROL <- inner_join(CONTROL, CONTROL_perm)  %>% filter(rank >=2)
INFECTED <- inner_join(INFECTED, INFECTED_perm)  %>% filter(rank >=2)

ALL$group <- "ALL"
CONTROL$group <- "CONTROL"
INFECTED$group <- "REACTOR"

ALL_conditional <- rbind(ALL, CONTROL, INFECTED)
colnames(ALL_conditional)
colnames(all_merged)
all_merged <- all_merged %>% select(-c("is_eGene", "pval_adj_BH", "qval", "pval_nominal_threshold"))
ALL_conditional <- ALL_conditional %>% select(-c("rank"))
all_merged$cohort <- "Permutation"
ALL_conditional$cohort <- "Conditional"

final_ALL <- rbind(all_merged, ALL_conditional)
ggplot(data  =final_ALL, aes(x=as.numeric(start_distance), fill = group)) + geom_histogram(bins = 40) +
scale_fill_manual(labels = c("All", "Control", "Reactor"), values = c("#542788","#2166ac", "#b2182b")) + 
  theme_bw() +
  scale_x_continuous(breaks = c(-1000000, -500000, 0, 500000, 1000000), labels = c("-1000", "-500", "0", "500", "1000")) +
  labs(y = "Count", x = "Cis-eQTL Distance to Transcriptional Start Site (TSS)") +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black")) + facet_wrap(~factor(cohort, levels = c("Permutation", "Conditional")))
ggsave("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/Symmetrical_cis_eQTL.pdf", width = 12, height = 12, dpi = 600)
