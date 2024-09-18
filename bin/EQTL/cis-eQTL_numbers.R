# Script to get eQTL numbers which are significant
# Also to get a violin plot of gene-wise FDR cutoffs


library(data.table)
library(tidyverse)
library(ggplot2)
library(UpSetR)
args = commandArgs(trailingOnly = T)
data_ALL <- fread(args[1])
data_CONTROL <- fread(args[2])
data_INFECTED <- fread(args[3])
ALL_nominal <- fread(args[4])
CONTROL_nominal <- fread(args[5])
INFECTED_nominal <- fread(args[6])

# Get rid of the replicated headings
ALL_nominal <- ALL_nominal %>% filter(phenotype_id != "phenotype_id")
CONTROL_nominal <- CONTROL_nominal %>% filter(phenotype_id != "phenotype_id")
INFECTED_nominal <- INFECTED_nominal %>% filter(phenotype_id != "phenotype_id")
# Plot gene wise threshold
threshold_plot_ALL <- data_ALL %>% select(1,19)
threshold_plot_CONTROL <- data_CONTROL %>% select(1,19)
threshold_plot_INFECTED <- data_INFECTED %>% select(1,19)

threshold_plot_ALL$Group <- "All"
threshold_plot_CONTROL$Group <- "Control"
threshold_plot_INFECTED$Group <- "Reactor"

threshold_plot <- rbind(threshold_plot_ALL, threshold_plot_CONTROL, threshold_plot_INFECTED)
ggplot(data = threshold_plot, aes(x = Group, y = -log10(pval_nominal_threshold), fill = Group)) + 
geom_violin(alpha = 0.5) + 
geom_boxplot(width = 0.07) + 
scale_fill_manual(values = c("#542788","#2166ac", "#b2182b")) + 
theme_bw() +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
    axis.title.y = element_text(size = 21, color = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15))
ggsave(args[7], width = 10, height = 12, dpi = 600)

# get pvalues for significance
all_control_pvalue <- wilcox.test(threshold_plot_ALL$pval_nominal_threshold, threshold_plot_CONTROL$pval_nominal_threshold)$p.value
control_reactor_pvalue <- wilcox.test(threshold_plot_CONTROL$pval_nominal_threshold, threshold_plot_INFECTED$pval_nominal_threshold)$p.value
all_reactor_pvalue <- wilcox.test(threshold_plot_ALL$pval_nominal_threshold, threshold_plot_INFECTED$pval_nominal_threshold)$p.value

# Plot upset plot of genes tested

listInput <- list(All = as.vector(threshold_plot_ALL$phenotype_id), Control = as.vector(threshold_plot_CONTROL$phenotype_id), Reactor = as.vector(threshold_plot_INFECTED$phenotype_id))

pdf(file=args[8], width = 12, height = 12, onefile=FALSE)
upset(fromList(listInput), order.by = "freq", sets.bar.color = c("#b2182b","#542788","#2166ac"),
sets.x.label = "cis-eGenes tested", point.size = 4, line.size = 2,
mainbar.y.label = "cis-eGenes tested intersections",
text.scale = 2.5, shade.alpha = 0.5)
dev.off()


ALL_nominal <- ALL_nominal %>% filter(phenotype_id != "phenotype_id")
CONTROL_nominal <- CONTROL_nominal %>% filter(phenotype_id != "phenotype_id")
INFECTED_nominal <- INFECTED_nominal %>% filter(phenotype_id != "phenotype_id")

data_ALL$is_eGene <- as.character(data_ALL$is_eGene)
data_CONTROL$is_eGene <- as.character(data_CONTROL$is_eGene)
data_INFECTED$is_eGene <- as.character(data_INFECTED$is_eGene)
dim(data_ALL)
dim(data_CONTROL)
dim(data_INFECTED)

data_ALL <- data_ALL %>% filter(is_eGene == "TRUE")
data_CONTROL <- data_CONTROL %>% filter(is_eGene == "TRUE")
data_INFECTED <- data_INFECTED %>% filter(is_eGene == "TRUE")

data_ALL <- left_join(data_ALL, ALL_nominal, by = c("phenotype_id" = "phenotype_id"))
data_CONTROL <- left_join(data_CONTROL, CONTROL_nominal, by = c("phenotype_id" = "phenotype_id"))
data_INFECTED <- left_join(data_INFECTED, INFECTED_nominal, by = c("phenotype_id" = "phenotype_id"))


dim(data_ALL)
dim(data_INFECTED)
dim(data_CONTROL)
data_ALL <- data_ALL %>% filter(as.numeric(pval_nominal.y) < as.numeric(pval_nominal_threshold))
data_CONTROL <- data_CONTROL %>% filter(as.numeric(pval_nominal.y) < as.numeric(pval_nominal_threshold))
data_INFECTED <- data_INFECTED %>% filter(as.numeric(pval_nominal.y) < as.numeric(pval_nominal_threshold))

rm(ALL_nominal, CONTROL_nominal, INFECTED_nominal)

number_ALL <- length(data_ALL$variant_id.y)
unique_ALL <- length(unique(data_ALL$variant_id.y))
number_CONTROL <- length(data_CONTROL$variant_id.y)
unique_CONTROL <- length(unique(data_CONTROL$variant_id.y))
number_INFECTED <- length(data_INFECTED$variant_id.y)
unique_INFECTED <- length(unique(data_INFECTED$variant_id.y))

write.table("Number of signicant cis-eQTLs in ALL group", file = args[9], sep = "\t", col.names = F )
write.table(number_ALL, file = args[9], sep = "\t", col.names = F, append = T)
write.table("Number of unique cis-eQTLs in ALL group", file = args[9], sep = "\t", col.names = F, append = T)
write.table(unique_ALL, file = args[9], sep = "\t", col.names = F, append = T)

write.table("Number of signicant cis-eQTLs in CONTROL group", file = args[9], sep = "\t", col.names = F, append = T )
write.table(number_CONTROL, file = args[9], sep = "\t", col.names = F, append = T)
write.table("Number of unique cis-eQTLs in CONTROL group", file = args[9], sep = "\t", col.names = F, append = T)
write.table(unique_CONTROL, file = args[9], sep = "\t", col.names = F, append = T)

write.table("Number of signicant cis-eQTLs in INFECTED group", file = args[9], sep = "\t", col.names = F, append = T )
write.table(number_INFECTED, file = args[9], sep = "\t", col.names = F, append = T)
write.table("Number of unique cis-eQTLs in INFECTED group", file = args[9], sep = "\t", col.names = F, append = T)
write.table(unique_INFECTED, file = args[9], sep = "\t", col.names = F, append = T)

# write files with cis-eQTL pairs for TWAS which are significant
write.table(data_ALL, file = args[10], sep = "\t", col.names = T, quote = F, row.names = F)
write.table(data_CONTROL, file = args[11],sep = "\t", col.names = T, quote = F, row.names = F)
write.table(data_INFECTED, file = args[12],sep = "\t", col.names = T, quote = F, row.names = F)

