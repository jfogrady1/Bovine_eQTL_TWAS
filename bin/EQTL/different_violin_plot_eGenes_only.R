library(data.table)
library(tidyverse)
library(ggplot2)
library(UpSetR)
args = commandArgs(trailingOnly = T)

data_ALL <- fread(args[1]) %>% filter(is_eGene == T)
data_CONTROL <- fread(args[2])  %>% filter(is_eGene == T)
data_INFECTED <- fread(args[3]) %>% filter(is_eGene == T)
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
geom_violin(alpha = 0.5, trim = TRUE) + 
geom_boxplot(width = 0.07, outlier.colour = NA) + 
scale_fill_manual(values = c("#542788","#2166ac", "#b2182b")) + 
theme_bw() +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
    axis.title.y = element_text(size = 21, color = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15))
ggsave(args[7], width = 10, height = 12, dpi = 600)

