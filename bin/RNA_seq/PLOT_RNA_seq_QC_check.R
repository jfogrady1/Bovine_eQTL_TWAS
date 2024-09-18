# Checking for quality control - that none of the libraries exhibited an abnormal distribution before or after filtering of lowly expressed genes.

# Libraries
library(tidyverse)
library(ggplot2)
library(tidyquant)
library(reshape2)
library(patchwork)

# Count matrix
# argument 1
counts <- read.table("count_matrix_clean.txt", header = T)
counts <- counts %>% select(-1)
counts <- log10(counts + 1)
meltdata <- melt(counts)
unique(meltdata$variable)

# All genes which were present
ggplot(meltdata,aes(value)) + geom_density(aes(group = variable)) +
  theme_tq(base_size = 15) +
  theme(axis.title  = element_text(colour = "black")) +
ylab("Density of raw gene counts per sample") +
xlab(expression(paste(log[10], "(counts + 1)")))


ggsave(filename = "Supp_Density_raw.tiff", width = 12, height = 8, dpi = 600)

# Filtering for genes which passed threshold for inclusion in eQTL study
# Will assess their ditribution
# Argument 2
counts <- read.table("../MatrixQTL_gene/Inputs/All_expression.txt", sep = ",", header = T)
counts <- counts %>% select(-1)
counts <- log10(counts + 1)
meltdata <- melt(counts)
infection <- rep("Control", (dim(counts)[1]*63))
infection <- c(infection, rep("bTB-reactor", (dim(counts)[1]*60)))
meltdata <- cbind(meltdata, infection)


# Plotting for all filtered genes in infected group
melt_data_infec = meltdata %>% filter(infection == "bTB-reactor")
melt_data_infec = melt_data_infec[order(melt_data_infec$variable,decreasing = F),]
melt_data_infec$variable <- factor(melt_data_infec$variable, levels = unique(melt_data_infec$variable))
head(melt_data_infec$variable)
Infec <- ggplot(melt_data_infec,aes(x = variable, y = value, fill = infection)) + geom_boxplot(alpha = 0.4, outlier.colour = "gray30")  + coord_flip() +
  theme_tq(base_size = 12) +
  theme(axis.title  = element_text(colour = "black")) +
  xlab("Sample") +
  ylab(expression(paste(log[10], "(counts + 1)"))) +
  ylim(0,5) +
  scale_fill_manual("", labels = "bTB-Reactor", values = c("#d7191c")) + 
  theme(legend.position = "bottom")


melt_data_cont = meltdata %>% filter(infection == "Control")
melt_data_cont = melt_data_cont[order(melt_data_cont$variable,decreasing = F),]
melt_data_cont$variable <- factor(melt_data_cont$variable, levels = unique(melt_data_cont$variable))
Cont <- ggplot(melt_data_cont,aes(x = variable, y = value, fill = infection)) + geom_boxplot(alpha = 0.4, outlier.colour = "gray30")  + coord_flip() +
  theme_tq(base_size = 12) +
  theme(axis.title  = element_text(colour = "black")) +
  xlab("Sample") +
  ylab(expression(paste(log[10], "(counts + 1)"))) +
  ylim(0,5) +
  scale_fill_manual("    Group", labels = "Control", values = c("#2c7bb6")) + 
  theme(legend.position = "bottom")

Cont + Infec + plot_layout(guides = "collect") & theme(legend.position = "right") + theme(legend.key.size = unit(2, 'cm'))
ggsave(filename = "Supp_Boxplots_RNA_filtered_genes.tiff", width = 15, height = 18, dpi = 600)
