# Checking for quality control - that none of the libraries exhibited an abnormal distribution before or after filtering of lowly expressed genes.

# Libraries
library(tidyverse)
library(ggplot2)
library(tidyquant)
library(reshape2)
library(patchwork)
args = commandArgs(trailingOnly = T)
# Count matrix
# argument 1
counts <- read.table(args[1], header = T)
counts <- counts %>% dplyr::select(-1)
counts <- log10(counts + 1)
meltdata <- melt(counts)
unique(meltdata$variable)

# All genes which were present
ggplot(meltdata,aes(value)) + geom_density(aes(group = variable)) +
theme_bw() +
  theme(axis.title  = element_text(colour = "black")) +
ylab("Density of raw gene counts per sample") +
xlab(expression(paste(log[10], "(counts + 1)"))) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 21, color = "black"),
        axis.title.x = element_text(size = 21, color = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "black"))


ggsave(filename = args[2], width = 12, height = 8, dpi = 600)