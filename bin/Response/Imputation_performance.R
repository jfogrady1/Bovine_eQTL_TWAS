library(data.table)
library(tidyverse)
args = commandArgs(trailingOnly = TRUE)
RNA_SNPs = fread(args[1], skip = 28)
head(RNA_SNPs)
Array_SNPs = fread(args[2], skip = 28)

RNA_SNPs$Group = "RNA+DNA"
Array_SNPs$Group = "DNA"

Combined = rbind(Array_SNPs, RNA_SNPs)

all(Array_SNPs$SNP == RNA_SNPs$SNP)
rm(RNA_SNPs)
rm(Array_SNPs)


Combined$af_categories <- cut(
  Combined$MAF,
  breaks = seq(0, 0.5, length.out = 11),
  labels = c("[0.00-0.05)","[0.05-0.10)","[0.10-0.15)","[0.15-0.20)","[0.20-0.25)","[0.25-0.30)","[0.30-0.35)", "[0.35-0.40)", "[0.40-0.45)", "[0.45-0.50]"),                       # Assigns numbers to categories
  include.lowest = TRUE      # Ensures 0 is included in the first category
)
Combined$af_categories <- factor(Combined$af_categories, levels = c("[0.00-0.05)","[0.05-0.10)","[0.10-0.15)","[0.15-0.20)","[0.20-0.25)","[0.25-0.30)","[0.30-0.35)", "[0.35-0.40)", "[0.40-0.45)", "[0.45-0.50]"))
ggplot(Combined, aes(x = af_categories, y = Rsq, fill = Group)) + geom_boxplot(outlier.colour = NA, width = 0.6, alpha = 0.8) + theme_bw() +
labs(
    x = "Allele Frequency Categories",
    y = "Rsq Value", # Improve y-axis label
    fill = "Target Set" # Improve legend title
  ) +
  scale_fill_brewer(palette = "Set2") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") + # Add reference line at y=0
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Rotate x-axis text
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"), # Adjust x-axis title
    axis.title.y = element_text(size = 16, face = "bold"), # Adjust y-axis title
    legend.title = element_text(size = 16, face = "bold"), # Bold legend title
    legend.text = element_text(size = 14) # Adjust legend text size
  )
ggsave(args[3], width = 12, height =12, dpi = 600)
