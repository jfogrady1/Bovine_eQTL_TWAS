library(tidyverse)
library(RColorBrewer)
library(ggsignif)
library(ggpubr)
library(cowplot)

args = commandArgs(trailingOnly = TRUE)
path = args[1]
hsq_files = list.files(path = path, pattern = "\\SNP_based\\.hsq$", full.names = TRUE)

SNP_heritability = data.frame(matrix(nrow = 0, ncol = 5))
colnames(SNP_heritability) <- c("Gene", "Vg", "h2", "h2_SE", "Pval")
for (file in hsq_files) {
    file_temp = gsub("//", "/", file)
    hsq = data.table::fread(file_temp, fill = TRUE)
    Vg = hsq$Variance[1] 
    h2 = hsq$Variance[4]
    SE = hsq$SE[4]
    P = hsq$Variance[9]
    file <- gsub(paste0(path, "/"), "", file)
    gene_name <- strsplit(file, "_")[[1]][1]
    row = c(gene_name, Vg, h2, SE, P)
    SNP_heritability = rbind(SNP_heritability, row)
}

colnames(SNP_heritability) <- c("Gene", "Vg", "h2", "h2_SE", "Pval")


SNP_heritability$h2 <- as.numeric(SNP_heritability$h2)
SNP_heritability$Pval <- as.numeric(SNP_heritability$Pval)
SNP_heritability$h2_SE <- as.numeric(SNP_heritability$h2_SE)
SNP_heritability$Padj <- p.adjust(SNP_heritability$Pval, method = "BH")



median(SNP_heritability$h2)
summary(SNP_heritability$h2)
sd(SNP_heritability$h2)
dim(SNP_heritability)

# Standard heritability
ALL_genes = ggplot(data = SNP_heritability, aes(x = "All Genes", y = h2)) + geom_boxplot(fill = "#1b9e77", outlier.colour = NA, alpha = 0.7) +
              labs(x = "",
              y = "Cis-heritability") +
              theme_bw() +
                theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14), # Rotate x-axis text
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16), # Adjust x-axis title
    axis.title.y = element_text(size = 16),# Adjust legend text size
  ) + guides(fill = "none") + stat_summary(fun=mean, geom="point", shape=15, size=4) + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.06))
ALL_genes
# Significant V non significant
SNP_heritability$significant_LRT <- if_else(SNP_heritability$Padj < 0.05, "Sig.", "Non-sig.")
table(SNP_heritability$significant_LRT)

ADD_Sig_V_nonsig = ggplot(data = SNP_heritability, aes(x = significant_LRT, y = h2, fill = significant_LRT)) + geom_boxplot(outlier.colour = NA, alpha = 0.7) +
theme_bw() +
              labs(x = "",
              y = "Cis-heritability") +
              scale_fill_manual(values = c("#d95f02","#7570b3")) +
                theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14), # Rotate x-axis text
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16), # Adjust x-axis title
    axis.title.y = element_text(size = 16),# Adjust legend text size
  ) + guides(fill = "none") + stat_summary(fun=mean, geom="point", shape=15, size=4) + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.06))


table(SNP_heritability$significant_LRT) # Non-sig.     Sig.    7806     6757 


SNP_heritability %>% filter(significant_LRT == "Sig.") %>% summarise(median(h2)) # 0.412539
SNP_heritability %>% filter(significant_LRT == "Non-sig.") %>% summarise(median(h2)) #0.0401925
SNP_heritability %>% summarize(mean(h2)) # 0.2416252



# Dominant heritability
hsq_files = list.files(path = path, pattern = "\\_add_domi\\.hsq$", full.names = TRUE)
SNP_heritability_dom = data.frame(matrix(nrow = 0, ncol = 7))
colnames(SNP_heritability_dom) <- c("Gene", "Vg_dom", "h2_dom", "h2_dom_add", "h2_together", "h2_dom_SE", "Pval_dom_component")
for (file in hsq_files) {
    file_temp = gsub("//", "/", file)
    hsq = data.table::fread(file_temp, fill = TRUE)
    Vg_dom = hsq$Variance[2]
    h2_dom_add = hsq$Variance[5]
    h2_together = hsq$Variance[8] 
    h2_dom = hsq$Variance[6]
    SE_dom = hsq$SE[6]
    P_dom = hsq$Variance[13]
    file <- gsub(paste0(path, "/"), "", file)
    gene_name <- strsplit(file, "_")[[1]][1]
    row = c(gene_name, Vg_dom, h2_dom, h2_dom_add, h2_together, SE_dom, P_dom)
    SNP_heritability_dom = rbind(SNP_heritability_dom, row)
}

colnames(SNP_heritability_dom) <- c("Gene", "Vg_dom", "h2_dom", "h2_dom_add", "h2_together", "h2_dom_SE", "Pval_dom_component")



SNP_heritability_dom$h2_dom <- as.numeric(SNP_heritability_dom$h2_dom)
SNP_heritability_dom$Pval_dom <- as.numeric(SNP_heritability_dom$Pval_dom)
SNP_heritability_dom$h2_dom_SE <- as.numeric(SNP_heritability_dom$h2_dom_SE)
SNP_heritability_dom$Padj_dom <- p.adjust(SNP_heritability_dom$Pval_dom, method = "BH")


SNP_heritability <- left_join(SNP_heritability, SNP_heritability_dom)
SNP_heritability_signif_add_both <- SNP_heritability %>% filter(Padj < 0.05)


SNP_heritability_signif_add_both <- SNP_heritability_signif_add_both %>% filter(!is.na(h2_dom))
dim(SNP_heritability_signif_add_both) # 5863 genes




# Reshape data for plotting
data_long <- reshape2::melt(SNP_heritability_signif_add_both, id.vars = "Gene", 
                            measure.vars = c("h2_dom_add", "h2_dom", "h2_together"))
data_long$value <- as.numeric(data_long$value)
data_long$variable <- factor(data_long$variable, levels = c("h2_together", "h2_dom_add", "h2_dom"), labels = c("Combined", "Additive", "Dominant"))


Add_V_dom = ggplot(data_long, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(outlier.colour =  NA, alpha = 0.7) +
  scale_fill_manual(values = c("#e7298a", "#66a61e", "#e6ab02")) +
  labs(x = "",
       y = "Cis_heritability") +
  theme_bw() +
              labs(x = "",
              y = "Cis-heritability") +
                theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14), # Rotate x-axis text
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16), # Adjust x-axis title
    axis.title.y = element_text(size = 16),# Adjust legend text size
  ) + guides(fill = "none") + stat_summary(fun=mean, geom="point", shape=15, size=4) + scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.06))
Add_V_dom

data_long %>% filter(variable == "Dominant") %>% summarise(mean(value))


SNP_heritability_signif_add_both %>% filter(Padj_dom < 0.05) %>% select(Gene, Padj_dom)
# Significant ones
#ENSBTAG00000015066 0.01009680
#ENSBTAG00000031825 0.02669096

SNP_heritability_signif_add_both %>% filter(Padj_dom < 0.05)

plot_grid(ALL_genes, ADD_Sig_V_nonsig,Add_V_dom, labels = c('a', 'b', 'c'), nrow =1, ncol = 3, label_size = 20)
ggsave(args[2], width = 12, height = 10, dpi = 600)
