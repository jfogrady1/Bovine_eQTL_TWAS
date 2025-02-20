######################################################
# Plot ieQTL results for O'Grady et al., 2024
#####################################################

library(data.table)
library(tidyverse)
library(vcfR)
args = commandArgs(trailingOnly = TRUE)
# First thing is to read in the data for top eQTLs
data <- fread(args[1]) %>% filter(pval_adj_bh < 0.25)
head(data)


vcf_all <- vcfR::read.vcfR()



vcf_all <- cbind(vcf_all@fix, vcf_all@gt)
vcf_all <- as.data.frame(vcf_all)
head(vcf_all)
vcf_all <- vcf_all %>% filter(ID %in% data$variant_id)


vcf_all[,colnames(vcf_all)[10:132]] <- lapply(vcf_all[,colnames(vcf_all)[10:132]], function(x) sub("0\\|0:.*", paste0("0"), x))
vcf_all[,colnames(vcf_all)[10:132]] <- lapply(vcf_all[,colnames(vcf_all)[10:132]], function(x) sub("0\\|1:.*", paste0("1"), x))
vcf_all[,colnames(vcf_all)[10:132]] <- lapply(vcf_all[,colnames(vcf_all)[10:132]], function(x) sub("1\\|0:.*", paste0("1"), x))
vcf_all[,colnames(vcf_all)[10:132]] <- lapply(vcf_all[,colnames(vcf_all)[10:132]], function(x) sub("1\\|1:.*", paste0("2"), x))

head(vcf_all)
rownames(vcf_all) <- vcf_all$ID
vcf_all <- vcf_all %>% select(-c(1,2,3,4,5,6,7,8,9))
head(vcf_all)
dim(vcf_all)
rownames(vcf_all)

# Read in expression data set
expression = fread(args[3]) %>% select(-c(1:4,6))
head(expression)
# read in expression data
rownames(expression) <- expression$gid
library(ggdist)

eQTL_plot <- function(gene_id, SNP_id, gene_name, vcf, counts, HOM, HET, HOM_ALT) {

  counts_gene <- counts %>% filter(gid == gene_id) %>% select(-c(gid)) %>% t() %>% as.data.frame() 
  
  colnames(counts_gene) <- "Expression"

  
  counts_gene$Sample <- rownames(counts_gene)
  counts_gene$Group <- c(rep("bTB-", 63), rep("bTB+", 60))
  counts_gene$Group <- as.factor(counts_gene$Group)
  
  
  print("HERE")
  
  vcf_temp <- vcf[SNP_id,]
  vcf_temp <- pivot_longer(vcf_temp, cols = colnames(vcf_temp), names_to = "Sample", values_to = "Genotype")
  print(head(vcf_temp))
  vcf_temp <- vcf_temp %>% mutate(Genotype = case_when(Genotype == "0" ~ HOM,
                                          Genotype == "1" ~ HET,
                                          Genotype == "2" ~ HOM_ALT))


  vcf_temp$Genotype <- factor(vcf_temp$Genotype, levels = c(HOM, HET, HOM_ALT))
  counts_gene <- left_join(counts_gene, vcf_temp, by = c("Sample" = "Sample"))
  print(colnames(counts_gene))
  counts_gene$density = paste0(counts_gene$Genotype, ":", counts_gene$Group)
  counts_gene$Group <- factor(counts_gene$Group, levels = c("bTB-", "bTB+"), labels = c("bTB-", "bTB+"))
  print(table(counts_gene$Genotype))
  plot_dom <- ggplot(counts_gene, aes(y = Expression, x = Genotype, shape = Group, colour = Group)) + 
    stat_halfeye(data = subset(counts_gene, Genotype == HOM & Group =="bTB-"), aes(x = Genotype, y = Expression), slab_colour =  "#2166ac", slab_fill = "#2166ac", adjust = .33, width = .33, alpha=0.5, position = position_nudge(x=0), side = "left") +
    stat_halfeye(data = subset(counts_gene, Genotype == HOM & Group =="bTB+"), aes(x = Genotype, y = Expression), slab_colour =  "#b2182b", slab_fill = "#b2182b", adjust = .33, width = .33, alpha=0.5, position = position_nudge(x=0), side = "right") +
    stat_halfeye(data = subset(counts_gene, Genotype == HET & Group =="bTB-"), aes(x = Genotype, y = Expression), slab_colour =  "#2166ac", slab_fill = "#2166ac", adjust = .33, width = .33, alpha=0.5, position = position_nudge(x=0), side = "left") +
    stat_halfeye(data = subset(counts_gene, Genotype == HET & Group =="bTB+"), aes(x = Genotype, y = Expression), slab_colour = "#b2182b", slab_fill = "#b2182b", adjust = .33, width = .33, alpha=0.5, position = position_nudge(x=0), side = "right") +
    stat_halfeye(data = subset(counts_gene, Genotype == HOM_ALT & Group =="bTB-"), aes(x = Genotype, y = Expression), slab_colour =  "#2166ac", slab_fill = "#2166ac", adjust = .33, width = .33, alpha=0.5, position = position_nudge(x=0), side = "left") +
    stat_halfeye(data = subset(counts_gene, Genotype == HOM_ALT & Group =="bTB+"), aes(x = Genotype, y = Expression), slab_colour = "#b2182b", slab_fill = "#b2182b", adjust = .33, width = .33, alpha=0.5, position = position_nudge(x=0), side = "right") +
    geom_point(size = 3) +
    theme_bw() + xlab(SNP_id) + ylab(paste0("Residualised expression of ", gene_name)) + labs(fill = SNP_id) +
    scale_x_discrete(labels =  c(HOM, HET, HOM_ALT)) +
    scale_colour_manual(values = c("#2166ac","#b2182b")) +
    geom_smooth(method = "lm", aes(group = Group, colour = Group), se = FALSE) +
    theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
    axis.title.y = element_text(size = 21, color = "black"),
    axis.title.x = element_text(size = 21, color = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15))

  return(plot_dom)


  
}


GBP4 = eQTL_plot(gene_id = "ENSBTAG00000014529", SNP_id = "3:53539162:G:T",  gene_name = "GBP4", vcf = vcf_all, counts = expression,  HOM = "GG", HET = "GT", HOM_ALT = "TT")
PCBP2 = eQTL_plot(gene_id = "ENSBTAG00000020757", SNP_id = "5:25861761:G:A",  gene_name = "PCBP2", vcf = vcf_all, counts = expression,  HOM = "GG", HET = "GA", HOM_ALT = "AA")
RELT = eQTL_plot(gene_id = "ENSBTAG00000016494", SNP_id = "15:53330033:T:C",  gene_name = "RELT", vcf = vcf_all, counts = expression,  HOM = "TT", HET = "TC", HOM_ALT = "CC")
PLD3 = eQTL_plot(gene_id = "ENSBTAG00000019150", SNP_id = "18:49353483:G:C",  gene_name = "PLD3", vcf = vcf_all, counts = expression,  HOM = "GG", HET = "GC", HOM_ALT = "CC")
PNPLA1 = eQTL_plot(gene_id = "ENSBTAG00000015118", SNP_id = "23:10034540:C:T",  gene_name = "PNPAL1", vcf = vcf_all, counts = expression,  HOM = "CC", HET = "CT", HOM_ALT = "TT")


GBP4
ggsave(args[4], width = 12, height = 12, dpi = 600)
PCBP2
ggsave(args[5], width = 12, height = 12, dpi = 600)
RELT
ggsave(args[6], width = 12, height = 12, dpi = 600)
PLD3
ggsave(args[7], width = 12, height = 12, dpi = 600)
PNPLA1
ggsave(args[8], width = 12, height = 12, dpi = 600)
