library(tidyverse)
library(data.table)
library(ggplot2)
library(ggpubr)
library(ggsignif)

args = commandArgs(trailing = T)
ALL_afc <- fread(args[1])
CON_afc <- fread(args[2])
INFEC_afc <- fread(args[3])

head(ALL_afc)
ALL_eQTL <- read.table(args[4], header = T) %>% filter(is_eGene == TRUE) 
CON_eQTL <- read.table(args[5], header = T) %>% filter(is_eGene == TRUE) 
INFEC_eQTL <- read.table(args[6], header = T) %>% filter(is_eGene == TRUE) 


ALL_afc <- ALL_afc %>% filter(ALL_afc$phenotype_id %in% ALL_eQTL$phenotype_id)
CON_afc <- CON_afc %>% filter(CON_afc$phenotype_id %in% CON_eQTL$phenotype_id)
INFEC_afc <- INFEC_afc %>% filter(INFEC_afc$phenotype_id %in% INFEC_eQTL$phenotype_id)


ALL_afc$group <- "ALL"
CON_afc$group <- "CONTROL"
INFEC_afc$group <- "INFECTED"




afc_plot <- rbind(ALL_afc, CON_afc, INFEC_afc)


wilcox.test(abs(CON_afc$log2_aFC), abs(INFEC_afc$log2_aFC))
comparisons = list(c("ALL", "CONTROL"), c("ALL", "INFECTED"), c("CONTROL", "INFECTED"))
ggplot(data = afc_plot, aes(x = group, y = abs(log2_aFC), fill = group)) + 
geom_boxplot() +
scale_fill_manual(breaks = c("ALL", "CONTROL", "INFECTED"), values = c("#542788","#2166ac", "#b2182b")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black")) +
        stat_signif(comparisons = comparisons,
                    test = "wilcox.test",
                    y_position = c(6.75,7,6.5), # Set custom y positions for bars
                    tip_length = 0.01)  # Adjust the tip length if needed

ggsave(args[7], width = 12, height = 12, dpi = 600)



ALL_af_plot <- ALL_afc %>% select(variant_id,af, slope) %>% mutate(group = "AAG")
CON_af_plot <- CON_afc %>% select(variant_id,af, slope) %>% mutate(group = "bTB-")
INFEC_af_plot <- INFEC_afc %>% select(variant_id,af, slope) %>% mutate(group = "bTB+")


af_afc_plot <- rbind(ALL_af_plot, CON_af_plot, INFEC_af_plot)
af_afc_plot$af2 <- 1 - af_afc_plot$af
af_afc_plot$maf <- with(af_afc_plot, pmin(af_afc_plot$af, af_afc_plot$af2))

af_afc_plot <- af_afc_plot %>% mutate(slope_group = ifelse(slope > 0, "Positive", "Negative"))

ggplot(af_afc_plot, aes(x = maf, y = slope, col = group, linetype = slope_group)) + 
geom_point(size = 2, alpha = 0.2) + 
geom_smooth(method = "loess", se = FALSE, linewidth = 2.5) +
scale_colour_manual(breaks = c("AAG", "bTB-", "bTB+"), values = c("#542788","#2166ac", "#b2182b")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black")) +coord_flip()
ggsave(args[8], width = 12, height = 12, dpi = 600)

median(abs(ALL_afc$log2_aFC)) # 0.2998068
median(abs(CON_afc$log2_aFC)) # 0.4294106
median(abs(INFEC_afc$log2_aFC)) # 0.5384912



ggplot(data = afc_plot, aes(x = slope,  y = log2_aFC, colour = group)) + 
geom_point(alpha = 0.25) +
scale_colour_manual(breaks = c("ALL", "CONTROL", "INFECTED"), values = c("#542788","#2166ac", "#b2182b")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black")) +
        geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 2)

ggsave(args[9], width = 12, height = 12, dpi = 600)

cor.test(ALL_afc$log2_aFC, ALL_afc$slope, method = "spearman", exact = FALSE) # 0.9320916
cor.test(CON_afc$log2_aFC, CON_afc$slope, method = "spearman", exact = FALSE) # 0.91
cor.test(INFEC_afc$log2_aFC, INFEC_afc$slope, method = "spearman", exact = FALSE) #0.90


#### Numbers of large eQTL
large_ALL <- ALL_afc %>% filter(abs(log2_aFC) >= 1) %>% dim() # 1098
large_ALL
large_ALL[1] / length(ALL_eQTL$phenotype_id) # 0.1644697


large_CON <- CON_afc %>% filter(abs(log2_aFC) >= 1) %>% dim() # 799
large_CON
large_CON[1] / length(CON_eQTL$phenotype_id) # 0.2404212

large_INFEC <- INFEC_afc %>% filter(abs(log2_aFC) >= 1) %>% dim() # 670
large_INFEC[1] / length(INFEC_eQTL$phenotype_id) # 0.2996422



CON_eQTL <- read.table(args[5], header = T) 
INFEC_eQTL <- read.table(args[6], header = T)

CON_eQTL_common <- CON_eQTL %>% filter(phenotype_id %in% INFEC_eQTL$phenotype_id)
INFEC_eQTL_common <- INFEC_eQTL %>% filter(phenotype_id %in% CON_eQTL_common$phenotype_id)


colnames(INFEC_eQTL_common) <- paste0("INFEC_", colnames(INFEC_eQTL_common))




# Need to get the most significant eQTL for a gene in either group and then get the corresponding eVariant in the other group

#1. join together - eQTLs
CON_INFEC <- left_join(CON_eQTL_common, INFEC_eQTL_common, by = c("phenotype_id" = "INFEC_phenotype_id"))


#2. Read in the nominal ones and get the eVariants
CONTROL_nominal <- fread(args[10])
INFECTED_nominal <- fread(args[11])
CONTROL_nominal <- CONTROL_nominal %>% filter(phenotype_id != "phenotype_id")
INFECTED_nominal <- INFECTED_nominal %>% filter(phenotype_id != "phenotype_id")


data_CONTROL_egene <- CON_eQTL_common #%>% filter(is_eGene == "TRUE")
data_INFECTED_egene <- INFEC_eQTL_common #%>% filter(is_eGene == "TRUE")

data_CONTROL <- left_join(data_CONTROL_egene, CONTROL_nominal, by = c("phenotype_id" = "phenotype_id"))
data_INFECTED <- left_join(data_INFECTED_egene, INFECTED_nominal, by = c("INFEC_phenotype_id" = "phenotype_id"))

# Here we determine whether a variant is an eVariant
print(head(data_CONTROL))
data_CONTROL <- data_CONTROL %>% mutate(CONTROL_is_eVariant = case_when(
    (as.numeric(pval_nominal.y) < as.numeric(pval_nominal_threshold)) & is_eGene == TRUE ~ "YES",
    (as.numeric(pval_nominal.y) > as.numeric(pval_nominal_threshold)) & is_eGene == TRUE ~ "NO",
    is_eGene == FALSE ~ "NO"))

data_INFECTED <- data_INFECTED %>%  mutate(INFEC_is_eVariant = case_when(
    as.numeric(pval_nominal) < as.numeric(INFEC_pval_nominal_threshold) & INFEC_is_eGene == TRUE ~ "YES",
    as.numeric(pval_nominal) > as.numeric(INFEC_pval_nominal_threshold) & INFEC_is_eGene == TRUE ~ "NO",
    INFEC_is_eGene == FALSE ~ "NO"))


data_CONTROL <- data_CONTROL %>% select(1,22,27,28,29,30)
colnames(data_CONTROL) <- c("CONTROL_phenotype_id", "CONTROL_variant_id", "CONTROL_pval_nominal", "CONTROL_slope", "CONTROL_slope_se", "CONTROL_is_eVariant")

data_INFECTED <- data_INFECTED %>% select(1,22,27,28,29,30)
colnames(data_INFECTED) <- c("INFECTED_phenotype_id", "INFECTED_variant_id", "INFECTED_pval_nominal", "INFECTED_slope", "INFECTED_slope_se", "INFECTED_is_eVariant")

data_CONTROL$CONTROL_variant_id <- paste0(data_CONTROL$CONTROL_variant_id, "-", data_CONTROL$CONTROL_phenotype_id, "-INFECTED") # will be used for indexing
data_INFECTED$INFECTED_variant_id <- paste0(data_INFECTED$INFECTED_variant_id, "-", data_INFECTED$INFECTED_phenotype_id, "-CONTROL") # note opposite as we will need to merge these to the other file


# Index eQTL file
head(CON_INFEC)
CON_INFEC$variant_id <- paste0(CON_INFEC$variant_id,"-",CON_INFEC$phenotype_id, "-CONTROL") # will be used for indexing
CON_INFEC$INFEC_variant_id <- paste0(CON_INFEC$INFEC_variant_id,"-",CON_INFEC$phenotype_id, "-INFECTED")

print(colnames(CON_INFEC))
# Get variant for gene with lowest P-value
CON_INFEC <- CON_INFEC %>% group_by(phenotype_id) %>% summarize(pval_min = min(INFEC_pval_adj_BH, pval_adj_BH),
                                                   variant_id = case_when(
                                                    INFEC_pval_adj_BH < pval_adj_BH ~ INFEC_variant_id,
                                                    INFEC_pval_adj_BH > pval_adj_BH ~ variant_id
                                                   ),
                                                   slope = case_when(
                                                    INFEC_pval_adj_BH < pval_adj_BH ~ INFEC_slope,
                                                    INFEC_pval_adj_BH > pval_adj_BH ~ slope
                                                   ),
                                                   top_group = case_when(
                                                    INFEC_pval_adj_BH < pval_adj_BH ~ "INFECTED",
                                                    INFEC_pval_adj_BH > pval_adj_BH ~ "CONTROL"
                                                   ),
                                                   CONTROL_is_eVariant_top = case_when(
                                                    is_eGene == TRUE ~ "YES",
                                                    is_eGene == FALSE ~ "NO"
                                                   ),
                                                   INFECTED_is_eVariant_top = case_when(
                                                    INFEC_is_eGene == TRUE ~ "YES",
                                                    INFEC_is_eGene == FALSE ~ "NO"
                                                   ))


CON_INFEC_merged <- left_join(CON_INFEC, data_CONTROL, by = c("variant_id" = "CONTROL_variant_id"))
CON_INFEC_merged <- left_join(CON_INFEC_merged, data_INFECTED, by = c("variant_id" = "INFECTED_variant_id"))




CON_INFEC_merged$INFECTED_is_eVariant <- if_else(is.na(CON_INFEC_merged$INFECTED_is_eVariant), CON_INFEC_merged$INFECTED_is_eVariant_top, CON_INFEC_merged$INFECTED_is_eVariant)
CON_INFEC_merged$CONTROL_is_eVariant <- if_else(is.na(CON_INFEC_merged$CONTROL_is_eVariant), CON_INFEC_merged$CONTROL_is_eVariant_top, CON_INFEC_merged$CONTROL_is_eVariant)



CON_INFEC_merged$CONTROL_slope <- as.numeric(CON_INFEC_merged$CONTROL_slope)
CON_INFEC_merged$INFECTED_slope <- as.numeric(CON_INFEC_merged$INFECTED_slope)


CON_INFEC_merged_final <- CON_INFEC_merged %>% mutate(CONTROL_slope_plot = case_when
                           (top_group == "CONTROL" ~ slope, 
                           top_group != "CONTROL" ~ CONTROL_slope))


CON_INFEC_merged_final <- CON_INFEC_merged_final %>% mutate(INFECTED_slope_plot = case_when
                           (top_group == "INFECTED" ~ as.numeric(slope), 
                           top_group != "INFECTED" ~ as.numeric(INFECTED_slope)))






CON_INFEC_merged_final <- CON_INFEC_merged_final %>% mutate(direction = case_when(
    CONTROL_slope_plot >= 0 & INFECTED_slope_plot >= 0 ~ "Consistent",
    CONTROL_slope_plot <= 0 & INFECTED_slope_plot <= 0 ~ "Consistent",
    TRUE ~ "Discordant"
))

table(CON_INFEC_merged_final$direction)
table(CON_INFEC_merged_final$INFECTED_is_eVariant)
dim(CON_INFEC_merged_final)

CON_INFEC_merged_final <- CON_INFEC_merged_final %>% mutate(grouping = case_when(
    CONTROL_is_eVariant == "YES" & INFECTED_is_eVariant == "YES" & direction == "Consistent" ~ "Shared",
    CONTROL_is_eVariant == "NO" & INFECTED_is_eVariant == "YES" ~ "Reactor_specific",
    CONTROL_is_eVariant == "YES" & INFECTED_is_eVariant == "NO" ~ "Control_specific",
    CONTROL_is_eVariant == "YES" & INFECTED_is_eVariant == "YES" & direction == "Discordant" ~ "Discordant",
    CONTROL_is_eVariant == "NO" & INFECTED_is_eVariant == "NO" ~ "not significant"

))



CON_INFEC_merged_final <- CON_INFEC_merged_final %>% mutate(size = case_when(
    grouping == "Shared" ~ 0.25,
    grouping == "Reactor_specific" ~ 0.25,
    grouping == "Control_specific" ~ 0.25,
    grouping == "Discordant" ~ 0.5,
    grouping == "not significant" ~ 0.2
))
colnames(CON_INFEC_merged_final)
CON_INFEC_merged_final$grouping <- factor(CON_INFEC_merged_final$grouping, levels = c("not significant", "Shared", "Control_specific", "Reactor_specific", "Discordant"),
labels = c("Not significant", "Shared", "Control magnified", "Reactor magnified", "Discordant"))

CON_INFEC_merged_final <- CON_INFEC_merged_final[order(CON_INFEC_merged_final$grouping),]
head(CON_INFEC_merged_final)

CON_INFEC_merged_final <- CON_INFEC_merged_final %>% filter(!is.na(CONTROL_slope_plot)) # 14479 variants tested in both
CON_INFEC_merged_final <- CON_INFEC_merged_final %>% filter(!is.na(INFECTED_slope_plot)) # 14479 variants tested in both
length(unique(CON_INFEC_merged_final$phenotype_id)) # 13871 variants tested in both
table(CON_INFEC_merged_final$grouping)
#  Not significant            Shared Control magnified Reactor magnified 
 #           10040              1143              2011               677 
  #     Discordant 
   #             0 




CON_INFEC_merged_final$size <- factor(CON_INFEC_merged_final$size)
CON_INFEC_merged_final$size

table(CON_INFEC_merged_final$grouping, CON_INFEC_merged_final$direction)
ggplot(data = CON_INFEC_merged_final, aes(x = CONTROL_slope_plot, y = INFECTED_slope_plot, col = grouping, size = size)) + geom_point() +
theme_bw() +
labs(colour = "eGene Category", y = "Effect size (slope) in Reactor cohort", x = "Effect size (slope) in Control cohort") +
scale_colour_manual(values = c("grey", "#9b6833","#2166ac", "#b2182b", "black")) +
scale_size_manual(values = c(0.5,2,6)) +
scale_y_continuous(limits = c(-3.5,3.5)) +
scale_x_continuous(limits = c(-3.5,3.5)) +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
    axis.title.y = element_text(size = 21, color = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15))
ggsave(args[12], width = 12, height = 12, dpi = 600)

rm(data_CONTROL)
rm(data_INFECTED)





# Check to see if they are overrepresented in genes

# DE data




DE_results_background <- fread(args[13]) %>% select(1) %>% as.vector()





CON_INFEC_merged_final_ORA <- CON_INFEC_merged_final %>% select(phenotype_id, grouping) 
dim(CON_INFEC_merged_final_ORA)
length(unique(CON_INFEC_merged_final_ORA$phenotype_id))


all_genes <- intersect(DE_results_background$V1, CON_INFEC_merged_final_ORA$phenotype_id)

DE <- fread(args[13]) %>% filter(padj < 0.05) %>% filter(V1 %in% all_genes) %>% select(1) %>% as.vector()


CON_INFEC_merged_final_ORA <- CON_INFEC_merged_final_ORA %>% filter(phenotype_id %in% all_genes)



table(CON_INFEC_merged_final_ORA$grouping)

Reactor_genes <- CON_INFEC_merged_final_ORA %>% filter(grouping == "Reactor magnified") %>% select(1) %>% filter(phenotype_id %in% all_genes)
Control_genes <- CON_INFEC_merged_final_ORA %>% filter(grouping == "Control magnified") %>% select(1) %>% filter(phenotype_id %in% all_genes)

head(Reactor_genes)

Reactor_genes



overlap <- length(intersect(Reactor_genes$phenotype_id, DE$V1))
only_regulated <- length(setdiff(Reactor_genes$phenotype_id, DE$V1))
only_de <- length(setdiff(DE$V1,Reactor_genes$phenotype_id))
neither <- length(setdiff(all_genes, union(Reactor_genes$phenotype_id, DE$V1)))

# Create a contingency table for reactor
# Count the number of overlaps and non-overlaps
contingency_table <- matrix(c(overlap, only_regulated, only_de, neither), 
                            nrow = 2, 
                            dimnames = list(c("Regulated", "Not_Regulated"),
                                            c("DE", "Not_DE")))

chisq.test_reactor <- chisq.test(contingency_table)
chisq.test_reactor



# Control
overlap <- length(intersect(Control_genes$phenotype_id, DE$V1))
only_regulated <- length(setdiff(Control_genes$phenotype_id, DE$V1))
only_de <- length(setdiff(DE$V1,Control_genes$phenotype_id))
neither <- length(setdiff(all_genes, union(Control_genes$phenotype_id, DE$V1)))

# Create a contingency table for reactor
# Count the number of overlaps and non-overlaps
contingency_table <- matrix(c(overlap, only_regulated, only_de, neither), 
                            nrow = 2, 
                            dimnames = list(c("Regulated", "Not_Regulated"),
                                            c("DE", "Not_DE")))

contingency_table
chisq.test_control <- chisq.test(contingency_table)
chisq.test_control





CON_INFEC_merged_final_ORA <- rbind(Reactor_genes, Control_genes)


overlap <- length(intersect(CON_INFEC_merged_final_ORA$phenotype_id, DE$V1))
only_regulated <- length(setdiff(CON_INFEC_merged_final_ORA$phenotype_id, DE$V1))
only_de <- length(setdiff(DE$V1,CON_INFEC_merged_final_ORA$phenotype_id))
neither <- length(setdiff(all_genes, union(CON_INFEC_merged_final_ORA$phenotype_id, DE$V1)))


# Create a contingency table for reactor
# Count the number of overlaps and non-overlaps
contingency_table <- matrix(c(overlap, only_regulated, only_de, neither), 
                            nrow = 2, 
                            dimnames = list(c("Regulated", "Not_Regulated"),
                                            c("DE", "Not_DE")))

contingency_table
chisq.test_all_magnified <- chisq.test(contingency_table)
chisq.test_all_magnified


p <- c(0.05, 0.5, 0.7)
p.adjust(p, method = "bonferroni")

chisq.test_reactor 
chisq.test_control
chisq.test_all_magnified


Reactor_de <- intersect(Reactor_genes$phenotype_id, DE$V1) %>% as.data.frame()
Control_de <- intersect(Control_genes$phenotype_id, DE$V1) %>% as.data.frame()

colnames(Reactor_de) <- "phenotype_id"
colnames(Control_de) <- "phenotype_id"










annotation <- read.csv(args[14]) %>% select(2,3)


annotation = fread("~/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf")
annotation <- annotation %>% filter(V3 == "gene")
head(annotation)
annotation <- annotation %>% separate(., V9, into = c("gene_id", "gene_version", "gene_name"), sep = ";")
annotation$gene_id <- gsub("^gene_id ", "", annotation$gene_id)
annotation$gene_id <- gsub('"', '', annotation$gene_id)

annotation$gene_name <- gsub("gene_name ", "", annotation$gene_name)
annotation$gene_name <- gsub("gene_source ", "", annotation$gene_name)
annotation$gene_name <- gsub('"', '', annotation$gene_name)
annotation$gene_name <- if_else(annotation$gene_name == " ensembl", annotation$gene_id, annotation$gene_name)
colnames(annotation)[1] <- "chr"
annotation <- annotation %>% dplyr::select(gene_id, gene_name)




Reactor_de <- left_join(Reactor_de, annotation, by = c("phenotype_id" = "gene_id"))
Control_de <- left_join(Control_de, annotation, by = c("phenotype_id" = "gene_id"))



Reactor_de$Symbol

CON_INFEC_merged_final_ORA <-left_join(CON_INFEC_merged_final_ORA, annotation, by = c("phenotype_id" = "gene_id"))



head(CON_INFEC_merged_final)
CON_INFEC_merged_final <- left_join(CON_INFEC_merged_final, annotation, by = c("phenotype_id" = "gene_id"))
colnames(CON_INFEC_merged_final)

CON_INFEC_merged_final %>% filter(grouping == "Reactor magnified")
head(CON_INFEC_merged_final)
library(vcfR)
vcf_control <- vcfR::read.vcfR(args[15])
vcf_infected <- vcfR::read.vcfR(args[16])

vcf_control <- cbind(vcf_control@fix, vcf_control@gt)
head(vcf_control)
vcf_control <- as.data.frame(vcf_control)

vcf_infected <- cbind(vcf_infected@fix, vcf_infected@gt)
head(vcf_infected)
vcf_infected <- as.data.frame(vcf_infected)

head(vcf_control)
vcf_control[,colnames(vcf_control)[10:72]] <- lapply(vcf_control[,colnames(vcf_control)[10:72]], function(x) sub("0\\|1:.*", paste0("1"), x))
vcf_control[,colnames(vcf_control)[10:72]] <- lapply(vcf_control[,colnames(vcf_control)[10:72]], function(x) sub("1\\|0:.*", paste0("1"), x))
vcf_control[,colnames(vcf_control)[10:72]] <- lapply(vcf_control[,colnames(vcf_control)[10:72]], function(x) sub("0\\|0:.*", paste0("0"), x))
vcf_control[,colnames(vcf_control)[10:72]] <- lapply(vcf_control[,colnames(vcf_control)[10:72]], function(x) sub("1\\|1:.*", paste0("2"), x))


vcf_infected[,colnames(vcf_infected)[10:69]] <- lapply(vcf_infected[,colnames(vcf_infected)[10:69]], function(x) sub("0\\|1:.*", paste0("1"), x))
vcf_infected[,colnames(vcf_infected)[10:69]] <- lapply(vcf_infected[,colnames(vcf_infected)[10:69]], function(x) sub("1\\|0:.*", paste0("1"), x))
vcf_infected[,colnames(vcf_infected)[10:69]] <- lapply(vcf_infected[,colnames(vcf_infected)[10:69]], function(x) sub("0\\|0:.*", paste0("0"), x))
vcf_infected[,colnames(vcf_infected)[10:69]] <- lapply(vcf_infected[,colnames(vcf_infected)[10:69]], function(x) sub("1\\|1:.*", paste0("2"), x))


# read in expression data
counts_control <- fread(args[17]) %>% select(-c(1:4,6)) #%>% as.matrix()
rownames(counts_control) <- counts_control$gid

counts_infected <- fread(args[18]) %>% select(-c(1:4,6)) #%>% as.matrix()

rownames(counts_infected) <- counts_infected$gid



# Get everything in order for regression
library(cowplot)
eQTL_plot <- function(gene_id, SNP_id, SNP_id2, gene_name, counts_control, counts_infected, vcf_control, vcf_infected, HOM, HET, HOM_ALT, HOM2, HET2, HOM_ALT2) {

  counts_control_gene <- counts_control %>% filter(gid == gene_id) %>% select(-c(gid)) %>% t() %>% as.data.frame() 
  counts_infected_gene <- counts_infected %>% filter(gid == gene_id) %>% select(-c(gid)) %>% t() %>% as.data.frame() 

  
  colnames(counts_control_gene) <- "Expression"
  colnames(counts_infected_gene) <- "Expression"

  
  counts_control_gene$Sample <- rownames(counts_control_gene)
  counts_infected_gene$Sample <- rownames(counts_infected_gene)
  
  
  
  print("HERE")
  
  vcf_temp <- vcf_control %>% filter(ID == SNP_id)
  vcf_temp <- pivot_longer(vcf_temp, cols = c(10:72), names_to = "Sample", values_to = "Genotype")
  vcf_temp_control <- vcf_temp %>% select(3,10,11)
  vcf_temp_control <- vcf_temp_control %>% mutate(Genotype = case_when(Genotype == "0" ~ HOM,
                                          Genotype == "1" ~ HET,
                                          Genotype == "2" ~ HOM_ALT))
  vcf_temp <- vcf_infected %>% filter(ID == SNP_id2)
  vcf_temp <- pivot_longer(vcf_temp, cols = c(10:69), names_to = "Sample", values_to = "Genotype")
  vcf_temp_infected <- vcf_temp %>% select(3,10,11)

  vcf_temp_infected <- vcf_temp_infected %>% mutate(Genotype = case_when(Genotype == "0" ~ HOM2,
                                          Genotype == "1" ~ HET2,
                                          Genotype == "2" ~ HOM_ALT2))

  vcf_temp_control$Genotype <- factor(vcf_temp_control$Genotype, levels = c(HOM, HET, HOM_ALT))
  vcf_temp_infected$Genotype <- factor(vcf_temp_infected$Genotype, levels = c(HOM2, HET, HOM_ALT2))
  counts_control_gene <- left_join(counts_control_gene, vcf_temp_control, by = c("Sample" = "Sample"))
  counts_infected_gene <- left_join(counts_infected_gene, vcf_temp_infected, by = c("Sample" = "Sample"))
  print(table(counts_control_gene$Genotype))
  print(table(counts_infected_gene$Genotype))
  p_control <- ggplot(counts_control_gene, aes(y = Expression, x = Genotype, fill = Genotype)) + 
    geom_boxplot(outlier.colour = NA, trim=FALSE, alpha = 0.5) +
    geom_jitter(shape=16, colour = "black", size = 3, position=position_jitter(0.2)) + 
    scale_fill_manual(values = c("#2166ac", "#2166ac", "#2166ac")) +
    theme_bw() + xlab(SNP_id) + ylab(paste0("Residualised expression of ", gene_name)) + labs(fill = SNP_id) +
    scale_x_discrete(labels =  c(HOM, HET, HOM_ALT))
  p_infected <- ggplot(counts_infected_gene, aes(y = Expression, x = Genotype, fill = Genotype)) + 
    geom_boxplot(outlier.colour = NA, trim=FALSE, alpha = 0.5) +
    geom_jitter(shape=16, colour = "black", position=position_jitter(0.2), size = 3) + 
    scale_fill_manual(values = c("#b2182b", "#b2182b", "#b2182b")) +
    theme_bw() + xlab(SNP_id2) + ylab(paste0("Residualised expression of ", gene_name)) + labs(fill = SNP_id) +
    scale_x_discrete(labels =  c(HOM2, HET2, HOM_ALT2))

  
  return(plot_grid(p_control, p_infected))
}



write.table(CON_INFEC_merged_final, args[19], sep = "\t", quote = F)

# IFITM3
IFITM3 <- eQTL_plot(gene_id = "ENSBTAG00000019015", SNP_id = "29:50294904:G:A", SNP_id2 = "29:50294904:G:A", gene_name = "IFITM3", counts_control = counts_control, counts_infected = counts_infected, vcf_control = vcf_control, vcf_infected = vcf_infected, HOM = "GG", HET = "GA", HOM_ALT = "AA", HOM2 = "GG", HET2 = "GA", HOM_ALT2 = "AA")
IFITM3
ggsave(args[20], width = 10, height = 10, dpi = 600)


IFI16 <- eQTL_plot(gene_id = "ENSBTAG00000011511", SNP_id = "3:10984726:G:A", SNP_id2 = "3:10984726:G:A", gene_name = "IFI16", counts_control = counts_control, counts_infected = counts_infected, vcf_control = vcf_control, vcf_infected = vcf_infected, HOM = "GG", HET = "GA", HOM_ALT = "AA", HOM2 = "GG", HET2 = "GA", HOM_ALT2 = "AA")
IFI16
ggsave(args[21], width = 10, height = 10, dpi = 600)






IL1R1 <- eQTL_plot(gene_id = "ENSBTAG00000005273", SNP_id = "11:6851248:C:G", SNP_id2 = "11:6851248:C:G", gene_name = "IL1R1", counts_control = counts_control, counts_infected = counts_infected, vcf_control = vcf_control, vcf_infected = vcf_infected, HOM = "CC", HET = "CG", HOM_ALT = "GG", HOM2 = "CC", HET2 = "CG", HOM_ALT2 = "GG")
IL1R1
ggsave(args[22], width = 10, height = 10, dpi = 600)

RGS10 <- eQTL_plot(gene_id = "ENSBTAG00000005273", SNP_id = "11:6851248:C:G", SNP_id2 = "11:6851248:C:G", gene_name = "IL1R1", counts_control = counts_control, counts_infected = counts_infected, vcf_control = vcf_control, vcf_infected = vcf_infected, HOM = "CC", HET = "CG", HOM_ALT = "GG", HOM2 = "CC", HET2 = "CG", HOM_ALT2 = "GG")
ggsave(args[23], width = 10, height = 10, dpi = 600)



