#Enrichment analysis of each Group and also DE genes
library(data.table)
library(qvalue)
library(tidyverse)
library(ggrepel)
library(gprofiler2)
args = commandArgs(trailingOnly = TRUE)
args[1] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt"
args[2] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_fdr0.05.txt"
args[3] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_fdr0.05.txt"
#args[1-3]
ALL_perm <- read.table(args[1], header = T) %>% 
filter(is_eGene == T) %>% select(1,20)
CONTROL_perm <- read.table(args[2], header = T) %>% 
filter(is_eGene == T) %>% select(1,20)
colnames(CONTROL_perm)[2] <- "Padj_Control"
INFECTED_perm <- read.table(args[3], header = T) %>% 
filter(is_eGene == T) %>% select(1,20)
colnames(INFECTED_perm)[2] <- "Padj_INFECTED"
listInput <- list(All = ALL_perm$phenotype_id, Control = CONTROL_perm$phenotype_id, Infected = INFECTED_perm$phenotype_id)
ALL_common <- intersect(intersect(listInput$All, listInput$Infected), listInput$Control)

# Get those eGenes common in bovine and ALL
ALL_infected <- intersect(listInput$All, listInput$Infected)
ALL_infected <- setdiff(ALL_infected, listInput$Control)

# Get eGenes common in Control and ALL
ALL_control <- intersect(listInput$All, listInput$Control)
ALL_control <- setdiff(ALL_control, listInput$Infected)
length(ALL_control)

# Now get everything into a data frame
# get p values in control group
control_all <- CONTROL_perm %>% filter(phenotype_id %in% ALL_common)

# Now get p-values in infected group
infected_all <- INFECTED_perm %>% filter(phenotype_id %in% ALL_common)

# Now join both together
# Here, now have p values in both for genes which were common and significant
plot_df <- left_join(control_all, infected_all)
plot_df$Type <- "All Groups"





# Now focus on the infected group
# read in dataframes again to get non-significant p values
CONTROL_perm <- read.table(args[2], header = T) %>% select(1,20)
colnames(CONTROL_perm)[2] <- "Padj_Control"
INFECTED_perm <- read.table(args[3], header = T) %>% select(1,20)
colnames(INFECTED_perm)[2] <- "Padj_INFECTED"
control_infected_all_subset <- CONTROL_perm %>% filter(phenotype_id %in% ALL_infected)
ALL_infected <- intersect(ALL_infected, control_infected_all_subset$phenotype_id)
infected_all_subset <- INFECTED_perm %>% filter(phenotype_id %in% ALL_infected)

# now join together
plot_df_infected <- left_join(control_infected_all_subset, infected_all_subset)
plot_df_infected$Type <- "Reactor only"

# Now focus on the control group 
infected_control_all_subset <- INFECTED_perm %>% filter(phenotype_id %in% ALL_control)
ALL_control <- intersect(ALL_control, infected_control_all_subset$phenotype_id)
control_all_subset <- CONTROL_perm %>% filter(phenotype_id %in% ALL_control)

# now join together
plot_df_control <- left_join(infected_control_all_subset, control_all_subset)
plot_df_control$Type <- "Control only"

dim(plot_df) # 1603 - will be used to bind eferything
dim(plot_df_infected) # 509
dim(plot_df_control) # 1593

plot_df_gg <- rbind(plot_df, plot_df_control, plot_df_infected)
ggplot(plot_df_gg, aes(x = -log10(Padj_INFECTED), y = -log10(Padj_Control), colour = Type)) + 
geom_point(alpha = 0.5) +
geom_vline(xintercept = -log10(0.05), col = "black", linetype="dashed", size = 0.5) + geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed", size = 0.5)  +
labs(x=expression(-log[10](italic(P)[adj])~"Reactor group"),
       y=expression(-log[10](italic(P)[adj])~"Control group"))+
scale_colour_manual(breaks = c("All Groups", "Control only", "Reactor only"), values = c("#542788","#2166ac", "#b2182b")) +
theme_bw() +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"))
ggsave("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/cis_comparison_reactor_control.pdf", width = 12, height = 12, dpi = 600)


# Get intersections
ALL_perm <- read.table(args[1], header = T) %>% select(1,20)
CONTROL_perm <- read.table(args[2], header = T) %>% select(1,20)
INFECTED_perm <- read.table(args[3], header = T) %>% select(1,20)

# INFECTED_ONLY
background_infected <- intersect(ALL_perm$phenotype_id, INFECTED_perm$phenotype_id)

print("Number of input genes for infected only")
print(length(background_infected)) # 14570

target_genes_infected <- plot_df_infected[order(plot_df_infected$Padj_INFECTED, decreasing = F),]



#target <- gorth(query = plot_df_infected$phenotype_id, source_organism = "btaurus",  # Note, rownames of the res_df are our target
#target_organism = "hsapiens", mthreshold = 1, filter_na = TRUE,
#numeric_ns = "ENTREZGENE_ACC")

print("Number of input genes for INFECTED only")
print(length(target_genes_infected$phenotype_id)) # 509

target_genes_infected <- list(target_genes_infected$phenotype_id)

results_infected <- gost(query = target_genes_infected,organism = "btaurus", correction_method = "fdr", ordered_query = T, domain_scope = "known", custom_bg = background_infected, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC"), evcodes = T)
results_infected$result$Group = "Reactor Only"
results_infected$result



# Infected_ALL  - Enrichment
background_infected <- intersect(ALL_perm$phenotype_id, INFECTED_perm$phenotype_id)

#background_infected <- gorth(query = background_infected, source_organism = "btaurus", 
 #       target_organism = "hsapiens", mthreshold = 1, filter_na = TRUE,
  #      numeric_ns = "ENTREZGENE_ACC")
print("Number of background genes for INFECTED ALL")
print(length(background_infected)) # 14570


plot_df_infected_all <- rbind(plot_df, plot_df_infected)
plot_df_infected_all <- plot_df_infected_all %>% filter(Padj_INFECTED < 0.05)
target_genes_infected <- plot_df_infected_all[order(plot_df_infected_all$Padj_INFECTED, decreasing = F),]
head(target_genes_infected)


#arget <- gorth(query = plot_df_infected_all$phenotype_id, source_organism = "btaurus",  # Note, rownames of the res_df are our target
#target_organism = "hsapiens", mthreshold = 1, filter_na = TRUE,
#numeric_ns = "ENTREZGENE_ACC")

print("Number of input genes for INFECTED ALL")
print(length(target_genes_infected$phenotype_id)) #2112

target_genes_infected <- list(target_genes_infected$phenotype_id)
results_infected_all <- gost(query = target_genes_infected,organism = "btaurus", correction_method = "fdr", ordered_query = T, domain_scope = "known", custom_bg = background_infected, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC", "WP"), evcodes = T)
results_infected_all$result$Group = "Reactor and All"
results_infected_all$result
# Control_ALL - ORA
background_control <- intersect(ALL_perm$phenotype_id, CONTROL_perm$phenotype_id)
print(length(background_control)) # 14530


plot_df_control_all <- rbind(plot_df, plot_df_control)
plot_df_control_all <- plot_df_control_all %>% filter(Padj_Control < 0.05)

target_genes_control_all <- plot_df_control_all[order(plot_df_control_all$Padj_Control, decreasing = F),]




print("Number of input genes for control ALL")
print(length(target_genes_control_all$phenotype_id)) # 3196

target_genes_control_all <- list(target_genes_control_all$phenotype_id)
results_control_all <- gost(query = target_genes_control_all,organism = "btaurus", correction_method = "fdr", ordered_query = T, domain_scope = "known", custom_bg = background_control, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC"), evcodes = T)
results_control_all$result$Group = "Control and all"
results_control_all$result

# Control_only
background_control <- intersect(ALL_perm$phenotype_id, CONTROL_perm$phenotype_id)


print("Number of background genes for control only")
print(length(background_control)) # 14530

target_genes_control <- plot_df_control[order(plot_df_control$Padj_Control, decreasing = F),]




print("Number of genes for control only")
print(length(target_genes_control$phenotype_id)) # 1593

target_genes_control <- list(target_genes_control)
results_control <- gost(query = target_genes_control,organism = "btaurus", correction_method = "fdr", ordered_query = T, domain_scope = "known", custom_bg = background_control, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC"), evcodes = T)
results_control$result$Group = "Control only"
results_control$result


# Now plot
all_results <- list(result = rbind(results_control$result, results_control_all$result, results_infected$result, results_infected_all$result),
                    meta = results_control$meta, results_control_all$meta, results_infected$meta, results_infected_all$meta)


results_file <- apply(all_results$result,2,as.character)
terms <- c("medium-chain fatty-acyl-CoA catabolic process",
        "succinyl-CoA metabolic process",
        "antigen processing and presentation of peptide antigen",
        "ribonucleoside bisphosphate metabolic process",
        "Endosomal/Vacuolar pathway",
        "ER-Phagosome pathway",
        "MHC class II protein complex",
        "Type I diabetes mellitus",
        "Autoimmune thyroid disease",
        "Antigen Presentation: Folding, assembly and peptide loading of class I MHC",
        "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell",
        "Leishmaniasis",
        "regulation of T cell apoptotic process",
        "PD-1 signaling",
        "Phosphorylation of CD3 and TCR zeta chains", 
        "negative regulation of macrophage inflammatory protein 1 alpha production",
        "negative regulation of chemokine (C-C motif) ligand 4 production",
        "negative regulation of activated CD8-positive, alpha-beta T cell apoptotic process",
        "negative regulation of chemokine (C-C motif) ligand 5 production",
        "negative regulation of interleukin-13 production" ,
        "Inflammatory bowel disease",
        "Th17 cell differentiation" ,
        "Th1 and Th2 cell differentiation",
        "Asthma",
        "Costimulation by the CD28 family")

all_results$result <- all_results$result %>%
mutate(label = ifelse(term_name %in% terms, term_name, NA))
colnames(result_lables)[11] <- "term_name"
typeof(results_file)
write.table(results_file, file ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/Gprofiler_cis_comparison_reactor_control_All_results.txt", sep = "\t", quote = F, row.names = F)
all_results$result <- all_results$result %>%
mutate(label = ifelse(term_name %in% terms, term_name, NA))
colnames(result_lables)[11] <- "term_name"
typeof(results_file)

plot <- gostplot(all_results, capped=FALSE, interactive=FALSE, pal=c(`GO:BP`="#ff9900", `GO:CC` ="purple", `KEGG`= "#dd4477",`REAC`="#3366cc")) + labs(y=bquote(~-log[10]~italic(P)[adj]), title = "", x = "Database") +
    theme_bw() +
    ylim(-1,-log10(0.0000001)) +
     geom_text_repel(aes(label = all_results$result$label), colour = "black", max.overlaps = 40,
                  nudge_x = 0, nudge_y = 1.2, size = 4) +
    theme(axis.text.x = element_text(angle = 90, size = 12, vjust = 0.5, hjust =0.2, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 12, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12)) + facet_wrap(~Group) 
plot

ggsave("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/Gprofiler_cis_eGenes_control_reactor.pdf", width = 12, height = 12, dpi = 600)


all_results$result[,c(11,17)]
#################
#### Part 2 ####

# Now read in the DE results



ALL_perm <- read.table(args[1], header = T)
length(ALL_perm$phenotype_id)

DE <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/ALL_results.txt", sep = "\t")
DE$phenotype_id <- rownames(DE)
common <- intersect(DE$phenotype_id, ALL_perm$phenotype_id) %>% as.vector()
length(common)
DE <- DE[common,]
ALL_perm <- ALL_perm %>% filter(ALL_perm$phenotype_id %in% common)

ALL <- left_join(ALL_perm, DE)
head(DE)
ALL$is_de <- if_else(ALL$padj < 0.05, TRUE, FALSE)
table(ALL$is_de)
table(ALL$is_eGene)

table(ALL$is_de)
table(ALL$is_eGene)
dim(DE)
ALL %>% filter(is_eGene == TRUE & is_de == TRUE) %>% dim()
#  data
total_genes <- length(rownames(ALL))
total_cis_eGenes <- 6676
total_DE_genes <- 2388
overlap_genes <- 1059

contingency_table <- matrix(c(overlap_genes, total_DE_genes - overlap_genes,
                               total_cis_eGenes - overlap_genes, total_genes - (total_DE_genes + total_cis_eGenes - overlap_genes)),
                             nrow = 2, byrow = TRUE)
chi_square_result <- chisq.test(contingency_table)
chi_square_result
p_value <- chi_square_result$p.value
X <- chi_square_result$statistic
p_value # 0.0.5852 not more than what would be expected by chance.
X # 0.2979
contingency_df <- as.data.frame(contingency_table)
head(contingency_df)
colnames(contingency_df) <- c("cis-eGene", "Non cis-eGene")
row.names(contingency_df) <- c("DE gene", "Non DE gene")
contingency_df$Gene_Group <- row.names(contingency_df)

# Reshape the data for plotting
library(reshape2)
contingency_melted <- melt(contingency_df, id.vars = "Gene_Group")

# Create the bar plot using ggplot2
ggplot(contingency_melted, aes(x = Gene_Group, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#7570b3", "#1b9e77")) +
  theme_bw() +
  labs(y = "Number of genes",
       x = "Differential expression status",
       fill = "cis-eGene status") +
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25, size = 5) +
  theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, hjust =0.2, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 12, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12))
1059 + 1329 + 5617 + 6607
ggsave("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/Overlap_cis_eGenes_DE_genes.pdf", width = 8, height = 10, dpi = 600)



print("Number of background genes for DE")
print(length(common))

# now have expression values
ALL_target_DE_genes <- ALL %>% filter(is_eGene == TRUE & is_de == TRUE) %>% select(phenotype_id, padj)


ALL_target_DE_genes <- ALL_target_DE_genes[order(ALL_target_DE_genes$padj, decreasing = F),]



print("Number of genes for control only")
print(length(common)) # 869

target_genes_DE <- list(ALL_target_DE_genes$phenotype_id)
results_DE <- gost(query = target_genes_DE,organism = "btaurus", correction_method = "fdr", ordered_query = T, domain_scope = "known", custom_bg = common, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC"), evcodes = T)

results_DE$result
lollipop_plot <- results_DE$result
results_DE$result
head(lollipop_plot)
lollipop_plot <- lollipop_plot[order(lollipop_plot$p_value, decreasing = F),]
head(lollipop_plot)
lollipop_plot %>% 
  mutate(term_name = fct_reorder(term_name, desc(p_value))) %>% ggplot(., aes(-log10(p_value), term_name, colour = source)) +
        geom_segment(aes(x = 0, y = term_name, xend = -log10(p_value), yend = term_name), color = "grey50") +
        geom_point(size = 4) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, hjust =0.2, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 12, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12)) +
        scale_colour_manual(breaks = c("GO:BP", "GO:CC", "KEGG"), values = c("#ff9900", "purple","#dd4477")) +
        xlim(0,-log10(0.0001))

ggsave("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/Enrichement_cis_eGenes_DE_genes_lollipop.pdf", width = 15, height = 10, dpi = 600)

results_file <- apply(results_DE$result,2,as.character)
typeof(results_file)

write.table(results_file, file ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/Gprofiler_cis_comparison_DE_genes.txt", sep = "\t", quote = F, row.names = F)
