# Script for getting source data

# Fig. 3A

# Script to get eQTL numbers which are significant
# Also to get a violin plot of gene-wise FDR cutoffs


library(data.table)
library(tidyverse)
library(ggplot2)
library(UpSetR)


ALL_perm <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt", header = T) %>% 
  filter(is_eGene == T) %>% select(7, 1,14) %>% mutate(Group = "AAG")
CONTROL_perm <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_fdr0.05.txt", header = T) %>% 
  filter(is_eGene == T) %>% select(7, 1,14) %>% mutate(Group = "bTB−")
INFECTED_perm <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_fdr0.05.txt", header = T) %>% 
  filter(is_eGene == T) %>% select(7, 1,14) %>% mutate(Group = "bTB+")


Fig_3A <- rbind(ALL_perm, CONTROL_perm, INFECTED_perm)

write.table(Fig_3A, "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_3A.txt", sep = "\t", quote = F, row.names = F)





# Fig 3B

ALL_eGene <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt", header = T) %>% 
  filter(is_eGene == T) %>% select(1) 
CONTROL_eGene <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_fdr0.05.txt", header = T) %>% 
  filter(is_eGene == T) %>% select(1)
INFECTED_eGene <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_fdr0.05.txt", header = T) %>% 
  filter(is_eGene == T) %>% select(1)

ALL <- read.table(gzfile("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_independent_qtl.txt.gz"), sep = "\t", header = T)
CONTROL <- read.table(gzfile("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_independent_qtl.txt.gz"), sep = "\t", header = T)
INFECTED <- read.table(gzfile("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_independent_qtl.txt.gz"), sep = "\t", header = T)

ALL <- inner_join(ALL, ALL_eGene)
CONTROL <- inner_join(CONTROL, CONTROL_eGene)
INFECTED <- inner_join(INFECTED, INFECTED_eGene)

dim(ALL)
dim(CONTROL)
dim(INFECTED)
# Percentage increase in eQTL numbers
all_increase <- ((nrow(ALL) - nrow(ALL_perm)) / nrow(ALL_perm)) * 100 
control_increase <- ((nrow(CONTROL) - nrow(CONTROL_perm)) / nrow(CONTROL_perm)) * 100 
infec_increase <- ((nrow(INFECTED) - nrow(INFECTED_perm)) / nrow(INFECTED_perm)) * 100 

ALL %>% filter(rank > 1) %>% dim()
CONTROL %>% filter(rank > 1) %>% dim()
INFECTED %>% filter(rank > 1) %>% dim()



all_increase
control_increase
infec_increase
mean(ALL$rank)
mean(CONTROL$rank)
mean(INFECTED$rank)

# Now plot number per group
ALL <- ALL[order(ALL$rank, decreasing = T),]
ALL <- ALL %>% group_by(phenotype_id) %>% summarize(conditional = max(rank),
                                                    Group = "ALL")
ALL_plot <- ALL %>% group_by(conditional) %>% summarize(Group = "ALL",
                                                        genes = n())

CONTROL <- CONTROL[order(CONTROL$rank, decreasing = T),]
CONTROL <- CONTROL %>% group_by(phenotype_id) %>% summarize(conditional = max(rank),
                                                            Group = "CONTROL")
CONTROL_plot <- CONTROL %>% group_by(conditional) %>% summarize(Group = "CONTROL",
                                                                genes = n())


INFECTED <- INFECTED[order(INFECTED$rank, decreasing = T),]
INFECTED <- INFECTED %>% group_by(phenotype_id) %>% summarize(conditional = max(rank),
                                                              Group = "INFECTED")
INFECTED_plot <- INFECTED %>% group_by(conditional) %>% summarize(Group = "INFECTED",
                                                                  genes = n())

conditional_plot <- rbind(ALL_plot, CONTROL_plot, INFECTED_plot)
conditional_plot$conditional <- as.character(conditional_plot$conditional)
conditional_plot$conditional <- factor(conditional_plot$conditional, levels = c("1","2","3","4","5","6","7","8","9","10","11"))

conditional_plot %>% group_by(Group) %>% summarize(sum(genes))

write.table(conditional_plot, "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_3B.txt", sep = "\t",  quote = F, row.names = F)











# Fig 3C



ALL <- read.table(gzfile("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_independent_qtl.txt.gz"), sep = "\t", header = T)
CONTROL <- read.table(gzfile("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_independent_qtl.txt.gz"), sep = "\t", header = T)
INFECTED <- read.table(gzfile("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_independent_qtl.txt.gz"), sep = "\t", header = T)
ALL <- inner_join(ALL, ALL_eGene)
CONTROL <- inner_join(CONTROL, CONTROL_eGene)
INFECTED <- inner_join(INFECTED, INFECTED_eGene)
# ALL
ALL_distance_perm <- ALL %>% filter(rank == 1) %>% select(start_distance) %>% abs()
ALL_distance_conditional <- ALL %>% filter(rank > 1) %>% select(start_distance) %>% abs()
ALL_distance_perm$Category <- "Top"
ALL_distance_conditional$Category <- "Conditional"
ALL_distance_perm$Group <- "ALL"
ALL_distance_conditional$Group <- "ALL"
ALL_distance_plot <- rbind(ALL_distance_perm, ALL_distance_conditional)
ALL_distance_plot$start_distance <- ALL_distance_plot$start_distance / 1000

# CONTROL
CONTROL_distance_perm <- CONTROL %>% filter(rank == 1) %>% select(start_distance) %>% abs()
CONTROL_distance_conditional <- CONTROL %>% filter(rank > 1) %>% select(start_distance) %>% abs()
CONTROL_distance_perm$Category <- "Top"
CONTROL_distance_conditional$Category <- "Conditional"
CONTROL_distance_perm$Group <- "CONTROL"
CONTROL_distance_conditional$Group <- "CONTROL"
CONTROL_distance_plot <- rbind(CONTROL_distance_perm, CONTROL_distance_conditional)
CONTROL_distance_plot$start_distance <- CONTROL_distance_plot$start_distance / 1000



# INFECTED
INFECTED_distance_perm <- INFECTED %>% filter(rank == 1) %>% select(start_distance) %>% abs()
INFECTED_distance_conditional <- INFECTED %>% filter(rank > 1) %>% select(start_distance) %>% abs() # filter for ranks >1 which is top cis-eQTL
INFECTED_distance_perm$Category <- "Top"
INFECTED_distance_conditional$Category <- "Conditional"
INFECTED_distance_perm$Group <- "INFECTED"
INFECTED_distance_conditional$Group <- "INFECTED"
INFECTED_distance_plot <- rbind(INFECTED_distance_perm, INFECTED_distance_conditional)
INFECTED_distance_plot$start_distance <- INFECTED_distance_plot$start_distance / 1000



merged_distance_plot <- rbind(ALL_distance_plot, CONTROL_distance_plot, INFECTED_distance_plot)
merged_distance_plot


wilcox.test(as.numeric(INFECTED_distance_plot$start_distance) ~ INFECTED_distance_plot$Category)
wilcox.test(as.numeric(CONTROL_distance_plot$start_distance) ~ CONTROL_distance_plot$Category)
wilcox.test(as.numeric(ALL_distance_plot$start_distance) ~ ALL_distance_plot$Category)


write.table(merged_distance_plot, "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_3c.txt", sep = "\t", quote = F, row.names = F)




# Fig. 3D

ALL <- read.table(gzfile("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_independent_qtl.txt.gz"), sep = "\t", header = T)
CONTROL <- read.table(gzfile("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_independent_qtl.txt.gz"), sep = "\t", header = T)
INFECTED <- read.table(gzfile("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_independent_qtl.txt.gz"), sep = "\t", header = T)
ALL <- inner_join(ALL, ALL_eGene)
CONTROL <- inner_join(CONTROL, CONTROL_eGene)
INFECTED <- inner_join(INFECTED, INFECTED_eGene)


# ALL
ALL_slope_distance <- ALL %>% select(start_distance, slope) 
ALL_slope_distance$start_distance <- abs(ALL_slope_distance$start_distance)
# CONTROL
CONTROL_slope_distance <- CONTROL %>% select(start_distance, slope)
CONTROL_slope_distance$start_distance <- abs(CONTROL_slope_distance$start_distance)

# INFECTED
INFECTED_slope_distance <- INFECTED %>% select(start_distance, slope) 
INFECTED_slope_distance$start_distance <- abs(INFECTED_slope_distance$start_distance)

ALL_slope_distance$Group <- "ALL"
CONTROL_slope_distance$Group <- "CONTROL"
INFECTED_slope_distance$Group <- "INFECTED"

final_point <- rbind(ALL_slope_distance, CONTROL_slope_distance, INFECTED_slope_distance)
final_point$slope <- abs(final_point$slope)
final_point$start_distance <- final_point$start_distance / 1000000

write.table(final_point, "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_3d.txt", sep = "\t", quote = F, row.names = F)

result_ALL <- cor.test(ALL_slope_distance[,1], abs(ALL_slope_distance[,2]), method = "spearman", exact = FALSE)
result_CONTROL <- cor.test(CONTROL_slope_distance[,1], abs(CONTROL_slope_distance[,2]), method = "spearman", exact = FALSE)
result_INFECTED <- cor.test(INFECTED_slope_distance[,1], abs(INFECTED_slope_distance[,2]), method = "spearman", exact = FALSE)

result_ALL$estimate
result_CONTROL$estimate
result_INFECTED$estimate


result_ALL$p.value
result_CONTROL$p.value
result_INFECTED$p.value
# Fig. 3E

########################################################################################
########################################################################################
## 1. AC methodology
########################################################################################
########################################################################################

#To calculate AC, we took the SNP-gene combination with the lowest p-value in the discovery cohort 
#for each gene if it was significant and matched these with the same SNP-gene pairs in the replication cohort if it was also significant. 
#We then determined the percentage of those significant SNP-gene pairs had the same allelic direction of effect compared to the discovery cohort. 
args = commandArgs(trailingOnly = TRUE)
# Note, significant in both so will need the two files, permutation (gene level significance and nominal p value)
args[1] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt"
args[2] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_fdr0.05.txt"
args[3] <-  "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_fdr0.05.txt"
args[4] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Blood.nominals.2rd.txt"
args[5] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Gtex_thresholds_determined.txt"
# Read in the data
# select top snp which were significant genome wide - note we are comparing the top eQTL in our dataset 
ALL <- read.table(args[1], header = T) %>% 
  filter(is_eGene == T) #%>% select(7, 1,13)
CONTROL <- read.table(args[2], header = T) %>% 
  filter(is_eGene == T) #%>% select(7, 1,13)
INFECTED <- read.table(args[3], header = T) %>% 
  filter(is_eGene == T) #%>% select(7, 1,13)
head(ALL)
ALL$association <- paste0(ALL$variant_id, "-", ALL$phenotype_id)
CONTROL$association <- paste0(CONTROL$variant_id, "-", CONTROL$phenotype_id)
INFECTED$association <- paste0(INFECTED$variant_id, "-", INFECTED$phenotype_id)

Cattle_gtex <- read.table(args[4], sep = " ", header =F)
Cattle_gtex_permute <- read.table(args[5], sep = "\t", header =T) # Note this file, we can determine signifant eQTLS which are < p-value nominal threshold (see methods)
Cattle_gtex_permute <- Cattle_gtex_permute %>% select(phenotype_id, pval_nominal_threshold)

# Modify the dataframe
colnames(Cattle_gtex) <- c("phenotype_id", "variant_id", "distance", "pval_nominal", "slope")
Cattle_gtex <- Cattle_gtex %>% select(2,1,4,5,3)
Cattle_gtex$variant_id <- gsub("_", ":", Cattle_gtex$variant_id)
Cattle_gtex$association <- paste0(Cattle_gtex$variant_id, "-", Cattle_gtex$phenotype_id)

# Join to get p_nominal_threshold for determining if an evariant is significant or not
Cattle_gtex <- left_join(Cattle_gtex, Cattle_gtex_permute)

# Now determine if an eVariant is significant or not
# Can then extract these for the AC comparison
Cattle_gtex$is_eVariant <- if_else(Cattle_gtex$pval_nominal < Cattle_gtex$pval_nominal_threshold, TRUE, FALSE)

Cattle_gtex_signif <- Cattle_gtex %>% filter(is_eVariant == TRUE)

head(Cattle_gtex_signif)
# Get common SNPs
gtex_ALL_signif <- Cattle_gtex_signif %>% filter(association %in% ALL$association)
gtex_CONTROL_signif <- Cattle_gtex_signif %>% filter(association %in% CONTROL$association)
gtex_INFECTED_signif <- Cattle_gtex_signif %>% filter(association %in% INFECTED$association)


ALL <- ALL %>% filter(association %in% gtex_ALL_signif$association)
CONTROL <- CONTROL %>% filter(association %in% gtex_CONTROL_signif$association)
INFECTED <- INFECTED %>% filter(association %in% gtex_INFECTED_signif$association)


# Now order based on association so they are in the same order
gtex_ALL_signif <- gtex_ALL_signif[order(gtex_ALL_signif$association),]
gtex_CONTROL_signif <- gtex_CONTROL_signif[order(gtex_CONTROL_signif$association),]
gtex_INFECTED_signif <- gtex_INFECTED_signif[order(gtex_INFECTED_signif$association),]
ALL <- ALL[order(ALL$association),]
CONTROL <- CONTROL[order(CONTROL$association),]
INFECTED <- INFECTED[order(INFECTED$association),]

head(ALL)
ALL_ac <- ALL %>% select(22,14)
colnames(ALL_ac)[2] <- "slope_discovery"
CONTROL_ac <- CONTROL %>% select(22,14)
colnames(CONTROL_ac)[2] <- "slope_discovery"
INFECTED_ac <- INFECTED %>% select(22,14)
colnames(INFECTED_ac)[2] <- "slope_discovery"
colnames(gtex_ALL_signif)
gtex_ALL_ac <- gtex_ALL_signif %>% select(association, slope)
gtex_CONTROL_ac <- gtex_CONTROL_signif %>% select(association, slope)
gtex_INFECTED_ac <- gtex_INFECTED_signif %>% select(association, slope)
head(CONTROL_ac)
ALL_ac <- left_join(ALL_ac, gtex_ALL_ac)
CONTROL_ac <- left_join(CONTROL_ac, gtex_CONTROL_ac)
INFECTED_ac <- left_join(INFECTED_ac, gtex_INFECTED_ac)

ALL_ac <- ALL_ac %>% mutate(Concordance = case_when(
  slope_discovery < 0 & slope < 0 ~ TRUE,
  slope_discovery > 0 & slope < 0 ~ FALSE,
  slope_discovery < 0 & slope > 0 ~ FALSE,
  slope_discovery > 0 & slope > 0 ~ TRUE))
ALL_ac_percentage <- table(ALL_ac$Concordance)[2] / nrow(ALL_ac) * 100

CONTROL_ac <- CONTROL_ac %>% mutate(Concordance = case_when(
  slope_discovery < 0 & slope < 0 ~ TRUE,
  slope_discovery > 0 & slope < 0 ~ FALSE,
  slope_discovery < 0 & slope > 0 ~ FALSE,
  slope_discovery > 0 & slope > 0 ~ TRUE))
CONTROL_ac_percentage <- table(CONTROL_ac$Concordance)[2] / nrow(CONTROL_ac) * 100

INFECTED_ac <- INFECTED_ac %>% mutate(Concordance = case_when(
  slope_discovery < 0 & slope < 0 ~ TRUE,
  slope_discovery > 0 & slope < 0 ~ FALSE,
  slope_discovery < 0 & slope > 0 ~ FALSE,
  slope_discovery > 0 & slope > 0 ~ TRUE))
INFECTED_ac_percentage <- table(INFECTED_ac$Concordance)[2] / nrow(INFECTED_ac) * 100

print("AC summary statistics")
ALL_ac_percentage # 98.79397
CONTROL_ac_percentage # 98.91304
INFECTED_ac_percentage # 99.06542

library(qvalue)


ALL <- read.table(args[1], header = TRUE) %>% filter(is_eGene == T) #%>% select(7, 1,13)
CONTROL <- read.table(args[2], header = T) %>% 
  filter(is_eGene == T) #%>% select(7, 1,13)
INFECTED <- read.table(args[3], header = T) %>% 
  filter(is_eGene == T) #%>% select(7, 1,13)
head(ALL)
ALL$association <- paste0(ALL$variant_id, "-", ALL$phenotype_id)
CONTROL$association <- paste0(CONTROL$variant_id, "-", CONTROL$phenotype_id)
INFECTED$association <- paste0(INFECTED$variant_id, "-", INFECTED$phenotype_id)



gtex_ALL <- Cattle_gtex %>% filter(association %in% ALL$association)
gtex_CONTROL <- Cattle_gtex %>% filter(association %in% CONTROL$association)
gtex_INFECTED <- Cattle_gtex %>% filter(association %in% INFECTED$association)


# Cross reference
gtex_ALL <- gtex_ALL[order(gtex_ALL$pval_nominal, decreasing = F),]
gtex_CONTROL <- gtex_CONTROL[order(gtex_CONTROL$pval_nominal, decreasing = F),]
gtex_INFECTED <- gtex_INFECTED[order(gtex_INFECTED$pval_nominal, decreasing = F),]


pi1_ALL <- 1 - pi0est(as.numeric(gtex_ALL$pval_nominal))$pi0
pi1_CONTROL <- 1- pi0est(as.numeric(gtex_CONTROL$pval_nominal))$pi0
pi1_INFECTED <- 1- pi0est(as.numeric(gtex_INFECTED$pval_nominal))$pi0

pi1_ALL #0.6138685
pi1_CONTROL #  0.6688508
pi1_INFECTED # 0.7981247

dim(gtex_ALL) #  4362 
dim(gtex_CONTROL) # 2160
dim(gtex_INFECTED) #  1497
# Number of bootstrap iterations
num_bootstraps <- 100

# Create an empty matrix to store the bootstrap π1 estimates for each group
bootstrap_ALL_estimates <- matrix(0, nrow = num_bootstraps, ncol = 1)  # Three groups: ALL, CONTROL, INFECTED
lambda_values <- seq(0.05,0.95,0.05)
num_bootstraps
# Bootstrap loop
set.seed(4587)
for (i in 1:num_bootstraps) {
  # Sample with replacement from the p-values
  sample_ps <- gtex_ALL$pval_nominal
  bootstrap_sample <- sample(sample_ps, replace = TRUE)
  bootstrap_sample <- sort(bootstrap_sample, decreasing = F)
  # Estimate π0 for the bootstrap sample
  bootstrap_pi0 <- pi0est(bootstrap_sample, lambda = lambda_values)
  
  # Calculate π1 for the bootstrap sample (1 - π0)
  bootstrap_pi1 <- 1 - bootstrap_pi0$pi0
  
  # Store the π1 estimate
  bootstrap_ALL_estimates[i,] <- bootstrap_pi1
}

set.seed(458798) # Only have to keep one seed
# Bootstrap loop
bootstrap_CONTROL_estimates <- matrix(0, nrow = num_bootstraps, ncol = 1)
for (i in 1:num_bootstraps) {
  # Sample with replacement from the p-values
  sample_ps <- gtex_CONTROL$pval_nominal
  bootstrap_sample <- sample(sample_ps, replace = TRUE)
  bootstrap_sample <- sort(bootstrap_sample, decreasing = F)
  
  # Estimate π0 for the bootstrap sample
  bootstrap_pi0 <- pi0est(bootstrap_sample, lambda = lambda_values)
  
  # Calculate π1 for the bootstrap sample (1 - π0)
  bootstrap_pi1 <- 1 - bootstrap_pi0$pi0
  
  # Store the π1 estimate
  bootstrap_CONTROL_estimates[i,] <- c(bootstrap_pi1)
}


# Bootstrap loop
bootstrap_INFECTED_estimates <- matrix(0, nrow = num_bootstraps, ncol = 1)
for (i in 1:num_bootstraps) {
  # Sample with replacement from the p-values
  sample_ps <- gtex_INFECTED$pval_nominal
  bootstrap_sample <- sample(sample_ps, replace = TRUE)
  bootstrap_sample <- sort(bootstrap_sample, decreasing = F)
  
  # Estimate π0 for the bootstrap sample
  bootstrap_pi0 <- pi0est(bootstrap_sample, lambda = lambda_values)
  
  # Calculate π1 for the bootstrap sample (1 - π0)
  bootstrap_pi1 <- 1 - bootstrap_pi0$pi0
  
  # Store the π1 estimate
  bootstrap_INFECTED_estimates[i,] <- c(bootstrap_pi1)
}
bootstrap_INFECTED_estimates
bootstrap_ALL_estimates
# Print the summary statistics for each group

colnames(bootstrap_ALL_estimates) <- c("Gtex_ALL")
colnames(bootstrap_CONTROL_estimates) <- c("Gtex_CONTROL")
colnames(bootstrap_INFECTED_estimates) <- c("Gtex_INFECTED")
bootstrap_ALL_estimates <- bootstrap_ALL_estimates %>% as.data.frame()
bootstrap_CONTROL_estimates <- bootstrap_CONTROL_estimates %>% as.data.frame()
bootstrap_INFECTED_estimates <- bootstrap_INFECTED_estimates %>% as.data.frame()


# ALL
mean_values_ALL <- colMeans(bootstrap_ALL_estimates)
sd_values_ALL <- apply(bootstrap_ALL_estimates, 2, sd)
sem_ALL <- sd_values_ALL / sqrt(100)
mean_values_ALL # 0.6247269
sem_ALL
sd_values_ALL #0.02478944

summary_data_ALL <- data.frame(
  Group = c('Gtex'),
  Category = c("ALL"),
  Mean = mean_values_ALL,
  SD = sd_values_ALL,
  SEM = sem_ALL
)
head(summary_data_ALL)

print("Summary Statistics Pioest")
# CONTROL
mean_values_CONTROL <- colMeans(bootstrap_CONTROL_estimates)
sd_values_CONTROL <- apply(bootstrap_CONTROL_estimates, 2, sd)
mean_values_CONTROL #0.6355099 
sem_CONTROL <- sd_values_CONTROL / sqrt(100)
sd_values_CONTROL # Gtex_CONTROL - 0.03327106 

summary_data_CONTROL <- data.frame(
  Group = c('Gtex'),
  Category = c("CONTROL"),
  Mean = mean_values_CONTROL,
  SD = sd_values_CONTROL,
  SEM = sem_CONTROL
)

# INFECTED
mean_values_INFECTED <- colMeans(bootstrap_INFECTED_estimates)
sd_values_INFECTED <- apply(bootstrap_INFECTED_estimates, 2, sd)
sem_INFECTED <- sd_values_INFECTED / sqrt(100)
mean_values_INFECTED # 0.7858296
sd_values_INFECTED  # 0.02731303

summary_data_INFECTED <- data.frame(
  Group = c('Gtex'),
  Category = c("INFECTED"),
  Mean = mean_values_INFECTED,
  SD = sd_values_INFECTED,
  SEM = sem_INFECTED
)



head(summary_data_ALL)
head(summary_data_CONTROL)
head(summary_data_INFECTED)

print("SUMMARY STATISTICS FOR TEXT")
mean_values_ALL
mean_values_CONTROL
mean_values_INFECTED
sem_ALL
sem_CONTROL
sem_INFECTED
summary_data <- rbind(summary_data_ALL, summary_data_CONTROL, summary_data_INFECTED)
summary_data

bootstrap_ALL_estimates$Group = "AAG"
bootstrap_CONTROL_estimates$Group = "bTB-"
bootstrap_INFECTED_estimates$Group = "bTB+"

colnames(bootstrap_ALL_estimates)[1] <- "Value"
colnames(bootstrap_CONTROL_estimates)[1] <- "Value"
colnames(bootstrap_INFECTED_estimates)[1] <- "Value"
pi1_df = rbind(bootstrap_ALL_estimates, bootstrap_CONTROL_estimates, bootstrap_INFECTED_estimates)



pi1_df
ggplot(pi1_df, aes(x = Group, y=Value, colour = Group)) + 
  geom_point(alpha = 0.1) +
  scale_colour_manual(breaks = c("AAG", "bTB-", "bTB+"), values = c("#542788","#2166ac", "#b2182b")) +
  stat_summary(fun=mean, geom="point", aes(color = Group), size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25, aes(color = Group)) + 
  theme_bw() +
  labs(x = "cis-eQTL cohort", y = "Pi1 statistic") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.position = "none")


ggsave("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/plotting/cis_replication_pi1_bootstrap_dotplot.pdf", width = 12, height = 12, dpi = 600)
#
write.table(pi1_df, file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_3e.txt", sep = "\t", quote = F, row.names = F)

ALL <- read.table(args[1], header = T) %>% 
  filter(is_eGene == T) #%>% select(7, 1,14)
CONTROL <- read.table(args[2], header = T) %>% 
  filter(is_eGene == T) #%>% select(7, 1,14)
INFECTED <- read.table(args[3], header = T) %>% 
  filter(is_eGene == T) #%>% select(7, 1,14)


# get association
ALL$association <- paste0(ALL$variant_id, "-", ALL$phenotype_id)
CONTROL$association <- paste0(CONTROL$variant_id, "-", CONTROL$phenotype_id)
INFECTED$association <- paste0(INFECTED$variant_id, "-", INFECTED$phenotype_id)

# filter for those which are common
gtex_ALL <- Cattle_gtex %>% filter(association %in% ALL$association)
gtex_CONTROL <- Cattle_gtex %>% filter(association %in% CONTROL$association)
gtex_INFECTED <- Cattle_gtex %>% filter(association %in% INFECTED$association)


# Croos reference
ALL <- ALL %>% filter(association %in% gtex_ALL$association)
CONTROL <- CONTROL %>% filter(association %in% gtex_CONTROL$association)
INFECTED <- INFECTED %>% filter(association %in% gtex_INFECTED$association)

# ------------------------------------------
# Get the error of each slope
# -----------------------------------------
sample_size <- 698 # Number of blood samples in supplementary information of cattle gtex
num_independent_vars <- 3 # 3 genotypes tested, check this
# 10 peer, 10 genotype PCs
covariates <- 20
# Calculate degrees of freedom
df <- sample_size - (num_independent_vars + covariates + 1)

alpha <- 0.05
t_critical <- qt(1 - alpha/2, df)
gtex_ALL$se_slope <- gtex_ALL$slope / t_critical # t = b - 0 / se(B) when you work back
gtex_CONTROL$se_slope <- gtex_CONTROL$slope / t_critical
gtex_INFECTED$se_slope <- gtex_INFECTED$slope / t_critical

# Setting up equation
ALL$Group <- "ALL"
CONTROL$Group <- "CONTROL"
INFECTED$Group <- "INFECTED"
gtex_ALL$Group <- "ALL"
gtex_CONTROL$Group <- "CONTROL"
gtex_INFECTED$Group <- "INFECTED"

colnames(ALL)[14] <- "slope_ref"
colnames(INFECTED)[14] <- "slope_ref"
colnames(CONTROL)[14] <- "slope_ref"

ALL <- left_join(ALL, gtex_ALL, by = c("association" = "association"))
CONTROL <- left_join(CONTROL, gtex_CONTROL, by = c("association" = "association"))
INFECTED <- left_join(INFECTED, gtex_INFECTED, by = c("association" = "association"))

final_point <- rbind(ALL, CONTROL, INFECTED)
final_point <- final_point %>% select(association, slope_ref, slope, Group.x)
final_point

final_point %>% group_by(Group.x) %>% summarize(Speaman = cor(slope_ref, slope, method = "spearman"),
                                                P = cor.test(slope_ref, slope, method = "spearman", exact = FALSE)[[3]])


cor.test(ALL$slope_ref, ALL$slope, method = "spearman", exact = F)$p.value


write.table(final_point, file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_3f.txt", sep = "\t", quote = F, row.names = F)







# Fig 4A
ALL_afc <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/ALL_AFC.txt")
CON_afc <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/CONTROL_AFC.txt")
INFEC_afc <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/afc/INFECTED_AFC.txt")

ALL_eQTL <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt", header = T) %>% filter(is_eGene == TRUE) 
CON_eQTL <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_fdr0.05.txt", header = T) %>% filter(is_eGene == TRUE) 
INFEC_eQTL <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_fdr0.05.txt", header = T) %>% filter(is_eGene == TRUE) 

ALL_afc <- ALL_afc %>% filter(ALL_afc$phenotype_id %in% ALL_eQTL$phenotype_id)
CON_afc <- CON_afc %>% filter(CON_afc$phenotype_id %in% CON_eQTL$phenotype_id)
INFEC_afc <- INFEC_afc %>% filter(INFEC_afc$phenotype_id %in% INFEC_eQTL$phenotype_id)


ALL_afc$group <- "ALL"
CON_afc$group <- "CONTROL"
INFEC_afc$group <- "INFECTED"




afc_plot <- rbind(ALL_afc, CON_afc, INFEC_afc)

afc_plot_source = afc_plot %>% select(phenotype_id, variant_id, group, log2_aFC)
afc_plot_source$log2_aFC <- abs(afc_plot_source$log2_aFC)


write.table(afc_plot_source, file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_4a.txt", sep = "\t", quote = F, row.names = F)





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
median(abs(ALL_afc$log2_aFC)) # 0.2998068
median(abs(CON_afc$log2_aFC)) # 0.4294106
median(abs(INFEC_afc$log2_aFC)) # 0.5384912


afc_plot_source <- afc_plot %>% select(phenotype_id, variant_id, slope, log2_aFC, group)

write.table(afc_plot_source, "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_4b.txt", sep = "\t", quote = F, row.names = F)

afc_plot_source %>% group_by(group) %>% summarize(cor.test(slope, log2_aFC, method = "spearman"))

afc_plot_source %>% filter(abs(log2_aFC) >=1) %>% group_by(group) %>% summarize(Group = n())

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





######################################################
# Plot ieQTL results for O'Grady et al., 2024
#####################################################

library(data.table)
library(tidyverse)
library(vcfR)

# First thing is to read in the data for top eQTLs
data <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/ieQTL/results/raw/ALL_TB_interaction.cis_qtl_top_assoc.txt.gz") %>% filter(pval_adj_bh < 0.25)
head(data)


vcf_all <- vcfR::read.vcfR("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz")



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
expression = fread('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_residualised_expression.txt') %>% select(-c(1:4,6))
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
  
  return(list(plot_dom, counts_gene))
  
  
  
}
PCBP2 = eQTL_plot(gene_id = "ENSBTAG00000020757", SNP_id = "5:25861761:G:A",  gene_name = "PCBP2", vcf = vcf_all, counts = expression,  HOM = "GG", HET = "GA", HOM_ALT = "AA")

write.table(PCBP2[[2]], file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_4c.txt", sep = "\t", quote = F, row.names = F)



# Fig 4 e - h
library(vcfR)
vcf_control <- vcfR::read.vcfR("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/CONTROL_IMPUTED_UPDATED.vcf.gz")
vcf_infected <- vcfR::read.vcfR("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/INFECTED_IMPUTED_UPDATED.vcf.gz")

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
counts_control <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/CONTROL_residualised_expression.txt") %>% select(-c(1:4,6)) #%>% as.matrix()
rownames(counts_control) <- counts_control$gid

counts_infected <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/INFECTED_residualised_expression.txt") %>% select(-c(1:4,6)) #%>% as.matrix()

rownames(counts_infected) <- counts_infected$gid

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
  
  
  return(list(plot_grid(p_control, p_infected), counts_control_gene, counts_infected_gene))
}
IFITM3 <- eQTL_plot(gene_id = "ENSBTAG00000019015", SNP_id = "29:50294904:G:A", SNP_id2 = "29:50294904:G:A", gene_name = "IFITM3", counts_control = counts_control, counts_infected = counts_infected, vcf_control = vcf_control, vcf_infected = vcf_infected, HOM = "GG", HET = "GA", HOM_ALT = "AA", HOM2 = "GG", HET2 = "GA", HOM_ALT2 = "AA")
write.table(IFITM3[[2]], file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_4f_control.txt", sep = "\t", quote = F, row.names = F)
write.table(IFITM3[[3]], file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_4f_infected.txt", sep = "\t", quote = F, row.names = F)

IFI16 <- eQTL_plot(gene_id = "ENSBTAG00000011511", SNP_id = "3:10984726:G:A", SNP_id2 = "3:10984726:G:A", gene_name = "IFI16", counts_control = counts_control, counts_infected = counts_infected, vcf_control = vcf_control, vcf_infected = vcf_infected, HOM = "GG", HET = "GA", HOM_ALT = "AA", HOM2 = "GG", HET2 = "GA", HOM_ALT2 = "AA")
write.table(IFI16[[2]], file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_4e_control.txt", sep = "\t", quote = F, row.names = F)
write.table(IFI16[[3]], file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_4e_infected.txt", sep = "\t", quote = F, row.names = F)


IL1R1 <- eQTL_plot(gene_id = "ENSBTAG00000005273", SNP_id = "11:6851248:C:G", SNP_id2 = "11:6851248:C:G", gene_name = "IL1R1", counts_control = counts_control, counts_infected = counts_infected, vcf_control = vcf_control, vcf_infected = vcf_infected, HOM = "CC", HET = "CG", HOM_ALT = "GG", HOM2 = "CC", HET2 = "CG", HOM_ALT2 = "GG")
write.table(IL1R1[[2]], file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_4g_control.txt", sep = "\t", quote = F, row.names = F)
write.table(IL1R1[[3]], file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_4g_infected.txt", sep = "\t", quote = F, row.names = F)


RGS10 <- eQTL_plot(gene_id = "ENSBTAG00000002647", SNP_id = "26:40304855:C:A", SNP_id2 = "26:40304855:C:A", gene_name = "RGS10", counts_control = counts_control, counts_infected = counts_infected, vcf_control = vcf_control, vcf_infected = vcf_infected, HOM = "CC", HET = "CG", HOM_ALT = "GG", HOM2 = "CC", HET2 = "CG", HOM_ALT2 = "GG")
RGS10
write.table(RGS10[[2]], file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_4h_control.txt", sep = "\t", quote = F, row.names = F)
write.table(RGS10[[3]], file = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/Source_data/Fig_4h_infected.txt", sep = "\t", quote = F, row.names = F)



