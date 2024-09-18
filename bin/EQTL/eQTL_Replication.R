# script to assess the replicaiton of EQTLs
library(data.table)
library(qvalue)
library(tidyverse)
library(ggsignif)
library(UpSetR)



args = commandArgs(trailingOnly = TRUE)
########################################################################################
########################################################################################
## 1. AC methodology
########################################################################################
########################################################################################

#To calculate AC, we took the SNP-gene combination with the lowest p-value in the discovery cohort 
#for each gene if it was significant and matched these with the same SNP-gene pairs in the replication cohort if it was also significant. 
#We then determined the percentage of those significant SNP-gene pairs had the same allelic direction of effect compared to the discovery cohort. 

# Note, significant in both so will need the two files, permutation (gene level significance and nominal p value)
#args[1] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt"
#args[2] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/CONTROL.cis_qtl_fdr0.05.txt"
#args[3] <-  "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/INFECTED.cis_qtl_fdr0.05.txt"
#args[4] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Blood.nominals.2rd.txt"
#args[5] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Gtex_thresholds_determined.txt"
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

########################################################################################
########################################################################################
## 2. Pioest methodology
########################################################################################
########################################################################################


#To calculate π1 we took the SNP-gene combination with the lowest p-value in the 
#discovery cohort for each gene if it was significant, 
#ordered the SNP-gene combinations from lowest to highest p-value, and matched these with the same SNP-gene combination in the replication cohort. 
#SNP-gene combinations not tested in the replication cohort were removed. We then calculated the proportion of true null p-values (π0) with the pi0est function of R, with
# parameter p = the p-values of the match SNP-gene combinations in the replication cohort, and calculated π1 with 1- π0. 

# Read in the data
# select top snp which were significant genome wide
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

# Bootstrap loop
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
set.seed(458798)
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
sem_ALL <- sd_values_ALL / sqrt(length(gtex_ALL$pval_nominal))
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
sem_CONTROL <- sd_values_CONTROL / sqrt(length(gtex_CONTROL$pval_nominal))
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
sem_INFECTED <- sd_values_INFECTED / sqrt(length(gtex_INFECTED$pval_nominal))
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


ggplot(summary_data, aes(x = Category, y = Mean, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM),
                position = position_dodge(width = 0.8), width = 0.4) +
  labs(
       x = "eQTL cohort", y = "Replication rate (Storey’s π1)", fill = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.position = "none") +
        scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.4), limits = c(0,1)) +
        scale_x_discrete(labels = c("All", "Control", "Infected")) +
        scale_fill_manual(breaks = c("ALL", "CONTROL", "INFECTED"), values = c("#542788","#2166ac", "#b2182b"))
ggsave(args[6], height = 12, width = 8, dpi = 600)
#ggsave("cis_replication_pi1_bootstrap.pdf", height = 12, width = 8, dpi = 600)
########################################################################################
########################################################################################
## 3. correlation methodology
########################################################################################
########################################################################################

# Here we take all matching in Gtex

# we calculated Rb by taking the SNP-gene combination with the lowest p-value in the discovery cohort 
# for each gene if it was significant, matching those to the replication dataset, 
# regardless of significance in the replication dataset.

# Now onto eQTL correlation effect size
# Read in the data
# select top snp which were significant genome wide
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
head(final_point)
ggplot(data = final_point, aes(x = slope_ref, y = slope, colour = Group.x)) + 
geom_point(alpha = 0.2) + scale_colour_discrete(labels = c("All", "Control", "Infected")) +
scale_colour_manual(labels = c("All", "Control", "Infected"), values = c("#542788","#2166ac", "#b2182b")) +
theme_bw() +
labs(colour = "eQTL cohort", x = "Effect size (slope) in discovery cohort", y = "Effect size (slope) in Gtex cohort") +
geom_smooth(method = "lm", se = FALSE) +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.position = "top")
ggsave(args[7], width = 12, height = 12, dpi = 600)
#ggsave("cis_replication_correlation.pdf", width = 12, height = 12, dpi = 600)
all_corr <- final_point %>% filter(Group.x == "ALL") 
all_corr <- cor.test(all_corr$slope_ref, all_corr$slope, method = "spearman", exact = FALSE)

con_corr <- final_point %>% filter(Group.x == "CONTROL") 
con_corr <- cor.test(con_corr$slope_ref, con_corr$slope, method = "spearman", exact = FALSE)

INFECTED_corr <- final_point %>% filter(Group.x == "INFECTED") 
INFECTED_corr <- cor.test(INFECTED_corr$slope_ref, INFECTED_corr$slope, method = "spearman", exact = FALSE)

all_corr$estimate # 0.7607341
con_corr$estimate # 0.7799759
INFECTED_corr$estimate # 0.7829873

all_corr$p.value # 0
con_corr$p.value # 0
INFECTED_corr$p.value # 0
