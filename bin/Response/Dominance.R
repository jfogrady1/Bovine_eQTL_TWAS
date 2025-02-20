######################################################
# Dominance and additive analysis for eQTLs
#####################################################

library(data.table)
library(tidyverse)
library(vcfR)
library(lmtest)

args = commandArgs(trailingOnly = TRUE)
# First thing is to read in the data for top eQTLs
data <- fread(args[1]) %>% filter(is_eGene == TRUE)
data_independent <-  fread(args[2]) %>% filter(phenotype_id %in% data$phenotype_id)


vcf_all <- vcfR::read.vcfR(args[3])



vcf_all <- cbind(vcf_all@fix, vcf_all@gt)
vcf_all <- as.data.frame(vcf_all)
head(vcf_all)
vcf_all <- vcf_all %>% filter(ID %in% data_independent$variant_id)

vcf_all[,colnames(vcf_all)[10:132]] <- lapply(vcf_all[,colnames(vcf_all)[10:132]], function(x) sub("0\\|0:.*", paste0("0"), x))
vcf_all[,colnames(vcf_all)[10:132]] <- lapply(vcf_all[,colnames(vcf_all)[10:132]], function(x) sub("0\\|1:.*", paste0("1"), x))
vcf_all[,colnames(vcf_all)[10:132]] <- lapply(vcf_all[,colnames(vcf_all)[10:132]], function(x) sub("1\\|0:.*", paste0("1"), x))
vcf_all[,colnames(vcf_all)[10:132]] <- lapply(vcf_all[,colnames(vcf_all)[10:132]], function(x) sub("1\\|1:.*", paste0("2"), x))

head(vcf_all)
vcf_all_dom <- vcf_all
head(vcf_all_dom)
vcf_all_dom[,colnames(vcf_all_dom)[10:132]] <- lapply(vcf_all_dom[,colnames(vcf_all_dom)[10:132]], function(x) sub("2", paste0("0"), x))

rownames(vcf_all) <- vcf_all$ID
rownames(vcf_all_dom) <- vcf_all_dom$ID
vcf_all <- vcf_all %>% select(-c(1,2,3,4,5,6,7,8,9))
vcf_all_dom <- vcf_all_dom  %>% select(-c(1,2,3,4,5,6,7,8,9))
# Read in expression data set
expression = fread(args[4])
expression <- as.data.frame(expression)
expression <- expression %>% select(4:length(colnames(expression)))
rownames(expression) <- expression$phenotype_id
expression <- expression %>% select(-1)
head(expression)



# Covariates
covariates_all = data.frame(t(read.table(args[5])))
covariate_names <- colnames(covariates_all)
covariate_names
formula <- as.formula(paste("Expression ~ Genotype +", paste(covariate_names, collapse = " + ")))
formula
covariate_names_dom <- c(covariate_names, "Genotype_dom")
formula_dom <- as.formula(paste("Expression ~ Genotype +", paste(covariate_names_dom, collapse = " + ")))
formula_dom


ALL_dom_results = data.frame(matrix(nrow = 1, ncol =8))
colnames(ALL_dom_results) = c("gene", "variant_id", "n_homo_ref", "n_het", "n_homo_alt", "additive_pval", "dominnat_pval", "LRT_P")
for (assoc in 1:length(rownames(data_independent))) {
    print(assoc)
    data_temp = data_independent[assoc,]
    gene = data_temp$phenotype_id
    variant_id = data_temp$variant_id
    Expression = as.numeric(expression[gene,])
    Genotype = as.numeric(vcf_all[variant_id,])
    Genotype_dom = as.numeric(vcf_all_dom[variant_id,])
    n_homo_ref = sum(Genotype == 0)
    n_het = sum(Genotype == 1)
    n_homo_alt = sum(Genotype == 2)

    lm_data = cbind(Expression, Genotype, covariates_all)
    lm_data_dom = cbind(lm_data, Genotype_dom)

    add = lm(formula, data = lm_data)
    dom = lm(formula_dom, data = lm_data_dom)
    
    if(n_homo_ref == 0 | n_homo_alt == 0) {
        print("Dominant model cannot be computed, not enough variation at locus")
        additive_pval = as.numeric(summary(add)$coefficients["Genotype",][4])
        dominant_pval = "Cannot_compute"
        LRT_P = "Cannot_compute"
    }
    else {
        additive_pval = as.numeric(summary(add)$coefficients["Genotype",][4])
        dominant_pval = as.numeric(summary(dom)$coefficients["Genotype_dom",][4])
        LRT_P =as.numeric(lrtest(add, dom)["Pr(>Chisq)"][2,])
    } 
    
    
    row = c(gene, variant_id, n_homo_ref, n_het, n_homo_alt, additive_pval, dominant_pval, LRT_P)
    ALL_dom_results = rbind(ALL_dom_results, row)
}

#rm(vcf_all_dom)

data <- fread(args[6]) %>% filter(is_eGene == TRUE)
data_independent <-  fread(args[7]) %>% filter(phenotype_id %in% data$phenotype_id)
head(data_independent)
# Read in SNP of interest
IFITM3 <- fread(args[8]) %>% filter(phenotype_id == "ENSBTAG00000019015") %>% filter(variant_id == "29:50294904:G:A") %>% as.data.frame() # For reviewer comment
IFITM3$beta_shape1 <- 0
IFITM3$beta_shape2 <- 0
IFITM3$true_df <- 0
IFITM3$pval_true_df <- 0
IFITM3$num_var <- 0
IFITM3$end_distance <- 0
IFITM3$pval_perm <- 0
IFITM3$pval_beta <- 0
IFITM3$rank <- 0

IFITM3 <- IFITM3[,colnames(data_independent)]

RGS10 <- fread(args[9]) %>% filter(phenotype_id == "ENSBTAG00000002647") %>% filter(variant_id == "26:40304855:C:A") %>% as.data.frame() # For reviewer comment
RGS10$beta_shape1 <- 0
RGS10$beta_shape2 <- 0
RGS10$true_df <- 0
RGS10$pval_true_df <- 0
RGS10$num_var <- 0
RGS10$end_distance <- 0
RGS10$pval_perm <- 0
RGS10$pval_beta <- 0
RGS10$rank <- 0

RGS10 <- RGS10[,colnames(data_independent)]


data_independent = rbind(data_independent, RGS10)

vcf_con <- vcfR::read.vcfR(args[10])


vcf_con <- cbind(vcf_con@fix, vcf_con@gt)
vcf_con <- as.data.frame(vcf_con)
vcf_con <- vcf_con %>% filter(ID %in% data_independent$variant_id)

vcf_con[,colnames(vcf_con)[10:length(colnames(vcf_con))]] <- lapply(vcf_con[,colnames(vcf_con)[10:length(colnames(vcf_con))]], function(x) sub("0\\|0:.*", paste0("0"), x))
vcf_con[,colnames(vcf_con)[10:length(colnames(vcf_con))]] <- lapply(vcf_con[,colnames(vcf_con)[10:length(colnames(vcf_con))]], function(x) sub("0\\|1:.*", paste0("1"), x))
vcf_con[,colnames(vcf_con)[10:length(colnames(vcf_con))]] <- lapply(vcf_con[,colnames(vcf_con)[10:length(colnames(vcf_con))]], function(x) sub("1\\|0:.*", paste0("1"), x))
vcf_con[,colnames(vcf_con)[10:length(colnames(vcf_con))]] <- lapply(vcf_con[,colnames(vcf_con)[10:length(colnames(vcf_con))]], function(x) sub("1\\|1:.*", paste0("2"), x))


vcf_con_dom <- vcf_con
vcf_con_dom[,colnames(vcf_con_dom)[10:length(colnames(vcf_con))]] <- lapply(vcf_con_dom[,colnames(vcf_con_dom)[10:length(colnames(vcf_con))]], function(x) sub("2", paste0("0"), x))

rownames(vcf_con) <- vcf_con$ID
rownames(vcf_con_dom) <- vcf_con_dom$ID
vcf_con <- vcf_con %>% select(-c(1,2,3,4,5,6,7,8,9))
vcf_con_dom <- vcf_con_dom  %>% select(-c(1,2,3,4,5,6,7,8,9))
# Read in expression data set
expression = fread(args[11])
expression <- as.data.frame(expression)
expression <- expression %>% select(4:length(colnames(expression)))
rownames(expression) <- expression$phenotype_id
expression <- expression %>% select(-1)
head(expression)



# Covariates
covariates_con = data.frame(t(read.table(args[12])))
covariate_names_con <- colnames(covariates_con)
covariate_names_con
formula <- as.formula(paste("Expression ~ Genotype +", paste(covariate_names_con, collapse = " + ")))
formula
covariate_names_con_dom <- c(covariate_names_con, "Genotype_dom")
formula_con_dom <- as.formula(paste("Expression ~ Genotype +", paste(covariate_names_con_dom, collapse = " + ")))
formula_con_dom


head(data_independent)


CON_dom_results = data.frame(matrix(nrow = 1, ncol =8))
colnames(CON_dom_results) = c("gene", "variant_id", "n_homo_ref", "n_het", "n_homo_alt", "additive_pval", "dominnat_pval", "LRT_P")
for (assoc in 1:length(rownames(data_independent))) {
    print(assoc)
    data_temp = data_independent[assoc,]
    gene = data_temp$phenotype_id
    variant_id = data_temp$variant_id
    Expression = as.numeric(expression[gene,])
    Genotype = as.numeric(vcf_con[variant_id,])
    Genotype_dom = as.numeric(vcf_con_dom[variant_id,])
    n_homo_ref = sum(Genotype == 0)
    n_het = sum(Genotype == 1)
    n_homo_alt = sum(Genotype == 2)

    lm_data = cbind(Expression, Genotype, covariates_con)
    lm_data_dom = cbind(lm_data, Genotype_dom)

    add = lm(formula, data = lm_data)
    dom = lm(formula_con_dom, data = lm_data_dom)
    
    if(n_homo_ref == 0 | n_homo_alt == 0) {
        print("Dominant model cannot be computed, not enough variation at locus")
        additive_pval = as.numeric(summary(add)$coefficients["Genotype",][4])
        dominant_pval = "Cannot_compute"
        LRT_P = "Cannot_compute"
    }
    else {
        additive_pval = as.numeric(summary(add)$coefficients["Genotype",][4])
        dominant_pval = as.numeric(summary(dom)$coefficients["Genotype_dom",][4])
        LRT_P =as.numeric(lrtest(add, dom)["Pr(>Chisq)"][2,])
    } 
    
    
    row = c(gene, variant_id, n_homo_ref, n_het, n_homo_alt, additive_pval, dominant_pval, LRT_P)
    CON_dom_results = rbind(CON_dom_results, row)
}


#rm(vcf_con_dom)

data <- fread(args[13]) %>% filter(is_eGene == TRUE)
data_independent <-  fread(args[14]) %>% filter(phenotype_id %in% data$phenotype_id)

IFITM3 <- fread(args[15]) %>% filter(phenotype_id == "ENSBTAG00000019015") %>% filter(variant_id == "29:50294904:G:A") %>% as.data.frame() # For reviewer comment
IFITM3$beta_shape1 <- 0
IFITM3$beta_shape2 <- 0
IFITM3$true_df <- 0
IFITM3$pval_true_df <- 0
IFITM3$num_var <- 0
IFITM3$end_distance <- 0
IFITM3$pval_perm <- 0
IFITM3$pval_beta <- 0
IFITM3$rank <- 0

IFITM3 <- IFITM3[,colnames(data_independent)]

RGS10 <- fread(args[16]) %>% filter(phenotype_id == "ENSBTAG00000002647") %>% filter(variant_id == "26:40304855:C:A") %>% as.data.frame() # For reviewer comment
RGS10$beta_shape1 <- 0
RGS10$beta_shape2 <- 0
RGS10$true_df <- 0
RGS10$pval_true_df <- 0
RGS10$num_var <- 0
RGS10$end_distance <- 0
RGS10$pval_perm <- 0
RGS10$pval_beta <- 0
RGS10$rank <- 0

RGS10 <- RGS10[,colnames(data_independent)]


data_independent = rbind(data_independent, RGS10)

vcf_reac <- vcfR::read.vcfR(args[17])
head(vcf_reac)
vcf_reac <- cbind(vcf_reac@fix, vcf_reac@gt)
vcf_reac <- as.data.frame(vcf_reac)
vcf_reac <- vcf_reac %>% filter(ID %in% data_independent$variant_id)
dim(vcf_reac)
vcf_reac[,colnames(vcf_reac)[10:length(colnames(vcf_reac))]] <- lapply(vcf_reac[,colnames(vcf_reac)[10:length(colnames(vcf_reac))]], function(x) sub("0\\|0:.*", paste0("0"), x))
vcf_reac[,colnames(vcf_reac)[10:length(colnames(vcf_reac))]] <- lapply(vcf_reac[,colnames(vcf_reac)[10:length(colnames(vcf_reac))]], function(x) sub("0\\|1:.*", paste0("1"), x))
vcf_reac[,colnames(vcf_reac)[10:length(colnames(vcf_reac))]] <- lapply(vcf_reac[,colnames(vcf_reac)[10:length(colnames(vcf_reac))]], function(x) sub("1\\|0:.*", paste0("1"), x))
vcf_reac[,colnames(vcf_reac)[10:length(colnames(vcf_reac))]] <- lapply(vcf_reac[,colnames(vcf_reac)[10:length(colnames(vcf_reac))]], function(x) sub("1\\|1:.*", paste0("2"), x))


vcf_reac_dom <- vcf_reac
head(vcf_reac)
vcf_reac_dom[,colnames(vcf_reac_dom)[10:length(colnames(vcf_reac))]] <- lapply(vcf_reac_dom[,colnames(vcf_reac_dom)[10:length(colnames(vcf_reac))]], function(x) sub("2", paste0("0"), x))
head(vcf_reac_dom)
rownames(vcf_reac) <- vcf_reac$ID
rownames(vcf_reac_dom) <- vcf_reac_dom$ID
vcf_reac <- vcf_reac %>% select(-c(1,2,3,4,5,6,7,8,9))
vcf_reac_dom <- vcf_reac_dom  %>% select(-c(1,2,3,4,5,6,7,8,9))
# Read in expression data set
expression = fread(args[18])
expression <- as.data.frame(expression)
expression <- expression %>% select(4:length(colnames(expression)))
rownames(expression) <- expression$phenotype_id
expression <- expression %>% select(-1)
head(expression)



# Covariates
covariates_reac = data.frame(t(read.table(args[19])))
covariate_names_reac <- colnames(covariates_reac)
covariate_names_reac
formula <- as.formula(paste("Expression ~ Genotype +", paste(covariate_names_reac, collapse = " + ")))
formula
covariate_names_reac_dom <- c(covariate_names_reac, "Genotype_dom")
formula_reac_dom <- as.formula(paste("Expression ~ Genotype +", paste(covariate_names_reac_dom, collapse = " + ")))
formula_reac_dom


head(data_independent)

REAC_dom_results = data.frame(matrix(nrow = 1, ncol =8))
colnames(REAC_dom_results) = c("gene", "variant_id", "n_homo_ref", "n_het", "n_homo_alt", "additive_pval", "dominnat_pval", "LRT_P")
for (assoc in 1:length(rownames(data_independent))) {
    print(assoc)
    data_temp = data_independent[assoc,]
    gene = data_temp$phenotype_id
    variant_id = data_temp$variant_id
    Expression = as.numeric(expression[gene,])
    Genotype = as.numeric(vcf_reac[variant_id,])
    Genotype_dom = as.numeric(vcf_reac_dom[variant_id,])
    n_homo_ref = sum(Genotype == 0)
    n_het = sum(Genotype == 1)
    n_homo_alt = sum(Genotype == 2)

    lm_data = cbind(Expression, Genotype, covariates_reac)
    lm_data_dom = cbind(lm_data, Genotype_dom)

    add = lm(formula, data = lm_data)
    dom = lm(formula_reac_dom, data = lm_data_dom)
    
    if(n_homo_ref == 0 | n_homo_alt == 0) {
        print("Dominant model cannot be computed, not enough variation at locus")
        additive_pval = as.numeric(summary(add)$coefficients["Genotype",][4])
        dominant_pval = "Cannot_compute"
        LRT_P = "Cannot_compute"
    }
    else {
        additive_pval = as.numeric(summary(add)$coefficients["Genotype",][4])
        dominant_pval = as.numeric(summary(dom)$coefficients["Genotype_dom",][4])
        LRT_P =as.numeric(lrtest(add, dom)["Pr(>Chisq)"][2,])
    } 
    
    
    row = c(gene, variant_id, n_homo_ref, n_het, n_homo_alt, additive_pval, dominant_pval, LRT_P)
    REAC_dom_results = rbind(REAC_dom_results, row)
}


## remove first rows and also the association with specific reviewer's SNPs as these are not eQTLs per say.

ALL_dom_results <- ALL_dom_results[-1,]
ALL_dom_results = ALL_dom_results %>% filter(dominnat_pval != "Cannot_compute") %>% mutate(LRT_P_adj = p.adjust(LRT_P, method = "BH"))

CON_dom_results <- CON_dom_results[-1,]
CON_dom_results = CON_dom_results %>% filter(dominnat_pval != "Cannot_compute") %>% mutate(LRT_P_adj = p.adjust(LRT_P, method = "BH"))


REAC_dom_results <- REAC_dom_results[-1,]
REAC_dom_results = REAC_dom_results %>% filter(dominnat_pval != "Cannot_compute") %>% mutate(LRT_P_adj = p.adjust(LRT_P, method = "BH"))


write.table(ALL_dom_results, file = args[20], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(CON_dom_results, file = args[21], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(REAC_dom_results, file = args[22], sep = "\t", quote = FALSE, row.names = FALSE)
