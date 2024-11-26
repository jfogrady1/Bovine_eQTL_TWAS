######################################################
# Dominance and additive analysis for eQTLs
#####################################################

library(data.table)
library(tidyverse)
library(vcfR)

# First thing is to read in the data for top eQTLs

data <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt") %>% filter(is_eGene == TRUE)
vcf_all <- vcfR::read.vcfR("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz")



vcf_all <- cbind(vcf_all@fix, vcf_all@gt)
head(vcf_all)
vcf_all <- as.data.frame(vcf_all)


vcf_all[,colnames(vcf_all)[10:132]] <- lapply(vcf_all[,colnames(vcf_all)[10:132]], function(x) sub("0\\|0:.*", paste0("0"), x))
vcf_all[,colnames(vcf_all)[10:132]] <- lapply(vcf_all[,colnames(vcf_all)[10:132]], function(x) sub("0\\|1:.*", paste0("1"), x))
vcf_all[,colnames(vcf_all)[10:132]] <- lapply(vcf_all[,colnames(vcf_all)[10:132]], function(x) sub("1\\|0:.*", paste0("1"), x))
vcf_all[,colnames(vcf_all)[10:132]] <- lapply(vcf_all[,colnames(vcf_all)[10:132]], function(x) sub("1\\|1:.*", paste0("2"), x))


rownames(vcf_all) <- vcf_all$ID

vcf_all <- vcf_all %>% select(-c(1,2,3,4,5,6,7,8,9))

# Read in expression data set
expression = fread('/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.expr_tmm_inv.bed.gz')
expression <- as.data.frame(expression)
expression <- expression %>% select(4:length(colnames(expression)))
rownames(expression) <- expression$phenotype_id
expression <- expression %>% select(-1)
head(expression)



# Covariates
covariates_all = data.frame(t(read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.covariates.txt")))
covariate_names <- colnames(covariates_all)
covariate_names
formula <- as.formula(paste("Expression ~ Genotype +", paste(covariate_names, collapse = " + ")))
formula
covariate_names_dom <- c(covariate_names, "Genotype_dom")
formula_dom <- as.formula(paste("Expression ~ Genotype +", paste(covariate_names_dom, collapse = " + ")))
formula_dom




ALL_dom_results = data.frame(matrix(nrow = 1, ncol =8))
colnames(ALL_dom_results) = c("gene", "variant_id", "n_homo_ref", "n_het", "n_homo_alt", "additive_pval", "dominnat_pval", "LRT_P")
for (assoc in 1:length(rownames(data))) {
    print(assoc)
    data_temp = data[assoc,]
    gene = data_temp$phenotype_id
    variant_id = data_temp$variant_id
    Expression = as.numeric(expression[gene,])
    Genotype = as.numeric(vcf_all[variant_id,])
    print(variant_id)
    Genotype_dom = as.numeric(vcf_all_dom[variant_id,])
    print(Genotype_dom)
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

View(ALL_dom_results)
ALL_dom_results = ALL_dom_results[-1,]
ALL_dom_results$LRT_adjust = p.adjust(ALL_dom_results$LRT_P, method = "BH")

n_homo_ref
n_homo_alt
n_het


alias(lm(formula = formula_dom, data = lm_data_dom))



head(pval)
data <- cbind(data, pval)
data %>% select(pval, pval_nominal) %>% View()
all(pval == data$pval_nominal)
data_temp
vcf_temp







head(covariates_all)

Expression = as.numeric(expression["ENSBTAG00000020035",])
Genotype = as.numeric(vcf_all["1:872341:T:C",])
Genotype_dom = as.numeric(vcf_all_dom["1:872341:T:C",])
Genotype_dom

covariates_all
covariates_all = cbind(Expression, Genotype, covariates_all)


test = lm(formula, data = covariates_all)
test
covariates_all = cbind(covariates_all, Genotype_dom)
covariates_all
formula
Genotype_dom
test_dominant = lm(Expression ~ Genotype + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
    PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + genotypePC1 + 
    genotypePC2 + genotypePC3 + genotypePC4 + genotypePC5 + Condition + 
    Age + Genotype_dom, data = covariates_all)
test_dominant


summary(test)$coefficients["Genotype",][4]


as.numeric(lrtest(test, test_dominant)["Pr(>Chisq)"][2,])

formula
 njup = as.vector(test$coefficients["Genotype",][4])
p

vcf_all_dom <- data.frame(lapply(vcf_all[,colnames(vcf_all)], function(x) sub("2", paste0("0"), x)))
rownames(vcf_all_dom) <- rownames(vcf_all)
head(data)

library(lmtest)

lrtest()

add_dom <- lm(as.numeric(expression["ENSBTAG00000020035",]) ~ as.numeric(vcf_all["1:872341:T:C",] ) + as.numeric(vcf_all_dom["1:872341:T:C",]))
add <- lm(as.numeric(expression["ENSBTAG00000020035",]) ~ as.numeric(vcf_all["1:872341:T:C",] ))






#lrtest(nested,complex)
result_temp <- lrtest(add,add_dom)
result_temp$chiq

as.numeric(expression["ENSBTAG00000020035",])

vcf_control <- vcfR::read.vcfR(args[15])
vcf_infected <- vcfR::read.vcfR(args[16])



