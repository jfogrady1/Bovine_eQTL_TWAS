# Script to match eQTL results, mediators QTL and GWAS data for TWAS

#libraries
library("data.table")
library(tidyverse)
args = commandArgs(trailingOnly = TRUE)
# Read in significant GWAS results
thresholds = fread(args[1])
thresholds <- thresholds %>% select(1,22,27,28,29)
colnames(thresholds) <- c("gene", "snp","p-value", "beta", "beta_se")
thresholds <- thresholds %>% select("snp", "gene", "beta", "p-value")
thresholds$FDR <- thresholds$`p-value` # dummy column for FDR
thresholds$`t-stat` <- 0 # dummy column
thresholds <- thresholds %>% select("gene", "snp","beta","t-stat", "p-value", "FDR")


# Now read in correlation results
correlations <- fread(args[2])

# Now bring in the GWAS
GWAS <- fread(args[3])
GWAS <- GWAS %>% select(-1)



# Now have SNPs in reference and GWAS panels which match
# Now filter to only include mediators with genes in reference panel
# Need to keep those with association in cis and which have GWAS data
correlations <- correlations %>% filter(mediator %in% thresholds$gene)
correlations <- correlations %>% filter(gene %in% thresholds$gene)
table(correlations$mediator %in% thresholds$gene)
table(correlations$gene %in% thresholds$gene)

# now get the snp positions
snp_pos <- data.frame(thresholds$snp)
snp_pos$snpid <- snp_pos$thresholds.snp
snp_pos <- snp_pos %>% separate(., snpid, into = c("chr", "pos"), sep = ":")
colnames(snp_pos) <- "snpid"

# Now get expression of mediators which are retained
exp = fread(args[4])
#TF = fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/Bos_taurus_TF.txt")
#coTF = fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/Bos_taurus_Cof.txt")

# Mediator locations
mediators_exp <- exp %>% filter(phenotype_id %in% correlations$mediator)

mediators_loc <- mediators_exp %>% select(4,1,2)
colnames(mediators_loc) <- c("snpid", "chr", "pos")

# Mediators intensity
mediators_exp <- mediators_exp %>% select(-c(1,2,3))
colnames(mediators_exp)[1] <- "Mediator"


# Finally read in the gene expression data 
exp = fread(args[4])
gene_exp <- exp %>% filter(phenotype_id %in% thresholds$gene)

gene_loc <- gene_exp %>% select(4,1,2)
colnames(gene_loc) <- c("snpid", "chr", "pos")

# gene intensity
gene_exp <- gene_exp %>% select(-c(1,2,3))
colnames(gene_exp)[1] <- "Mediator"

# write to files
# cis-eQTLs
# snp pos
# mediator-qtls
# GWAS
# mediator intensities
# mediator positions
write.table(thresholds, file = args[5], sep = "\t", col.names = T, row.names = F, quote = F)
write.table(snp_pos, file = args[6], sep = "\t", col.names = T, row.names = F, quote = F)
write.table(correlations, file = args[7], sep = "\t", col.names = T, row.names = F, quote = F)
write.table(GWAS, file = args[8], sep = "\t", col.names = T, row.names =F, quote = F)
write.table(mediators_loc, file = args[9], sep = "\t", col.names = T, row.names = F, quote = F)
write.table(mediators_exp, file = args[10], sep = "\t", col.names = T, row.names =F, quote = F)
write.table(gene_loc, file = args[11], sep = "\t", col.names = T, row.names =F, quote = F )
write.table(gene_exp, file = args[12], sep = "\t", col.names = T, row.names =F, quote = F )

# final SNP file to extract from main vcf before TWAS
write.table(snp_pos$snpid, file = args[13], sep = "\t", col.names = F, row.names =F, quote = F) # no colname for plink
dim(thresholds)
