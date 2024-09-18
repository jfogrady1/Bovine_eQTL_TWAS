# Assoications between mediators and genes
# Critical for MOSTWAS power - also very poorly described in the protocol
# will use correlation analysis and correct for multiple testing using BH method
# Will order by p value then and only inlcude top 5 in MOSTWAS predictive models.



library(tidyverse)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
# Expression data 
expression_raw <- fread(args[1])
head(expression_raw)

# Expression location of the all data 
# Note will need to subset out mediators and genes locations
gene_loc <- expression_raw %>% select(4,1,2)
TFs <- read.csv(args[2], sep = "\t", header = T) %>% dplyr::select(3)
CoF <- read.csv(args[3], sep = "\t", header = T) %>% dplyr::select(3)

# specify colnames
colnames(gene_loc)[1] <- "id"
colnames(gene_loc)[2] <- "Chromosome"



# Final set of mediators
# need to find those which are expressed
mediators <- rbind(TFs, CoF)
mediators <- left_join(mediators, gene_loc, by = c("Ensembl" = "id")) #%>% drop_na()
mediators <- mediators %>% filter(!is.na(mediators$Chromosome)) # exlcude those not on autosomes

mediators <- mediators[!duplicated(mediators$Ensembl),] # do not consider duplicates

# Isolate those genes which are not TFs / CoTFs
expression <- expression_raw %>% select(-c(1:3)) 
colnames(expression)[1] <- "id"
mediators <- mediators[mediators$Ensembl %in% expression$id,]
expression <- expression %>% filter(!(expression$id %in% mediators$Ensembl))


# Get their locations
gene_loc <- gene_loc %>% filter(!(gene_loc$id %in% mediators$Ensembl))
colnames(mediators)[1] <- "id"
mediators <- mediators %>% dplyr::select(1,2,3)
colnames(mediators)[3] <- "pos"

# write location
# args
write.table(mediators, file = args[4], sep = "\t", row.names = T, quote = F, col.names = T)

# get expression data (dosage for mediators)
# Getting their intensities
mediators <- mediators %>% dplyr::select(1)
mediators <- left_join(mediators, expression_raw[,-c(1:3)], by = c("id" = "phenotype_id"))

write.table(mediators, file = args[5], sep = "\t", row.names = T, quote = F, col.names = T)
# Convert to matrices
rownames(mediators) <- mediators[,1]
mediators <- mediators[,-1]
mediators <- as.matrix(mediators)
expression <- as.matrix(expression)
head(expression)
rownames(expression) <- expression[,1]
expression <- expression[,-1]
expression <- as.matrix(expression)
expression_mat <- apply(expression, 2, as.numeric)
rownames(expression_mat) <- rownames(expression)

# Correlation
cors <- t(cor(t(mediators), t(expression_mat)))
dim(cors)
df <- ncol(expression_mat) - 2
t_vals <- cors * sqrt(df) / sqrt(1 - cors ^ 2)

p_vals <- 2 * pt(-abs(t_vals), df = 123 - 2)
p_vals <- as.vector(p_vals)
cors_vec <- as.vector(cors)
bh_vals <- p.adjust(p_vals, method = "BH")
rows <- dim(mediators)[1] * dim(expression_mat)[1]
rows
final_df <- data.frame(matrix(nrow = rows , ncol = 6))
colnames(final_df) <- c("mediator", "gene", "t", "r", "p", "fdr")
final_df$mediator <- rep(rownames(mediators), each = dim(expression_mat)[1])
final_df$gene <- rep(rownames(expression_mat), dim(mediators)[1])
final_df$t <- as.vector(t_vals)
final_df$p <- as.vector(p_vals)
final_df$r <- cors_vec

final_df_corrected <- final_df %>% mutate(fdr = p.adjust(p, method = "BH"))
final_df_corrected <- final_df_corrected %>% filter(fdr < 0.01)
final_df_corrected <- final_df_corrected[order(final_df_corrected$fdr, decreasing = F),]
head(final_df_corrected)
write.table(final_df_corrected, args[6], row.names = F, col.names = T, quote = F, sep = "\t")

