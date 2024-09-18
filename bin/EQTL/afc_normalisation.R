library(DESeq2)
library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
#args[1] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/Quantification/ALL_matrix_clean.txt"
#args[2] <- ""
#args[3] <- "ALL"
#args[4] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl.txt.gz"
#args[5] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.covariates.txt"
data <- fread(args[1]) %>% select(-1)

if (args[3] == "ALL") {
    rows <-data$Geneid
    head(data)
    data <- as.matrix(data[,-c(1,125,126,127,128)])
} else {
    rows <-data$Geneid
    head(data)
    data <- as.matrix(data[,-c(1)])
}

rownames(data) <- rows
head(data)

colData <- as.data.frame(matrix(ncol = 1, nrow = length(colnames(data))))
colnames(colData) <- "Test"
colData$Test <- "1"
colData <- as.matrix(colData)


dds <- DESeqDataSetFromMatrix(countData =data, 
                              colData = colData, 
                              design = ~ 1)

dds <- estimateSizeFactors(dds)

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)


log2_normalized_counts <- log2(normalized_counts + 1)


qtl <- fread(args[4])

qtl <- qtl %>% select(1) # get genes

# filter log2 for these genes
log2_normalised_counts <- log2_normalized_counts[qtl$phenotype_id,]
log2_nomalised_counts <- as.data.frame(log2_normalised_counts)
dim(log2_normalised_counts)
head(log2_normalised_counts)
all(rownames(log2_normalised_counts) == qtl$phenotype_id)
names <- colnames(log2_normalised_counts)
log2_normalised_counts <- as.data.frame(log2_normalised_counts)
log2_normalised_counts$phenotype_id <- qtl$phenotype_id
log2_normalised_counts <- log2_normalised_counts %>% select(phenotype_id, names)


#write.table(log2_normalized_counts, file = args[2], quote = FALSE, row.names = FALSE)



######################
######################
# covariates
######################
######################
data <- fread(args[5])
colnames(data)[1] <- c("id")

rownames(data) <- 1:nrow(data)
write.table(data, file = args[6], row.names = F, col.names = T, sep = "\t", quote = F)


# Get in bedfile information
bed = fread(args[7])
bed <- bed[,1:4]

log2_normalised_bed <- left_join(bed, log2_normalised_counts)

write.table(log2_normalised_bed, file = args[2], row.names = F, col.names = T, sep = "\t", quote = F)
