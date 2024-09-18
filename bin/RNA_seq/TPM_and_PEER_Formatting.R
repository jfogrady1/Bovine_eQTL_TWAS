# Script for TPM normalization and getting 

## 1: Read in the data

## 2. Perform batch correction, save as adjusted counts

## 3. Perform TPM on RAW counts and select genes which pass threshold

## 4. Perform DESEQ2 and vst on adjusted (batch corrected) counts controlling for condition

## 5. Merge PCA with covariate data to get covariate matrix for PEER

# 3.  Output DESEQ2 normalised counts (control, infected, all) separetly for PEER.



# Necessary packages
library(tidyverse)
library(tibble)
library(sva)
library(DESeq2)
library(vsn)
library(RNOmni)
args = commandArgs(trailingOnly=TRUE)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Input data format ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# arg 1 = counts
counts_raw <- read.table(args[1], header = T)
counts_raw$Geneid <- gsub("GeneID:", "", counts_raw$Geneid) %>% as.numeric() # will be numeric from now on
new_order = sort(colnames(counts_raw[2:length(colnames(counts_raw))]))
counts_raw <- counts_raw[,c("Geneid", new_order)]

Geneid <- counts_raw[,1]

# get samples
all_samples <- colnames(counts_raw)[2:length(colnames(counts_raw))]
TB_samples <- subset(all_samples, grepl(glob2rx("T0*"), all_samples))
Con_samples <- subset(all_samples, grepl(glob2rx("C0*"), all_samples))
all_samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Covariate data formatting & batch correction  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# arg 2 - covariate data base
covariate_data <- read.table(args[2], header = T)
covariate_data$Age <- as.numeric(covariate_data$Age)

# CombatSeq set up
group = covariate_data$Condition
batch = covariate_data$Batch
counts <- counts_raw %>% select(-1)
counts <- as.matrix(counts)

cols = colnames(counts)
cols = sort(cols)

# order the samples 
counts <- counts[,cols]
counts_raw <- counts_raw[,cols]

covariate_data <- covariate_data[order(covariate_data$Sample),]


counts <- ComBat_seq(counts, batch=batch, group=group)




Con_covariate <- covariate_data[covariate_data[,1] %in% Con_samples,]
TB_covariate <- covariate_data[covariate_data[,1] %in% TB_samples,]

# PCA data = Arg 3 - all PCA = list sorted and pasted
All_PCA <- read.table(args[3], header = T)
# ARG4 = Control -PCA
Con_PCA <- read.table(args[4], header = T)
# ARG5 = infected PCA
TB_PCA <- read.table(args[5], header = T)


# Now merge together
covariate_data <- left_join(covariate_data, All_PCA) # join by sample, PCA file has this id
Con_covariate <- left_join(Con_covariate, Con_PCA)
TB_covariate <- left_join(TB_covariate, TB_PCA)

# Can now take out batch variable as this has been corrected for using COMBAT SEQ
covariate_data <- covariate_data %>% select(-c("Batch"))
Con_covariate <- Con_covariate %>% select(-c("Batch"))
TB_covariate <- TB_covariate %>% select(-c("Batch"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DESEQ2 Normalisation  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ALL
## Coldata formatting
coldata <- covariate_data
coldata <- as.data.frame(coldata)
rownames(coldata) <- covariate_data$Sample
coldata$Condition <- as.factor(coldata$Condition)
coldata <- coldata %>% select(-1) # Get rid of sample ID column - needs to be a matrix
counts <- counts[,colnames(counts) %in% rownames(coldata)] # ensure correct smaples
all(rownames(coldata) == colnames(counts))

dds_ALL <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = coldata,
                                  design = ~ 1)
dds_ALL <- DESeq(dds_ALL) # normalising
dds_counts <- counts(dds_ALL, normalized=TRUE)
# Control
CONTROL_counts <- counts[,colnames(counts) %in% Con_samples]
coldata <- Con_covariate
rownames(coldata) <- Con_covariate$Sample
coldata$Condition <- as.factor(coldata$Condition)
coldata <- coldata %>% select(-1) # Get rid of sample ID column - needs to be a matrix
all(rownames(coldata) == colnames(CONTROL_counts))
dds_CONTROL <- DESeqDataSetFromMatrix(countData = CONTROL_counts,
                                      colData = coldata,
                                      design = ~ 1) # Need something in the design formula
dds_CONTROL <- DESeq(dds_CONTROL) # normalising
dds_con_counts <- counts(dds_CONTROL, normalized=TRUE)

# INFECTED
INFECTED_counts <- counts[,colnames(counts) %in% TB_samples]
coldata <- TB_covariate
rownames(coldata) <- TB_covariate$Sample
coldata$Condition <- as.factor(coldata$Condition)
coldata <- coldata %>% select(-1) # Get rid of sample ID column - needs to be a matrix
all(rownames(coldata) == colnames(INFECTED_counts))
dds_INFECTED <- DESeqDataSetFromMatrix(countData = INFECTED_counts,
                                       colData = coldata,
                                       design = ~ 1) # Need something in the design formula
dds_INFECTED <- DESeq(dds_INFECTED) # normalising
dds_tb_counts <- counts(dds_INFECTED, normalized=TRUE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Remove genes with no location info  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# arg6 = annotation file
annotation = read.csv(args[6]) %>% select(1,3,6,7,8)
annotation$length = abs(annotation$Start_location - annotation$End_location)
annotation$Gene_ID <- as.numeric(annotation$Gene_ID)
colnames(annotation)[1] <- "Geneid"
counts_raw <- cbind(Geneid, counts_raw)
counts_raw <- left_join(counts_raw, annotation) # join to raw counts as we will need to perform TPM on these
#dim(counts_raw) # [1] 30574     8







# remove genes with no genomic location
counts_raw <- counts_raw[!is.na(counts_raw$length),]
length <- as.data.frame(counts_raw$length)
rownames(length) <- counts_raw$Geneid
#dim(counts_raw) # [1] 30073     9


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Gtex raw filter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Filtering with raw exression
threshold = as.integer(length(all_samples)) * 0.2
rownames(counts_raw) <- counts_raw$Geneid
counts <- counts_raw %>% select(all_of(all_samples))
rownames(counts) <- rownames(counts_raw)
rows <- apply(counts,1,function(x) sum(x >= 6))
counts <-cbind(counts, rows)
counts <- counts %>% filter(rows > threshold) %>% select(-c(rows))
length <- length %>% filter(rownames(length) %in% rownames(counts))
counts <- cbind(counts, length)

# ensure rows are equal
table(rownames(length) == rownames(counts))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TPM Normalisation and filtering ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

counts_rpk <-  counts[, 1:(as.integer(length(all_samples)) + 1)] / counts$`counts_raw$length`
tpm.mat <- as.data.frame(t( t(counts_rpk) * 1e6 / colSums(counts_rpk) ) )
tpm.mat <- tpm.mat %>% select(-c("counts_raw$length"))


# Filtering for TPM normalisation
rows <- apply(tpm.mat,1,function(x) sum(x >= 0.1))
tpm.mat <- cbind(tpm.mat, rows)
tpm.mat <- tpm.mat %>% filter(rows > threshold) %>% select(-c(rows))






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Write outputs to respective files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



# PEER DATAFRAMES = Raw expression
rownames(dds_counts) <- Geneid
rownames(dds_con_counts) <- Geneid
rownames(dds_tb_counts) <- Geneid
head(dds_tb_counts)
All_peer <-  dds_counts %>% as.data.frame(.) %>% select(all_of(all_samples)) %>% filter(rownames(.) %in% rownames(tpm.mat)) %>% as.data.frame()
Con_peer <- dds_con_counts %>% as.data.frame(.) %>% select(all_of(Con_samples)) %>% filter(rownames(.) %in% rownames(tpm.mat)) %>% as.data.frame()
TB_peer <- dds_tb_counts %>% as.data.frame(.) %>% select(all_of(TB_samples)) %>% filter(rownames(.)%in% rownames(tpm.mat)) %>% as.data.frame()

# inverse normal transformation
INT <- function(x) {
  rows <- rownames(x)
  X_INT = data.frame(matrix(nrow = 0, ncol = ncol(x)))
  colnames(X_INT) <- colnames(x) # empty data frame with same size
  for (row in 1:nrow(x)) {
    norm <- RNOmni::RankNorm(as.vector(as.numeric(x[row,])))
    X_INT <- rbind(X_INT, norm)
  }
  print(dim(X_INT))
  rownames(X_INT) <- rows
  colnames(X_INT) <- colnames(x)
  return(X_INT)
}

All_Tran <- INT(All_peer)
TB_Tran <- INT(TB_peer)
Con_Tran <- INT(Con_peer)


All_peer <- t(All_Tran) %>% as.data.frame()
Con_peer <- t(Con_Tran) %>% as.data.frame()
TB_peer <- t(TB_Tran) %>% as.data.frame()


counts_DESEQ <- as.data.frame(t(All_peer))
Con_counts_DESEQ <- as.data.frame(t(Con_peer))
TB_counts_DESEQ <- as.data.frame(t(TB_peer))


# format DESEQ datasets
# need to remove row number from these
counts_DESEQ$Geneid <- rownames(counts_DESEQ)
counts_DESEQ <- counts_DESEQ %>% select(all_of(c("Geneid", all_samples)))
nrow(counts_DESEQ)
rownames(counts_DESEQ) <- 1:nrow(counts_DESEQ)

Con_counts_DESEQ$Geneid <- rownames(Con_counts_DESEQ)
Con_counts_DESEQ <- Con_counts_DESEQ %>% select(all_of(c("Geneid", Con_samples)))
rownames(Con_counts_DESEQ) <- 1:nrow(Con_counts_DESEQ)

TB_counts_DESEQ$Geneid <- rownames(TB_counts_DESEQ)
TB_counts_DESEQ <- TB_counts_DESEQ %>% select(all_of(c("Geneid", TB_samples)))
rownames(TB_counts_DESEQ) <- 1:nrow(TB_counts_DESEQ)


# Gene Location dataset
location <- counts_raw %>% select(Geneid,Chromosome,Start_location,End_location) %>% filter(Geneid %in% counts_DESEQ$Geneid)
rownames(location) <- 1:nrow(location)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Final Outputs  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Important to keep files with same nomenclature so they can be grouped as tuples in nextflow


# PEER FILES
write.table(All_peer, file = "ALL_PEER.txt", quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Con_peer, file = "CONTROL_PEER.txt", quote = F, row.names = T, col.names = T, sep = "\t")
write.table(TB_peer, file = "INFECTED_PEER.txt", quote = F, row.names = T, col.names = T, sep = "\t")


# DESEQ files
write.table(counts_DESEQ, file = "ALL_DESEQ.txt", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(Con_counts_DESEQ, file = "CONTROL_DESEQ.txt", quote = F, row.names = F, col.names = T, sep = "\t")
write.table(TB_counts_DESEQ, file = "INFECTED_DESEQ.txt", quote = F, row.names = F, col.names = T, sep = "\t")

# covariate files
rownames(covariate_data) <- covariate_data$Sample
covariate_data <- covariate_data %>% select(-c("Sample"))
covariate_data$Condition[covariate_data$Condition =="Control"] <- 0
covariate_data$Condition[covariate_data$Condition == "Infected"] <- 1
# remove Condition from each inidividual group
Con_covariate <- Con_covariate %>% select(-c("Condition"))
TB_covariate <- TB_covariate %>% select(-c("Condition"))
rownames(Con_covariate) <- Con_covariate$Sample
rownames(TB_covariate) <- TB_covariate$Sample
Con_covariate <- Con_covariate %>% select(-c("Sample"))
TB_covariate <- TB_covariate %>% select(-c("Sample"))
write.table(covariate_data, file = "ALL_COV.txt", quote = F, row.names = T, col.names = T, sep = "\t")
write.table(Con_covariate, file = "CONTROL_COV.txt", quote = F, row.names = T, col.names = T, sep = "\t")
write.table(TB_covariate, file = "INFECTED_COV.txt", quote = F, row.names = T, col.names = T, sep = "\t")

# Gene location files
write.table(location, file = "gene_location.txt", quote = F, sep = "\t", row.names = T, col.names = T)