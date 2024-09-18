# Load in important packages
# Necessary packages
library(tidyverse)
library(tibble)
library(DESeq2)
library(vsn)
args = commandArgs(trailingOnly=TRUE)

# Read in the counts
counts_raw <- read.table(args[1], header = T) # nolint
Geneid <- counts_raw[,1]

# get samples
all_samples <- colnames(counts_raw)[2:length(colnames(counts_raw))]
TB_samples <- subset(all_samples, grepl(glob2rx("T0*"), all_samples)) #note this will need to be changed if different nomenclature used
Con_samples <- subset(all_samples, grepl(glob2rx("C0*"), all_samples))
TB_samples

# arg6 = annotation file
args[2]
annotation = read.csv(args[2], sep = "\t", header =F, skip = 5) #%>% select(1,3,6,7,8)
annotation <- annotation %>% separate(., V9, into = paste0("V", 9:18), sep = "; ") %>% filter(V3 == "gene") %>% select(V9, V1, V4, V5)
annotation$V9 <- gsub("gene_id ", "", annotation$V9)
colnames(annotation) <- c("Geneid", "chr", "Start_location", "End_location")
annotation$length = abs(annotation$Start_location - annotation$End_location)                         
colnames(annotation)[1] <- "Geneid"
counts_raw <- left_join(counts_raw, annotation) # join to raw counts as we will need to perform TPM on these


###############################
# TPM normalisation
###############################
counts_rpk <-  counts_raw[, 2:124] / (counts_raw$length) # divide by length
rownames(counts_rpk) <- counts_raw$Geneid # save gene names
counts_rpk <- drop_na(counts_rpk) # remove those with no start/end position
tpm.mat <- as.data.frame(t( t(counts_rpk) * 1e6 / colSums(counts_rpk) ) ) # normlaise by TPM

Con_samples_raw <- c("Geneid", Con_samples)
TB_samples_raw <- c("Geneid", TB_samples)
control_raw <- counts_raw %>% select(all_of(Con_samples_raw))
TB_raw <- counts_raw %>% select(all_of(TB_samples_raw))


control_tpm <- tpm.mat %>% select(all_of(Con_samples))
TB_tpm <- tpm.mat %>% select(all_of(TB_samples))
head(control_tpm)
write.table(counts_raw, file = args[3], sep = "\t", row.names = T, col.names = T)
write.table(control_raw, file = args[4], sep = "\t", row.names = T, col.names = T)
write.table(TB_raw, file = args[5], sep = "\t", row.names = T, col.names = T)
write.table(tpm.mat, file = args[6], sep = "\t", row.names = T, col.names = T)
write.table(control_tpm, file = args[7], sep = "\t", row.names = T, col.names = T)
write.table(TB_tpm, file = args[8], sep = "\t", row.names = T, col.names = T)
