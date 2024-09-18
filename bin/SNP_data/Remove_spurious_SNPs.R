# R script to determine what variants are correctly remapped and which variants need
# to be strand flipped.

# input files
# ThermoFisher meta data file - variants mapped to UMD3.1 and also pairs
# Updated MAP file: variants remapped to ARSUCD 1.2 using liftover
# Rob Schnables Coordiante files and ref/alt files and just the ref files

# Rob reference allele determines whether we need to flip the alleles
# THe Thermofisher file provides an insight into the alleles present
# WIll do a final check with robs allele pairs and will match up any discrepancies
# The ref/alt file from rob schnabel will determine the final correct alternative allele

library(tidyverse)
library(Biostrings)
library(seqinr)
library(dplyr)

# ------------------------------------------------------------------------------
# File input specification

args = commandArgs(trailingOnly=TRUE)
#1 ThermoFisher file

args[1]
TF <- read.csv(args[1], skip = 16, header = T)
TF <- TF %>% select(1,2,5,6,8,12,13)
colnames(TF) <- c("AX_id","TF_affy_ID", "TF_chr", "TF_pos", "TF_strand", "TF_ref", "TF_alt" )


# 2. Updated map file
args[2]
LOver <- as.data.frame(read.csv(args[2], sep = "\t", header = F))
colnames(LOver) <- c("Lover_ID", "Lover_chr", "Lover_pos")


# 3. Rob Schnabel_coordiantes
args[3]
Rob <- as.data.frame(read.csv(args[3], sep = "\t", header = F))
Rob <- Rob %>% select(2,1,4)
colnames(Rob) <- c("Rob_ID", "Rob_chr", "Rob_pos")


# 4. Rob Schnabel ref allele
args[4]
Rob_ref <- as.data.frame(read.csv(args[4], sep = "\t", header = F))
colnames(Rob_ref) <- c("Rob_ref_id", "Rob_ref_allele")




# Perform the merging
master <- left_join(TF, LOver, by = c("TF_affy_ID" = "Lover_ID")) %>% 
  left_join(., Rob, by = c("TF_affy_ID" = "Rob_ID")) %>% 
  left_join(., Rob_ref, by = c("TF_affy_ID" = "Rob_ref_id")) %>% drop_na()




# -------------------------------------------------------------------------------
# Perform filtering 

master$Lover_chr[master$Lover_chr == "X"] <- "31" # change X
master$Lover_chr <- as.numeric(master$Lover_chr) # change to numeric


# Filter for instances where LiftOver Chr and Schnabels Chr and Pos Match - key point
master_filtered <-  master %>% filter((Lover_chr == Rob_chr) &
                                        (Lover_pos == Rob_pos)) # filter to match chr and positions in both methods




# Removing poorly annotated SNPS - impossible to determine which strand they are on when remapping
# or if they require themselves to be flipped
# There are only about 20K of these so they can be discarded
master_remove_IDs <- master_filtered %>% filter((TF_ref == "A" & TF_alt == "T") |
                                                  (TF_ref == "T" & TF_alt == "A") |
                                                  (TF_ref == "C" & TF_alt == "G") |
                                                  (TF_ref == "G" & TF_alt == "C")) %>% select(2)




# Remove the above IDs from the dataset
master_final <- master_filtered %>% anti_join(., master_remove_IDs, by = c("TF_affy_ID"))

dim(master_remove_IDs)
write.table(master_remove_IDs, file = args[5], row.names = F, col.names = T)
# ------------------------------------------------------------------------------
# File to update coordinates
# We know the variants to retain and their updated coordinates, will now write these
# to a file which will be used in PLINK to update their respective coordinates
# 1. final snps
ids <-  master_final %>% select(1)

write.table(ids, file = args[6], col.names = F, quote = F, sep = "\t", row.names = F)

# 2. Updating coordinates
new_map <- master_final %>% select(1, 8, 9)
write.table(new_map, file = args[7], col.names = F, quote = F, sep = "\t", row.names = F)

#3. write master_file - required for next process

write.table(master_final, file = args[8], sep = "\t", quote = F, row.names = F, col.names = T)