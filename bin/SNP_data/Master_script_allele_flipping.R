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

# Packages
setwd("C:/Users/crtuser/OneDrive/OneDrive - University College Dublin/Analysis/Imputation")
library(tidyverse)
#BiocManager::install("Biostrings")
#install.packages("seqinr")
library(Biostrings)
library(seqinr)
library(dplyr)

# ------------------------------------------------------------------------------
# File input specification


#1 ThermoFisher file

TF <- read.csv("Axiom_GW_Bos_SNP_1.na35.annot.csv", skip = 16, header = T)
TF <- TF %>% select(1,2,5,6,8,12,13)
head(TF)
dim(TF)
colnames(TF) <- c("AX_id","TF_affy_ID", "TF_chr", "TF_pos", "TF_strand", "TF_ref", "TF_alt" )
# 2. Updated map file

LOver <- as.data.frame(read.csv("../Remapping/ARS_liftOver_updated.map", sep = "\t", header = F))
colnames(LOver) <- c("Lover_ID", "Lover_chr", "Lover_pos")
dim(LOver)
head(LOver)

# 3. Rob Schnabel_coordiantes

Rob <- as.data.frame(read.csv("Rob_Schnabel/9913_ARS1.2_648875_BOS1_marker_name_180910.map", sep = "\t", header = F))
Rob <- Rob %>% select(2,1,4)
colnames(Rob) <- c("Rob_ID", "Rob_chr", "Rob_pos")
dim(Rob)

# 4. Rob Schnabel ref allele

Rob_ref <- as.data.frame(read.csv("Rob_Schnabel/9913_ARS1.2_648875_BOS1_marker_name_180910.REF_ALLELE", sep = "\t", header = F))
head(Rob_ref)                       
colnames(Rob_ref) <- c("Rob_ref_id", "Rob_ref_allele")
head(LOver)



# Perform the merging
master <- left_join(TF, LOver, by = c("TF_affy_ID" = "Lover_ID")) %>% 
  left_join(., Rob, by = c("TF_affy_ID" = "Rob_ID")) %>% 
  left_join(., Rob_ref, by = c("TF_affy_ID" = "Rob_ref_id")) %>% drop_na()




# -------------------------------------------------------------------------------
# Perform filtering 



master$Lover_chr[master$Lover_chr == "X"] <- "31" # change X
master$Lover_chr <- as.numeric(master$Lover_chr) # change to numeric


# Filter for instances where LiftOver Chr and Schnabels Chr and Pos Match
master_filtered <-  master %>% filter((Lover_chr == Rob_chr) &
                                     (Lover_pos == Rob_pos)) # filter to match chr and positions in both methods

head(master_filtered,100)

dim(master_filtered)
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
write.table(master_remove_IDs, file = "Palindromic_SNPs_removed.txt", row.names = F, col.names = T)
# ------------------------------------------------------------------------------
# File to update coordinates
# We know the variants to retain and their updated coordinates, will now write these
# to a file which will be used in PLINK to update their respective coordinates
dim(master_final)
# 1. final snps
ids <-  master_final %>% select(1)
head(ids)
write.table(ids, file = "Final_raw_IDs_4_analysis", col.names = F, quote = F, sep = "\t", row.names = F)

# 2. Updating coordinates
colnames(master_final)
new_map <- master_final %>% select(1, 8, 9)
head(new_map)
write.table(new_map, file = "ARS_LOver_Rob_final.map", col.names = F, quote = F, sep = "\t", row.names = F)

rm(TF)


# ----------------------------------------------------------------------------------------------
# VCF modification, extremely important

# Extract the header to write to the file after
vcf_header <- read.csv("../PLINK_files/VCF_4_last_imputation_prep.vcf", sep = "\t", header = F)
vcf_header <- vcf_header[1:33,1]

vcf_header

# Specify the vcf data which will be edited
vcf <- as.data.frame(read.csv("../PLINK_files/VCF_4_last_imputation_prep.vcf", sep = "\t", skip = 34)) # make sure to say colnames = T when writing
head(vcf)
colnames(vcf)[1] <- "#CHROM"
colnames(vcf)[4] <- "VCF_REF"
colnames(vcf)[5] <- "VCF_ALT"


# join vcf to master ID file
master_final <- master_final %>% left_join(., vcf, by = c("AX_id" = "ID")) %>% drop_na() #drop_na here as SNPS on X need ot be gone



dim(master_final)
# Determining which SNPs need to be flipped and respective genotypes
master_flip <- master_final %>% mutate(Flip = case_when(
  Rob_ref_allele == VCF_REF | Rob_ref_allele == comp(seq = VCF_REF, forceToLower = F) ~ "N", # do nothing as essentially the same genotype on both assemblied
  Rob_ref_allele == VCF_ALT | Rob_ref_allele == comp(seq = VCF_ALT, forceToLower = F) ~ "Y", # these genotypes will need to be flipped 0 - 1 and 1 - 0 
  VCF_ALT == "." & Rob_ref_allele == VCF_REF | Rob_ref_allele  == comp(seq = VCF_REF, forceToLower = F) ~ "N", # Need to include this as PLINK messes everything up
  VCF_ALT == "." & Rob_ref_allele != VCF_REF | Rob_ref_allele != comp(seq = VCF_REF, forceToLower = F) ~ "Y")) # Don't know if these are on the same strand - missing


sum(table(master_flip$Flip)) == dim(master_flip)[1] # check that all conditions satisfied
table(master_flip$Flip) # see how many we need to FLIP
write.table(master_flip, file = "Supplementary_Table_Flipping.csv", col.names = T, row.names = F, sep = ",")

# Perform the flipping
# Many if_else statemnts, reason for so many is because need to account for "." allele which is an issue with PLINK
master_flip_IDs <- master_flip %>%  
  
  # 1. change the VCF_alt allele
  mutate(VCF_ALT = case_when(Flip == "Y" & VCF_ALT != "." ~ VCF_REF, # Just swapping the VCF REF for the alt
                             Flip == "Y" & VCF_ALT == "." ~ Rob_ref_allele,# If flip Yes and do not know alt == "., t"hen ALT = REF
                             Flip == "N" & VCF_ALT == "." & VCF_REF == Rob_ref_allele & VCF_REF == TF_ref ~ TF_alt, # If Flip N and ALT
                             Flip == "N" & VCF_ALT == "." & VCF_REF == Rob_ref_allele & VCF_REF != TF_ref ~ TF_ref,
                             Flip == "N" & VCF_ALT != "." & VCF_REF == Rob_ref_allele & VCF_REF == TF_ref ~ TF_alt,
                             Flip == "N" & VCF_ALT != "." & VCF_REF == Rob_ref_allele & VCF_REF != TF_ref ~ TF_ref,
                             Flip == "N" & VCF_ALT == TF_alt & VCF_REF == comp(Rob_ref_allele, forceToLower = F) & VCF_REF == TF_ref ~ comp(TF_alt, forceToLower =  F),
                             Flip == "N" & VCF_ALT == TF_ref & VCF_REF == comp(Rob_ref_allele, forceToLower = F) ~ comp(TF_ref, forceToLower = F),
                             Flip == "N" & VCF_ALT == "." & VCF_REF == TF_ref & Rob_ref_allele == comp(VCF_REF, forceToLower = F) ~ comp(TF_alt, forceToLower = F),
                             Flip == "N" & VCF_ALT == "." & VCF_REF == comp(Rob_ref_allele, forceToLower = F) ~ comp(TF_ref, forceToLower = F))) %>%
                             
  # 2. Change the reference allele
  mutate(VCF_REF = case_when(Flip == "Y" & VCF_ALT == VCF_REF ~ Rob_ref_allele,
                             Flip == "N" ~ Rob_ref_allele,
                             Flip == "N" & VCF_REF == comp(Rob_ref_allele, forceToLower = F) & VCF_REF == TF_ref & VCF_ALT == TF_alt ~ TF_ref,
                             Flip == "Y" & VCF_ALT == Rob_ref_allele & VCF_ALT == TF_ref ~ TF_ref,
                             Flip == "Y" & VCF_ALT == Rob_ref_allele & Rob_ref_allele != TF_ref ~ Rob_ref_allele)) %>%
  # 4. change the Alt and Ref allele again
  mutate(VCF_ALT = case_when(Flip == "N" ~ VCF_ALT,
                             Flip == "Y" & VCF_ALT == Rob_ref_allele & VCF_ALT == TF_ref ~ TF_alt,
                             Flip == "Y" & VCF_REF == Rob_ref_allele & VCF_REF != TF_ref ~ TF_ref,
                             Flip == "Y" ~ VCF_ALT),
         VCF_REF = case_when(Flip == "Y" ~ VCF_REF,
                             Flip == "N" ~ VCF_REF)) %>%
  # 5, Final VCF alt again to account for instances where alleles are the same
  mutate(VCF_ALT = case_when(
                             Flip == "Y" ~ VCF_ALT,
                             Flip == "N" ~ VCF_ALT,
                             VCF_ALT == VCF_REF & VCF_REF == TF_ref ~ TF_alt,
                             Flip == "N" & VCF_REF == Rob_ref_allele & Rob_ref_allele == comp(TF_ref, forceToLower = F) ~ comp(TF_alt, forceToLower = F))) %>%
  
  # 6. Actual flipping of genotypes
  mutate(across(21:147, ~ case_when(
                      Flip == "Y" & . == "./." ~ "./.", # Account for missing data
                      Flip == "Y" & . == "0/0" ~ "1/1",
                      Flip == "Y" & . == "0/1" ~ "1/0",
                      Flip == "Y" & . == "1/1" ~ "0/0",
                      Flip == "N" & . == "0/0" ~ "0/0",
                      Flip == "N" & . == "0/1" ~ "0/1",
                      Flip == "N" & . == "1/1" ~ "1/1",
                      Flip == "N" & . == "./." ~ "./."))) # mutate flip ID again


# Last step of script is to ensure the proper alt allele

# Actual ref and alternative allele pair
ref_alt_pair <- read.csv("Rob_Schnabel/9913_ARS1.2_648875_BOS1_marker_name_180910.REF", sep = "\t", header = F) %>% select(1,4,5)

# change column names
colnames(ref_alt_pair) <- c("ID", "ref", "alt")

# filter for the SNPs we have
ref_alt_pair <- ref_alt_pair %>% filter(ID  %in% master_flip_IDs$TF_affy_ID)


# Select the Reference allele
# Merge to our master file with flipped genotypes
real_ref <- master_flip_IDs %>% select(TF_affy_ID, Rob_ref_allele)
ref_alt_pair <- left_join(ref_alt_pair, real_ref, by = c("ID" = "TF_affy_ID"))

# Determine the correct alternative allele for our VCF file
ref_alt_pair <- ref_alt_pair %>% mutate(Rob_alt_allele = case_when(
  Rob_ref_allele == ref ~ alt,
  Rob_ref_allele == alt ~ ref
))

# Select the alternative allele
real_alt <- ref_alt_pair %>% select(1,5)

# Merge the alternative allele
# Edit the alternative allele in the VCF file
# Note, the incorrect VCF_alt is a consequence of not knowing if it is a T or a C.
master_flip_IDs <- right_join(master_flip_IDs, real_alt, by = c("TF_affy_ID" = "ID")) 
master_flip_IDs <- master_flip_IDs %>% mutate(VCF_ALT = case_when(VCF_ALT == Rob_alt_allele ~ VCF_ALT,
                             VCF_ALT != Rob_alt_allele ~ Rob_alt_allele))

# select appropriate columns
master_flip_IDs_ID_and_FLIP_only <- master_flip_IDs %>% select(TF_affy_ID, Flip)
head(master_flip_IDs)
master_flip_IDs <- master_flip_IDs %>% select(2,13:147)
colnames(master_flip_IDs)
master_flip_IDs <- master_flip_IDs %>% select(2,3,1,4:136)

colnames(master_flip_IDs)[1] <- "Chrom"
colnames(master_flip_IDs)[3] <- "ID"
colnames(master_flip_IDs)[4] <- "REF"
colnames(master_flip_IDs)[5] <- "ALT"
head(master_flip_IDs)
master_flip_IDs <- master_flip_IDs[order(master_flip_IDs$Chrom, master_flip_IDs$POS),]
colnames(master_flip_IDs)[1] <- "#CHROM"


View(master_flip)
# --------------------------------------------------------------------------------
# Write to a file
# header
head(master_flip_IDs)
write.table(master_flip_IDs_ID_and_FLIP_only, file = "ID_Flip_decision.txt", sep = "\t")
write.table(vcf_header, file = "VCF_header", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(master_flip_IDs, file = "Axiom_Remapped_Flipped_FINAL_A1.vcf",  sep = "\t", col.names = T, quote = F, row.names = F) 
snp_subset <- master_flip_IDs %>% select (1,2)
write.table(snp_subset, file = "regions_file_subset.txt", sep = "\t", col.names = T, quote =F, row.names = F)





# Looking at allele frequenceis

master_flip_IDs[master_flip_IDs$POS==83066775,]


final_check <- left_join(master_flip_IDs, master_flip, by = c("ID" = "TF_affy_ID")) %>% filter(Flip == "Y") %>% select(1,2,3,4,5,10,11,12,13,148,149,150,151,156,157,158,159,283)
colnames(final_check)
head(final_check)



# Generating file with no flipping of alleles
# Perform the flipping
# Many if_else statemnts, reason for so many is because need to account for "." allele which is an issue with PLINK

master_NOT_flip_IDs <- master_flip %>%  
  
  # 1. change the VCF_alt allele
  mutate(VCF_ALT = case_when(Flip == "Y" & VCF_ALT != "." ~ VCF_REF, # Just swapping the VCF REF for the alt
                             Flip == "Y" & VCF_ALT == "." ~ Rob_ref_allele,# If flip Yes and do not know alt == "., t"hen ALT = REF
                             Flip == "N" & VCF_ALT == "." & VCF_REF == Rob_ref_allele & VCF_REF == TF_ref ~ TF_alt, # If Flip N and ALT
                             Flip == "N" & VCF_ALT == "." & VCF_REF == Rob_ref_allele & VCF_REF != TF_ref ~ TF_ref,
                             Flip == "N" & VCF_ALT != "." & VCF_REF == Rob_ref_allele & VCF_REF == TF_ref ~ TF_alt,
                             Flip == "N" & VCF_ALT != "." & VCF_REF == Rob_ref_allele & VCF_REF != TF_ref ~ TF_ref,
                             Flip == "N" & VCF_ALT == TF_alt & VCF_REF == comp(Rob_ref_allele, forceToLower = F) & VCF_REF == TF_ref ~ comp(TF_alt, forceToLower =  F),
                             Flip == "N" & VCF_ALT == TF_ref & VCF_REF == comp(Rob_ref_allele, forceToLower = F) ~ comp(TF_ref, forceToLower = F),
                             Flip == "N" & VCF_ALT == "." & VCF_REF == TF_ref & Rob_ref_allele == comp(VCF_REF, forceToLower = F) ~ comp(TF_alt, forceToLower = F),
                             Flip == "N" & VCF_ALT == "." & VCF_REF == comp(Rob_ref_allele, forceToLower = F) ~ comp(TF_ref, forceToLower = F))) %>%
  
  # 2. Change the reference allele
  mutate(VCF_REF = case_when(Flip == "Y" & VCF_ALT == VCF_REF ~ Rob_ref_allele,
                             Flip == "N" ~ Rob_ref_allele,
                             Flip == "N" & VCF_REF == comp(Rob_ref_allele, forceToLower = F) & VCF_REF == TF_ref & VCF_ALT == TF_alt ~ TF_ref,
                             Flip == "Y" & VCF_ALT == Rob_ref_allele & VCF_ALT == TF_ref ~ TF_ref,
                             Flip == "Y" & VCF_ALT == Rob_ref_allele & Rob_ref_allele != TF_ref ~ Rob_ref_allele)) %>%
  # 4. change the Alt and Ref allele again
  mutate(VCF_ALT = case_when(Flip == "N" ~ VCF_ALT,
                             Flip == "Y" & VCF_ALT == Rob_ref_allele & VCF_ALT == TF_ref ~ TF_alt,
                             Flip == "Y" & VCF_REF == Rob_ref_allele & VCF_REF != TF_ref ~ TF_ref,
                             Flip == "Y" ~ VCF_ALT),
         VCF_REF = case_when(Flip == "Y" ~ VCF_REF,
                             Flip == "N" ~ VCF_REF)) %>%
  # 5, Final VCF alt again to account for instances where alleles are the same
  mutate(VCF_ALT = case_when(
    Flip == "Y" ~ VCF_ALT,
    Flip == "N" ~ VCF_ALT,
    VCF_ALT == VCF_REF & VCF_REF == TF_ref ~ TF_alt,
    Flip == "N" & VCF_REF == Rob_ref_allele & Rob_ref_allele == comp(TF_ref, forceToLower = F) ~ comp(TF_alt, forceToLower = F)))
  

# Last step of script is to ensure the proper alt allele

# Actual ref and alternative allele pair
ref_alt_pair <- read.csv("Rob_Schnabel/9913_ARS1.2_648875_BOS1_marker_name_180910.REF", sep = "\t", header = F) %>% select(1,4,5)

# change column names
colnames(ref_alt_pair) <- c("ID", "ref", "alt")

# filter for the SNPs we have
ref_alt_pair <- ref_alt_pair %>% filter(ID  %in% master_NOT_flip_IDs$TF_affy_ID)


# Select the Reference allele
# Merge to our master file with flipped genotypes
real_ref <- master_NOT_flip_IDs %>% select(TF_affy_ID, Rob_ref_allele)
ref_alt_pair <- left_join(ref_alt_pair, real_ref, by = c("ID" = "TF_affy_ID"))

# Determine the correct alternative allele for our VCF file
ref_alt_pair <- ref_alt_pair %>% mutate(Rob_alt_allele = case_when(
  Rob_ref_allele == ref ~ alt,
  Rob_ref_allele == alt ~ ref
))

# Select the alternative allele
real_alt <- ref_alt_pair %>% select(1,5)

# Merge the alternative allele
# Edit the alternative allele in the VCF file
# Note, the incorrect VCF_alt is a consequence of not knowing if it is a T or a C.
master_NOT_flip_IDs <- right_join(master_NOT_flip_IDs, real_alt, by = c("TF_affy_ID" = "ID")) 
master_NOT_flip_IDs <- master_NOT_flip_IDs %>% mutate(VCF_ALT = case_when(VCF_ALT == Rob_alt_allele ~ VCF_ALT,
                                                                  VCF_ALT != Rob_alt_allele ~ Rob_alt_allele))

# select appropriate columns

master_NOT_flip_IDs <- master_NOT_flip_IDs %>% select(2,13:149)
colnames(master_NOT_flip_IDs)
master_NOT_flip_IDs <- master_NOT_flip_IDs %>% select(2,3,1,4:138)

colnames(master_NOT_flip_IDs)[1] <- "Chrom"
colnames(master_NOT_flip_IDs)[3] <- "ID"
colnames(master_NOT_flip_IDs)[4] <- "REF"
colnames(master_NOT_flip_IDs)[5] <- "ALT"
head(master_NOT_flip_IDs)
master_NOT_flip_IDs <- master_NOT_flip_IDs[order(master_NOT_flip_IDs$Chrom, master_NOT_flip_IDs$POS),]
colnames(master_NOT_flip_IDs)[1] <- "#CHROM"


head(master_NOT_flip_IDs)
# --------------------------------------------------------------------------------
# Write to a file
# header

write.table(vcf_header, file = "VCF_header", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(master_NOT_flip_IDs, file = "Axiom_Remapped_Unflipped_FINAL.vcf",  sep = "\t", col.names = T, quote = F, row.names = F) 
snp_subset <- master_NOT_flip_IDs %>% select (1,2)
write.table(snp_subset, file = "regions_file_subset.txt", sep = "\t", col.names = T, quote =F, row.names = F)