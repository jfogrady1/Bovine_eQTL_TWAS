library(tidyverse)
library(Biostrings)
library(seqinr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)


# Gather the master file
master_final = read.table(args[1], sep = "\t", header = T)


# VCF modification, extremely important


# Extract the header to write to the file after
vcf_header <- read.csv(args[2], sep = "\t", header = F)
vcf_header <- vcf_header[1:33,1]


# Specify the vcf data which will be edited
vcf <- as.data.frame(read.csv(args[2], sep = "\t", skip = 34)) # make sure to say colnames = T when writing

colnames(vcf)[1] <- "#CHROM"
colnames(vcf)[4] <- "VCF_REF"
colnames(vcf)[5] <- "VCF_ALT"


# join vcf to master ID file
master_final <- master_final %>% left_join(., vcf, by = c("AX_id" = "ID")) %>% drop_na() #drop_na here as SNPS on X need ot be gone




# Determining which SNPs need to be flipped and respective genotypes
master_flip <- master_final %>% mutate(Flip = case_when(
  Rob_ref_allele == VCF_REF | Rob_ref_allele == comp(seq = VCF_REF, forceToLower = F) ~ "N", # do nothing as essentially the same genotype on both assemblied
  Rob_ref_allele == VCF_ALT | Rob_ref_allele == comp(seq = VCF_ALT, forceToLower = F) ~ "Y", # these genotypes will need to be flipped 0 - 1 and 1 - 0 
  VCF_ALT == "." & Rob_ref_allele == VCF_REF | Rob_ref_allele  == comp(seq = VCF_REF, forceToLower = F) ~ "N", # Need to include this as PLINK messes everything up
  VCF_ALT == "." & Rob_ref_allele != VCF_REF | Rob_ref_allele != comp(seq = VCF_REF, forceToLower = F) ~ "Y")) # Don't know if these are on the same strand - missing


sum(table(master_flip$Flip)) == dim(master_flip)[1] # check that all conditions satisfied
table(master_flip$Flip) # see how many we need to FLIP
write.table(master_flip, file = "Supplementary_Table_Flipping.txt", col.names = T, row.names = F, sep = "\t")

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
  mutate(across(21:142, ~ case_when(
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
args[3]
ref_alt_pair <- read.csv(args[3], sep = "\t", header = F) %>% select(1,4,5)

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

master_flip_IDs <- master_flip_IDs %>% select(2,13:145)

master_flip_IDs <- master_flip_IDs %>% select(2,3,1,4:132)

colnames(master_flip_IDs)[1] <- "Chrom"
colnames(master_flip_IDs)[3] <- "ID"
colnames(master_flip_IDs)[4] <- "REF"
colnames(master_flip_IDs)[5] <- "ALT"
head(master_flip_IDs)
master_flip_IDs <- master_flip_IDs[order(master_flip_IDs$Chrom, master_flip_IDs$POS),]
colnames(master_flip_IDs)[1] <- "#CHROM"



# --------------------------------------------------------------------------------
# Write to a file
# header

write.table(master_flip_IDs_ID_and_FLIP_only, file = args[4], sep = "\t")
write.table(vcf_header, file = args[5], sep = "\t", quote = F, row.names = F, col.names = F)
write.table(master_flip_IDs, file = args[6],  sep = "\t", col.names = T, quote = F, row.names = F) 
snp_subset <- master_flip_IDs %>% select (1,2)
write.table(snp_subset, file = args[7], sep = "\t", col.names = T, quote = F, row.names = F)




  