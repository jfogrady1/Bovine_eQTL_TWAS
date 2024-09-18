library(tidyverse)
library(data.table)
# Data wrangling for GWAS data to ensure that only consider eQTLs 
# and that REF/alternative alleles are the same
args = commandArgs(trailingOnly=TRUE)
CH <- args[1] # ARS remapped
CH_beta <- args[2] # UMD GWAS
Run6 <- args[3] # RUN 6 1000 BULLs
ARS <- args[4] # ARS BQSR REF/ALT alleles
gwas <- args[5]
# read in remapped
# note ID column is the old UMD position
# Toms GWAS

CH <- fread(CH)

# Note, REF and ALT alleles refer to ARS SNPs, not UMD, hence need to determine if effect size is flipped


# Read in siobhans data to get the positions of the variants, will use these to merge with 1000 BUll
CH_beta <- fread(CH_beta)
colnames(CH_beta) <- c("UMD_coord", "CHR_UMD", "POS_UMD", "P_value", "Beta")
CH_beta <- CH_beta %>% dplyr::select(-c("P_value"))

#1000 bull run 6
# Some samples have rsids tied to beta values but we need to get chr and pos for these
Bull <- fread(Run6)
colnames(Bull)[1] <- "Chr"
Bull$Chr <- gsub("^", "Chr", Bull$Chr)
Bull$UMD_final <- paste0(Bull$Chr, ":", Bull$POS)
Bull$Ref_SNPs <- if_else(Bull$Ref_SNPs == ".", paste0(Bull$Chr, ":", Bull$POS), Bull$Ref_SNPs)

# Identify chromosome position for SNPs on UMD3.1 and merge these with ref and alt alleles
CH <- left_join(CH, Bull, by = c("ID" = "Ref_SNPs"))
CH <- left_join(CH, CH_beta, by = c("UMD_final" = "UMD_coord"))
rm(Bull)
rm(CH_beta)
colnames(CH) <- c("ARS_CHR", "ARS_POS", "GWAS_ID", "P", "FDR", "Ref", "Alt", "UMD_CHR", "UMD_POS",
                  "UMD_REF", "UMD_ALT", "Strand", "UMD_final", "CHR_UMD", "POS_UMD", "Beta")


CH <- CH %>% select(1,2,3,4,5,10,11,12,13,14,15,16)



# Now have beta values for all ARS SNPs with matching UMD IDs
# Also have REF/ALT alleles for UMD3.1 based on 1000 BULLS which was used to impute Siobhan's data

# Now need to get ARS_1.2 REF/ALT allele pairs
ARS_reference <- fread(ARS)




# Get ID for identification
CH$ARS_ID <- paste0(CH$ARS_CHR, ":", CH$ARS_POS)
colnames(ARS_reference)[1] <- "Chr"

ARS_reference$ARS_ID <- paste0(ARS_reference$Chr, ":", ARS_reference$POS)
ARS_reference <- ARS_reference %>% select(9,4,5)



CH_final <- left_join(CH, ARS_reference, by = c("ARS_ID" = "ARS_ID"))
rm(ARS_reference)
CH_final <- CH_final %>% mutate(flip = case_when(
                           UMD_REF == "A" & UMD_ALT == "T" & REF == "A" & ALT == "T" ~ "No",
                           UMD_REF == "T" & UMD_ALT == "A" & REF == "T" & ALT == "A" ~ "No",
                           UMD_REF == "C" & UMD_ALT == "G" & REF == "C" & ALT == "G" ~ "No",
                           UMD_REF == "G" & UMD_ALT == "C" & REF == "G" & ALT == "C" ~ "No", # Perfect agreement and palindromic
                           UMD_REF == "A" & UMD_ALT == "G" & REF == "A" & ALT == "G" ~  "No",
                           UMD_REF == "A" & UMD_ALT == "C" & REF == "A" & ALT == "C" ~  "No",
                           UMD_REF == "G" & UMD_ALT == "T" & REF == "G" & ALT == "T" ~  "No",
                           UMD_REF == "G" & UMD_ALT == "A" & REF == "G" & ALT == "A" ~  "No",
                           UMD_REF == "C" & UMD_ALT == "A" & REF == "C" & ALT == "A" ~  "No",
                           UMD_REF == "C" & UMD_ALT == "T" & REF == "C" & ALT == "T" ~  "No",
                           UMD_REF == "T" & UMD_ALT == "C" & REF == "T" & ALT == "C" ~  "No",
                           UMD_REF == "T" & UMD_ALT == "G" & REF == "T" & ALT == "G" ~  "No", # pefect variants
                           UMD_REF == "A" & UMD_ALT == "C" & REF == "T" & ALT == "G" ~ "Yes",
                           UMD_REF == "C" & UMD_ALT == "A" & REF == "G" & ALT == "T" ~ "Yes",
                           UMD_REF == "T" & UMD_ALT == "G" & REF == "A" & ALT == "C" ~ "Yes",
                           UMD_REF == "G" & UMD_ALT == "T" & REF == "C" & ALT == "A" ~  "Yes", # strand flip 
                           UMD_REF == "A" & UMD_ALT == "T" & REF == "T" & ALT == "A" ~  "palindrome_flip",
                           UMD_REF == "C" & UMD_ALT == "G" & REF == "G" & ALT == "C" ~  "palindrome_flip",
                           UMD_REF == "T" & UMD_ALT == "A" & REF == "A" & ALT == "T" ~  "palindrome_flip",
                           UMD_REF == "G" & UMD_ALT == "C" & REF == "C" & ALT == "G" ~  "palindrome_flip"))

CH_final <- CH_final %>% filter(flip != "palindrome_flip") # Remove palindromic SNPs
CH_final$Beta_corrected <- if_else(CH_final$flip == "Yes", (CH_final$Beta * -1), CH_final$Beta)

# Select necessary columns

CH_final <- CH_final %>% select("GWAS_ID", "ARS_ID", "ARS_CHR", "ARS_POS", "REF", "ALT", "P", "CHR_UMD", "POS_UMD", "Beta_corrected")


# Now onto estimating the Z scores - need to do this for each GWAS
if (gwas == "CH") {
    # https://stats.stackexchange.com/questions/337070/compute-standard-error-from-beta-p-value-sample-size-and-the-number-of-regres
    N = 2039 
    cov = 4 
    df = (N - 1 ) - (cov - 1)
    CH_final$t_val <- qt(CH_final$P/2, df = df)
    CH_final$B_se <- abs(CH_final$Beta_corrected/abs(CH_final$t_val)) # Calculating standard error
    CH_final$Z <- CH_final$Beta_corrected / abs(CH_final$B_se)
} else if (gwas == "LM") {
    N = 1964 
    cov = 4 
    df = (N - 1 ) - (cov - 1)
    CH_final$t_val <- qt(CH_final$P/2, df = df)
    CH_final$B_se <- abs(CH_final$Beta_corrected/abs(CH_final$t_val)) # Calculating standard error
    CH_final$Z <- CH_final$Beta_corrected / abs(CH_final$B_se)

} else if (gwas == "HOFR") {
    N = 1502 
    cov = 4 
    df = (N - 1 ) - (cov - 1)
    CH_final$t_val <- qt(CH_final$P/2, df = df)
    CH_final$B_se <- abs(CH_final$Beta_corrected/abs(CH_final$t_val)) # Calculating standard error
    CH_final$Z <- CH_final$Beta_corrected / abs(CH_final$B_se)
} else if (gwas == "ALL") {
    N = 7346 
    cov = 5  # extra breed component
    df = (N - 1) - (cov - 1)
    df
    CH_final$t_val <- qt(CH_final$P/2, df = df)
    CH_final$B_se <- abs(CH_final$Beta_corrected/abs(CH_final$t_val)) # Calculating standard error
    CH_final$Z <- CH_final$Beta_corrected / abs(CH_final$B_se) # wald statistic
}

CH_final$ARS_SNP <- paste0(CH_final$ARS_ID, ":", CH_final$REF, ":", CH_final$ALT)
CH_final <- CH_final %>% select(ARS_SNP, ARS_CHR, ARS_POS, REF, ALT, P, Beta_corrected, B_se, Z)

write.table(CH_final, args[6], sep = "\t", row.names = T, col.names = T, quote = F)
