#############################################################################################################################
################################################# Variant Remapping #########################################################
#############################################################################################################################


# Script #1 in SNP data processing
# This script remaps varaints from UMD3.1 (bosTau6) to ARSUCD1.2 using LIFTOVER from UCSC
# Much more acccurate than Robert Schaefer's method

args = commandArgs(trailingOnly=TRUE)


# Install packages if required
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") # 

# load in libraries
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(tidyverse))


# Download chain file and unzip it and specfiy it to chain
# https://genome.ucsc.edu/goldenPath/help/chain.html
# can use wget() on command line
chain <- import.chain(args[1])

# Read in our file with IDs and old positions
# Dome data wrangling like convert to numeric and remove chr etc
# This numeric format is important for rtracklayer.
# skip the first 16 lines to get to the proper data 
annotation_file <- as.data.frame(read.csv(args[2], skip = 16))
annotation_file$Chromosome <- sub("^", "chr", annotation_file$Chromosome)
annotation_file$Physical.Position <- as.numeric(annotation_file$Physical.Position)
annotation_file$Position.End <- as.numeric(annotation_file$Position.End)

# Check for NAs and remove them (GPos object can't have NAs)
#which(is.na(annotation_file$Physical.Position))
#which(is.na(annotation_file$Position.End))

# Remove NAs
annotation_file <- annotation_file[complete.cases(annotation_file),]


# Set up the 'Granges' object with names (chr), positions, strand and names as metadata
gr <- GPos(seqnames = Rle(annotation_file$Chromosome),
           pos = IRanges(start = annotation_file$Physical.Position, end = annotation_file$Position.End),
           strand = Rle(strand(annotation_file$Strand)),
           names = annotation_file$Affy.SNP.ID)


# Perform the liftOver analysis 
UCDARS <- liftOver(gr, chain)

# Convert to a dataframe
UCDARS_df <- as.data.frame(UCDARS)

# Select the columns which are applicable to us
# Change the names
UCDARS_df <- UCDARS_df %>% select(names, seqnames, start)
colnames(UCDARS_df) <- c("Probe", "chr", "pos")

# Sub out chr from chromosmoe column (PLINK doesn't like this)
UCDARS_df$chr <- sub("chr", "", UCDARS_df$chr)

# Write to a file whcih will be used by PLINK in analysis.
write.table(UCDARS_df, file = args[3], col.names = F, quote = F, sep = "\t", row.names = F)
