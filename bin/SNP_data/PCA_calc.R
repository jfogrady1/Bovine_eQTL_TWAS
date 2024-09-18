# PCA of genomic data from VCF format

# Install packages
#BiocManager::install("SNPRelate")
#BiocManager::install("gdsfmt")
library(SNPRelate)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# Reformat VCF files
# Create a GDS file

#arg[1]
vcf.fn <- args[1]
snpgdsVCF2GDS(vcf.fn, "VCF.gds", method="biallelic.only")
snpgdsSummary("VCF.gds")


# open the genofile
genofile <- snpgdsOpen("VCF.gds")
head(read.gdsn(index.gdsn(genofile, "snp.rs.id")))


# PCA
# Perform the PCA - set snp.id = NULL to ensure TBSNPs are analysed
pca <- snpgdsPCA(genofile, snp.id=NULL, num.thread=2)

# convert to a percentage
pc.percent <- pca$varprop*100
head(pc.percent)

# convert to a table for writing
tab_ALL <- data.frame(sample.id = pca$sample.id,
                      EV1 = pca$eigenvect[,1],    # the first eigenvector
                      EV2 = pca$eigenvect[,2],    # the second eigenvector
                      EV3 = pca$eigenvect[,3],    # the third eigenvector.
                      stringsAsFactors = FALSE)
pc.percent
ggplot(data = tab_ALL, aes(x = EV1, y = EV2)) +
  geom_point() +
  xlab(paste("PC1(",round(pc.percent[1],2),"%)")) +
  ylab(paste("PC2(",round(pc.percent[2],2),"%)")) +
  geom_text(aes(label=sample.id), nudge_x = 0.0005, nudge_y = 0.005) 

# arg 2
ggsave(args[2], dpi = 600, width = 12, height = 8)

colnames(tab_ALL)[1] <- "Sample"

#arg3
write.table(tab_ALL, file = args[3], sep = "\t", row.names = F, quote = F)       
