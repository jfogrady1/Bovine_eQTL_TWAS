#!/usr/bin/env Rscript

##############################################################################################################################################################
##############################################################################################################################################################
###################################################### R Script for differential Expression Analysis #########################################################
##############################################################################################################################################################
##############################################################################################################################################################


args = commandArgs(trailingOnly=TRUE)

##############################################################################################################################################################
##################################################### 1. Package installation and loading libraries ##########################################################
##############################################################################################################################################################

suppressPackageStartupMessages(library(tidyverse, quietly = T))
suppressPackageStartupMessages(library(dplyr, quietly = T))
suppressPackageStartupMessages(library(DESeq2, quietly = T))
suppressPackageStartupMessages(library(ggplot2, quietly = T))
suppressPackageStartupMessages(library(vsn, quietly = T))
suppressPackageStartupMessages(library(pheatmap, quietly = T))
suppressPackageStartupMessages(library(apeglm, quietly = T))
suppressPackageStartupMessages(library(hexbin, quietly = T))
suppressPackageStartupMessages(library(tidyquant, quietly = T))
suppressPackageStartupMessages(library(sva, quietly = T))
library(ggrepel)
# Set the proper directory where the files are located

##############################################################################################################################################################
##################################################### 2. Sorting input data into the correct format ##########################################################
##############################################################################################################################################################

# DESeq2 requires 

  # 1. Count matrix, 
  # 2. Metadata file  for input into the DESeq2 object (Generated below)

# Specify the count data as a matrix
# Specify the metadata file as a dataframe - annoying but that is how it is
counts <- as.matrix(read.table(args[1], skip = 0, header = T, sep="\t",row.names=1)) # matrix file
coldata <- as.data.frame(read.table(args[2], sep = '', skip = 0, header = T, row.names=1)) # metadata file

cols = colnames(counts)
cols = sort(cols)

# order the samples 
counts <- counts[,cols]
coldate <- coldata[cols,]


# Select conditions we want
# Convert to factor type for R for certain variables
# Can do it for continuous variables but in R it is just complicated so leave
### batch correction
group = coldata$Condition
batch = coldata$Batch
Geneid <- rownames(counts)
counts <- ComBat_seq(counts, batch=batch, group=group)
rownames(counts) <- Geneid



coldata %>% select(-c("Batch"))

coldata$Condition <- factor(coldata$Condition, levels = c("Control", "Infected"))


# Check to make sure that row names of metadata and column names of count data are the same 
# This is required by DeSeq2
# If not, please do some data wrangling to change them (see below what I did. Had to remove some samples when Genotype covariates were included)



#all(rownames(coldata) == colnames(counts)) # Final check

#####################################
### Generate the DESeq object########
#####################################

#When you add a covariate with a large mean and small SD, you are creating a near collinearity with the intercept. 
#This is similar to a nearly confounded design: the coefficients can trade off and give the same likelihood, 
#so you are impairing the estimation of the coefficients. Collinearity is part of a typical linear models course.

#It is simple to resolve: don't include covariates with a large mean value, just center and scale the continuous-valued covariates before adding to the design. Yes it is recommended.

#This is why DESeq2 prints the message to center and scale covariates when the mean is large (Mike Love).

admix <- read.table(args[4], header = F)
colnames(admix)[1] <- "Holstein" 
admix <- admix %>% select(1)
coldata <- cbind(coldata, admix)
coldata$Age <- scale(coldata$Age, center = T)
coldata$Holstein <- scale(coldata$Holstein, center = T)



# Making object: This class is used for downstream analysis in DESEQ2
dds <- DESeqDataSetFromMatrix(countData = counts, # count matrix
                              colData = coldata, # metadata 
                              design = ~ Age + Holstein + Condition) # Columns from metadata: Here we are measuring the effect of condition controlling for


dds$Condition <- factor(dds$Condition, levels = c("Control","Infected"))
# Filter out lowly expressed genes
# This is not required but will speed up computational flow. 
# Can put in genes which pass eQTL filtering here

eQTL_genes <- read.table(args[5]) %>% select(1)
eQTL_genes$Geneid <- sub("^", "GeneID:", eQTL_genes$Geneid) %>% as.character()
genesTokeep <- which(rownames(dds) %in% eQTL_genes$Geneid)
dds <- dds[genesTokeep, ]





###############################################################################################################################################################
################################################## 3. QC and Data Examination of Raw Counts ###################################################################
###############################################################################################################################################################

# Compare counts across/between samples, see how similar they are to each other and identify sources of variation

# This is done on the *****RAW COUNTS*******
# I.e., do not normalize before generating the DESeq2 object

# First log transform the counts to improve visualization of the clustering
# Method: Variance stabilizing tranformation (VST) - similar to putting counts on log2 scale. The purpose of this is to;
  # - Remove the dependence of the variance on the mean - especially high variance of log data when mean is low.
    # - The genes with the same mean do not have the same standard deviations but experimental wide trend of variance is flattened. Genes with high row variance are interesting ones
    # - Idea is to flatten the experimental wide trend of variance
    # - Can use rlog (similar) but it is slower and not recommended
# The Blind argument as to whether vst/rlog should be blind to the sample information specified in the design formula of DESeq2 object is set to FALSE and is recommended by Michael Love who created DESeq2
# Calculate the across sample variability (blind = FALSE)


# Performing variance stabilization
# dds_vst is a class DESeqTransform
dds_vst <- vst(dds, blind = FALSE)

# Assay is used to extract the matrix ofnormalized values
#head(assay(dds_vst), 3)

# Plot of the variance on the mean. Should be relatively small
#meanSdPlot(assay(dds)) # Before transformation, lots of variability
#meanSdPlot(assay(dds_vst)) # small after transformation. difference is tiny 


## Extract data to view in ggplot as PCA
# Note this data is not controled for counfounders. I don't know how to visualise transcriptomic data while
# Taking account of covaraites
pcaData <- plotPCA(dds_vst, intgroup=c("Condition"), ntop = 1500, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  scale_colour_discrete(labels=c('Control', 'Reactor')) +
  geom_text(aes(label=rownames(pcaData)), nudge_x = 0.0005, nudge_y = 0.5)

ggsave("PCA.tiff", width=10, height=8, dpi=600)

####################################################################################################################################################################
#################################################### 4.  DESeq Analysis ############################################################################################
####################################################################################################################################################################


# command to perform the steps detailed above
dds <- DESeq(dds, quiet = T)

# plotting the dispersion and shrinkage
plotDispEsts(dds)

dds$Condition
# Extract results and specify required p-value
res <- results(dds, alpha = 0.05, contrast=c("Condition","Infected","Control"))

#summary(res)

# MA plot
# in DESeq2, the function plotMA shows the log2 fold changes attributable to a given variable over the 
# mean of normalized counts for all the samples in the DESeqDataSet.

plotMA(res, ylim=c(-8,8)) # can see noise from log fold change in genes with low counts (expected)

###############################################################################################################################################################################
######################################################### LFC shrinkage #######################################################################################################
###############################################################################################################################################################################
# DO NOT USE LFCSHRINK METHOD (HERE)- USE APGLEM METHOD BELOW
# Many differentially expressed genes are those which have low counts - because variability is higher in those genes with low mean counts
# This is known as heteroskedasicity
# As sample size increases, a gene can be found to be significnatly different and the null (LFC > | < 0) and rejected even if size of difference is small
# Also, If counts have a high coefficient of variation, estimated LFCs will have high variance = false positives
# Overcome this by shrinking estimates towards zero when available information for a gene is low
# This is done with an empirical Bayes procedure as done above

# First do GLM to obtain MLEs for the LFCs and fit a normal distribution for observed LFCs. This normal distribution is a prior
# Prior information used in second round of GLM and the maximum a posteriori (MAP) are kept as final estimates of LFC and used in P-value estimates
# Standard error calculated fromt he posterior's curvature at its maximum
# Overall, the MAP LFCs are baiased towards 0

#resultsNames(dds) # list coefficient names
#res <- lfcShrink(dds, coef= "Condition_infected_vs_control", type="normal") # OBSOLETE

#summary(res, alpha = 0.05)

#plotMA(res, ylim=c(-2,2))

###################
## apeglm method###
###################

# Above method (normal applies an adaptive normally distributed prior)
# This can have be overly stringent and lead to a loss of genes  or overly aggressive shrinkage of true LFCs

# This method uses an empirical Bayes procedure
# This approach uses a heavy tailed Cauchy distribution (t distribution with 1 degree of freedom) as a prior instead of the normal distribution
# Then the approximate Posterior Estimation for the GLM (apeglm) is calculated
# Overall, variance is reduced but preserves true large effect sizes

#resultsNames(dds) # list coefficient names
res <- lfcShrink(dds, coef= "Condition_Infected_vs_Control", type="apeglm", quiet = T) # USE Condition as this is what we are testing for
resultsNames(dds)
#summary(res, alpha = 0.05) # Can see that the number of DE genes has fallen

# plot MA with the shrunken values
#plotMA(res, ylim=c(-5,5)) # can see that very few outliers 


res_df <- as.data.frame(res)

res_final <- res_df %>%
  filter(padj <= 0.05)



res_final <- res_final[order(-res_final$log2FoldChange),  ]

################################################################################################################################################################################
################################################### Hypothesis testing for differential expression #############################################################################
################################################################################################################################################################################

# DESEq2 uses a WALD test: shrunken estimate of LFC / se (from MAP) to generate a Z statistic compared to standard normal distribution
# Resulting p-values are generated from genes which pass an independent filtering step are adjusted for multiple testing

# Genes which have little or no chance of being dec=tected as differentially expressed are omitted from testing (criteria needs to be independent of test statistic)
# The average expressions terength of each gene across all samples as its filter criterion and it omits all genes with mean normalized counts < filtering threshold from multiple testing adjustment

# Default hypothesis is that there is no log fold change. Eventually, every gene will be detected if power is there (sample size)
# Shift to detect a biologically significant magnitude change
# May be wise to do a composite null hypothesis of |Bir| < Theta where Bir is the shrunken lfc estimation and theta is a value 


###############################################
## Detecting outliers with COOKs distance #####
###############################################

# Cook's distance defined within each gene for each sample as the scaled distance that the coefficient vector Bi of a linear model would move if the sample were removed and model refitted
# If a gene has a Cook's distance > 0.99 quantile of the F(p, m-p) distribution where p is the number of model parameters and m is the number of samples
# Motivation that removing a single sample should not remove Bi outside of a 99% confidence region around Bi fit using all the samples
# These outlines should be examined and action depends on sample size
# DESeq2 replaces outlier counts with an imputed value (trimmed mean over all samples scaled by size factor) which is more conservative then leaving value out
# This is for samples with > 7 replicates


##########################################
###### Extracting necessary files ########
##########################################

normalized.counts <- as.data.frame(counts(dds, normalized=TRUE))
write.table(normalized.counts, file="DESEQ2_norm_counts.conVinfec.txt", sep = "\t", row.names = T, col.names = T, quote = F)


#####################################################################################################################
###### Plotting #####################################################################################################
#####################################################################################################################


##### Gene Names ############
annotation <- read.csv(args[3], sep = ",") # annotation file
res_df$Gene_ID <- rownames(res_df) # Getting row names their own vector
res_df$Gene_ID <- sub("GeneID:", "", res_df$Gene_ID) %>% as.numeric()
res_df <- left_join(res_df, annotation, by = c("Gene_ID" = "Gene_ID")) %>% select(1:8)
res_df <- res_df %>% select(8,7,6,1,2,3,4,5)


###### Colours ##############
# add a column of NAs
res_df$diffexpressed <- "Not DE"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res_df$diffexpressed[res_df$log2FoldChange > 0 & res_df$padj < 0.05] <- "DE Up"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_df$diffexpressed[res_df$log2FoldChange < 0 & res_df$padj < 0.05] <- "DE Down"

res_df <- res_df %>% mutate(PLOT_Symbol =
                              case_when(
                                log2FoldChange < -.1 & padj < 0.00005 ~ Symbol,
                                log2FoldChange > 1 & padj < 0.00005 ~ Symbol,
                                log2FoldChange > 0 & padj < 0.0000005 ~ Symbol,
                                log2FoldChange < -0.5 & padj < 0.01 ~ Symbol,
                                log2FoldChange > 1.5 & padj < 0.01 ~ Symbol,
                                FALSE ~ ""))


# Plotting
ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=PLOT_Symbol)) +
  geom_point(size = 1) + 
  scale_color_manual("Comparison", values=c("#2166ac", "#b2182b", "grey")) +
  labs(x=bquote(~log[2]~'fold change'),
       y=bquote(~-log[10]~italic(P)[adj]))+
  scale_x_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_y_continuous(limits = c(0,7)) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
  theme_tq() +
  geom_text_repel(colour = "black", max.overlaps = 40, size = 3) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 21, color = "black", face = "bold"),
        axis.title.x = element_text(size = 21, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "black", face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 5)))


ggsave("Diff_Expression.png", width=10, height=8, dpi=600)


library(gprofiler2)
res_df
genes <- res_df %>% filter(padj < 0.01) %>% dplyr::select(2) %>% as.list()
genes

results <- gost(query = genes,organism = "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = res_df$Symbol, domain_scope = "known", user_threshold = 0.05, sources = c("GO:BP", "KEGG", "REAC"), evcodes = T)
results$result$query <- ""
result_df <- results$result
result_df$term_id
terms <- c("interferon-mediated signaling pathway", "innate immune response",
           "defense response to virus",
           "regulation of cytokine-mediated signaling pathway",
           "type II interferon-mediated signaling pathway",
           "regulation of chromatin binding",
           "NOD-like receptor signaling pathway",
           "RIG-I-like receptor signaling pathway",
           "ISG15 antiviral mechanism",
           "Viral life cycle - HIV-1",
           "receptor signaling pathway via STAT",
           "Antiviral mechanism by IFN-stimulated genes")
results$result <- results$result %>%
  mutate(label = ifelse(term_name %in% terms, term_name, NA))
colnames(result_lables)[11] <- "term_name"
plot <- gostplot(results, capped=F, interactive=F, pal=c(`GO:BP`="#ff9900",KEGG= "#dd4477",REAC="#3366cc",WP="#0099c6")) + labs(y=bquote(~-log[10]~italic(P)[adj]), title = "", x = "Database") +
  theme_bigstatsr() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 21, color = "black", face = "bold"),
        axis.title.x = element_text(size = 21, color = "black", face = "bold")) +
  geom_text_repel(aes(label = results$result$label), colour = "black", max.overlaps = 40,
                  nudge_x = 0, nudge_y = 1.2, size = 4)


publish_gostplot(plot, highlight_terms = NULL, filename="Enrichment_FDR01.png", width = 15, height = 8)
write.table(res_df, file="DE_Genes_Con_V_Infected.txt", sep = "\t", quote = F, col.names = T, row.names = F)
Enrich_res <- as.data.frame(do.call(cbind, results$result))
Enrich_res <- apply(Enrich_res,2,as.character)
write.table(Enrich_res, file = "Enrichment_FDR1.txt", sep = "\t", quote = F, col.names = T, row.names = F)









