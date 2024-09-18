### by John O'Grady on September 29, 2020
options(stringsAsFactors = FALSE)
library(edgeR)
library(preprocessCore)
library(RNOmni)
library(BiocManager)
library(data.table)
library(R.utils)
library(SNPRelate)
library(tidyverse)

#----------------------------------------------------------------------------
### functions
# Transform rows to a standard normal distribution
inverse_normal_transform = function(x) {
    qnorm(rank(x) / (length(x)+1))
}
#----------------------------------------------------------------------------
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
### data input <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
args <- commandArgs(trailingOnly = TRUE)
file_counts = args[1]# Counts file. Row is gene, column is sample; rowname is gene id, colname is sample id
file_tpm  = args[2] # TPM file. Row is gene, column is sample; rowname is gene id, colname is sample id
annot_file = args[3]  # annotation file
vcf.fn = args[4]# Input data for genotype PCA. genotype data from imputation (VCF format)
group = args[5]

if (!file.exists(file_counts)) { stop("Can not find the file_counts") }
if (!file.exists(file_tpm)) { stop("Can not find the file_tpm") }
if (!file.exists(annot_file)) { stop("Can not find the annot_file") }
if (!file.exists(vcf.fn)) { stop("Can not find the vcf.fn") }


#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
### main program
#----------------------------------------------------------------------------
# Input data for TMM calculation
# Read counts matrix. Row is gene, column is sample; rowname is gene id, colname is sample id
Counts = read.table(file_counts, row.names = 2)
# Remove columns we do not need which are only in the large file
length_cols <- length(colnames(Counts))
length_cols
head(Counts)
if (length_cols > 75) {
    Counts <- Counts %>% select(-c(row.names, chr, Start_location, End_location, length))
} else {
    Counts <- Counts %>% select(-c(row.names))
}


# TPM matrix. Row is gene, column is sample; rowname is gene id, colname is sample id
TPM = read.table(file_tpm, row.names = 1)


## Note, different dimensions of counts as for some we could not calculate gene length, these will be excluded but will be retained for normalisation by TMM

## 1. prepare TMM
samids = colnames(Counts) # sample id
expr_counts = Counts
expr = DGEList(counts=expr_counts) # counts
nsamples = length(samids) # sample number
ngenes = nrow(expr_counts) # gene number

# calculate TMM
y = calcNormFactors(expr, method="TMM")
TMM = cpm(y,normalized.lib.sizes=T)

# expression thresholds
count_threshold = 6
tpm_threshold = 0.1
sample_frac_threshold = 0.2

#keep the genes with >=0.1 tpm and >=6 read counts in >=20% samples.
expr_tpm = TPM
#expr_counts <- expr_counts[rownames(TPM),]
tpm_th = rowSums(expr_tpm >= tpm_threshold)
count_th = rowSums(expr_counts >= count_threshold)
ctrl1 = tpm_th >= (sample_frac_threshold * nsamples)
ctrl2 = count_th >= (sample_frac_threshold * nsamples)
mask = ctrl1 & ctrl2
table(mask)
TMM_pass = TMM[mask,] ##row is gene; column is sample


###expression values (TMM) were inverse normal transformed across samples.
TMM_inv = t(apply(TMM_pass, MARGIN = 1, FUN = inverse_normal_transform)) #apply to each row, each row represents one gene, observed values for all the samples. scale across samples.


#----------------------------------------------------------------------------
### 2. prepare bed file
#dir.create("bed",showWarnings=F)

region_annot = read.csv(annot_file, header = F, skip = 5, sep = "\t") #%>% select(1,6,7,8) # load gtf file
region_annot <- region_annot %>% select(V10, V1, V4, V5)
colnames(region_annot) <- c("Geneid", "chr", "start", "end")
region_annot$Geneid <- gsub(";", "", region_annot$Geneid)
geneid = region_annot$Geneid
head(region_annot)
region_annot$phenotype_id <- geneid
expr_matrix = TMM_inv[rownames(TMM_inv) %in% geneid,] # expr_matrix TMM_inv

# prepare bed file for tensorQTL
bed_annot = region_annot[region_annot$Geneid %in% rownames(expr_matrix),]

rownames(bed_annot) <- bed_annot$Geneid
rownames(bed_annot) == rownames(expr_matrix)
expr_matrix <- expr_matrix[as.character(bed_annot$Geneid),]

bed = cbind(bed_annot,expr_matrix[as.character(bed_annot$Geneid),])
table(rownames(bed) == rownames(expr_matrix))

bed <- bed[,-1]
bed = bed[bed[,1] %in% as.character(1:29),]
head(bed)
bed[,1] = as.numeric(bed[,1])
bed = bed[order(bed[,1],bed[,2]),]
colnames(bed)[1] = "#Chr"
head(bed)
# output bed file
fwrite(bed, file = args[6], sep = "\t") # temp bed file


#----------------------------------------------------------------------------
### 3. Expression PCA estimation

# Generate PCA for expresion matrix 
expr_pca <-t(bed[,-(1:4)]) # keep just expression data
# Centering or scaling does not have an effect
prcompResult<-prcomp(expr_pca, scale = T, center = T) #This should take less than a minute.

resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)

print(resultRunElbow)

# Upper K limit
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
resultRunBE<-PCAForQTL::runBE(expr_pca,B=20,alpha=0.05)
print(resultRunBE$numOfPCsChosen)

K_elbow<-resultRunElbow #13.
K_BE<-resultRunBE$numOfPCsChosen #15.

PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow","BE"),values=c(K_elbow,K_BE),
                         titleText=paste0(args[5], "- Peripheral Blood"))


ggsave(args[8], width = 12, height = 8, dpi = 600)
#---------------------------------------------
### 4. Genotype PCA
snpgdsVCF2GDS(vcf.fn, args[9], method = "biallelic.only")
genofile <-snpgdsOpen(args[9])
ccm_pca <- snpgdsPCA(genofile,num.thread=8)

pca_genotype <- ccm_pca$eigenvect[, 1:30]
colnames(pca_genotype) <- paste0("pc", 1:30)
rownames(pca_genotype) <- ccm_pca$sample.id
pca_genotype0 = data.frame(SampleID=ccm_pca$sample.id,pca_genotype)
pca_var0 = data.frame(pc=1:30,eigenval=ccm_pca$eigenval[1:30],varprop=ccm_pca$varprop[1:30])
# output
fwrite(pca_genotype0, args[10], sep = "\t", row.names = F, quote = FALSE)
fwrite(pca_var0, args[11], sep = "\t", row.names = F, quote = FALSE)
#----------------------------------------------------------------------------


### 5. merge all covariates
# Genotype
pca_vect = pca_genotype0
pca_vect = pca_vect[,2:6]

# PCA expression
pca_qtl_vect =prcompResult$x
cov_exp = pca_qtl_vect[,1:K_elbow]
n_samples = nrow(pca_vect)
rownames(pca_vect)
# known covariates other than genotype PCs
covariates = read.table(args[12], row.names = 1, header = T)
rows <- rownames(pca_vect)
covariates = covariates[rows,]
covariates
if (args[5] != "ALL") {
    covariates = covariates %>% select(-Condition)
}
knownCovariates <- cbind(pca_vect, covariates)
identical(rownames(knownCovariates),rownames(expr_pca))
colnames(knownCovariates)[1:5]<-paste0("genotypePC",1:5)
head(knownCovariates)
if (args[5] == "ALL") {
    knownCovariates$Condition = if_else(knownCovariates$Condition == "Control", 1, 2)    
}

knownCovariatesFiltered<-PCAForQTL::filterKnownCovariates(knownCovariates,cov_exp,unadjustedR2_cutoff=0.9)

PCsTop<-scale(cov_exp)

cov_output = cbind(PCsTop, knownCovariatesFiltered)
cov_output = t(cov_output)
head(cov_output)
# output
write.table(cov_output, args[13], sep="\t",quote=F, row.names = T, col.names = T)
#----------------------------------------------------------------------------
cat("done.\n")
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------