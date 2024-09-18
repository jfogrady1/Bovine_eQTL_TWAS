library(MOSTWAS)
library(tidyverse)
library(dplyr)
library(vcfR)
burdenTest_local <- function(wgt,
                       snps,
                       sumStats,
                       snpAnnot = NULL,
                       onlyCis = F,
                       Z = NULL,
                       beta = NULL,
                       se = NULL,
                       chr,
                       pos,
                       ref,
                       pval,
                       R2cutoff,
                       alpha,
                       locChrom,
                       nperms = 1e3){
  
  load(wgt)
  
  pieces = strsplit(wgt,'/')
  fn = pieces[[1]][length(pieces[[1]])]
  geneInt = strsplit(fn,'[.]')[[1]][1]
  colnames(sumStats)[which(colnames(sumStats) == "chr")] = 'Chromosome'
  colnames(sumStats)[which(colnames(sumStats) == "pos")] = 'Position'
head(sumStats)  
  if (R2 <= 0.01){
    return(paste0(geneInt,
                  ' is not predicted at R2 > ',
                  R2cutoff))
  }
  
  if (!'GenPos' %in% colnames(sumStats)){
    sumStats$GenPos = paste(sumStats$Chromosome,sumStats$Position,sep = ':')
  }
dim(sumStats)  
  
  if (onlyCis){
    
    Model = subset(Model,Chromosome == locChrom)
    if (R2 <= R2cutoff){
      return(paste0(geneInt,
                    ' is not locally predicted at R2 > ',
                    R2cutoff))
    }
    
  }
  
  require(dplyr)
  if ('Mediator' %in% colnames(Model)){
    Model =
      as.data.frame(Model %>%
                      dplyr::group_by(SNP,Chromosome,Position) %>%
                      dplyr::summarize(sum(Effect)))
    colnames(Model) = c('SNP','Chromosome','Position','Effect')
  } else {
    Model =
      as.data.frame(Model %>%
                      dplyr::group_by(SNP,Chromosome,Position) %>%
                      dplyr::summarize(sum(Effect)))
    colnames(Model) = c('SNP','Chromosome','Position','Effect')}
  Model$GenPos = paste(Model$Chromosome,Model$Position,sep = ':')
head(Model)  
  sumS = subset(sumStats,GenPos %in% Model$GenPos)
  if (!is.null(Z)){ colnames(sumS)[which(colnames(sumS) == Z)] = 'Z' }
  if (is.null(beta) & is.null(se)) {
    sumS$Beta = sumS$Z
    sumS$SE = 1
    beta = 'Beta'
    se = 'SE'
  }
  dim(Model)
  dim(sumS)
  colnames(sumS)[which(colnames(sumS) == beta)] = 'Beta'
  colnames(sumS)[which(colnames(sumS) == se)] = 'SE'
  colnames(sumS)[which(colnames(sumS) == chr)] = 'Chromosome'
  colnames(sumS)[which(colnames(sumS) == pos)] = 'Position'
  colnames(sumS)[which(colnames(sumS) == ref)] = 'REF'
  colnames(sumS)[which(colnames(sumS) == pval)] = 'P'
  if (sum(!Model$GenPos %in% sumS$GenPos) > 0){
    ModelOut = subset(Model,!GenPos %in% sumS$GenPos)
  }
  else {
    print(paste0("SNPs not available to compute SNP-trait LD covariate matrix"))
    ModelOut = list(geneInt,"Cannot compute LD covariance matrix","Cannot compute LD covariance matrix","Cannot compute LD covariance matrix","Cannot compute LD covariance matrix","Cannot compute LD covariance matrix") 
    return(ModelOut)# Will not work otherwise
  }
  Model = subset(Model,GenPos %in% sumS$GenPos)
  
  sumS = sumS[match(Model$GenPos,sumS$GenPos),]
  
  if (nrow(sumS) == 0){
    return(list(geneInt,'No Overlapping SNPs.','No Overlapping SNPs.','No Overlapping SNPs.','No Overlapping SNPs.','No Overlapping SNPs.'))
  }
  
  if (is.null(snpAnnot)){
    
    snpAnnot = as.data.frame(matrix(nrow = nrow(Model),
                                    ncol = 3))
    colnames(snpAnnot) = c('SNP','REF','ALT')
    snpAnnot$SNP = Model$SNP
    snpAnnot$REF = sumS$REF
    snpAnnot$ALT = sapply(strsplit(as.character(snpAnnot$SNP),':'),
                          function(x) x[4])
  }
  
  annot = subset(snpAnnot, SNP %in% Model$SNP)
  Model = subset(Model, SNP %in% annot$SNP)
  sumS = subset(sumS, GenPos %in% Model$GenPos)
  
  annot = merge(annot,
                Model[,c('SNP','GenPos')],
                by = 'SNP')
  
  flip = function(i,
                  snpAnnot,
                  sumS){
    
    search = sumS[i,]
    ref = snpAnnot[snpAnnot$GenPos == search$GenPos,]
    if (toupper(ref$REF[1]) != toupper(search$REF[1])){
      return(-1*search$Beta[1])
    } else {return(search$Beta)}
    
  }
  
  sumS$Flip = sapply(1:nrow(sumS),
                     flip,
                     snpAnnot = annot,
                     sumS = sumS)
  
  calculateTWAS <- function(effects,
                            Z,
                            LD,
                            indices){
    effects = effects[indices]
    twasZ = as.numeric(effects %*% Z)
    twasr2pred = as.numeric(effects %*% LD %*% effects)
    if (twasr2pred > 0){
      twas = as.numeric(twasZ/sqrt(twasr2pred))
    } else {
      twas = 0
    }
    return(twas)
  }
  
  if (exists('ModelOut')){
    snpCur = subset(snps, SNP %in% c(as.character(Model$SNP),
                                     as.character(ModelOut$SNP)))}
  else {
    snpCur = subset(snps,SNP %in% as.character(Model$SNP))
  }
    
  if (wgt == "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/models/CONTROL_MeTWAS_models_10_mediators/ENSBTAG00000017276.wgt.med.RData") {
    print(paste0("SNPs not available to compute SNP-trait LD covariate matrix"))
    ModelOut = list(geneInt,"Cannot compute LD covariance matrix","Cannot compute LD covariance matrix","Cannot compute LD covariance matrix","Cannot compute LD covariance matrix","Cannot compute LD covariance matrix") 
    return(ModelOut)
  }
  snpCur = snpCur[match(c(as.character(Model$SNP),
                          as.character(ModelOut$SNP)),snpCur$SNP),]
  genos = as.matrix(snpCur[,-1])
  LD = (genos %*% t(genos)) / (ncol(genos)-1)
  miss.LD = LD[1:nrow(Model),
               (nrow(Model)+1):nrow(LD)] %*%
    solve(LD[(nrow(Model)+1):nrow(LD),
             (nrow(Model)+1):nrow(LD)] + .1 * diag(nrow(ModelOut)))
  Z = as.numeric(sumS$Flip)/as.numeric(sumS$SE)
  impz = t(miss.LD) %*% as.vector(Z)
  r2pred = diag((t(miss.LD) %*% LD[1:nrow(Model),
                                   1:nrow(Model)]) %*%
                  miss.LD)
  newZ = as.numeric(impz) / sqrt(r2pred)
  Z.tot = c(Z,newZ)
  Model.tot = rbind(Model,ModelOut)
  
  
  twasLD = as.numeric(Model.tot$Effect %*% Z.tot) /
    sqrt(as.numeric(Model.tot$Effect %*% LD %*% Model.tot$Effect))
  P = 2*pnorm(-abs(twasLD))
  
  if (P <= alpha){
    permutationLD = boot::boot(data = Model$Effect,
                               statistic = calculateTWAS,
                               R = nperms,
                               sim = 'permutation',
                               Z = Z,
                               LD = LD[1:nrow(Model),1:nrow(Model)])
    permute.p = (nperms * mean(abs(permutationLD$t) >
                                 abs(permutationLD$t0)) + 1)/(nperms+1)
  } else {
    permute.p = 1}
  
  return(list(Gene = geneInt,
              Z = twasLD,
              P = 2*pnorm(-abs(twasLD)),
              permute.P = permute.p,
              topSNP = sumS$GenPos[which.min(sumS$P)],
              topSNP.P = min(sumS$P)))
}



args = commandArgs(trailingOnly=TRUE)
args1 <- args[1]
arg <- unlist(strsplit(args1, "_"))[1]


# "CH_GWAS_FINAL.txt" args[1]
sumStats <- as.data.frame(read.table(args[1], sep = "\t", header = T)) #%>% select(-X) 
# Calculating standard error
head(sumStats)



# args [2] "../eqtl_study/eqtl_nextflow/results/SNP_data/Imputation_Performance/ALL_GT.RECODED.txt"
snps <- vcf_data <- read.vcfR(args[2])
head(snps)
snps_geno <- snps@gt
head(snps@fix)
snps_id <- snps@fix %>% as.data.frame() %>% select(ID, CHROM, POS)
head(snps_id)
colnames(snps_geno) <- gsub("_.*", "", colnames(snps_geno))
snps_geno <- snps_geno[,-1]
snps_geno <- as.data.frame(snps_geno)

snps_geno <- as.data.frame(snps_geno)
snps_geno <- apply(snps_geno, 2, function(x) {gsub("0/1", 1, x)})
snps_geno <- apply(snps_geno, 2, function(x) {gsub("1/0", 1, x)})
snps_geno <- apply(snps_geno, 2, function(x) {gsub("0/0", 0, x)})
snps_geno <- apply(snps_geno, 2, function(x) {gsub("1/1", 2, x)})
snps_geno <- apply(snps_geno[,1:length(colnames(snps_geno))], 2, as.numeric) # convert to numeric

snps_id <- as.data.frame(snps_id)
snps <- cbind(snps_id, snps_geno)
colnames(snps)[1] <- "SNP"
snps <- as.data.frame(snps)
head(snps)
typeof(snps$C001)
snps <- snps[,-c(2:3)]
sumStats$Chr <- sumStats$ARS_CHR
sumStats$Pos <- sumStats$ARS_POS
sumStats$REF <- sumStats$REF
sumStats$ALT <- sumStats$ALT

head(snps)


#sumStats <- sumStats[sumStats$GWAS_ID %in% snps$SNP,] # check that all SNPs get through




# renaming


sumStats <- sumStats %>% select(2,3,1,6,4,5,7,9,8)
colnames(sumStats) <- c("chr", "pos", "GWAS_ID", "pval", "ref", "alt", "beta", "Z", "se")
head(sumStats)
#sumStats <- sumStats %>% dplyr::select(- c(10, 11))
typeof(sumStats$pval)
#args[3] = path to MeTWAS model
model_list <- list.files(path = paste0(args[3]))
print(length(model_list))
TWAS_final <- list("Gene","Z","P","permute.P","TOP_SNP","topSNP.P") # set up empty list


colnames(snps)[1] <- "SNP"
colnames(snps)[2] <- "chr"
colnames(snps)[3] <- "pos"
typeof(snps$pos)
counter = 0
for (model in model_list){
  counter = counter + 1
  print(counter)
  model = paste0(args[3], "/",model)
  print(model)
  TWAS <- burdenTest_local(wgt = paste0(model), ### .RData file name for a given gene's MOSTWAS model
                       snps = snps, ### data.frame/data.table for SNP dosages
                       sumStats = sumStats, ### data.frame/data.table for GWAS summary statistics
                       snpAnnot = NULL,
                       Z = "Z",
                       onlyCis = F, ### toggle to include only local SNPs
                       beta = "beta", ### column name for effect sizes in GWAS summary stats
                       se = "se", ### column name for standard errors for betas in GWAS summary stats
                       chr = "chr", ### column name for SNP chromosome locations in GWAS summary stats
                       pos = "pos", ### column name for SNP basepair position in GWAS summary stats
                       ref = "ref", ### column name for SNP reference allele in GWAS summary stats
                       pval = "pval", ### column name for P-value in GWAS summary stats
                       R2cutoff = 0.01, ### cross-validation R2 cutoff to test a gene
                       alpha = 0.05/length(model_list), ### P-value cutoff in TWAS association to conduct permuation test (0.05 / 4082 = number of tests)
                       nperms = 1e3)### number of permutations in permuation test)
   TWAS_final <- purrr::map2(TWAS_final, TWAS, rbind)
}
TWAS_df <- as.data.frame(TWAS_final)
colnames(TWAS_df) <- TWAS_df[1,]
TWAS_df <- TWAS_df[-1,]
write.table(TWAS_df, args[4], sep = "\t", quote = F, row.names = F)
