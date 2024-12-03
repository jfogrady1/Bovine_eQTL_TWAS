library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)



ALL_cis_eVariants <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_raw_cis_eQTL_associations.txt.gz")
gtf <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.expr_tmm_inv.bed.gz")
plink_files = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_genotypes_tensor"
bim = read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_genotypes_tensor.bim")
fam = read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_genotypes_tensor.fam")
covariates = data.frame(t(read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.covariates.txt")))
genelist = (unique(ALL_cis_eVariants$phenotype_id))
temp_direc = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/heritability/ALL_files/"



# Get covariates in right format
discrete_covariates = covariates$Condition
discrete_covariates <- cbind(fam$V1, fam$V2, discrete_covariates)
write.table(discrete_covariates, file=paste(temp_direc,"discrete_covar.covar", sep = ""), quote=F, col.names=F, row.names=F)

continuous <- covariates %>% select(-c(Condition))
continuous <- cbind(fam$V1, fam$V2, continuous)
write.table(continuous, file=paste(temp_direc,"continuous_covar.qcovar", sep = ""), quote=F, col.names=F, row.names=F)

#gcta = "/home/workspace/jogrady/my-bin/gcta64/gcta64"

num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)

foreach(i = 1:length(genelist), .packages = c("dplyr")) %dopar% {
    cat(i, "/", length(genelist), "\n")
    
    # 2. Get gene info
    gene <- genelist[i]
    gene
    geneinfo <- gtf[match(gene, gtf$phenotype_id),]
    colnames(geneinfo)[1] <- "Chr"
    chr <- geneinfo[1] 
    c <- chr$Chr
    start <- geneinfo$start - 1e6  # 1Mb upstream of the gene
    end <- geneinfo$end + 1e6      # 1Mb downstream of the gene  
    start <- if_else(start < 0, 0, start)
    end <- if_else(end < 0, 0, end)

    # 3. Filter SNPs in the same chromosome
    chrsnps <- subset(bim, bim$V1 == c)
    # 4. Filter for cis-SNPs
    cissnps <- subset(chrsnps, chrsnps$V4 >= start & chrsnps$V4 <= end)
    # 5. Create SNP list
    snplist <- cissnps[,2]
    write.table(snplist, file=paste(temp_direc,"tmp.cis.SNPlist", sep = ""), quote=F, col.names=F, row.names=F)

    # write phenotype
    pheno <- as.numeric(geneinfo[1,5:length(colnames(geneinfo))])
    pheno
    pheno = cbind(fam$V1, fam$V2, pheno)    
    write.table(pheno, file=paste(temp_direc,gene,".phen", sep = ""), quote=F, col.names=F, row.names=F)
    
    # GCTA command to generate GR<
    runGCTAgrm <- paste0(gcta, " --bfile ", plink_files, " --make-grm --autosome-num 29 --extract ", temp_direc, "tmp.cis.SNPlist", "  --thread-num 1 --out ", temp_direc, gene)
    runGCTAgrmd <- paste0(gcta, " --bfile ", plink_files, " --make-grm-d --autosome-num 29 --extract ", temp_direc, "tmp.cis.SNPlist", "  --thread-num 1 --out ", temp_direc, gene,"_domi")
    write.table(paste0(temp_direc, gene,"\n",temp_direc, gene,"_domi.d","\n"), file = paste0(temp_direc,"add_domi_grm.txt"), quote=F, col.names=F, row.names=F)

    runGREML <- paste0(gcta, " --mgrm ", temp_direc, "add_domi_grm.txt --reml --reml-no-constrain --reml-bendV --reml-lrt 1 2 --pheno ", temp_direc, gene,".phen --covar ",  temp_direc, "discrete_covar.covar --qcovar ", temp_direc, "continuous_covar.qcovar --out ", temp_direc, gene,"_add_domi")
    
    system(runGCTAgrm)
    system(runGCTAgrmd)
    system(runGREML)
}

#/ENSBTAG00000044172.grm.id

stopCluster(cl)
