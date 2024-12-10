library(tidyverse)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

ALL_cis_eVariants <- fread(args[1])
print("HERE")
gtf <- fread(args[2])
plink_files = args[3]
bim = read.table(args[4])
fam = read.table(args[5])
gene = args[6]
temp_direc = args[7]
temp_direc2 = args[8]
txt_files = gsub("add_h2", "", temp_direc)


gcta = args[9]




# 2. Get gene info]
print(gene)
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
write.table(snplist, file=paste(temp_direc, gene, ".cis.SNPlist", sep = ""), quote=F, col.names=F, row.names=F)
write.table(snplist, file=paste(temp_direc2, gene, ".cis.SNPlist", sep = ""), quote=F, col.names=F, row.names=F)

# write phenotype
pheno <- as.numeric(geneinfo[1,5:length(colnames(geneinfo))])
pheno
pheno = cbind(fam$V1, fam$V2, pheno)    
write.table(pheno, file=paste(temp_direc,gene,".phen", sep = ""), quote=F, col.names=F, row.names=F)

# GCTA command to generate GR<
runGCTAgrm <- paste0(gcta, " --bfile ", plink_files, " --make-grm --autosome-num 29 --extract ", temp_direc, gene, ".cis.SNPlist", "  --thread-num 1 --out ", temp_direc, gene)
runGCTAgrmd <- paste0(gcta, " --bfile ", plink_files, " --make-grm-d --autosome-num 29 --extract ", temp_direc, gene, ".cis.SNPlist", "  --thread-num 1 --out ", temp_direc2, gene,"_domi")
write.table(paste0(temp_direc, gene,"\n",temp_direc2, gene,"_domi.d","\n"), file = paste0(txt_files,gene,"_add_domi_grm.txt"), quote=F, col.names=F, row.names=F)
runGREML_narrow <- paste0(gcta, " --reml --grm ", temp_direc, gene, " --pheno ", temp_direc, gene, ".phen --covar ", txt_files, "discrete_covar.covar --qcovar ", txt_files, "continuous_covar.qcovar --out ", temp_direc, gene,"_SNP_based") # Calculate SNP based heritability
runGREML_comparison <- paste0(gcta, " --mgrm ", txt_files,gene,"_add_domi_grm.txt --reml --reml-lrt 2 --pheno ", temp_direc, gene,".phen --covar ",  txt_files, "discrete_covar.covar --qcovar ", txt_files, "continuous_covar.qcovar --out ", temp_direc, gene,"_add_domi") # Comparison of the two

system(runGCTAgrm)
system(runGCTAgrmd)
system(runGREML_narrow)
system(runGREML_comparison)
