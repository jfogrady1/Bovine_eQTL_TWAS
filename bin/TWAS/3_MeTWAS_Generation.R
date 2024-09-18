# Get everything set up
# Getting everything set up
# changed the gcta64 to gcta and also commented in some key components
# local model 2 to increase distance to 1Mb, comments in
# MeTWAS model generation
# Also note the heritability cutoff




###############################################################################################################################################
################################################### Packages ##################################################################################
###############################################################################################################################################

args = commandArgs(trailingOnly=TRUE)
library("bigsnpr")
library("MOSTWAS")
library("MatrixEQTL")
library("dplyr")
library("stringr")
library("data.table")
trainLocalModel_local <- function(phenoInt,
                           midSNP,
                           mediator,
                           medLocs,
                           covariates,
                           cisDist = 1e6,
                           nfolds = 5,
                           seed = 1218,
                           verbose = T){

  if (verbose){
    print(phenoInt)
  }
  colnames(mediator)[1] = 'Mediator'
  pheno = as.numeric(mediator[mediator$Mediator == phenoInt,-1])
  ml = subset(medLocs,geneid == phenoInt)
  w = which(midSNP$map$chromosome == ml$chr[1] &
              midSNP$map$physical.pos < ml$right[1] + cisDist &
              midSNP$map$physical.pos > ml$left[1] - cisDist)
  if (length(w) == 0 | length(w) == 1){
    return(list(Model = data.frame(SNP = NA,
                                   Chromosome = NA,
                                   Position = NA,
                                   Effect = NA),
                Predicted = rep(NA,length = length(pheno)),
                CVR2 = 0))
  }
  midSNP.cur = bigsnpr::snp_attach(subset(midSNP,ind.col = w))

  set.seed(seed)
  print("Performing train test split")
  train = caret::createFolds(y = pheno,
                             k=nfolds,
                             returnTrain = T)
  set.seed(seed)
  test = caret::createFolds(y = pheno,
                            k = nfolds,
                            returnTrain = F)

  pred.blup = pred.enet = vector(mode = 'numeric',
                                 length = length(pheno))
  df = cbind(pheno,t(as.matrix(covariates[,-1])))
  colnames(df)[1] = 'pheno'
  pheno = as.numeric(resid(lm(pheno ~ .,data = as.data.frame(df))))

  mod.enet = glmnet::cv.glmnet(x = midSNP.cur$genotypes[],
                               y = pheno,
                               nfolds = nfolds,
                               keep = T)
  pred.enet = mod.enet$fit.preval[,which.min(mod.enet$cvm)]

  for (i in 1:5){
    mod.blup = rrBLUP::mixed.solve(y = pheno[train[[i]]],
                               Z = midSNP.cur$genotypes[train[[i]],])
    pred.blup[test[[i]]] = as.numeric(midSNP.cur$genotypes[test[[i]],] %*%
                                        mod.blup$u)
  }
  enet.R2 = adjR2(pheno,pred.enet)
  blup.R2 = adjR2(pheno,pred.blup)
  model = ifelse(blup.R2 >= enet.R2,
                 'LMM','Elastic net')

  if (model == 'Elastic net'){

    mod.df.enet = data.frame(SNP = midSNP.cur$map$marker.ID,
                             Chromosome = midSNP.cur$map$chromosome,
                             Position = midSNP.cur$map$physical.pos,
                             Effect = coef(mod.enet,s='lambda.min')[-1])
    colnames(mod.df.enet) = c('SNP','Chromosome','Position','Effect')
    mod.df.enet = subset(mod.df.enet,Effect != 0)
    if (nrow(mod.df.enet) != 0){
      return(list(Model = mod.df.enet,
                  Predicted = pred.enet,
                  CVR2 = max(enet.R2,blup.R2)))
    } else { model = 'LMM' }
    }
    if (model == 'LMM'){
      mod.blup = rrBLUP::mixed.solve(y = pheno,
                                     Z = midSNP.cur$genotypes[])
      mod.df.blup = data.frame(SNP = midSNP.cur$map$marker.ID,
                               Chromosome = midSNP.cur$map$chromosome,
                               Position = midSNP.cur$map$physical.pos,
                               Effect = mod.blup$u)
      colnames(mod.df.blup) = c('SNP','Chromosome','Position','Effect')
      return(list(Model = mod.df.blup,
                  Predicted = pred.blup,
                  CVR2 = max(enet.R2,blup.R2)))

    }
  }


MeTWAS_local <- function(geneInt,
                   snpObj,
                   mediator,
                   medLocs,
                   covariates,
                   dimNumeric,
                   qtlFull,
                   h2Pcutoff = 0.05,
                   numMed = 10,
                   seed = 1218,
                   k = 5,
                   cisDist = 1e6,
                   parallel = F,
                   prune = F,
                   ldThresh = .5,
                   cores,
                   verbose = T,
                   R2Cutoff = 0.01,
                   modelDir,
                   tempFolder,
                   gctaFolder = NULL){
  
  set.seed(seed)
  if (!dir.exists(modelDir)){
    dir.create(modelDir)
  }
  
  if (is.null(gctaFolder)){
    gctaFolder = ''
  }
  
  print('GATHERING MEDIATORS')
  pheno = as.numeric(mediator[mediator$Mediator == geneInt,-1])
  pheno = (pheno - mean(pheno))/sd(pheno)
  medList = gatherMediators(geneInt,qtlFull,numMed)
  
  print('ESTIMATING HERITABILITY')
  lll = c(geneInt,medList)
  lll = lll[!is.na(lll)]
  w = c()
  for (i in 1:length(lll)){
    
    ml = subset(medLocs,geneid == lll[i])
    w = c(w,which(snpObj$map$chromosome == ml$chr[1] &
                    snpObj$map$physical.pos < ml$right[1] + cisDist &
                    snpObj$map$physical.pos > ml$left[1] - cisDist))
    
  }
  if (length(w) == 0){
    return('SNPs not found')
  }
  midSNPfile = subset(snpObj,ind.col = w)
  midSNP = bigsnpr::snp_attach(midSNPfile)
  midSNP$fam$affection = pheno
  
  if (prune){
    print('RUNNING LD CLUMPING')
    keep = bigsnpr::snp_clumping(
      midSNP$genotypes,
      infos.chr = midSNP$map$chromosome,
      ind.row = rows_along(midSNP$genotypes),
      S = NULL,
      thr.r2 = ldThresh,
      size = 100/ldThresh,
      exclude = NULL,
      ncores = ifelse(parallel,cores,1)
    )
    midSNP = bigsnpr::snp_attach(subset(midSNP,ind.col=keep))
    
  }
  
  tmpBed = paste0(tempFolder,"MeTWAS_",geneInt,".bed")
  bigsnpr::snp_writeBed(midSNP,tmpBed)
  
  system(paste(paste0(gctaFolder,'gcta64'),
               '--bfile',strsplit(tmpBed,'.bed')[[1]][1],
               '--autosome-num 29',
               '--autosome --make-grm',
               '--out',strsplit(tmpBed,'.bed')[[1]][1]),
         intern = !verbose)
  
  phenFile = paste0(strsplit(tmpBed,'.bed')[[1]][1],'.phen')
  covarFile = paste0(strsplit(tmpBed,'.bed')[[1]][1],'.qcovar')
  write.table(midSNP$fam[,c(1,2,6)],phenFile,row.names = F,
              col.names = F, quote = F)
  write.table(cbind(midSNP$fam[,1:2],
                    t(covariates[1:dimNumeric,-1])),
              covarFile,row.names = F,col.names = F,quote=F)
  system(paste(paste0(gctaFolder,'gcta64'),
               '--reml --reml-no-constrain --reml-bendV',
               '--grm',strsplit(tmpBed,'.bed')[[1]][1],
               '--pheno',phenFile,
               '--qcovar',covarFile,
               '--autosome-num 29',
               '--out',paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi')),
         intern = !verbose)
  
  if (file.exists(paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi.hsq'))){
    a = data.table::fread(paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi.hsq'),fill=T)
    herit = list(h2 = a$Variance[a$Source == 'V(G)/Vp'],
                 P = a$Variance[a$Source == 'Pval'])
  } else {
    system(paste(paste0(gctaFolder,'gcta64'),
                 '--reml',
                 '--grm',strsplit(tmpBed,'.bed')[[1]][1],
                 '--pheno',phenFile,
                 '--qcovar',covarFile,
                 '--autosome-num 29',
                 '--out',paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi')),
           intern = !verbose)
    if (file.exists(paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi.hsq'))){
      a = data.table::fread(paste0(strsplit(tmpBed,'.bed')[[1]][1],'_multi.hsq'),fill=T)
      herit = list(h2 = a$Variance[a$Source == 'V(G)/Vp'],
                   P = a$Variance[a$Source == 'Pval'])
    } else {
      herit = list(h2 = 0,P=1)
    }}
  
  print("HERITABILITY ESTIMATION COMPLETE")
  if (herit$P > h2Pcutoff) {
    return(paste(geneInt,
                 'is not germline heritable at P <',
                 h2Pcutoff))
  }
  
  
  print('TRAINING MEDIATORS')
  
  if (parallel) {
    medTrainList = parallel::mclapply(medList,
                                      trainLocalModel_local,
                                      midSNP = midSNP,
                                      mediator = mediator,
                                      medLocs = medLocs,
                                      covariates = covariates,
                                      cisDist = cisDist,
                                      seed = seed,
                                      nfolds = k,
                                      mc.cores = cores)
  }
  if (!parallel){
    medTrainList = lapply(medList,
                          trainLocalModel_local,
                          midSNP = midSNP,
                          mediator = mediator,
                          medLocs = medLocs,
                          covariates = covariates,
                          cisDist = cisDist,
                          seed = seed,
                          nfolds = k)
  }
  print("Completed Training of Mediators")
  names(medTrainList) = medList
  medTrainList = medTrainList[sapply(medTrainList,length) > 0]
  
  print('FITTING MEDIATORS')
  fe.R2 = 0
  if (length(medTrainList) > 0){
    medTrainList = medTrainList[as.numeric(which(sapply(medTrainList,
                                                        function(x) x[3]) >= .01))]
    if (length(medTrainList) > 0){
      fixedEffects = as.data.frame(matrix(ncol = length(medTrainList),
                                          nrow = ncol(mediator)-1))
      colnames(fixedEffects) = names(medTrainList)
      for (i in 1:ncol(fixedEffects)){
        fixedEffects[,i] = medTrainList[[i]][2]
      }
      fixedEffects = as.data.frame(apply(fixedEffects,2,scale))
      fixedEffects$pheno = pheno
      lmCVFit <- lm(pheno~.,
                    data = fixedEffects)
      
      fe.R2 = fe.R2 + adjR2(as.numeric(predict(lmCVFit)),pheno)
      
      trans.mod.df = as.data.frame(abind::abind(lapply(1:length(medTrainList),
                                                       amplifyTrans,
                                                       medTrainList = medTrainList,
                                                       lmCaretObj = lmCVFit),
                                                along = 1))
      trans.mod.df$Effect = as.numeric(as.character(trans.mod.df$Effect))
      trans.mod.df = subset(trans.mod.df,SNP != 'Intercept')
      rownames(trans.mod.df) = NULL
      pheno = pheno - as.numeric(predict(lmCVFit))
    } else {
      pheno = pheno
      fixedEffects = NULL}
  }
  
  
  print('FITTING LOCAL-GENOTYPES')
  
  cisGenoMod = trainLocalModel_local(phenoInt = geneInt,
                               midSNP = midSNP,
                               mediator = mediator,
                               medLocs = medLocs,
                               covariates = covariates,
                               seed = seed,
                               nfolds = k,
                               cisDist = cisDist)
  if (is.null(cisGenoMod$Model)){
    cisGenoMod$Model = as.data.frame(matrix(nrow = 1,ncol = 4))
    colnames(cisGenoMod$Model) = c('SNP',
                                   'Chromosome',
                                   'Position',
                                   'Effect')
    cisGenoMod$Model[1,] = 0
  }
  
  
  cisGenoMod$Model$Mediator = 'Cis'
  if (exists('trans.mod.df')){
    cisGenoMod$Model = rbind(cisGenoMod$Model,trans.mod.df)
  }
  cisGenoMod$Model = subset(cisGenoMod$Model,Effect!=0)
  cisGenoMod$CVR2 = cisGenoMod$CVR2 + fe.R2
  cisGenoMod$CVR2.cis = cisGenoMod$CVR2 - fe.R2
  cisGenoMod$h2 = herit$h2
  cisGenoMod$h2.P = herit$P
  
  
  
  Model = cisGenoMod$Model
  R2 = cisGenoMod$CVR2
  Predicted = cisGenoMod$Predicted
  Mediators = cisGenoMod$medlist
  CisR2 = cisGenoMod$CVR2.cis
  h2 = abs(herit$h2)
  h2.Pvalue = herit$P
  print('****************************')
  print(paste0('R2 = ',R2))
  if (R2 < R2Cutoff){ return('CV R2 < 0.01.') }
  ## REMOVE THE NEXT LINE
  CorMat = cbind(Predicted,fixedEffects)
  save(Model,R2,Predicted,Mediators,CisR2,h2,h2.Pvalue,CorMat,
       file = paste0(modelDir,geneInt,'.wgt.med.RData'))
  rm(midSNP)
  file.remove(midSNPfile)
  
}

bedfile <- args[1]
rds <- snp_readBed(bedfile, backingfile = sub_bed(bedfile))
bigsnp_file <- snp_attach(rds)

# args[1] ALL_EXP_Mediators.txt'
mediator = fread(args[2])
colnames(mediator)
# 'ALL_EXP_NoMediators.txt'
exp = fread(args[3])

colnames(mediator)[1] <- "id"
colnames(exp)[1] <- "id"




colnames(exp) = colnames(mediator)
colnames(exp)
mediator = rbind(mediator,exp)
mediator = mediator[!duplicated(mediator$id),]
dim(mediator)

colnames(mediator)[1] <- "Mediator"

#args[3] 'ALL_EXP_MediatorsLoc.txt'
medLocs = fread(args[4])
colnames(medLocs) <- c("geneid", "chr", "left")
medLocs$right = medLocs$left + 1
#args[4] "gene_locnoMediators.txt"
geneLocs = fread(args[5])

colnames(geneLocs)[1] <- "id"
colnames(geneLocs)[3] <- "left"
geneLocs$right <- geneLocs$left + 1
colnames(medLocs) = colnames(geneLocs)
length(unique(medLocs$id))
medLocs = rbind(geneLocs,medLocs)#
head(medLocs)
medLocs = medLocs[!duplicated(medLocs$id),]
colnames(medLocs)[2] <- "chr"
dim(medLocs)
colnames(medLocs)[3] <- "left"
colnames(medLocs)[4] <- "right"
colnames(medLocs)



#"../eqtl_study/eqtl_nextflow/work/be/c5b14521fd3d75bbe08861fafaf163/ALL_COV_MatrixQTL_format.txt"
covariates = fread(args[6])
covariates = covariates %>% filter(id %in% c("genotypePC1","genotypePC2","genotypePC3","genotypePC4","genotypePC5"))
# "CH_Cis_FINAL.txt"
qtl_local <- fread(args[7]) 
#'CH_Trans_FINAL.txt'
#qtl_distal <- fread(args[8])

# "ALL_MEDIATOR_EQTLSFDR1.txt"
qtl_mediator <- fread(args[8]) 
qtl_mediator <- qtl_mediator %>% select(1,2,3,5,6)#%>% dplyr::select(snps, gene, t, beta, statistic, pvalue, FDR)
qtl_mediator$beta = 0 # keeping this as a dummy variable

head(qtl_mediator)

# final set up

qtlFull <- qtl_local %>% dplyr::select("snp", "gene", "beta", "t-stat", "p-value", "FDR")
qtl_mediator <- qtl_mediator %>% select(mediator, gene, beta, t, p, fdr)
colnames(qtl_mediator)[1] <- "SNP"
colnames(qtl_mediator)[4] <- "t-stat"
colnames(qtl_mediator)[5] <- "p-value"
colnames(qtl_mediator)[6] <- "FDR"
colnames(qtlFull)[1] <- "SNP"
colnames(qtlFull)[4] <- "t-stat"
colnames(qtlFull)[5] <- "p-value"
colnames(qtlFull)
colnames(qtl_mediator)
qtlFull <- rbind(qtlFull, qtl_mediator)
tail(qtlFull)
head(mediator)
head(medLocs)
colnames(mediator)[1] <- "Mediator"
colnames(medLocs) <- c("geneid", "chr", "left", "right")
mediator$Mediator <- as.character(mediator$Mediator)
medLocs$geneid <- as.character(medLocs$geneid)
geneList <- mediator$Mediator
length(unique(geneList))
counter = 0
parts <- unlist(strsplit(args[1], "/"))
part_of_interest = parts[10]
split_string = unlist(strsplit(part_of_interest, "_"))
group = split_string[1]
for (g in geneList){
  counter = counter + 1
  print(counter)
  print(g)
  if (!paste0(g,'.wgt.med.RData') %in% list.files(paste0("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/models/",group,"_MeTWAS_models_10_mediators/"))){
    MeTWAS_local(geneInt = g,
             snpObj = bigsnp_file,
             mediator = mediator,
             medLocs = medLocs,
             covariates = covariates,
             dimNumeric = 4,
             qtlFull = qtlFull,
             h2Pcutoff = 0.05,
             numMed = 10,
             seed = 1218,
             k = 5,
             cisDist = 1e6,
             parallel = T,
             prune = F,
             ldThresh = .5,
             cores = 1,
             verbose = T,
             R2Cutoff = 0.01,
             modelDir = paste0("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/models/",group,"_MeTWAS_models_10_mediators/"),
             tempFolder = paste0("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/models/",group,"_temp_folder/"),
             gctaFolder = paste0("/home/workspace/jogrady/my-bin/gcta64/"))
  }
  #system("rm *_sub*")
}

