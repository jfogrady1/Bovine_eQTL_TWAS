library(tidyverse)
library(data.table)
args = commandArgs(trailingOnly = TRUE)
bim = read.table(args[1])
fam = read.table(args[2])
covariates = data.frame(t(read.table(args[3])))


discrete_covariates = covariates$Condition
discrete_covariates <- cbind(fam$V1, fam$V2, discrete_covariates)
write.table(discrete_covariates, file=paste(args[4], sep = ""), quote=F, col.names=F, row.names=F)

continuous <- covariates %>% select(-c(Condition))
continuous <- cbind(fam$V1, fam$V2, continuous)
write.table(continuous, file=paste(args[5], sep = ""), quote=F, col.names=F, row.names=F)