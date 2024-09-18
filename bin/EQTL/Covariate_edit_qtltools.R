#  Need to modify covariates as qtltools requires the first column to be called "id"
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
data <- fread(args[1])
colnames(data)[1] <- c("id")
rownames(data) <- 1:nrow(data)
write.table(data, file = args[2], row.names = F, col.names = T, sep = "\t", quote = F)
