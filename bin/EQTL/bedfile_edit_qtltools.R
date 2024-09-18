#  Need to modify covariates as qtltools requires the first column to be called "id"
library(data.table)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
data <- fread(args[1])
group = args[3]
if (group == "ALL") {
    ids <- colnames(data)[5:127] 
}

if (group == "CONTROL") {
    ids <- colnames(data)[5:67] 
}

if (group == "ALL") {
    ids <- colnames(data)[5:64] 
}

64 -5 

colnames(data)[1:4] <- c("#chr", "start", "end", "gid")

data$pid <- NA
data$strand <- NA
length(colnames(data))
# chr, start, end, pid, gid, strand, sample
data <- data %>% select(1,2,3,(length(colnames(data)) - 1), 4, (length(colnames(data))),5:(length(colnames(data))-2))
colnames(data)[1] <- "chr"
data$chr <- gsub("chr", "", data$chr)
rownames(data) <- 1:nrow(data)
colnames(data)[1] <- "#chr"
write.table(data, file = args[4], row.names = F, col.names = T, sep = "\t", quote = F)
