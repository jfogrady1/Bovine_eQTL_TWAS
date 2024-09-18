library(data.table)
library(tidyverse)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library("ggsci")
library(UpSetR)
args = commandArgs(trailingOnly = T)

#args[1] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_CH_TWAS_results.txt"
#args[2] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_LM_TWAS_results.txt"
#args[3] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_HOFR_TWAS_results.txt"
#args[4] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_ALL_TWAS_results.txt"
#args[5] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_CH_TWAS_results.txt"
#args[6] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_LM_TWAS_results.txt"
#args[7] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_HOFR_TWAS_results.txt"
#args[8] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/CONTROL_ALL_TWAS_results.txt"#
#args[9] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_CH_TWAS_results.txt"
#args[10] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_LM_TWAS_results.txt"
#args[11] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_HOFR_TWAS_results.txt"
#args[12] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/INFECTED_ALL_TWAS_results.txt"

ALL_CH <- fread(args[1]) %>% mutate(Group = "AA")
ALL_LM <- fread(args[2])  %>% mutate(Group = "AA")
ALL_HF <- fread(args[3])  %>% mutate(Group = "AA")
ALL_Multi <- fread(args[4])  %>% mutate(Group = "AA")
CON_CH <- fread(args[5])  %>% mutate(Group = "CT")
CON_LM <- fread(args[6])  %>% mutate(Group = "CT")
CON_HF <- fread(args[7])  %>% mutate(Group = "CT")
CON_Multi <- fread(args[8])  %>% mutate(Group = "CT")
REA_CH <- fread(args[9])  %>% mutate(Group = "RT")
REA_LM <- fread(args[10])  %>% mutate(Group = "RT")
REA_HF <- fread(args[11])  %>% mutate(Group = "RT")
REA_Multi <- fread(args[12])  %>% mutate(Group = "RT")

ALL_LM$GWAS <- "LM"
CON_LM$GWAS <- "LM"
REA_LM$GWAS <- "LM"

ALL_HF$GWAS <- "HF"
CON_HF$GWAS <- "HF"
REA_HF$GWAS <- "HF"

ALL_CH$GWAS <- "CH"
CON_CH$GWAS <- "CH"
REA_CH$GWAS <- "CH"

ALL_Multi$GWAS <- "Multi"
CON_Multi$GWAS <- "Multi"
REA_Multi$GWAS <- "Multi"

ALL_HF <- ALL_HF %>% filter(as.numeric(permute.P) < 0.05 & as.numeric(P) < 0.05/dim(ALL_HF)[1])
ALL_LM <- ALL_LM %>% filter(as.numeric(permute.P) < 0.05 & as.numeric(P) < 0.05/dim(ALL_LM)[1])
ALL_CH <- ALL_CH %>% filter(as.numeric(permute.P) < 0.05 & as.numeric(P) < 0.05/dim(ALL_CH)[1])
ALL_Multi <- ALL_Multi %>% filter(as.numeric(permute.P) < 0.05 & as.numeric(P) < 0.05/dim(ALL_Multi)[1])



CON_HF <- CON_HF %>% filter(as.numeric(permute.P) < 0.05 & as.numeric(P) < 0.05/dim(CON_HF)[1])
CON_LM <- CON_LM %>% filter(as.numeric(permute.P) < 0.05 & as.numeric(P) < 0.05/dim(CON_LM)[1])
CON_CH <- CON_CH %>% filter(as.numeric(permute.P) < 0.05 & as.numeric(P) < 0.05/dim(CON_CH)[1])
CON_Multi <- CON_Multi %>% filter(as.numeric(permute.P) < 0.05 & as.numeric(P) < 0.05/dim(CON_Multi)[1])


REA_HF <- REA_HF %>% filter(as.numeric(permute.P) < 0.05 & as.numeric(P) < 0.05/dim(REA_HF)[1])
REA_LM <- REA_LM %>% filter(as.numeric(permute.P) < 0.05 & as.numeric(P) < 0.05/dim(REA_LM)[1])
REA_CH <- REA_CH %>% filter(as.numeric(permute.P) < 0.05 & as.numeric(P) < 0.05/dim(REA_CH)[1])
REA_Multi <- REA_Multi %>% filter(as.numeric(permute.P) < 0.05 & as.numeric(P) < 0.05/dim(REA_Multi)[1])



REA_HF




listInput <- list(ALL_CH = ALL_CH$Gene, 
                  ALL_LM = ALL_LM$Gene,
                  ALL_HF = ALL_HF$Gene,
                  ALL_Multi = ALL_Multi$Gene,
                  CON_CH = CON_CH$Gene, 
                  CON_LM = CON_LM$Gene,
                  CON_HF = CON_HF$Gene,
                  CON_Multi = CON_Multi$Gene,
                  REA_CH = REA_CH$Gene, 
                  REA_LM = REA_LM$Gene,
                  REA_HF = REA_HF$Gene,
                  REA_Multi = REA_Multi$Gene
                  )
listInput
upset(fromList(listInput), order.by = "freq", nsets = 12,
sets.x.label = "TWAS genes", point.size = 4, line.size = 2,
mainbar.y.label = "TWAS genes intersections", sets.bar.color = brewer.pal(n=12, "Paired"),
text.scale = 2, shade.alpha = 0.5)


pdf(args[13], width = 15, height = 12)
upset(fromList(listInput), order.by = "freq", nsets = 12,
sets.x.label = "TWAS genes", point.size = 4, line.size = 2,
mainbar.y.label = "TWAS genes intersections", sets.bar.color = brewer.pal(n=12, "Paired"),
text.scale = 2, shade.alpha = 0.5)
dev.off()

annotation <- fread(args[14]) %>% select(3,2)
annotation <- annotation[!(duplicated(annotation$Ensembl)),]
dim(annotation)

# Get intersections
df2 <- data.frame(gene=unique(unlist(listInput)))
head(df2)

df1 <- lapply(listInput,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()
head(df_int)
df_int %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))
length(unique(df2$gene))

df_int$int <- gsub("_", ":", df_int$int)
head(df_int)
df_int <- left_join(df_int, annotation, by = c("gene" = "Ensembl"))
write.table(df_int, args[15], sep = "\t", row.names = F, col.names = T, quote = F)
