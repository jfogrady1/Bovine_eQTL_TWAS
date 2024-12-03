library(tidyverse)
library(gprofiler2)
library(data.table)

DE_results <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/ALL_results.txt")

DE_results <- DE_results[order(DE_results$padj, decreasing = F),]
head(DE_results)
DEG_0.05 <- DE_results %>% filter(padj < 0.05) %>% select(Symbol) %>% as.vector()

background <- DE_results %>% select(Symbol) %>% as.vector()

length(DEG_0.05)
DEG_0.05 <- DEG_0.05$Symbol
background <- background$Symbol


results_0.05 <- gost(query = DEG_0.05,organism = "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = background, domain_scope = "known", user_threshold = 0.05, sources = c("GO:BP", "KEGG", "REAC"), evcodes = T)
results_0.05$result$query <- ""
result_df <- results_0.05$result

View(result_df)


results_0.01 <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/RNA-seq/DESEQ2/Gprofiler_results.txt")

library(VennDiagram)

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")



library(ggVennDiagram)
ggVennDiagram(list(results_0.01$term_name, result_df$term_name), color = "black", lwd = 0.8, lty = 1,
category.names = c("FDR < 0.01", "FDR < 0.05")) +
scale_fill_distiller(palette = "RdBu") +
         scale_x_continuous(expand = expansion(mult = .2))

ggsave("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/gprofiler/Venn_diagram.pdf", width = 12, height = 12, dpi = 600)


result_0.01_private <- results_0.01 %>% filter(!(term_name %in% result_df$term_name))

dim(result_0.01_private)

result_0.01_private$term_name

View(result_0.01_private)

result_0.05_private <- result_df %>% filter(!(term_name %in% results_0.01$term_name))

dim(result_0.05_private)
View(result_0.05_private)
