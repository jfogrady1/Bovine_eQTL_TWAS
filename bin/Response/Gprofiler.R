library(tidyverse)
library(gprofiler2)
library(data.table)
library(ggplot2)
library(ggrepel)

args = commandArgs(trailingOnly = TRUE)
DE_results <- fread(args[1])

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




results_0.01 <- fread()

library(VennDiagram)

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")




result_0.01_private <- results_0.01 %>% filter(!(term_name %in% result_df$term_name))

dim(result_0.01_private)

result_0.01_private$term_name


result_0.05_private <- result_df %>% filter(!(term_name %in% results_0.01$term_name))



results <- gost(query = DEG_0.05,organism = "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = background, domain_scope = "known", user_threshold = 0.05, sources = c("GO:BP", "KEGG", "REAC"), evcodes = T)

terms <- c("MDA-5 signaling pathway",
        "defense response to virus",
        "Interferon Signaling",
        "cytokine production",
        "type I interferon-mediated signaling pathway",
        "Human immunodeficiency virus 1 infection",
        "interleukin-27-mediated signaling pathway",
        "interferon-mediated signaling pathway", "innate immune response",
        "MyD88-independent toll-like receptor signaling pathway",
        "regulation of cytokine-mediated signaling pathway",
        "type II interferon-mediated signaling pathway",
        "programmed necrotic cell death",
        "NOD-like receptor signaling pathway",
        "RIG-I-like receptor signaling pathway",
        "ISG15 antiviral mechanism",
        "Yersinia infection",
        "regulation of transformation of host cell by virus",
        "Antiviral mechanism by IFN-stimulated genes",
        "ISG15 antiviral mechanism ",
        "complex of collagen trimers",
        "NF-kB activation through FADD/RIP-1 pathway mediated by caspase-8 and -10",
        "Host-pathogen interaction of human coronaviruses - interferon induction",
        "Toll Like Receptor 4 (TLR4) Cascade",
        "Toll Like Receptor 3 (TLR3) Cascade",
        "Integrin cell surface interactions"
        )
results$result <- results$result %>%
mutate(label = ifelse(term_name %in% terms, term_name, NA))
result_df <- results$result
    result_df <- result_df %>% mutate(alpha_value = ifelse(term_name %in% terms, 1, 0.5))
    color_palette <- viridis(5)
    pos <- position_jitter(width = 0.3, seed = 3)
    my_plot <- ggplot(result_df, aes(x = source, y = -log10(p_value))) +
    geom_jitter(aes(color = source), alpha = result_df$alpha_value, position = pos, size = 3) +  # Use shape instead of color
    scale_color_brewer(palette = "Dark2" , name = "Source") +  # Apply the color palette
    scale_y_continuous(limits = c(0, 13), breaks = seq(0:13)) +
    labs(y = expression(-log[10](italic(P)[adj])), x = "") +
    theme_bw() +
    labs(size = "Intersection\nSize") +
    theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black", face = "bold"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 21, color = "black", face = "bold"),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15)) +
    geom_label_repel(aes(x=source, y=-log10(p_value)), label = results$result$label, max.overlaps = 20,
                   size=3.0, color='black', fontface = "bold",
                   fill='#FFFFFF33',
                   position = pos,
                   box.padding = 1,
                   point.padding = 0,
                   force = 4
    ) +
    geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
    guides(color = guide_legend(override.aes = list(size = 5)))
my_plot
View(results$result)
ggsave(args[3], width = 12, height = 12, dpi = 600)
