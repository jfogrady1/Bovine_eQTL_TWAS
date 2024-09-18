##############################################################################################################################################################
##################################################### 1. Package installation and loading libraries ##########################################################
##############################################################################################################################################################

suppressPackageStartupMessages(library(tidyverse, quietly = T))
suppressPackageStartupMessages(library(dplyr, quietly = T))
suppressPackageStartupMessages(library(DESeq2, quietly = T))
suppressPackageStartupMessages(library(ggplot2, quietly = T))
suppressPackageStartupMessages(library(vsn, quietly = T))
suppressPackageStartupMessages(library(pheatmap, quietly = T))
suppressPackageStartupMessages(library(apeglm, quietly = T))
suppressPackageStartupMessages(library(hexbin, quietly = T))
suppressPackageStartupMessages(library(tidyquant, quietly = T))
suppressPackageStartupMessages(library(sva, quietly = T))
library(data.table)
library(gprofiler2)
library(ggrepel)
library(viridis)
library("RColorBrewer")
# Set the proper directory where the files are located

args = commandArgs(trailingOnly = T)
##############################################################################################################################################################
##################################################### 2. Sorting input data into the correct format ##########################################################
##############################################################################################################################################################

# DESeq2 requires 

  # 1. Count matrix, 
  # 2. Metadata file  for input into the DESeq2 object (Generated below)

counts <- as.matrix(read.table(args[1], skip = 0, header = T, sep="\t",row.names=1)) # matrix file
coldata <- as.data.frame(read.table(args[2], sep = '', skip = 0, header = T, row.names=1)) # metadata file

# Specify the count data as a matrix
# Specify the metadata file as a dataframe - annoying but that is how it is

res_no_sva <- function(counts_in, coldata_in) {
     output <- list()

    all(rownames(coldata) == colnames(counts))


    # Read in the eigen vector data from plink
    genotype_PCs <- read.table(args[3])
    genotype_PCs <- genotype_PCs %>% select(-1)
    genotype_PCs$V2 <- gsub("_.*", "", genotype_PCs$V2)
    rownames(genotype_PCs) <- genotype_PCs$V2
    genotype_PCs <- genotype_PCs %>% select(-1)
    genotype_PCs <- genotype_PCs %>% select(1,2)
    colnames(genotype_PCs) <- paste0("Genotype_PC", 1:2)


    # Read in other known covaraites
    coldata <- cbind(coldata_in, genotype_PCs)
    coldata$Batch <- factor(coldata$Batch, levels = c("1", "2"))
    coldata$Condition <- factor(coldata$Condition, levels = c("Control", "Infected"), labels = c("Control", "Reactor"))
    coldata$Age <- scale(coldata$Age, center = TRUE)
    ddsMat <- DESeqDataSetFromMatrix(countData = counts_in,
                                    colData = coldata,
                                    design = ~ Batch + Genotype_PC1 + Genotype_PC2 + Age + Condition)

    # Filter for low gene counts
    keep <- rowSums(counts(ddsMat) >= 6) >= (ncol(counts) * 0.2) # remove low count genes
    ddsMat <- ddsMat[keep,]

    # Bring in Admixture analysis for PCA
    file <-args[4]


    cat <- c(rep("Control", 63), rep("Infected", 60)) %>% as.data.frame()
    rows <- c(
  "C001", "C002", "C003", "C004", "C005", "C006", "C007", "C008", "C009", "C010",
  "C011", "C012", "C013", "C014", "C015", "C016", "C017", "C018", "C019", "C020",
  "C021", "C022", "C023", "C024", "C025", "C026", "C027", "C028", "C029", "C030",
  "C031", "C033", "C034", "C035", "C036", "C037", "C038", "C039", "C040", "C041",
  "C042", "C043", "C044", "C045", "C046", "C047", "C048", "C049", "C050", "C051",
  "C052", "C053", "C054", "C055", "C056", "C057", "C058", "C059", "C060", "C061",
  "C062", "C063", "C064", "T001", "T002", "T003", "T004", "T005", "T006", "T007",
  "T008", "T009", "T010", "T011", "T013", "T014", "T015", "T016", "T017", "T018",
  "T019", "T020", "T022", "T023", "T024", "T026", "T028", "T029", "T030", "T031",
  "T032", "T033", "T034", "T035", "T036", "T037", "T038", "T039", "T040", "T041",
  "T042", "T043", "T044", "T045", "T046", "T047", "T048", "T049", "T050", "T051",
  "T052", "T053", "T054", "T055", "T056", "T057", "T058", "T059", "T060", "T061",
  "T062", "T063", "T064")
    rows <- as.data.frame(rows)
    

    alist <- pophelper::readQ(files = file)
    attributes(alist)
    attributes(alist[[1]])
    rownames(alist[[1]]) <- rows$rows
    rownames(alist$SNP_Pruned.2.Q) <- rows$rows
    alist$SNP_Pruned.2.Q <- alist$SNP_Pruned.2.Q %>%  mutate(ADMIX_bin = cut((Cluster1 * 100), right = F, breaks=c(0.000,1,20, 40, 60, 80, 100), 
    labels = c("0", "1-20", "20-40", "40-60", "60-80", "80-100"))) %>% select(ADMIX_bin)

    ddsvsd <- vst(ddsMat)
    pcaData <- plotPCA(ddsvsd, intgroup=c("Condition"), ntop = 1500, returnData=TRUE) # top 1500 variable genes
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    pcaData <- cbind(pcaData, alist$SNP_Pruned.2.Q$ADMIX_bin)
    colnames(pcaData)[6] <- "%Holstein"
    pca <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape = `%Holstein`)) +
    geom_point(size=4, alpha = 0.8) +
    scale_color_manual(values = c("#2166ac", "#b2182b")) +
    scale_shape_manual(values = c(15,16,17,18,11,13)) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() + theme_bw() +
    theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12))

    # Perform the differential expression
    dds <- DESeq(ddsMat)
    res <- results(dds)
    res <- results(dds, contrast=c("Condition","Control","Reactor"))
    # Perform the LF shrink method
    res <- lfcShrink(dds, coef="Condition_Reactor_vs_Control", type="apeglm")




    res_df <- as.data.frame(res)
    head(res_df)
    
    res_df$padj <- as.numeric(res_df$padj)
    
    all_symbols <- gconvert(query = rownames(res_df), organism = "btaurus", 
         target="ENSG", mthreshold = Inf, filter_na = TRUE)
    

    res_df$Symbol <- all_symbols$name 
    

    res_df$diffexpressed <- "Not DE"
    # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
    res_df$diffexpressed[res_df$log2FoldChange > 0 & res_df$padj < 0.05] <- "DE Up"
    # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
    res_df$diffexpressed[res_df$log2FoldChange < 0 & res_df$padj < 0.05] <- "DE Down"

    res_df <- res_df %>% mutate(PLOT_Symbol =
                                case_when(
                                    log2FoldChange < -.1 & padj < 0.00005 ~ Symbol,
                                    log2FoldChange > 1 & padj < 0.00005 ~ Symbol,
                                    log2FoldChange > 0 & padj < 0.0000005 ~ Symbol,
                                    log2FoldChange < -0.5 & padj < 0.01 ~ Symbol,
                                    log2FoldChange > 1.5 & padj < 0.001 ~ Symbol,
                                    padj < 0.00005 ~ Symbol,
                                    FALSE ~ ""))

    tab <- table(res_df$diffexpressed)
    tab
    DE_down <- tab[1]
    DE_up <- tab[2]
    Volcano <- ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=PLOT_Symbol)) +
    geom_point(size = 1, alpha = 0.5) + # increased point size and added alpha for transparency
    scale_color_manual("Comparison", values=c("#2166ac", "#b2182b", "grey")) +
    labs(x=expression(log[2]("fold change")),
       y=expression(-log[10](italic(P)[adj])))+
    scale_x_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3)) +
    scale_y_continuous(limits = c(0,8.2), breaks = c(0,1,2,3,4,5,6,7,8)) +
    geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
    geom_text_repel(colour = "black", fontface = 4, max.overlaps = 40, size = 3.5) +
    theme_bw() + # changed to a classic theme for a clean look
    theme(axis.text.x = element_text(size = 15, colour = "black"),
        axis.text.y = element_text(size = 15, colour = "black"),
        axis.title.y = element_text(size = 21, color = "black"),
        axis.title.x = element_text(size = 21, color = "black"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "black", face = "bold"),
        panel.grid.minor = element_blank()) +
    annotate("text", x=-3, y=-log10(1e-8), size = 5, label=paste(DE_down), col="#2166ac", fontface = "bold", alpha = 0.8) +
    annotate("text", x=-2.66, y=-log10(1e-8), label= "↓", col="#2166ac", size = 10, fontface = "bold", alpha = 0.8) +
    annotate("text", x=3, y=-log10(1e-8), size = 5, label=paste(DE_up), col="#b2182b", fontface = "bold", alpha = 0.8) +
    annotate("text", x=2.66, y=-log10(1e-8), label="\u2191", col="#b2182b", size = 10, fontface = "bold", alpha = 0.8) +
    guides(color = guide_legend(override.aes = list(size = 4))) + coord_flip()
  
    # removed panel background
    res_final <- res_df %>%
    filter(padj <= 0.05)
    res_final <- res_final %>% select(-c("diffexpressed", "PLOT_Symbol"))
    res_df <- res_df %>% filter(!is.na(padj)) # Only include genes tested
 
    res_final <- res_final[order(res_final$padj, decreasing = F),]

    ############################################
    ############################################
    ############################################
    # Enrichment analysis ######################
    ############################################

    res_filtered <- res_final %>% filter(padj < 0.01) # highly significant results
    genes <- rownames(res_filtered)

   
    genes <- list(genes)
    results <- gost(query = genes,organism = "btaurus", correction_method = "fdr", ordered_query = T, domain_scope = "known", custom_bg = rownames(res_df), user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC"), evcodes = T)
  
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

    # Basic plot from gprofiler
    plot <- gostplot(results, capped=FALSE, interactive=FALSE, pal=c(`GO:BP`="#ff9900", `GO:CC` ="purple", `KEGG`= "#dd4477",`REAC`="#3366cc")) + labs(y=bquote(~-log[10]~italic(P)[adj]), title = "", x = "Database") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12)) +
         geom_text_repel(aes(label = results$result$label), colour = "black", max.overlaps = 40,
                  nudge_x = 0, nudge_y = 1.2, size = 3, fontface = "bold")
    
    #############################################
    #############################################
    # Generate my own jitter plot ###############
    #############################################
    
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

    results$meta
    # Store outputs
    output$volcano <- Volcano # Volcano plot
    output$all_results <- res_final # All significant genes (bovine IDS + gene_symbol )
    output$signif_results <- res_df # All tested genes (bovine IDs + gene symbol)
    output$gprofiler_results <- results$result # Terms for Gene Ontology analysis
    output$gprofiler_plot <- plot # GOST plot
    output$my_gprofiler <- my_plot
    output$pca <- pca
    return(output)

}
result_normal <- res_no_sva(counts,coldata)

Kirsten <- c("ENSBTAG00000000306",
             "ENSBTAG00000000507",
             "ENSBTAG00000001060",
             "ENSBTAG00000002758",
             "ENSBTAG00000003553",
             "ENSBTAG00000003650",
             "ENSBTAG00000004305",
             "ENSBTAG00000006806",
             "ENSBTAG00000008182",
             "ENSBTAG00000008353",
             "ENSBTAG00000009354",
             "ENSBTAG00000013125",
             "ENSBTAG00000016163",
             "ENSBTAG00000019716",
             "ENSBTAG00000021766",
             "ENSBTAG00000031707",
             "ENSBTAG00000035224",
             "ENSBTAG00000037608",
             "ENSBTAG00000039037")

res_kirsten <- result_normal$signif_results[Kirsten,]
total_genes_FDR1 <- dim(result_normal$all_results %>% filter(padj < 0.01))[1]
total_genes <- dim(result_normal$signif_results)[1]
result_df <- result_normal$signif_results %>% filter(padj < 0.01) %>% select(Symbol)
res_kirsten
res_kirsten <- res_kirsten %>% filter(res_kirsten$Symbol %in% result_df$Symbol)
head(res_kirsten)

total_genes <- dim(result_normal$signif_results)[1]
genes_of_interest = 16 # number of genes tested in our study
total_genes_FDR1 <- dim(result_normal$all_results %>% filter(padj < 0.01))[1]
overlap_genes <- dim(res_kirsten)[1]

# Hyper geometric test for enrichemnt of kirsten's set in ours
p_value <- phyper(overlap_genes -1, total_genes_FDR1, total_genes - total_genes_FDR1, genes_of_interest, lower.tail = F)
p_value

result_normal$hypergeom_kirsten <- p_value


# Write files

# First plots

result_normal$volcano
ggsave(args[5], width = 12, height = 12, dpi = 600)
result_normal$my_gprofiler
ggsave(args[6], width = 12, height = 12, dpi = 600)
result_normal$pca
ggsave(args[7], width = 12, height = 12, dpi = 600)

# DE results
write.table(result_normal$signif_results, file = args[8], sep = "\t", row.names = T, col.names = T, quote = F)
write.table(result_normal$all_results, file = args[9], sep = "\t", row.names = T, col.names = T, quote = F)

gem <- result_normal$gprofiler_results[,c("query", "term_id", "query_size", "intersection_size", "term_name", "p_value", "intersection")]
colnames(gem)[6] <- "FDR"

write.table(gem, file = args[10],sep = "\t", row.names = T, col.names = T, quote = F)

write.table(p_value, file = args[11], sep = "\t", row.names = T, col.names = F, quote = F)
