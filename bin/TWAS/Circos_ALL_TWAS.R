# Plotting circos of TWAS results
library(data.table)
library(tidyverse)
library(circlize)
args = commandArgs(trailingOnly = TRUE)
Reac_CH <- fread(args[1])
Reac_HF <- fread(args[2])
Reac_LM <- fread(args[3])
Reac_ALL <- fread(args[4])

ensemble <- fread(args[5])
ensemble <- ensemble %>% filter(V3 == "gene")
head(ensemble)
ensemble <- ensemble %>% separate(., V9, into = c("gene_id", "gene_version", "gene_name"), sep = ";")
ensemble$gene_id <- gsub("^gene_id ", "", ensemble$gene_id)
ensemble$gene_id <- gsub('"', '', ensemble$gene_id)

ensemble$gene_name <- gsub("gene_name ", "", ensemble$gene_name)
ensemble$gene_name <- gsub("gene_source ", "", ensemble$gene_name)
ensemble$gene_name <- gsub('"', '', ensemble$gene_name)
ensemble$gene_name <- if_else(ensemble$gene_name == " ensembl", ensemble$gene_id, ensemble$gene_name)
colnames(ensemble)[1] <- "chr"
ensemble <- ensemble %>% dplyr::select(gene_id, gene_name, chr, V4)
colnames(ensemble)[4] <- "pos"
#Reac_CH <- left_join(Reac_CH, ensemble, by = c("Gene" = "gene_id"))
head(Reac_CH)
Reac_CH$type <- "Charolais"
Reac_HF$type <- "Holstein-Friesian"
Reac_LM$type <- "Limousin"
Reac_ALL$type <- "Multi"
Reac_CH$Z <- as.numeric(Reac_CH$Z)
Reac_HF$Z <- as.numeric(Reac_HF$Z)
Reac_LM$Z <- as.numeric(Reac_LM$Z)
Reac_ALL$Z <- as.numeric(Reac_ALL$Z)


# get bonferroni 
cut_off_CH <- 0.05 / length(Reac_CH$Z)
cut_off_LM <- 0.05 / length(Reac_LM$Z)
cut_off_HF <- 0.05 / length(Reac_HF$Z)
cut_off_ALL<- 0.05 / length(Reac_ALL$Z)

cut_off_CH = -log10(cut_off_CH)
cut_off_CH = cut_off_CH * -1


cut_off_LM = -log10(cut_off_LM)
cut_off_LM = cut_off_LM * -1


cut_off_HF = -log10(cut_off_HF)
cut_off_HF = cut_off_HF * -1


cut_off_ALL = -log10(cut_off_ALL)
cut_off_ALL = cut_off_ALL * -1

Reac_CH$Z <- if_else(Reac_CH$Z == "Cannot compute LD covariance matrix" | Reac_CH$Z == "No Overlapping SNPs.", 0,  as.numeric(Reac_CH$Z))
Reac_HF$Z <- if_else(Reac_HF$Z == "Cannot compute LD covariance matrix" | Reac_HF$Z == "No Overlapping SNPs.", 0,  as.numeric(Reac_HF$Z))
Reac_LM$Z <- if_else(Reac_LM$Z == "Cannot compute LD covariance matrix" | Reac_LM$Z == "No Overlapping SNPs.", 0,  as.numeric(Reac_LM$Z))
Reac_ALL$Z <- if_else(Reac_ALL$Z == "Cannot compute LD covariance matrix" | Reac_ALL$Z == "No Overlapping SNPs.", 0,  as.numeric(Reac_ALL$Z))

Reac_CH$P <- if_else(Reac_CH$P == "Cannot compute LD covariance matrix" | Reac_CH$P == "No Overlapping SNPs.", 1,  as.numeric(Reac_CH$P))
Reac_HF$P <- if_else(Reac_HF$P == "Cannot compute LD covariance matrix" | Reac_HF$P == "No Overlapping SNPs.", 1,  as.numeric(Reac_HF$P))
Reac_LM$P <- if_else(Reac_LM$P == "Cannot compute LD covariance matrix" | Reac_LM$P == "No Overlapping SNPs.", 1,  as.numeric(Reac_LM$P))
Reac_ALL$P <- if_else(Reac_ALL$P == "Cannot compute LD covariance matrix" | Reac_ALL$P == "No Overlapping SNPs.", 1,  as.numeric(Reac_ALL$P))


Reac_CH$permute.P <- if_else(Reac_CH$permute.P == "Cannot compute LD covariance matrix" | Reac_CH$permute.P == "No Overlapping SNpermute.Ps.", 1,  as.numeric(Reac_CH$permute.P))
Reac_HF$permute.P <- if_else(Reac_HF$permute.P == "Cannot compute LD covariance matrix" | Reac_HF$permute.P == "No Overlapping SNpermute.Ps.", 1,  as.numeric(Reac_HF$permute.P))
Reac_LM$permute.P <- if_else(Reac_LM$permute.P == "Cannot compute LD covariance matrix" | Reac_LM$permute.P == "No Overlapping SNpermute.Ps.", 1,  as.numeric(Reac_LM$permute.P))
Reac_ALL$permute.P <- if_else(Reac_ALL$permute.P == "Cannot compute LD covariance matrix" | Reac_ALL$permute.P == "No Overlapping SNpermute.Ps.", 1,  as.numeric(Reac_ALL$permute.P))

Reac_CH <- left_join(Reac_CH, ensemble, by = c("Gene" = "gene_id"))
Reac_HF <- left_join(Reac_HF, ensemble, by = c("Gene" = "gene_id"))
Reac_LM <- left_join(Reac_LM, ensemble, by = c("Gene" = "gene_id"))
Reac_ALL <- left_join(Reac_ALL, ensemble, by = c("Gene" = "gene_id"))
head(Reac_CH)
Reac_CH$end <- Reac_CH$pos +1
Reac_HF$end <- Reac_HF$pos +1
Reac_LM$end <- Reac_LM$pos +1
Reac_ALL$end <- Reac_ALL$pos +1

Reac_CH$end <- Reac_CH$pos +1
dim(Reac_CH)
Reac_CH$signif <- if_else(as.numeric(Reac_CH$P) < (0.05/length(Reac_CH$P)), "Sig", "NotSig")
Reac_HF$signif <- if_else(as.numeric(Reac_HF$P) < (0.05/length(Reac_HF$P)), "Sig", "NotSig")
Reac_LM$signif <- if_else(as.numeric(Reac_LM$P) < (0.05/length(Reac_LM$P)), "Sig", "NotSig")
Reac_ALL$signif <- if_else(as.numeric(Reac_ALL$P) < (0.05/length(Reac_ALL$P)), "Sig", "NotSig")
Reac_CH$permute_sig <- if_else(as.numeric(Reac_CH$permute.P) < 0.05, "Sig", "NotSig")
Reac_HF$permute_sig <- if_else(as.numeric(Reac_HF$permute.P) < 0.05, "Sig", "NotSig")
Reac_LM$permute_sig <- if_else(as.numeric(Reac_LM$permute.P) < 0.05, "Sig", "NotSig")
Reac_ALL$permute_sig <- if_else(as.numeric(Reac_ALL$permute.P) < 0.05, "Sig", "NotSig")

head(Reac_CH)

Reac_CH <- Reac_CH %>% dplyr::select(9,10,11,3,8,7,12,13)
Reac_HF <- Reac_HF %>% dplyr::select(9,10,11,3,8,7,12,13)
Reac_LM <- Reac_LM %>% dplyr::select(9,10,11,3,8,7,12,13)
Reac_ALL <- Reac_ALL %>% dplyr::select(9,10,11,3,8,7,12,13)
head(Reac_HF)
head(Reac_LM)
Reac_CH$chr <- gsub("^", "chr", Reac_CH$chr)
Reac_HF$chr <- gsub("^", "chr", Reac_HF$chr)
Reac_LM$chr <- gsub("^", "chr", Reac_LM$chr)
Reac_ALL$chr <- gsub("^", "chr", Reac_ALL$chr)
head(Reac_CH)
colnames(Reac_CH) <- c("chr", "start", "end", "P", "gene", "group", "nominal", "permute")
colnames(Reac_HF) <- c("chr", "start", "end", "P", "gene", "group", "nominal", "permute")
colnames(Reac_LM) <- c("chr", "start", "end", "P", "gene", "group", "nominal", "permute")
colnames(Reac_ALL) <- c("chr", "start", "end", "P", "gene", "group", "nominal", "permute")
head(Reac_CH)
head(Reac_LM)




col = colorRamp2(c(-1.3,100), c("steelblue", "darkred"))

head(Reac_CH)

Reac_CH$P <- -log10(as.numeric(Reac_CH$P))

Reac_CH$P <- Reac_CH$P * -1
head(Reac_CH)


# Charolais plotting


circos.clear()
pdf(file = args[6],   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 15)
col_text <- "grey40"
#circos.par("track.height"=0.3, gap.degree=3, start.degree = 90, cell.padding=c(0, 0, 0, 0))
circos.par("track.height" = 3, start.degree = 90, gap.degree = 5)
circos.initializeWithIdeogram(species = "bosTau9", chromosome.index = paste0("chr", c(1:29)))


circos.text(sector.index="chr1",track.index = 1,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
            get.cell.meta.data("cell.ylim")-max(get.cell.meta.data("cell.ylim"))/2, labels = "Charolais",facing = "clockwise", 
            niceFacing = TRUE, adj = c(0,0),cex = 0.5)
circos.genomicTrackPlotRegion(Reac_CH, stack = F, numeric.column = "P", ylim = range(Reac_CH$P), panel.fun = function(region, value, ...) {
    circos.lines(CELL_META$cell.xlim, c(cut_off_CH, cut_off_CH), lty = 2, col = "black")
    y = Reac_CH[chr == CELL_META$sector.index,]
    colours = c()
    for (row in 1:nrow(y)) {
        curr = y[row,]
        head(curr)
        colours = if(curr$P < cut_off_CH & curr$permute == "Sig") {
            colours = c(colours,"darkred")
        } else if (curr$P < cut_off_CH) {
        colours = c(colours,"darkorange")
        } else {
            colours = c(colours,"steelblue")
        }
    }
    circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = colours)
    colours = c()
}, bg.border=T, track.height=0.06)

ch_signif = Reac_CH %>% filter(P < cut_off_CH & permute == "Sig")
Reac_CH
circos.yaxis("right", sector.index="chr1", labels.cex = 0.4)
circos.genomicLabels(ch_signif, labels.column=5, cex=0.25, col="black", line_lwd=0, line_col="black", 
side="inside", connection_height=0.01, labels_height=0.002, padding = 0.0001)

# Limousin
circos.text(sector.index="chr1",track.index = 4,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
            get.cell.meta.data("cell.ylim")-max(get.cell.meta.data("cell.ylim"))/2, labels = "Limousin",facing = "clockwise", 
            niceFacing = TRUE, adj = c(0,0),cex = 0.5)
Reac_LM$P <- as.numeric(Reac_LM$P)
Reac_LM$P <- -log10(Reac_LM$P)
Reac_LM$P <- Reac_LM$P * -1
head(Reac_LM)

circos.genomicTrackPlotRegion(Reac_LM, stack = F, numeric.column = "P", ylim = range(Reac_LM$P), panel.fun = function(region, value, ...) {
    circos.lines(CELL_META$cell.xlim, c(cut_off_LM, cut_off_LM), lty = 2, col = "black")
    y = Reac_LM[chr == CELL_META$sector.index,]
    colours = c()
    for (row in 1:nrow(y)) {
        curr = y[row,]
        head(curr)
        colours = if(curr$P < cut_off_LM & curr$permute == "Sig") {
            colours = c(colours,"darkred")
        } else if (curr$P < cut_off_LM) {
        colours = c(colours,"darkorange")
        } else {
            colours = c(colours,"steelblue")
        }
    }
    circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = colours)
    colours = c()
}, bg.border=T, track.height=0.06)
head(Reac_CH)
circos.yaxis("right", sector.index="chr1", labels.cex = 0.4)
lm_signif = Reac_LM %>% filter(P < cut_off_LM & permute == "Sig")
circos.genomicLabels(lm_signif, labels.column=5, cex=0.25, col="black", line_lwd=0.5, line_col="black", 
side="inside", connection_height=0.01, labels_height=0.02)

# HOFR
circos.text(sector.index="chr1",track.index = 7,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
            get.cell.meta.data("cell.ylim")-max(get.cell.meta.data("cell.ylim"))/2, labels = "Holstein-Friesian",facing = "clockwise", 
            niceFacing = TRUE, adj = c(0,0),cex = 0.5)
head(Reac_LM)
Reac_HF$P <- as.numeric(Reac_HF$P)
Reac_HF$P <- -log10(Reac_HF$P)
Reac_HF$P <- Reac_HF$P * -1
head(Reac_HF)

circos.genomicTrackPlotRegion(Reac_HF, stack = F, numeric.column = "P", ylim = range(Reac_HF$P), panel.fun = function(region, value, ...) {
    circos.lines(CELL_META$cell.xlim, c(cut_off_HF, cut_off_HF), lty = 2, col = "black")
    y = Reac_HF[chr == CELL_META$sector.index,]
    colours = c()
    for (row in 1:nrow(y)) {
        curr = y[row,]
        head(curr)
        colours = if(curr$P < cut_off_HF & curr$permute == "Sig") {
            colours = c(colours,"darkred")
        } else if (curr$P < cut_off_HF) {
        colours = c(colours,"darkorange")
        } else {
            colours = c(colours,"steelblue")
        }
    }
    circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = colours)
    colours = c()
}, bg.border=T, track.height=0.06)
hf_signif = Reac_HF %>% filter(P < cut_off_HF & permute == "Sig")
circos.yaxis("right", sector.index="chr1", labels.cex = 0.4)
circos.genomicLabels(hf_signif, labels.column=5, cex=0.25, col="black", line_lwd=0.5, line_col="black", 
side="inside", connection_height=0.01, labels_height=0.02)


# ALL
head(Reac_LM)
circos.text(sector.index="chr1",track.index = 10,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
            get.cell.meta.data("cell.ylim")-max(get.cell.meta.data("cell.ylim"))/2, labels = "Multi-breed",facing = "clockwise", 
            niceFacing = TRUE, adj = c(0,0),cex = 0.5)
Reac_ALL$P <- as.numeric(Reac_ALL$P)
Reac_ALL$P <- -log10(Reac_ALL$P)
Reac_ALL$P <- Reac_ALL$P * -1
head(Reac_LM)

circos.genomicTrackPlotRegion(Reac_ALL, stack = F, numeric.column = "P", ylim = range(Reac_ALL$P), panel.fun = function(region, value, ...) {
    circos.lines(CELL_META$cell.xlim, c(cut_off_ALL, cut_off_ALL), lty = 2, col = "black")
    y = Reac_ALL[chr == CELL_META$sector.index,]
    colours = c()
    for (row in 1:nrow(y)) {
        curr = y[row,]
        head(curr)
        colours = if(curr$P < cut_off_ALL & curr$permute == "Sig") {
            colours = c(colours,"darkred")
        } else if (curr$P < cut_off_ALL) {
        colours = c(colours,"darkorange")
        } else {
            colours = c(colours,"steelblue")
        }
    }
    circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = colours)
    colours = c()
}, bg.border=T, track.height=0.06)
all_signif = Reac_ALL %>% filter(P < cut_off_ALL & permute == "Sig")
circos.yaxis("right", sector.index="chr1", labels.cex = 0.4)
circos.genomicLabels(all_signif, labels.column=5, cex=0.25, col="black", line_lwd=0.5, line_col="black", 
side="inside", connection_height=0.01, labels_height=0.02)
dev.off()


