library(tidyverse)
library(data.table)
eQTL <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL.expr_tmm_inv.bed.gz") %>% select(1,2,4)
data <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_Supplementary.txt")
colnames(eQTL)[1] <- "chr"
head(data)
head(eQTL)

plotting <- left_join(data, eQTL, by = c("Gene" = "phenotype_id"))

plotting$chr <- as.numeric(plotting$chr)
plotting$P <- as.numeric(plotting$P)
plotting$permute.P <- as.numeric(plotting$permute.P)
plotting$permute.P <- if_else(is.na(plotting$permute.P), 1, plotting$permute.P)
plotting$P <- if_else(is.na(plotting$P), 1, plotting$P)
plotting <- transform(plotting, cum_position = ave(start, chr, FUN = cumsum))
(0.05/(15828/4))
nrow(plotting)
thresh = 0.05/(15828/4)
-log10(0.05)
snpsOfInterest <- plotting %>% filter(as.numeric(P) < thresh) %>% filter(permute.P < 0.05) %>% select(Symbol) %>% as.vector()
snpsOfInterest
library(ggrepel)
typeof(plotting$start)

# Prepare the dataset
don <- plotting %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(start)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(plotting, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, start) %>%
  mutate( BPcum=start+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(-log10(as.numeric(P))>(-log10(thresh)), "yes", "no")) %>%
  mutate( is_annotate=ifelse(Symbol %in% snpsOfInterest$Symbol & permute.P < 0.05, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
-log10(thresh)
# Make the plot
ggplot(don, aes(x=BPcum, y=-log10(as.numeric(P)))) +
    
    # Show all points
    geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis

    # Add highlighted points
    geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2.5) +

    # Add highlighted points
    geom_point(data=subset(don, is_annotate=="yes"), color="darkred", size=4) +
  
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=Symbol), size=2, max.overlaps = 20) +

    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) + facet_wrap(vars(GWAS), nrow = 4, scales = "free")

ggsave("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/TWAS/results/ALL_TWAS_Manhattan.pdf", width = 15, height = 12, dpi = 600)

