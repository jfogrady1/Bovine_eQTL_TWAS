library(tidyverse)
library(ggpubr)
library(tidyquant)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
results <- read.csv("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/ALL_imputation_summary_statistics_sorted.info", sep = "\t", header = T, skip = 28)
colnames(results) <- c("SNP", "REF", "ALT","AF","MAF", "AvgCall", "R2", "Genotyped", "LooRsq", "ER", "ER2", "Dose0","Dose1")

head(results, 40)
# Calculating R2 at typed genetic varaints

mybins <- c(-0.000000001,0.00001*100,0.0001*100,0.0005*100,0.001*100,0.005*100,0.01*100,0.05*100,0.10*100,0.20*100,0.4*100,0.50*100)
test_ER <- results %>% filter(ER2 != "-")  %>% dplyr::select(ER2, MAF, R2) %>% as.data.frame()

head(results)

test_ER$R2 <- as.numeric(test_ER$R2)
test_ER$ER2 <- as.numeric(test_ER$ER2)
test_ER$MAF <- as.numeric(test_ER$MAF * 100)
head(test_ER)
test_ER$bins <- cut(test_ER$MAF,breaks = mybins)
dim(results)
#labels = c(0, '0.02','0.05','0.1', '0.2', '0.5', '1.0','2.0','5.0','10.0','15.0','20.0','30','40','50')

test_ER <- test_ER %>% filter(bins != "NA")
plotting <- test_ER %>% group_by(bins) %>% dplyr::summarise(ER2_mean = mean(ER2), R2_mean = mean(R2))
plotting$bins <- c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 20, 40, 50)

plotting <- as.data.frame(plotting)

head(plotting)

# Empirical R2
ggplot(plotting, aes(x = bins, y = ER2_mean)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  xlab("Minor allele frequency (%)") +
  ylab(bquote('Mean' ~ ER^"2")) +
  scale_y_continuous(limits=c(0,1),breaks = seq(0, 1.1, by = .1)) +
  coord_trans(x = "log10", xlim=c(0.001,50), expand = T, clip = "off") + annotation_logticks(scaled = F, outside = F, short = unit(0, "cm")) +
  scale_x_continuous(breaks = c(0.001,0.01,0.05,0.1,0.5,1,5,10,50), labels = c(0.001,0.01,0.05,0.1,0.5,1,5,10,50), limits = c(0.001, 50)) +
  theme_bw(base_size = 15) +
  theme(axis.title  = element_text(colour = "black"))
ggsave(args[2], width=10, height=8, dpi=600)



ggplot(plotting, aes(x = bins, y = R2_mean)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  xlab("Minor allele frequency (%)") +
  ylab(bquote('Mean Dosage' ~ R^"2")) +
  scale_y_continuous(limits=c(0,1),breaks = seq(0, 1.1, by = .1)) +
  coord_trans(x = "log10", xlim=c(0.001,50), expand = T, clip = "off") + annotation_logticks(scaled = F, outside = F, short = unit(0, "cm")) +
  scale_x_continuous(breaks = c(0.001,0.01,0.05,0.1,0.5,1,5,10,50), labels = c(0.001,0.01,0.05,0.1,0.5,1,5,10,50), limits = c(0.001, 50)) +
  theme_tq(base_size = 15) +
  theme(axis.title  = element_text(colour = "black"))

ggsave(args[3], width=10, height=8, dpi=600)

# ALL
# R2
# R2 for imputed variants
ALLtest_R <- results %>% dplyr::select(MAF, R2) %>% as.data.frame()
ALLtest_R$R2 <- as.numeric(ALLtest_R$R2)
ALLtest_R$MAF <- as.numeric(ALLtest_R$MAF * 100)
ALLtest_R$bins <- cut(ALLtest_R$MAF,breaks = mybins)

#labels = c(0, '0.02','0.05','0.1', '0.2', '0.5', '1.0','2.0','5.0','10.0','15.0','20.0','30','40','50')
ALLtest_R <- ALLtest_R %>% filter(bins != "NA")
plotting2 <- ALLtest_R %>% group_by(bins) %>% dplyr::summarise(R2_mean = mean(R2))
plotting2$bins <- c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 20, 40, 50)

ggplot(plotting2, aes(x = bins, y = R2_mean)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  xlab("Minor allele frequency (%)") +
  ylab(bquote('Mean Dosage' ~ R^"2")) +
  scale_y_continuous(limits=c(0,1),breaks = seq(0, 1, by = .1)) +
  coord_trans(x = "log10", xlim=c(0.001,50), expand = T, clip = "off") + annotation_logticks(scaled = F, outside = F, short = unit(0, "cm")) +
  scale_x_continuous(breaks = c(0.001,0.01,0.05,0.1,0.5,1,5,10,50), labels = c(0.001,0.01,0.05,0.1,0.5,1,5,10,50), limits = c(0.001, 50)) +
  theme_tq(base_size = 15) +
  theme(axis.title  = element_text(colour = "black"))

ggsave(args[4], width=10, height=8, dpi=600)


R2_All <- ggplot(plotting2, aes(x = bins, y = R2_mean)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  xlab("Minor allele frequency (%)") +
  ylab(bquote('Mean Dosage' ~ R^"2")) +
  scale_y_continuous(limits=c(0,1),breaks = seq(0, 1, by = .1)) +
  coord_trans(x = "log10", xlim=c(0.001,50), expand = T, clip = "off") + annotation_logticks(scaled = F, outside = F, short = unit(0, "cm")) +
  scale_x_continuous(breaks = c(0.001,0.01,0.05,0.1,0.5,1,5,10,50), labels = c(0.001,0.01,0.05,0.1,0.5,1,5,10,50), limits = c(0.001, 50)) +
  theme_bw(base_size = 15) +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15, angle = 45, vjust = 0.8),
        axis.text.y = element_text(colour = "black", size = 15))



ER2 <- ggplot(plotting, aes(x = bins, y = ER2_mean)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  xlab("Minor allele frequency (%)") +
  ylab(bquote('Mean' ~ ER^"2")) +
  scale_y_continuous(limits=c(0,1),breaks = seq(0, 1.1, by = .1)) +
  coord_trans(x = "log10", xlim=c(0.001,50), expand = T, clip = "off") + annotation_logticks(scaled = F, outside = F, short = unit(0, "cm")) +
  scale_x_continuous(breaks = c(0.001,0.01,0.05,0.1,0.5,1,5,10,50), labels = c(0.001,0.01,0.05,0.1,0.5,1,5,10,50), limits = c(0.001, 50)) +
  theme_bw(base_size = 15) +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15, angle = 45, vjust = 0.85),
        axis.text.y = element_text(colour = "black", size = 15))

ggarrange(R2_All, ER2, labels = "AUTO", ncol = 2, align = "v", font.label = list(size = 14, color = "black", face = "bold", family = NULL))
ggsave(filename = args[5], width=14, height=10, dpi=600)

