# load tidyverse package
library(tidyverse)
library(pophelper)

args = commandArgs(trailingOnly = TRUE)


###################### Read in PCA data


# read in data
#args[1] and #args[2]
args[1] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.eigenvec"
args[2] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.eigenval"
args[3] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/UMD_Genotype_pca_pve.pdf"
args[4] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/EQTL/Master_Sample_Sheet_textfile.txt"
args[5] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.2.Q"
args[7] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/Correlation_Admixture1_Holstein_pedigree.pdf"
args[7] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/UMD_Genotype_pca.pdf"
pca <- read_table(args[1], col_names = FALSE)
eigenval <- scan(args[2])

# sort out the pca data
# remove nuisance column = family ID
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
pca$ind <- gsub("_.*", "", pca$ind)
pca$ind
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_bw() +
theme(axis.title  = element_text(colour = "black"),
  axis.title.x = element_text(colour = "black", size = 18),
  axis.title.y = element_text(colour = "black", size = 18),
  axis.text.x = element_text(colour = "black", size = 15),
  legend.title = element_text(size = 15, color = "black", face = "bold"),
  legend.text = element_text(size = 15),
  axis.text.y = element_text(colour = "black", size = 15)) 
pve

ggsave(args[3], width = 15, height = 10, dpi = 600)

################## Read in pedigree data for plot
library(data.table)
Sample_database_2022_master <- fread(args[4], fill = T)

colnames(Sample_database_2022_master)

Sample_database_2022_master <- Sample_database_2022_master %>% dplyr::select(Sample_ID, HO_percent, FR_percent) %>% as.data.frame()
Sample_database_2022_master
pca$ind
Sample_database_2022_master <- Sample_database_2022_master[Sample_database_2022_master$Sample_ID %in% pca$ind, ]
dim(Sample_database_2022_master)

################################# Read in admixture data
file <- args[5]


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
head(rows)

alist <- pophelper::readQ(files = file)
alist
alist$SNP_Pruned.2.Q
attributes(alist)
attributes(alist[[1]])
rownames(alist[[1]]) <- rows$rows
rownames(alist$SNP_Pruned.2.Q) <- rows$rows

if(length(unique(sapply(alist,nrow)))==1) alist <- lapply(alist,"rownames<-",rows$rows)
lapply(alist, rownames)[1:2]
alist$SNP_Pruned.2.Q
attributes(alist[[1]])


Sample_database_2022_master$HO_percent
admix <- alist$SNP_Pruned.2.Q %>% dplyr::select(1)
admix$sample.id <- rownames(admix)
admix$sample.id
PCA_plot <- left_join(pca, Sample_database_2022_master, by = c("ind" = "Sample_ID"))
PCA_plot <- left_join(PCA_plot, admix, by = c("ind" = "sample.id"))
PCA_plot$Cluster1

correlation_plot <- PCA_plot %>% filter(!(HO_percent == ""))
dim(correlation_plot)
correlation_plot <- correlation_plot %>% filter(!(FR_percent == ""))
dim(correlation_plot)
correlation_plot$HF_percent <- (as.numeric(correlation_plot$HO_percent) + as.numeric(correlation_plot$FR_percent))
plot(correlation_plot$Cluster1, correlation_plot$HO_percent)
cor.test(as.numeric(correlation_plot$HO_percent), as.numeric(correlation_plot$Cluster1), method = "pearson")

ggplot(data = correlation_plot, aes(x = as.numeric(Cluster1), y = as.numeric(HO_percent))) + geom_point(size = 2, alpha = 0.5) + theme_bw() + theme(axis.title  = element_text(colour = "black"),
                                                                                  axis.title.x = element_text(colour = "black", size = 18),
                                                                                  axis.title.y = element_text(colour = "black", size = 18),
                                                                                  axis.text.x = element_text(colour = "black", size = 15),
                                                                                  legend.title = element_text(size = 15, color = "black", face = "bold"),
                                                                                  legend.text = element_text(size = 15),
                                                                                  axis.text.y = element_text(colour = "black", size = 15)) + geom_smooth(method = "lm", se = TRUE) +
                                                                                  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 100), col = "darkred") +
                                                                                  labs(x = "Admixture Component 1", y = "Reported Pedigree Holstein %")
ggsave(args[6], width = 15, height = 10, dpi = 600)
hist(correlation_plot$Cluster1)
PCA_plot$HO_percent[is.na(PCA_plot$HO_percent)] <- 0
PCA_plot <- cbind(PCA_plot, cat)
colnames(PCA_plot)[25] <- "Category"

PCA_plot <- PCA_plot %>% mutate(HOFR_bin = cut(HO_percent, right = F, breaks=c(0.000,1,20, 40, 60, 80, 100), labels = c("0-Missing", "1-20", "20-40", "40-60", "60-80", "80-100")))
PCA_plot <- PCA_plot %>% mutate(ADMIX_bin = cut((Cluster1 * 100), right = F, breaks=c(0.000,0.001,20, 40, 60, 80, 100), labels = c("0", "1-20", "20-40", "40-60", "60-80", "80-100")))
PCA_plot$HOFR_bin <- factor(PCA_plot$HOFR_bin, levels = c("0-Missing", "1-20", "20-40", "40-60", "60-80", "80-100"))
PCA_plot$ADMIX_bin <- factor(PCA_plot$ADMIX_bin, levels = c("0", "1-20", "20-40", "40-60", "60-80", "80-100"))

ggplot(data = PCA_plot, aes(x = PC1, y = PC2, colour = Cluster1, size = HOFR_bin, shape = Category)) +
  geom_point() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  labs(colour = "Admixture Component 1", shape = "Group", size = "Pedigree Holstein %") +
  scale_colour_gradientn(colours = c("darkred","red","orange", "grey","steelblue","blue","darkblue"),
                         values = c(1,0.8,0.6, 0.5, 0.4,0.2,0)) + theme_bw() + theme(axis.title  = element_text(colour = "black"),
                                                                                  axis.title.x = element_text(colour = "black", size = 18),
                                                                                  axis.title.y = element_text(colour = "black", size = 18),
                                                                                  axis.text.x = element_text(colour = "black", size = 15),
                                                                                  legend.title = element_text(size = 15, color = "black", face = "bold"),
                                                                                  legend.text = element_text(size = 15),
                                                                                  axis.text.y = element_text(colour = "black", size = 15)) +
  guides(shape = guide_legend(override.aes = list(size = 4)))
ggsave(args[7], width = 15, height = 10, dpi = 600)
