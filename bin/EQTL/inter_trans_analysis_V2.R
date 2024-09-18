########################
# Load in the libraries
########################

# Load in library
library(SNPRelate)
library(data.table)
library(tidyverse)
library(devtools)
library(vcfR)
library(ggplot2)
library(ggridges)
args = commandArgs(trailingOnly = T)
set.seed(42)


args[1] =  "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz"
args[2] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL.cis_qtl_fdr0.05.txt"
args[3] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_trans_FDR_corrected.txt"
args[4] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/inter_chrom_trans_cis_SNPs.txt"
args[5] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/inter_chrom_trans_cis_LD.ld"
args[6] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/inter_chrom_trans_cis_SNPs_NULL.txt"
args[7] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/inter_chrom_trans_cis_LD_NULL.ld"
args[8] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/inter_trans_V_mean_median_permuted.pdf"
args[9] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/inter_trans_V_5_permuted.pdf"
args[10] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/inter_trans_cis_permutation_raw.txt"
args[11] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/results/ALL_interchromosomal_LD_corrected.txt"
args[12] = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/software/circos/current/transqtl/data/inter_chrom_links.txt"


vcf = args[1]
cis = args[2]
trans = args[3]

vcf_data <- read.vcfR(vcf)


# TOP cis-eQTLs
cis = fread(cis)
cis = cis %>% filter(is_eGene == TRUE)
# Trans-eQTL results
trans = fread(trans, header = F) %>% filter(V10 < 0.05)
dim(trans)
trans$intra = if_else(trans$V2 == trans$V5, TRUE, FALSE)
trans <- trans %>% filter(intra == FALSE)
length(unique(trans$V1)) # 81 genes on different chromosomes
dim(trans)
#Not using this as we are using all significant SNPs in trans
trans <- trans %>% group_by(V1) %>% filter(V7 == min(V7))
trans <- trans[!duplicated(trans$V7),] # remove ties

# Filter cis results for common genes
cis <- cis %>% filter(phenotype_id %in% trans$V1)
trans <- trans %>% filter(V1 %in% cis$phenotype_id)

cis <- cis[order(cis$phenotype_id),]
trans <- trans[order(trans$V1),]
# get snps write to a file and calculate interchromosomal LD between these SNPs
snps <- c(cis$variant_id, trans$V4)
length(snps) # trans-eQTLs which are also cis-eQTLs
dim(trans) # 23 trans-eGenes that are cis-eGenes
write.table(snps, args[4], row.names = FALSE, col.names = FALSE, quote = FALSE,sep = "\t")
if (!file.exists(args[5])) {
system(paste0("plink --cow --vcf /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz ",
              "--r2 inter-chr --ld-window-r2 0 --extract /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/inter_chrom_trans_cis_SNPs.txt ",
              "--keep-allele-order --out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/inter_chrom_trans_cis_LD"))
}
# Read in the data
LD_data <- fread(args[5])
LD_data$pair <- paste0(LD_data$SNP_A,"-", LD_data$SNP_B)
LD_data$pair_reverse <- paste0(LD_data$SNP_B, "-",LD_data$SNP_A)
LD_data$same <- if_else(LD_data$SNP_A == LD_data$SNP_B, TRUE, FALSE)
LD_data$pair
LD_data$R2 <- as.numeric(LD_data$R2)
head(cis)
original_pairs <- cbind(cis$variant_id, trans$V4)
colnames(original_pairs) <- c("cis", "trans")
original_pairs <- as.data.frame(original_pairs)
original_pairs$pairs <- paste0(original_pairs$cis, "-", original_pairs$trans)

LD_data_pair1 <- LD_data %>% filter(pair %in% original_pairs$pairs) %>% select(7,8)
LD_data_pair2 <- LD_data %>% filter(pair_reverse %in% original_pairs$pairs) %>% select(7,9)
dim(cis)

#Now join both by pair
head(LD_data_pair2)
original_pairs <- left_join(original_pairs, LD_data_pair1,by = c("pairs" = "pair") )
original_pairs <- left_join(original_pairs, LD_data_pair2, by = c("pairs" = "pair_reverse"))

original_pairs <- original_pairs %>% mutate(R2_final = case_when(is.na(R2.x) ~ R2.y,
                                                                 is.na(R2.y) ~ R2.x)) %>% select(-c(R2.x, R2.y))

original_pairs$R <- sqrt(original_pairs$R2_final)
boxplot(original_pairs$R)
mean(original_pairs$R)

# Now we will sample the variants from the same chromosome and same allele frequency as putative trans-eQTLs


# Get VCF data ready
vcf_data <- vcf_data@fix %>% as.data.frame() %>% select(1,3,8)
vcf_data <- vcf_data %>% separate(., INFO, into = c("AF"), sep = ";") %>% select(1,2,3)
vcf_data$AF <- gsub("AF=", "", vcf_data$AF) %>% as.numeric()

null_snps <- as.data.frame(matrix(ncol = 4, nrow = 0))
colnames(null_snps) <- c("CHROM", "trans", "AF", "cis")
dim(trans)
for (i in 1:nrow(trans)) {
    trans_temp <- trans[i,]
    cis_temp <- cis[i,] # need this to keep track of pairs
    cis_var <- cis_temp %>% select(variant_id)
    vcf_data_temp <- vcf_data %>% as.data.frame() %>% filter(CHROM == trans_temp$V5) # filter for same chromosome as trans
    af <- vcf_data_temp %>% filter(vcf_data_temp$ID == trans_temp$V4) %>% select(AF) # get allele frequency
    vcf_data_temp <- vcf_data_temp %>% filter(as.numeric(AF) == as.numeric(af)) # filter for same allele frequency
    ran_snp <- sample_n(vcf_data_temp,1000, replace = T) 
    ran_snp <- as.data.frame(ran_snp)
    ran_snp$cis <- cis_var$variant_id
    null_snps <- rbind(null_snps,ran_snp)
}

null_snps$pairs <- paste0(null_snps$ID, "-", null_snps$cis)
null_snps_only <- null_snps$ID
null_snps_only <- c(null_snps_only, cis$cis_id)

# now we have SNPs which are null
# compute the LD between these and the cis eQTL SNPs
# Note, this file will be much larger
write.table(null_snps, args[6], row.names = F, col.names = F, quote = F,sep = "\t")

if (!file.exists(args[7])) {
system(paste0("plink --cow --vcf /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz ",
              "--r2 inter-chr --ld-window-r2 0  --keep-allele-order --extract /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/inter_chrom_trans_cis_SNPs_NULL.txt ",
              "--out /home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/transQTL/inter_chrom_trans_cis_LD_NULL"))
}
LD_null <- fread(args[7])

# filter for pairs which we want to compare
LD_null$pairs <- paste0(LD_null$SNP_A, "-", LD_null$SNP_B)
LD_null$pairs_reverse <- paste0(LD_null$SNP_B, "-", LD_null$SNP_A)
LD_null <- LD_null %>% filter(LD_null$pairs %in% null_snps$pairs | LD_null$pairs_reverse %in% null_snps$pairs)
dim(original_pairs)
# calculate position of variants
#LD_null$POS_trans <- paste0(LD_null$CHR1, ":", LD_null$POS1)
#LD_null$POS_cis <- paste0(LD_null$CHR2, ":", LD_null$POS2)
dim(LD_null)

# filter for associations of interest
LD_null <- LD_null %>% filter(pairs %in% null_snps$pairs | pairs_reverse %in% null_snps$pairs)



final_stats <- as.data.frame(matrix(ncol = 3, nrow = 0))
temp_stats <- as.data.frame(matrix(ncol = 3, nrow = 0))

colnames(final_stats) <- c("R2", "pairs","permutation")
colnames(temp_stats) <- c("R2", "pairs", "permutation")
# Now need to sample exactly the same number of variants
for (i in 1:10000) {
    print(i)
    temp_stats <- as.data.frame(matrix(ncol = 3, nrow = 0))
    colnames(temp_stats) <- c("R2", "pairs", "permutation")
    for(n in 1:length(cis$variant_id)) { # note this is the same order as trans
        trans_temp <- trans[n,] # select row
        trans_chr <- trans_temp$V5 # select chr
        cis_temp <- cis[n,] # select top cis eQTL
        cis_var <- cis_temp$variant_id

        # Get possible comparisons
        LD_null_temp <- LD_null %>% filter(SNP_A %in% cis_var | SNP_B %in% cis_var) 
        LD_null_temp <- LD_null_temp %>% filter(CHR_A %in% trans_chr | CHR_B %in% trans_chr) # null snps with same AF as trans-eQTL
        # Sample
        sample_R2_temp <- LD_null_temp[sample(nrow(LD_null_temp), 1),] 
        # Get columns of interest
        sample_R2_temp <- sample_R2_temp %>% select(R2, pairs)
        sample_R2_temp$permutation <- i
       temp_stats <- rbind(temp_stats, sample_R2_temp)
    }
    final_stats <- rbind(final_stats, temp_stats)
}
head(final_stats)
final_stats$R <- sqrt(final_stats$R2)
stats_df_Mean <- final_stats %>% group_by(permutation) %>% summarise(Value = "Mean", R = mean(R)) %>% select(Value, R)
stats_df_Median <- final_stats %>% group_by(permutation) %>% summarise(Value = "Median", R = median(R)) %>% select(Value, R)
head(original_pairs)
original_pairs_plotting<- original_pairs %>% select(R) %>% mutate(Value = "Observed", R = R) %>% select(Value, R) 
head(stats_df_Mean)
intra_plotting <- rbind(stats_df_Mean, stats_df_Median, original_pairs_plotting)
head(LD_trans_plotting)
library(viridis)
ggplot(intra_plotting, aes(x = Value, y = R, fill = Value)) + 
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    ## custom bandwidth
    adjust = .5, 
    ## adjust height
    width = .8, 
    ## move geom to the right
    justification = -.2, 
    ## remove slab interval
    .width = 0, 
    point_colour = NA
  ) + 
  geom_boxplot(
    width = .15) + 
    stat_boxplot(geom = "errorbar",width = 0.15) +
    stat_summary(fun.y=mean, geom="point", size = 5
    ) +
        scale_fill_viridis(discrete = TRUE, option = "D", alpha = 0.8) +
        theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12))
ggsave(args[8], width = 12, height = 12, dpi = 600)
# Test
num = sample(seq(1:10000), 10)
data_plot <- final_stats %>% filter(final_stats$permutation %in% num)
data_plot$Type = "null"
data_plot$Type = "null"
observed <- original_pairs %>% select(R2_final, pairs,R)
observed$permutation <- "observed"
observed <- observed %>% select(1,2,4,3)
colnames(observed)[1] <- "R2"
observed$Type = "Observed"
min_plot <- min_plot %>% select(4,5,3)
head(min_plot)
head(data_plot)
data_plot <- data_plot %>% select(4,5,3)
observed <- observed %>% select(4,5,3)
data_plot <- rbind(data_plot, observed)
data_plot$permutation <- as.factor(data_plot$permutation)
data_plot$Type <- as.factor(data_plot$Type)

ggplot(data = data_plot, aes(y = permutation, x = R, fill = Type)) + 
geom_density_ridges(scale = 1, jittered_points = TRUE, quantile_lines = T, alpha = 0.5) + theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12))
ggsave(args[9], width = 12, height = 12, dpi = 600)


head(final_stats)
write.table(final_stats, file = args[10], row.names = F, col.names = T, quote = F, sep ="\t")


original_pairs
# First generate links
original_pairs <- cbind(cis$phenotype_id, original_pairs)
links <- original_pairs
links
links$cis_var <- links$cis
links$trans_var <- links$trans
head(cis_var)
links <- links %>% separate(., cis_var, into = c("cis_chr", "cis_pos"), sep = ":")
links <- links %>% separate(., trans_var, into = c("trans_chr", "trans_pos"), sep = ":")
links <- links %>% select(2,7,8,3,9,10,6)
head(links)
# get colour
links <- links %>% mutate(linetype = case_when(R > 0 & R < .1 ~ paste0("color=vvdpyellow_a1,thickness=2p"), 
                          R >= .1 & R < .2 ~ paste0("color=vvdporange_a1,thickness=2p"),
                          R >= .2 & R < .4 ~ paste0("color=vvdpred_a1,thickness=2p"),
                          R > .4 ~ paste0("color=black_a1,thickness=2p")))
links$cis_pos2 <- as.numeric(links$cis_pos) + 1
links$trans_pos2 <- as.numeric(links$trans_pos) + 1

links<- links %>% select(cis_chr, cis_pos, cis_pos2,trans_chr, trans_pos,trans_pos2,linetype)
links$cis_chr <- gsub("^", "btau", links$cis_chr)
links$trans_chr <- gsub("^", "btau", links$trans_chr)
write.table(links, args[12], sep = " ", row.names = F, col.names = F, quote = F)


#Filter trans-eQTLs now
# Load in the data
trans = args[3]
trans = fread(trans, header = F) %>% filter(V10 < 0.05)
trans$intra = if_else(trans$V2 == trans$V5, TRUE, FALSE)
trans <- trans %>% filter(intra == FALSE)
length(unique(trans$V1))
trans_no_cis <- trans %>% filter(!(trans$V1 %in% original_pairs[,1]))  # keep these
length(unique(trans_no_cis$V1))

# Filter our set for genes
original_filtered <- original_pairs %>% filter(R^2 < 0.01)
trans_with_cis_LD_filtered <- trans %>% filter(V1 %in% original_filtered[,1])

# Now bind the two together
trans_LD_corrected <- rbind(trans_no_cis, trans_with_cis_LD_filtered)

# Now filter for FDR < 0.01
hist(trans_LD_corrected$V10)
trans_LD_corrected <- trans_LD_corrected %>% filter(V10 < 0.01)
unique(trans_LD_corrected$V1)
write.table(trans_LD_corrected, args[11], col.names = F, row.names = F, quote = F, sep = "\t")
