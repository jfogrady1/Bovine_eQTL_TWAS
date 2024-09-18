library(data.table)
library(ggplot2)
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)
mbvFindBestMatch <- function(mbv_df){
  res = dplyr::transmute(mbv_df, mbv_genotype_id = SampleID,
                         het_consistent_frac = n_het_consistent/n_het_covered,
                         hom_consistent_frac = n_hom_consistent/n_hom_covered)

  #Identify best het
  best_het = dplyr::arrange(res, -het_consistent_frac) %>% dplyr::filter(dplyr::row_number() == 1)
  other_het = dplyr::arrange(res, -het_consistent_frac) %>% dplyr::filter(dplyr::row_number() > 1)
  best_row = dplyr::mutate(best_het, het_min_dist = min(best_het$het_consistent_frac - other_het$het_consistent_frac),
                           hom_min_dist = min(best_het$hom_consistent_frac - other_het$hom_consistent_frac),
                           distance = sqrt(het_min_dist^2 + hom_min_dist^2))

  #Compare against best hom
  best_hom = dplyr::arrange(res, -hom_consistent_frac) %>% dplyr::filter(dplyr::row_number() == 1)
  if(best_row$mbv_genotype_id != best_hom$mbv_genotype_id){
    best_row = dplyr::mutate(best_row, het_consistent_frac = as.numeric(NA), hom_consistent_frac = as.numeric(NA),
                             het_min_dist = as.numeric(NA), hom_min_dist = as.numeric(NA))
  }
  return(best_row)
}

string_vector <- c(
  'C001', 'C002', 'C003', 'C004', 'C005', 'C006', 'C007', 'C008', 'C009', 'C010',
  'C011', 'C012', 'C013', 'C014', 'C015', 'C016', 'C017', 'C018', 'C019', 'C020',
  'C021', 'C022', 'C023', 'C024', 'C025', 'C026', 'C027', 'C028', 'C029', 'C030',
  'C031', 'C033', 'C034', 'C035', 'C036', 'C037', 'C038', 'C039', 'C040', 'C041',
  'C042', 'C043', 'C044', 'C045', 'C046', 'C047', 'C048', 'C049', 'C050', 'C051',
  'C052', 'C053', 'C054', 'C055', 'C056', 'C057', 'C058', 'C059', 'C060', 'C061',
  'C062', 'C063', 'C064', 'T001', 'T002', 'T003', 'T004', 'T005', 'T006', 'T007',
  'T008', 'T009', 'T010', 'T011', 'T013', 'T014', 'T015', 'T016', 'T017', 'T018',
  'T019', 'T020', 'T022', 'T023', 'T024', 'T026', 'T028', 'T029', 'T030', 'T031',
  'T032', 'T033', 'T034', 'T035', 'T036', 'T037', 'T038', 'T039', 'T040', 'T041',
  'T042', 'T043', 'T044', 'T045', 'T046', 'T047', 'T048', 'T049', 'T050', 'T051',
  'T052', 'T053', 'T054', 'T055', 'T056', 'T057', 'T058', 'T059', 'T060', 'T061',
  'T062', 'T063', 'T064'
)

final_df_match <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(final_df_match) <- c("sample_cur", "mbv_genotype_id", "het_consistent_frac", "hom_consistent_frac", "het_min_dist", "hom_min_dist", "distance")


for (sample in string_vector) {
    sample_cur = as.character(sample)
    data <- fread(paste0("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/sample_check/",sample_cur,".mbv_output.txt"))
    best_match <- mbvFindBestMatch(data)
    data <- cbind(sample_cur, best_match)
    final_df_match <- rbind(final_df_match, data)
}


sample_cur = as.character("T064")
data <- fread(paste0("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/sample_check/",sample_cur,".mbv_output.txt"))
head(data)
colnames(data)
data <- data %>% select(1,9,10)
data$mbv_genotype_id <- "T064"
data$match <- "FALSE"
data$label <- NA
colnames(data)[1] <- "sample_cur"
data$distance <- NA
data$het_min_dist <- NA
data$hom_min_dist <- NA
colnames(data)
data <- data %>% select(1,4,2,3,8,9,7,5,6)
colnames(final_df_match)
colnames(final_df_match)
library(ggrepel)
data <- data %>% filter(mbv_genotype_id != sample_cur)
head(final_df_match)
all(final_df_match$sample_cur == final_df_match$mbv_genotype_id)
final_df_match$match = if_else(final_df_match$sample_cur == final_df_match$mbv_genotype_id, "True", "False")
final_df_match$label = if_else(final_df_match$het_consistent_frac < 0.85 | final_df_match$hom_consistent_frac < 0.92, final_df_match$sample_cur, NA)
colnames(data) <- colnames(final_df_match)
final_df <- rbind(final_df_match, data)
ggplot(data = final_df, aes(x = het_consistent_frac, hom_consistent_frac, colour = match)) + geom_point(size = 2) + theme_bw() +
scale_colour_manual(values = c("steelblue", "darkred")) + labs(x = "Concordance at HET", y = "Concordance at HOM", colour = "DNA/RNA\nmatch") +
ylim(0,1) + xlim(0,1) +
theme(axis.title  = element_text(colour = "black"),
  axis.title.x = element_text(colour = "black", size = 18),
  axis.title.y = element_text(colour = "black", size = 18),
  axis.text.x = element_text(colour = "black", size = 15),
  legend.title = element_text(size = 15, color = "black"),
  legend.text = element_text(size = 15),
  axis.text.y = element_text(colour = "black", size = 15)) 
ggsave(args[1], width = 12, height = 12, dpi = 600)
