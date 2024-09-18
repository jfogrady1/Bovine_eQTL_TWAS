
# Plot number of conditional cis-eGenes
library(UpSetR)
library(tidyverse)
library(ggridges)
library(ggplot2)


args = commandArgs(trailing = T)

########################################################################################
########################################################################################
## 1. Upset plot of cis-eGenes
########################################################################################
########################################################################################

# Now onto upset plot of eGenes in reference groups
#args[1-3]
ALL_perm <- read.table(args[1], header = T) %>% 
filter(is_eGene == T) %>% select(7, 1,14)
CONTROL_perm <- read.table(args[2], header = T) %>% 
filter(is_eGene == T) %>% select(7, 1,14)
INFECTED_perm <- read.table(args[3], header = T) %>% 
filter(is_eGene == T) %>% select(7, 1,14)

listInput <- list(All = ALL_perm$phenotype_id, Control = CONTROL_perm$phenotype_id, Infected = INFECTED_perm$phenotype_id)

pdf(file=args[4], width = 12, height = 12, onefile=FALSE)
upset(fromList(listInput), order.by = "freq", sets.bar.color = c("#542788","#2166ac", "#b2182b"),
sets.x.label = "cis-eGenes", point.size = 4, line.size = 2,
mainbar.y.label = "cis-eGene intersections",
text.scale = 2.5, shade.alpha = 0.5)
dev.off()

# reinstall to get the eGenes (FDR < 0.05) significant at genome wide level
ALL_eGene <- read.table(args[1], header = T) %>% 
filter(is_eGene == T) %>% select(1)
CONTROL_eGene <- read.table(args[2], header = T) %>% 
filter(is_eGene == T) %>% select(1)
INFECTED_eGene <- read.table(args[3], header = T) %>% 
filter(is_eGene == T) %>% select(1)

dim(ALL_eGene)
########################################################################################
########################################################################################
## 2. Number of genes with number of eQTLs
########################################################################################
########################################################################################

#args[5-7]
# Now onto number of conditional eQTLs
ALL <- read.table(gzfile(args[5]), sep = "\t", header = T)
CONTROL <- read.table(gzfile(args[6]), sep = "\t", header = T)
INFECTED <- read.table(gzfile(args[7]), sep = "\t", header = T)


ALL <- inner_join(ALL, ALL_eGene)
CONTROL <- inner_join(CONTROL, CONTROL_eGene)
INFECTED <- inner_join(INFECTED, INFECTED_eGene)
# Percentage increase in eQTL numbers
all_increase <- ((nrow(ALL) - nrow(ALL_perm)) / nrow(ALL_perm)) * 100 
control_increase <- ((nrow(CONTROL) - nrow(CONTROL_perm)) / nrow(CONTROL_perm)) * 100 
infec_increase <- ((nrow(INFECTED) - nrow(INFECTED_perm)) / nrow(INFECTED_perm)) * 100 

ALL %>% filter(rank > 1) %>% dim()
CONTROL %>% filter(rank > 1) %>% dim()
INFECTED %>% filter(rank > 1) %>% dim()



all_increase
control_increase
infec_increase

# 78.83949
# 22.94136
# 5.631159


mean(ALL$rank)
mean(CONTROL$rank)
mean(INFECTED$rank)

# Now plot number per group
ALL <- ALL[order(ALL$rank, decreasing = T),]
ALL <- ALL %>% group_by(phenotype_id) %>% summarize(conditional = max(rank),
                                                    Group = "ALL")
ALL_plot <- ALL %>% group_by(conditional) %>% summarize(Group = "ALL",
                                                  genes = n())

CONTROL <- CONTROL[order(CONTROL$rank, decreasing = T),]
CONTROL <- CONTROL %>% group_by(phenotype_id) %>% summarize(conditional = max(rank),
                                                    Group = "CONTROL")
CONTROL_plot <- CONTROL %>% group_by(conditional) %>% summarize(Group = "CONTROL",
                                                  genes = n())


INFECTED <- INFECTED[order(INFECTED$rank, decreasing = T),]
INFECTED <- INFECTED %>% group_by(phenotype_id) %>% summarize(conditional = max(rank),
                                                    Group = "INFECTED")
INFECTED_plot <- INFECTED %>% group_by(conditional) %>% summarize(Group = "INFECTED",
                                                  genes = n())

conditional_plot <- rbind(ALL_plot, CONTROL_plot, INFECTED_plot)
conditional_plot$conditional <- as.character(conditional_plot$conditional)
conditional_plot$conditional <- factor(conditional_plot$conditional, levels = c("1","2","3","4","5","6","7","8","9","10","11"))



# Main plot
ggplot(data = conditional_plot, aes(x = conditional, y = genes, fill = Group)) + 
geom_bar(alpha = 0.8, stat="identity", position = "dodge") +
scale_fill_manual(breaks = c("ALL", "CONTROL", "INFECTED"), values = c("#542788","#2166ac", "#b2182b")) + 
scale_y_continuous(breaks = c(0,500,1000,1500,2000,2500,3000,3500,4000), limits = c(0,3500)) +
scale_x_discrete(breaks = c("1","2","3","4","5","6","7","8","9","10","11")) +
  theme_bw() +
  labs(y = "Genes with eQTLs", x = "# of conditional eQTLs") +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"))
ggsave(args[8], height = 12, width = 10, dpi = 600)

conditional_plot$filter <- as.numeric(conditional_plot$conditional)

# sub plot
conditional_plot %>% filter(conditional_plot$filter >= 4) %>% 
ggplot(., aes(x = conditional, y = genes, fill = Group)) + 
geom_bar(alpha = 0.8, stat="identity", position = "dodge") +
scale_fill_manual(labels = c("All", "Control", "Infected"), values = c("#542788","#2166ac", "#b2182b")) + 
scale_y_continuous(breaks = c(0,10,20,30,40,50,100,150,200,250,300,350), limits = c(0,360)) +
scale_x_discrete(breaks = c("4","5","6","7","8","9","10","11")) +
  theme_bw() +
  labs(y = "Genes with eQTLs", x = "# of conditional eQTLs") +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"))

ggsave(args[9], height = 12, width = 10, dpi = 600)



########################################################################################
########################################################################################
## 3. Are top eQTLs closer to TSS than conditional ones?
########################################################################################
########################################################################################


# Now onto estimating how close top and conditional eQTLs are to TSS
#args[5-7]
# Now onto number of conditional eQTLs
ALL <- read.table(gzfile(args[5]), sep = "\t", header = T)
CONTROL <- read.table(gzfile(args[6]), sep = "\t", header = T)
INFECTED <- read.table(gzfile(args[7]), sep = "\t", header = T)
ALL <- inner_join(ALL, ALL_eGene)
CONTROL <- inner_join(CONTROL, CONTROL_eGene)
INFECTED <- inner_join(INFECTED, INFECTED_eGene)

# ALL
ALL_distance_perm <- ALL %>% filter(rank == 1) %>% select(start_distance) %>% abs()
ALL_distance_conditional <- ALL %>% filter(rank > 1) %>% select(start_distance) %>% abs()
ALL_distance_perm$Category <- "Top"
ALL_distance_conditional$Category <- "Conditional"
ALL_distance_perm$Group <- "ALL"
ALL_distance_conditional$Group <- "ALL"
ALL_distance_plot <- rbind(ALL_distance_perm, ALL_distance_conditional)
ALL_distance_plot$start_distance <- ALL_distance_plot$start_distance / 1000

# CONTROL
CONTROL_distance_perm <- CONTROL %>% filter(rank == 1) %>% select(start_distance) %>% abs()
CONTROL_distance_conditional <- CONTROL %>% filter(rank > 1) %>% select(start_distance) %>% abs()
CONTROL_distance_perm$Category <- "Top"
CONTROL_distance_conditional$Category <- "Conditional"
CONTROL_distance_perm$Group <- "CONTROL"
CONTROL_distance_conditional$Group <- "CONTROL"
CONTROL_distance_plot <- rbind(CONTROL_distance_perm, CONTROL_distance_conditional)
CONTROL_distance_plot$start_distance <- CONTROL_distance_plot$start_distance / 1000



# INFECTED
INFECTED_distance_perm <- INFECTED %>% filter(rank == 1) %>% select(start_distance) %>% abs()
INFECTED_distance_conditional <- INFECTED %>% filter(rank > 1) %>% select(start_distance) %>% abs() # filter for ranks >1 which is top cis-eQTL
INFECTED_distance_perm$Category <- "Top"
INFECTED_distance_conditional$Category <- "Conditional"
INFECTED_distance_perm$Group <- "INFECTED"
INFECTED_distance_conditional$Group <- "INFECTED"
INFECTED_distance_plot <- rbind(INFECTED_distance_perm, INFECTED_distance_conditional)
INFECTED_distance_plot$start_distance <- INFECTED_distance_plot$start_distance / 1000

merged_distance_plot <- rbind(ALL_distance_plot, CONTROL_distance_plot, INFECTED_distance_plot)

ggplot(merged_distance_plot, aes(y=Category, x=start_distance, fill=Group)) + 
geom_density_ridges(scale = 0.8) +
  theme_bw() +
  labs(y = "Category of eQTL", x = "Distance to transcriptional start site (kb)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        strip.background = element_blank(), strip.text.x = element_blank()) + facet_wrap(~Group) +
        scale_fill_manual(labels = c("All", "Control", "Infected"), values = c("#542788","#2166ac", "#b2182b")) +
        scale_x_continuous(expand = c(0, 0), breaks = seq(0,1000,by = 250))

ggsave(args[10], height = 12, width = 10, dpi = 600)


# Actually test this to see if it is significantly different
all_result <- wilcox.test(abs(ALL_distance_perm$start_distance), abs(ALL_distance_conditional$start_distance))
infected_result <- wilcox.test(abs(INFECTED_distance_perm$start_distance), abs(INFECTED_distance_conditional$start_distance))
control_result <- wilcox.test(abs(CONTROL_distance_perm$start_distance), abs(CONTROL_distance_conditional$start_distance))

all_result$p.value
control_result$p.value
infected_result$p.value


#0
#4.00794e-75
#9.773674e-20

########################################################################################
########################################################################################
## 4. Is there a relationship between eQTL effect size and distance to TSS
########################################################################################
########################################################################################

# Here will use top and conditional eQTLs
# See if correlation in effect size estimates exists between close and far way()

ALL <- read.table(gzfile(args[5]), sep = "\t", header = T)
CONTROL <- read.table(gzfile(args[6]), sep = "\t", header = T)
INFECTED <- read.table(gzfile(args[7]), sep = "\t", header = T)
ALL <- inner_join(ALL, ALL_eGene)
CONTROL <- inner_join(CONTROL, CONTROL_eGene)
INFECTED <- inner_join(INFECTED, INFECTED_eGene)


# ALL
ALL_slope_distance <- ALL %>% select(start_distance, slope) 
ALL_slope_distance$start_distance <- abs(ALL_slope_distance$start_distance)
# CONTROL
CONTROL_slope_distance <- CONTROL %>% select(start_distance, slope)
CONTROL_slope_distance$start_distance <- abs(CONTROL_slope_distance$start_distance)

# INFECTED
INFECTED_slope_distance <- INFECTED %>% select(start_distance, slope) 
INFECTED_slope_distance$start_distance <- abs(INFECTED_slope_distance$start_distance)

ALL_slope_distance$Group <- "ALL"
CONTROL_slope_distance$Group <- "CONTROL"
INFECTED_slope_distance$Group <- "INFECTED"

final_point <- rbind(ALL_slope_distance, CONTROL_slope_distance, INFECTED_slope_distance)
final_point$slope <- abs(final_point$slope)
final_point$start_distance <- final_point$start_distance / 1000000

ggplot(data = final_point, aes(y = slope, x = start_distance, colour = Group)) + 
geom_point(alpha = 0.2) + scale_colour_discrete(labels = c("All", "Control", "Infected")) +
scale_colour_manual(labels = c("All", "Control", "Infected"), values = c("#542788","#2166ac", "#b2182b")) +
theme_bw() + labs(colour = "eQTL cohort", y = "Absolute effect size (slope)", x = "Distance to TSS (Mbs)") +
geom_smooth(method = "lm", se = FALSE, colour = "black") +
theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.position = "none", strip.background = element_blank(), strip.text.x = element_blank()) + facet_wrap(~Group)

ggsave(args[11], height = 12, width = 10, dpi = 600)

# Actually test this to see if correlation signficiant


result_ALL <- cor.test(ALL_slope_distance[,1], abs(ALL_slope_distance[,2]), method = "spearman", exact = FALSE)
result_CONTROL <- cor.test(CONTROL_slope_distance[,1], abs(CONTROL_slope_distance[,2]), method = "spearman", exact = FALSE)
result_INFECTED <- cor.test(INFECTED_slope_distance[,1], abs(INFECTED_slope_distance[,2]), method = "spearman", exact = FALSE)

result_ALL$p.value #5.781677e-235
result_CONTROL$p.value # 1.646611e-38
result_INFECTED$p.value # 1.734306e-09

result_ALL$estimate # -0.2938419  --> rho
result_CONTROL$estimate # -0.2019752 --> rho 
result_INFECTED$estimate # -0.1264643 -->rho