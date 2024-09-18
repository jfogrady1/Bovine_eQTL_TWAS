library("tidyverse")
library("data.table")

args = commandArgs(trailingOnly = TRUE)
data <- fread(args[1])

ggplot(data = data, aes(x = DST)) + geom_histogram(bins = 100) +
theme_bw() +
labs(x = "Identity by state (IBS)", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.position = "none") #+
       # scale_x_continuous(breaks = c(0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9), limits = c(0.6,1))
ggsave(args[2], width = 12, height = 12, dpi = 600)
