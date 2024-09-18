# Admixture plotting results
args = commandArgs(trailingOnly=TRUE)
#remotes::install_github('royfrancis/pophelper')
library(pophelper)
library(tidyverse)
# Args 1 is the Q file
file <- args[1]
head(data)
alist <- readQ(files = file)

attributes(alist[[1]])
cat <- c(rep("Control", 63), rep("Reactor", 60)) %>% as.data.frame()
colnames(cat)[1] <- "Group"
plotQ(alist,exportpath=args[2], imgtype="pdf",
      clustercol=c("steelblue","darkred"),
      showlegend=T, legendpos="right",
      legendlab = c("Holstein  ", "Non-Holstein"),
      showindlab=T,
      showyaxis=T,showticks=T,indlabcol="black",
      indlabangle=90,indlabvjust=1, basesize=4,
      legendkeysize = 5,
      height = 8, width = 17, dpi=600,
      splab = "K=2", splabsize = 8, grplab=cat, grplabpos = 0.5, grplabsize=2.5,linesize=0.8,pointsize=3,indlabwithgrplab=T,
      sortind="Cluster2", grplabheight = -4,
      outputfilename = args[3])

