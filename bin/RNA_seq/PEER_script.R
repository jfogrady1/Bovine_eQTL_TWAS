############################################################################
############### Peer analysis ##############################################
############################################################################

## This needs to be done in the conda rpeer_env environment created on command line
library(peer)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
expression <- read.csv(args[1], header=F, row.names = 1,sep = "\t", skip = 1) #1,4,7,5
covariates <- read.csv(args[2], header=F, row.names = 1, sep = "\t", skip = 1)
args1 <- args[1]
# Get output in correct format
out <- unlist(strsplit(args1, "_"))[1]
print(out)
# need to format covariates
if (out == "ALL") {
  params = 14
  
} else if (out == "CONTROL") {
  params = 5
} else {
  params = 10
}

print("Expression data and covariate data read in")

expression <- apply(expression, c(1,2), as.double)
covariates <- apply(covariates, c(1,2), as.double)


print("Dimensions of Transcriptomic data:") 
dim(expression)

print("Dimensions of Covariates data ")
dim(covariates)


print("Initialising the model")
model=PEER()


print("Setting the expression data to be inferred")
PEER_setPhenoMean(model, expression)

# add the mean estimation as default since it is recommended in the tutorial
# of peer. will return Nk+1 factors
# does not add anything
#PEER_setAdd_mean(model, TRUE)


print("Set number of hidden factors to identify") # n/4 as recommended
PEER_setNk(model, params)


print("Dimensions of the model:")
dim(PEER_getPhenoMean(model))


print("Introducing covariates into the model")
PEER_setCovariates(model, as.matrix(covariates))

print("Number of factors to be inferred:") 
PEER_getNk(model) # should include extra covaraite above


print("Setting paramaters")
PEER_setPriorAlpha(model,0.001,0.1)
PEER_setPriorEps(model,0.1,10) # all default gamma distributed (precision and noise)

print("Updating and running the model")
PEER_update(model)


factors = PEER_getX(model)

print("Dimensions of factors inferred")
dim(factors)

print("Dimensions of weights inferred")
weights = PEER_getW(model)

dim(weights)

print("Dimensions of precision inferred")
precision = PEER_getAlpha(model)
dim(precision)

print("Dimensions of residuals")
residuals = PEER_getResiduals(model)
dim(residuals)
colnames(residuals) <- colnames(expression)
rownames(residuals) <- rownames(expression)


jpeg(paste0(out,"_MODEL.jpg"))
PEER_plotModel(model)

jpeg(paste0(out,"_FACTORS.jpg"),)
plot(1.0 / precision, xlab="Factors", ylab="Factor relevance", main="")


write.csv(t(residuals), file = paste0(out,"_RESIDUALS.txt"), sep = "\t")
write.csv(t(factors), file = paste0(out,"_FACTORS.txt"), sep = "\t")
write.csv(t(weights), file = paste0(out,"_WEIGHTS.txt"), sep = "\t")
write.csv(t(precision),file = paste0(out,"_PRECISION.txt"), sep = "\t")
dev.off()

