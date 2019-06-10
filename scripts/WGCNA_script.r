# Rscript WGCNA_script.r --matrix data.csv --power 8 --module_size 5 --o out

library(WGCNA)
library(dplyr)
library(flashClust)
library(reshape2)

args <- commandArgs( trailingOnly=TRUE )
args

if( "--matrix" %in% args ){
  matrix.idx <- grep("--matrix", args)
  matrix <- args[ matrix.idx + 1 ]
  matrix <- read.csv(matrix, header = TRUE)
  datExpr=matrix[,2:ncol(matrix)]
} else {
  stop("please enter the path to the matrix with prefix '--matrix'")
}

if( "--power" %in% args ){
  power.idx <- grep("--power", args)
  power <- args[ power.idx + 1 ]
  power <- as.numeric(power)
} else {
  stop("please enter the value for power with prefix '--power'")
}

if( "--module_size" %in% args ){
  size.idx <- grep("--module_size", args)
  size <- args[ size.idx + 1 ]
  size <- as.numeric(size)
} else {
  stop("please enter the value for size with prefix '--size'")
}

if( "--o" %in% args ){
  out.idx <- grep("--o", args)
  out <- args[ out.idx + 1 ]
} else {
  stop("please enter the path to out with prefix '--o'")
}

# Set up
set.seed(0)
par(mfrow = c(2,1))
enableWGCNAThreads()

# Calculate soft thresholds
sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,powerVector = c(c(1:10), seq(from = 12, to=30, by=2)),corFnc = cor,corOptions = list(use = 'p'),networkType = "signed hybrid")

#turn adjacency matrix into a topological overlap matrix (TOM) to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(datExpr,networkType = "signed hybrid", TOMType = "signed", power = power)
colnames(TOM) =rownames(TOM) =colnames(matrix[,2:ncol(matrix)])
dissTOM=1-TOM

#hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average")

#set the diagonal of the dissimilarity to NA 
diag(dissTOM) = NA

# find modules
dynamicColors = labels2colors(cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = size))

# create and save module colour data frame
module.df <- data.frame("variables" = colnames(TOM), "module" = dynamicColors)
write.csv(module.df, paste0(out,".csv"), row.names = FALSE)

tmp <- melt(dissTOM^power)

write.csv(tmp, paste0(out,"_dissTOM.csv"), row.names = FALSE)
# create a heatmap
svg(paste0(out,".svg"), width=12,height=9)
TOMplot(dissTOM^power, geneTree, as.character(dynamicColors))
dev.off()