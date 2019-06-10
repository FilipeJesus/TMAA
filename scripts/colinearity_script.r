# Rscript colinearity_script.r --matrix matrix.csv --o out.csv

library(dplyr)

args <- commandArgs( trailingOnly=TRUE )
args

if( "--matrix" %in% args ){
  matrix.idx <- grep("--matrix", args)
  matrix <- args[ matrix.idx + 1 ]
  matrix <- read.csv(matrix, header = TRUE)
} else {
  stop("please enter the path to the matrix matrix with prefix '--matrix'")
}

if( "--o" %in% args ){
  out.idx <- grep("--o", args)
  out <- args[ out.idx + 1 ]
} else {
  stop("please enter the path to out with prefix '--o'")
}

# correlation matric using spearmans rank
matrix_dist <- cor(matrix[2:ncol(matrix)], method = "spearman")
# remove matrix_dist metabolite which has a correlation of greater than 0.9 with any other metabolite
matrix_non_colin <- matrix[,2:ncol(matrix)]
tmp <- apply(matrix_dist,2,function(x) colnames(matrix_dist)[x > 0.9][1])
matrix_non_colin <- matrix_non_colin[,unique(tmp)]
# Build final dataframe and save
results <- cbind("id" = matrix[,1],matrix_non_colin)
write.csv(results, paste0(out), row.names = FALSE)