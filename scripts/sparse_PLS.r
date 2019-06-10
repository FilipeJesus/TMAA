# Rscript colinearity_script.r --matrix matrix.csv

library("dplyr")
library("spls")

args <- commandArgs( trailingOnly=TRUE )
args

if( "--response" %in% args ){
  response.idx <- grep("--response", args)
  response <- args[ response.idx + 1 ]
  response <- read.csv(response, header = TRUE)
  response <- response[order(response[,1]),]
} else {
  stop("please enter the path to the response matrix with prefix '--response'")
}

if( "--predictor" %in% args ){
  predictor.idx <- grep("--predictor", args)
  predictor <- args[ predictor.idx + 1 ]
  predictor <- read.csv(predictor, header = TRUE)
  predictor <- predictor[order(predictor[,1]),]
} else {
  stop("please enter the path to the predictor matrix with prefix '--predictor'")
}

if( "--o" %in% args ){
  out.idx <- grep("--o", args)
  out <- args[ out.idx + 1 ]
} else {
  stop("please enter the path to out with prefix '--o'")
}

response <- response[response[,1] %in% predictor[,1],]
predictor <- predictor[predictor[,1] %in% response[,1],]

rownames(response) <- response[,1]
rownames(predictor) <- predictor[,1]

predictor <- predictor[rownames(response),]

predictor[,2:ncol(predictor)] <- scale(x = predictor[,2:ncol(predictor)],center = TRUE,scale = FALSE)
response[,2:ncol(response)] <- scale(x = response[,2:ncol(response)],center = TRUE,scale = FALSE)

start.time <- Sys.time()

#cv <- cv.spls(x = predictor[,2:ncol(predictor)], y = response[,2:ncol(response)], eta = seq(0.85,0.97,0.03), fold=3, K = 5:10, plot.it=FALSE )

fit1 <- spls(x = predictor[,2:ncol(predictor)], y = response[,2:ncol(response)], eta = 0.85, K = 5)

end.time <- Sys.time()
time.taken <- start.time - end.time

saveRDS(fit1, paste0(out,"_sparse_PLS_matrix.rds"))

print(time.taken)
