# Rscript O2_PLS.r --response response.csv --predictor predictor.csv --o out

library("dplyr")
library("OmicsPLS")

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
cv <- crossval_o2m_adjR2(predictor[,2:ncol(predictor)], response[,2:ncol(response)], 1:10, 0:10, 0:10, nr_folds = 3, nr_cores = 6)
print(cv)

cv <- cv[cv$MSE == min(cv$MSE),]
fit0 =o2m(predictor[,2:ncol(predictor)], response[,2:ncol(response)], n = cv$n, nx = cv$nx, ny = cv$ny)

end.time <- Sys.time()
time.taken <- start.time - end.time

print(summary(fit0))


saveRDS(fit0, paste0(out, "_O2PLS.rds"))

print(time.taken)
