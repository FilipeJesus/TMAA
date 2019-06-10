# Rscript diff_test.r --response response.csv --expression expression.csv --logic_variable adhd_16_case --id_variable id

library(dplyr)
library(mdscore)

args <- commandArgs( trailingOnly=TRUE )
args

if( "--response" %in% args ){
  response.idx <- grep("--response", args)
  response <- args[ response.idx + 1 ]
  response <- read.csv(response, header = TRUE)
} else {
  stop("please enter the path to the response matrix with prefix '--response'")
}

if( "--expression" %in% args ){
  expression.idx <- grep("--expression", args)
  expression <- args[ expression.idx + 1 ]
  expression <- read.csv(expression, header = TRUE)
} else {
  stop("please enter the path to the expression matrix with prefix '--expression'")
}

if( "--o" %in% args ){
  out.idx <- grep("--o", args)
  out <- args[ out.idx + 1 ]
} else {
  stop("please enter the path to out with prefix '--o'")
}

if( "--logic_variable" %in% args ){
  variable.idx <- grep("--logic_variable", args)
  variable <- args[ variable.idx + 1 ]
} else {
  stop("please enter the path to the variable matrix with prefix '--variable'")
}

if( "--id_variable" %in% args ){
  id.idx <- grep("--id_variable", args)
  id <- args[ id.idx + 1 ]
} else {
  stop("please enter the path to the id variable with prefix '--id_variable'")
}

cat("start\n")

# Clean data
response <- response[!is.na(response[,variable]),]

# Combine the two data frames
combined_df <- left_join(response[,c(id,variable)],expression, by = id, copy = TRUE)

cat("running logistic regression\n")
# run logistic regression and calculate p values
results <- apply(combined_df[,3:ncol(combined_df)], 2, function(x) {
  tmp <- glm(formula = case ~ data,data = data.frame("case"=combined_df[,variable], "data"=x),family = binomial(link = 'logit'), control = list(maxit = 50))
  tmp <- wald.test(model = tmp, terms = 2)
  return(c(tmp$W,tmp$pvalue))
})

# Clean results
results <- as.data.frame(t(results))
colnames(results) <- c("W","p_val")
results$id <- rownames(results)
results$p_adj <- p.adjust(results$p_val,method = "bonferroni")

cat("saving results\n")
# Create results file
write.csv(results, out, row.names = FALSE)
