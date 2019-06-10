enrichment_test <- function(x,x.id,groups,population){
  results <- data.frame(matrix(ncol = 4),stringsAsFactors = FALSE)
  colnames(results) <- c("association","variable","probability","number seen")
  k <- length(x)
  for(i in unique(groups)){
    if(length(groups[groups==i]) > 2){
      m <- length(groups[groups==i])
      n <- length(groups) - m
      q <- length(groups[groups==i & population %in% x])
      results <- rbind(results, c(x.id,i,1-phyper(q=q,m=m,n=n,k=k),q),stringsAsFactors = FALSE)
    }
  }
  results$probability <- as.numeric(results$probability)
  results$p_adj <- p.adjust(results$probability, method = "bonferroni")
  return(results[-1,])
}