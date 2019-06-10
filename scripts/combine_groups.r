library(dplyr)

args <- commandArgs( trailingOnly=TRUE )
args

if( "--biological" %in% args ){
  biological.idx <- grep("--biological", args)
  biological <- args[ biological.idx + 1 ]
  biological <- read.csv(biological, header = TRUE)
  biological <- biological[,grepl("PC.$", colnames(biological))]
} else {
  stop("please enter the path to the biological matrix with prefix '--biological'")
}

if( "--WGCNA" %in% args ){
  WGCNA.idx <- grep("--WGCNA", args)
  WGCNA <- args[ WGCNA.idx + 1 ]
  WGCNA <- read.csv(WGCNA, header = TRUE)
  WGCNA <- WGCNA[,grepl("PC.$", colnames(WGCNA))]
} else {
  stop("please enter the path to the WGCNA matrix with prefix '--WGCNA'")
}

if( "--o" %in% args ){
  out.idx <- grep("--o", args)
  out <- args[ out.idx + 1 ]
} else {
  stop("please enter the path to out with prefix '--o'")
}

biological <- strsplit(x = colnames(biological),split = "_", fixed = TRUE) %>% as.data.frame(.) %>% t(.)
biological <- sapply(unique(biological[,1]),function(x) return(c(x,length(biological[as.character(biological[,1])==x,1]),"biological"))) %>% t(.)
WGCNA <- strsplit(x = colnames(WGCNA),split = "_", fixed = TRUE) %>% as.data.frame(.) %>% t(.)
WGCNA <- sapply(unique(WGCNA[,1]),function(x) return(c(x,length(WGCNA[as.character(WGCNA[,1])==x,1]),"WGCNA"))) %>% t(.)

print(head(biological))
print(head(WGCNA))

joint <- rbind(biological,WGCNA)
colnames(joint) <- c("variable","PC number", "type")

write.csv(joint, paste0(out), row.names = FALSE)