# Rscript enrichment_PCA.r --matrix data.csv --module module.csv --metadata metadata.csv --variable group --id id --o out

library(dplyr)
source('../functions/enrichment.r')

args <- commandArgs( trailingOnly=TRUE )
args

if( "--matrix" %in% args ){
  matrix.idx <- grep("--matrix", args)
  matrix <- args[ matrix.idx + 1 ]
  matrix <- read.csv(matrix, header = TRUE)
} else {
  stop("please enter the path to the matrix with prefix '--matrix'")
}

if( "--module" %in% args ){
  module.idx <- grep("--module", args)
  module <- args[ module.idx + 1 ]
  module <- read.csv(module, header = TRUE)
} else {
  stop("please enter the path to the module with prefix '--module'")
}

if( "--metadata" %in% args ){
  metadata.idx <- grep("--metadata", args)
  metadata <- args[ metadata.idx + 1 ]
  metadata <- read.csv(metadata, header = TRUE)
} else {
  stop("please enter the path to the metadata with prefix '--metadata'")
}

if( "--variable" %in% args ){
  variable.idx <- grep("--variable", args)
  variable <- args[ variable.idx + 1 ]
} else {
  stop("please enter the path to the variable with prefix '--variable'")
}

if( "--id" %in% args ){
  id.idx <- grep("--id", args)
  id <- args[ id.idx + 1 ]
} else {
  stop("please enter the path to the id with prefix '--id'")
}

if( "--o" %in% args ){
  out.idx <- grep("--o", args)
  out <- args[ out.idx + 1 ]
} else {
  stop("please enter the path to out with prefix '--o'")
}

cat("check 1\n")
metadata <- metadata[metadata$Variable.Name %in% module$variables,]
dim(metadata)
dim(module)
modules_enrich <- data.frame("colour"=NA,"group"=NA,"variable"=NA,"num"=NA,"p_adj"=NA)
for(x in unique(module$module)){
  tmp = enrichment_test(x=module[module$module ==x,"variables"],x.id = paste0(x,"_Name"),groups = metadata[,variable],population = metadata[,1])
  modules_enrich <- rbind(modules_enrich,data.frame("colour"=x,"group"=variable,"variable"=tmp$variable,"num"=tmp$`number seen`,"p_adj"=tmp$p_adj))
}
modules_enrich <- modules_enrich[modules_enrich$colour != "grey",] %>% na.omit()
modules_enrich <- modules_enrich[modules_enrich$p_adj <= 0.05, ]

print(head(modules_enrich))
write.csv(modules_enrich, paste0(out,"_module_enrich.csv"), row.names = FALSE)

cat("check 2\n")

# PCA analysis
pc <- data.frame(matrix(nrow = nrow(matrix),ncol = 0))
pc_meta <- list()
for(x in unique(modules_enrich[modules_enrich$group == variable,"colour"])){
  tmp <- module[module$module == x,1]
  tmp <- prcomp(matrix[,as.character(tmp)], center = TRUE,scale. = TRUE)
  tmp_x <- tmp$x[,1:which(summary(tmp)$importance["Cumulative Proportion",] > 0.95)[1]]
  tmp_rotate <- tmp$rotation[,1:which(summary(tmp)$importance["Cumulative Proportion",] > 0.95)[1]]
  if(class(tmp_x) == "numeric"){
    tmp_x <- as.data.frame(tmp_x)
    colnames(tmp_x) <- paste0(x,"_PC1")
    tmp_rotate <- as.data.frame(tmp_rotate)
    colnames(tmp_rotate) <- paste0(x,"_PC1")
  } else {
    colnames(tmp_x) <- paste0(x,"_",colnames(tmp_x))
    colnames(tmp_rotate) <- paste0(x,"_",colnames(tmp_rotate))
  }
  pc <- cbind(pc,tmp_x)
  pc_meta <- append(list(tmp_rotate),pc_meta,0)
}
names(pc_meta) <- unique(modules_enrich[modules_enrich$threshold == 5 & modules_enrich$group == "Name","colour"])
saveRDS(pc_meta,paste0(out,"_module_res.rds"))

head(as.character(module[module$module=="grey","variables"]))

cat("check 3\n")

pc_f <- cbind("id" = matrix$id, pc[,grep("PC1",colnames(pc))])
pc_f <- cbind(pc_f, matrix[,as.character(module[module$module=="grey","variables"])])
write.csv(pc_f, paste0(out,"_module_PC1_res.csv"), row.names = FALSE)

cat("check 4\n")

pc <- cbind("id" = matrix$id, pc)
pc <- cbind(pc, matrix[,as.character(module[module$module=="grey","variables"])])
write.csv(pc, paste0(out,"_module_res.csv"), row.names = FALSE)
