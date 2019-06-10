
args <- commandArgs( trailingOnly=TRUE )
args

if( "--file" %in% args ){
  file.idx <- grep("--file", args)
  file <- args[ file.idx + 1 ]
  file <- readRDS(file)
} else {
  stop("please enter the path to the association file with prefix '--file'")
}

if( "--o" %in% args ){
  out.idx <- grep("--o", args)
  out <- args[ out.idx + 1 ]
} else {
  stop("please enter the path to out with prefix '--o'")
}
node.color = c("red","purple","green")
length(file)
nodes <- data.frame("id"=NA, "lab"=NA,"axis"=NA,"radius"=NA, "size"=NA, "color"=NA)
for(x in 1:length(file)){
  if(nrow(nodes)==1){
    nodes[1,] <- c("id"=1,"lab"=file[[x]]$name,"axis"=2,"radius"=1,"size"=0.1,"color"=node.color[x])
  } else {
    nodes <- rbind("id"=nodes,"lab"=c(nrow(nodes)+1,file[[x]]$name,"axis"=2,"radius"=1,"size"=0.1,"color"=node.color[x]))
  }
  print(file[[x]])
  tmp <- file[[x]]$geneNames
  nodes <- rbind(nodes,data.frame("id"=(nrow(nodes)+1):(nrow(nodes)+length(tmp)), "lab"=tmp,"axis"=rep_len(3, length.out = length(tmp)),"radius"=(nrow(nodes[nodes$axis==3,])+1):(nrow(nodes[nodes$axis==3,])+length(tmp)), "size"=rep_len(0.1,length.out = length(tmp)), "color"=rep_len("blue",length.out = length(tmp))))
  tmp <- file[[x]]$responseNames
  nodes <- rbind(nodes,data.frame("id"=(nrow(nodes)+1):(nrow(nodes)+length(tmp)), "lab"=tmp,"axis"=rep_len(1, length.out = length(tmp)),"radius"=(nrow(nodes[nodes$axis==1,])+1):(nrow(nodes[nodes$axis==1,])+length(tmp)), "size"=rep_len(0.1,length.out = length(tmp)), "color"=rep_len("orange",length.out = length(tmp))))
}
nodes <- nodes[!duplicated(nodes$lab),]
nodes$id <- 1:nrow(nodes)

edges <- data.frame("id1"=NA,"id2"=NA,"weight"=NA,"color"=NA) 
for(x in 1:length(file)){
  tmp <- file[[x]]$geneNames
  for(i in tmp){
    node.id <- nodes[nodes$lab == i,"id"]
    association.id <- nodes[nodes$lab == file[[x]]$name,"id"]
    edges <- rbind(edges, data.frame("id1"=association.id,"id2"=node.id,"weight"=1,"color"="#9ECAE1"))
  }
  tmp <- file[[x]]$responseNames
  for(i in tmp){
    node.id <- nodes[nodes$lab == i,"id"]
    association.id <- nodes[nodes$lab == file[[x]]$name,"id"]
    edges <- rbind(edges, data.frame("id1"=node.id,"id2"= association.id,"weight"=1,"color"="#BAE4B3"))
  }
}

edges <- na.omit(edges)
nodes <- na.omit(nodes)

write.csv(edges, paste0(out, "_edges.csv"))
write.csv(nodes, paste0(out, "_nodes.csv"))