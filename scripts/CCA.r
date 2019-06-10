# Rscript CCA.r --response ../results/metab_pct_5__module_res.csv --predictor ../data/methyl_data.csv --type many_single --reps 5 --p 4 --o ../results/ --it 10

library("dplyr")
library("genalg")
library("doMC")
library("parallel")
source("../functions/my_CCA.r")
library("e1071")
source("../functions/genealg3.R")
library("foreach")
library("CCA")

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

if( "--type" %in% args ){
  type.idx <- grep("--type", args)
  type <- args[ type.idx + 1 ]
  if(!type %in% c('single_single', 'single_many', 'many_single', 'many_many')){
    stop("please enter the path to the type variable with prefix '--type', options are 'single_single', 'single_many', 'many_single', 'many_many'")
  }
} else {
  stop("please enter the path to the type variable with prefix '--type', options are 'single_single', 'single_many', 'many_single', 'many_many'")
}

if( "--reps" %in% args ){
  reps.idx <- grep("--reps", args)
  reps <- args[ reps.idx + 1 ]
  reps <- as.numeric(reps)
} else {
  stop("please enter the path to the reps variable with prefix '--reps'")
}

if( "--it" %in% args ){
  it.idx <- grep("--it", args)
  it <- args[ it.idx + 1 ]
  it <- as.numeric(it)
} else {
  stop("please enter the path to the iteration variable with prefix '--it'")
}

if( "--p" %in% args ){
  p.idx <- grep("--p", args)
  p <- args[ p.idx + 1 ]
  p <- as.numeric(it)
} else {
  stop("please enter the path to the cores variable with prefix '--p'")
}

if( "--o" %in% args ){
  out.idx <- grep("--o", args)
  out <- args[ out.idx + 1 ]
} else {
  stop("please enter the path to out with prefix '--o'")
}

registerDoMC()
options(cores=p)
getDoParWorkers()

l.predictor = ncol(predictor)-1
l.response = ncol(response)-1

predictor.names <- colnames(predictor)[2:ncol(predictor)]
response.names <- colnames(response)[2:ncol(response)]

response <- response[response[,1] %in% predictor[,1],]
predictor <- predictor[predictor[,1] %in% response[,1],]

rownames(response) <- response[,1]
rownames(predictor) <- predictor[,1]

predictor <- predictor[rownames(response),]

start.time <- Sys.time()

if(type == "single_single"){
  source("../functions/genealg3.R")
  source("../functions/gene-based.cca3.2.R")
  source("../functions/genealgTwoCromosom.R")
  results = foreach(i=1:l.predictor, .combine = rbind) %dopar%{
    dna  = as.matrix(predictor[,i+1])
    assoc = apply(response[,2:ncol(response)],2,function(x) gene.assoc(dna,as.matrix(x)))
  }
  colnames(results)=colnames(response)[2:ncol(response)]
  rownames(results)=colnames(predictor)[2:ncol(predictor)]
} else if(type == "single_many"){
  source("../functions/genealg3.R")
  source("../functions/gene-based.cca3.2.R")
  library("foreach")
  library("CCA")
  # iterate for each gene
  results = foreach(i=1:l.predictor) %dopar%{
    predictor.pruned = as.matrix(predictor[,i+1])
    predictor.name = colnames(predictor)[i+1]

    verbose = 0
    nind=ncol(response)
    nsnps=1

    # define fitness function (if not use forech, can be defined outside)
    # calculate the gene based association value for phenotype selection
    fitness3 = function(ids){
      if(all(ids==FALSE))
        return(30)
      else{
        ix = as.logical(ids)
        return(gene.assoc(predictor.pruned,response[,ix]))
      }
    }

    # calculate all single gene association for all the response.names
    pvalues = apply(response,2,function(f) gene.assoc(predictor.pruned,f))
    min.pvalue = min(pvalues)

    cat(paste(predictor.name," min uni pvalue: ",min.pvalue,"\n"))
    list2 = list()

    j=1 # index for repetitions (max rules extracted)
    z=1 # index for, in case of non-converge, not infinite loop
    maxIt = 20 # max value for z

    while (j<=reps && z<maxIt){
      rbga.results = genalg3(size=l.response,  zeroToOneRatio=50,iters=it,elitism= 20,
                            evalFunc=fitness3,popSize =100 , verbose=F)

      cat(paste(predictor.name," ",rbga.results$best[it]))
      cat("\n")
      if(rbga.results$best[it]<=min.pvalue){ #only use results higher than single max single point association
        list1 = list()
        list1[["pvalue"]]= rbga.results$best[it]
        list1[["responseInd"]] = rbga.results$population[which.min(rbga.results$evaluations),]
        list1[["responseNam"]] = colnames(response[2:ncol(response)])[as.logical(rbga.results$population[which.min(rbga.results$evaluations),])]
        #list1[["corY.Yscore"]]=cc(predictor.pruned,response[,as.logical(rbga.results$population[which.min(rbga.results$evaluations),])])$scores
        list2[[j]]=list1
        j=j+1
      }
      z = z+1
    }

    list2[["name"]]=predictor.name
    list2[["minpvalue"]]=min.pvalue
    list2

  }
} else if(type =="many_single"){
  library("doMC")
  library(genalg)
  library(flexclust)
  source("../functions/gene-based.cca3.2.R")
  source("../functions/genealg3.R")
  results = for(i in 1:l.response){
    response.name = colnames(response)[i+1]

    # define fitness function (if not use forech, can be defined outside)
    # calculate the gene based association value for gene selection
    fitness3 = function(ids){
      verbose = F
      if(all(ids==FALSE))
        return(30)
      else{
        ids = as.logical(ids)
        selpredictor = colnames(predictor[,2:ncol(predictor)])[ids]
        predictorIdx = which(colnames(predictor)%in%selpredictor)
          #ix = as.logical(ids)
        return(gene.assoc(as.matrix(predictor[,(predictorIdx)]),response[,(i+1)]))
      }
    }

    list2 = list()
    # population size for the gene selection approach
    populationSize = 600
    for(j in 1:reps){
      rbga.results = genalg3(size=l.predictor,  zeroToOneRatio=500,iters=it, elitism= 30,
                            evalFunc=fitness3,popSize =populationSize , verbose=F) #popsize=300
      cat(paste(response.name," ",rbga.results$best[it]))
      cat("\n")
      list1 = list()
      list1[["pvalue"]]= rbga.results$best[it]
      list1[["predictorInd"]] = rbga.results$population[which.min(rbga.results$evaluations),]
      list1[["predictorNames"]] = colnames(predictor[,as.logical(rbga.results$population[which.min(rbga.results$evaluations),])])

      list2[[j]]=list1
    }
    list2[["name"]]=response.name
    print(paste("rep",i,"done"))
    list2
  }
} else if(type == "many_many"){
  source("../functions/genealgTwoCromosom.R")
  # calculate the gene based association value for gene /phenotype selection
  results = list()
  require("flexclust")
  fitness3 = function(ids){
    dimX = l.predictor
    dimY = l.response
    verbose = F
    ids = as.logical(ids)
    idsX = ids[1:dimX]
    idsY = ids[(dimX+1):(dimX+dimY)]
    if(all(idsX==FALSE) | all(idsY==FALSE)){
      return(30)
    } else{
      selGenes = predictor.names[idsX]
      #ix = as.logical(ids)
      return(assoc_alg(as.matrix(predictor[,(selGenes)]),as.matrix(response[,(which(idsY))]),verbose))
    }
  }

  # population size for the gene/phenotype selection approach
  populationSize = 1000
  mutationChanceX = 1/l.predictor  #number of total genes dependent mutation chance for gene selection
  mutationChanceY= 1/l.response #number of total phenotype dependent mutation chance for phen selection

  for(j in 1:reps){
    rbga.results = genalg3(sizeX=l.predictor, sizeY=l.response,  zeroToOneRatioX=1000,zeroToOneRatioY=5 ,iters=it, elitism= 10, evalFunc=fitness3,popSize =populationSize , verbose=F)
    #cat(paste(phenName," ",rbga.results$best[it]))
    #cat("\n")
    list1 = list()
    list1[["pvalue"]]= rbga.results$best[it]
    bestInd = rbga.results$population[which.min(rbga.results$evaluations),]

    predictorInd = bestInd[1:l.predictor]
    list1[["predictorInd"]] = predictorInd
    list1[["geneNames"]] = predictor.names[as.logical(predictorInd)]

    responseInd = bestInd[(l.predictor+1):(l.predictor+l.response)]
    list1[["responseInd"]] = responseInd
    list1[["responseNames"]] = response.names[as.logical(responseInd)]
    results[[j]]=list1
  }
}

end.time <- Sys.time()
time.taken <- start.time-end.time

if(type == "single_single"){
  write.csv(results, paste0(out,"_single_single_association.csv"))
} else {
  saveRDS(results, paste0(out,"_", type, "_association.rds"))
}

print(time.taken)
