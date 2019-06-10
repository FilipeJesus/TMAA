assoc_alg <- function(dna,phen,verbose){
  nsnps=max(ncol(dna),1,na.rm=T); ntraits=max(ncol(phen),1,na.rm=T);
  if (nsnps == 1){
    nind = length(dna)
  } else {
    nind=nrow(dna)
  }
  all.can = as.numeric(stats::cancor(dna,phen)$cor)
  #all.can = as.numeric(kcca(dna,phen)$cor)
  all.can[all.can>0.99]=0.99
  # F-test
  fapp = function(nsnps,ntraits,nind,all.can)
  {
    wilks = prod(1-all.can^2)
    w = nind - 1 - 0.5*(nsnps+ntraits+1)
    t = sqrt( (nsnps^2 * ntraits^2 - 4) / (nsnps^2 + ntraits^2 - 5) )
    if (nsnps*ntraits == 2) t=1
    df1 = nsnps * ntraits
    df2 = w*t - 0.5*nsnps*ntraits + 1
    f = ( (1-wilks^(1/t)) / (wilks^(1/t)) ) * df2/df1
    p = pf(f,df1,df2,lower.tail=F)
    p
  }
  p = fapp(nsnps,ntraits,nind,all.can)
  p
}

# define fitness function 
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

monitor <- function(obj) {
      minEval = min(obj$evaluations, na.rm = TRUE);
      #plot(obj, type="hist");
      cat(sum(hamming.distance(obj$population)))
      cat("\n")
      cat(minEval)
      cat("\n")
    }