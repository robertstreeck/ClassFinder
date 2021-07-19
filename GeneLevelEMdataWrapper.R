EnrichCalc = function(ChIP, Inp){
  ChIP = ChIP/sum(ChIP)
  Inp = Inp/sum(Inp)
  ChIP/(ChIP + Inp)
}

GeneLevelClusteringFromFile = function(CountFile, path, ChIPcolumns = c(1,2), Inputcolumns = c(3,4), 
                               ncores = parallel::detectCores()/2, nmodel = 3){
  require(Rfast)
  source("/Users/streeck/Desktop/EM-Project/EMcalc.R")
  
  GeneCounts = read.delim(paste(path, CountFile, sep = ""), skip = 1)
  GeneCounts[,2:6] = NULL
  
  S = as.matrix(GeneCounts[,ChIPcolumns+1])
  R = as.matrix(GeneCounts[,Inputcolumns+1])
  N = R + S
  keep = rowsums(N>0) == length(ChIPcolumns)
  R = R[keep,]
  S = S[keep,]
  N = N[keep,]
  
  
  fit = BinomEMwrapperParallel(N, S, nmodel, ncores = ncores)
  
  fit$Input = R
  fit$InputFiles = colnames(GeneCounts)[Inputcolumns+1]
  fit$ChIP = S
  fit$ChIPFiles = colnames(GeneCounts)[ChIPcolumns+1]
  fit$GeneIDs = GeneCounts[,1]
  fit$Enrichment = EnrichCalc(S, R)
  
  return(fit)
}


GeneLevelClusteringFromDataFrame = function(counts, 
                                         ncores = parallel::detectCores()/2, nmodel = 3){
  require(Rfast)
  source("/Users/streeck/Desktop/EM-Project/EMcalc.R")
  
  ncol = dim(counts)[2]
  div = ncol/2
  S = as.matrix(counts[,1:div])
  R = as.matrix(counts[,(div+1):ncol])
  N = R + S
  keep = rowsums(N>0) == div
  R = R[keep,]
  S = S[keep,]
  N = N[keep,]
  
  
  fit = BinomEMwrapperParallel(N, S, nmodel, ncores = ncores)
  
  fit$Input = as.matrix(counts[,1:div])
  fit$InputFiles = colnames(counts)[(div+1):ncol]
  fit$ChIP = as.matrix(counts[,(div+1):ncol])
  fit$included = keep
  fit$ChIPFiles = colnames(counts)[1:div]
  fit$GeneIDs = rownames(counts)
  fit$Enrichment = EnrichCalc(S, R)
  
  return(fit)
}
