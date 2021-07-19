EnrichCalc = function(ChIP, Inp){
  require(Rfast)
  nrow = dim(ChIP)[1]
  ncol = dim(ChIP)[2]
  ChIP = ChIP/matrix(colsums(ChIP),nrow = nrow ,ncol = ncol,byrow = T)
  Inp = Inp/matrix(colsums(Inp),nrow = nrow ,ncol = ncol,byrow = T)
  ChIP/(ChIP + Inp)
}



BinEMClusterFromBAM = function(ChIPBamlist, InputBamList, gr, path, target, nmodel = 7, binsize = 200, saveBed = "none", 
                               ncores = parallel::detectCores()/2, paired.end = "midpoint", tlenFilter = c(70, 500), filteredFlag = 1024){
  # ChIPBamlist is a list of ChIP Bam files to be used for processing
  # InputBamList is the list of matching input bam files (repeat input files where necessary)
  # gr is the GeneomicRanges object of the chromosome sizes. e.g.
      #genome = read.delim("/Users/streeck/Genomes/DmelBDGP6.91/chrNameLength.txt", header = F, stringsAsFactors = F)
      #genome = genome[1:7,]
      #gr = GRanges(genome[,1], IRanges(1, as.integer(genome[,2])))
  # path is the path to the BAM file directory
  source("/Users/streeck/Desktop/EM-Project/EMcalc.R")
  require(bamsignals)
  
  ### read in Bam files
  print("Reading ChIP files")
  ChIPBamCounts = parallel::mcmapply(bamsignals::bamProfile, bampath=paste(path, ChIPBamlist, sep = ""), MoreArgs =
                                            list(gr=gr,binsize = binsize, paired.end=paired.end, tlenFilter = tlenFilter, 
                                                 filteredFlag = filteredFlag, verbose = F),
                                          mc.cores = ncores, SIMPLIFY = F)
  
  print("Reading Input files")
  InputBamCounts = parallel::mcmapply(bamsignals::bamProfile, bampath=paste(path, InputBamList, sep = ""), MoreArgs =
                                       list(gr=gr,binsize = binsize, paired.end=paired.end, tlenFilter = tlenFilter, 
                                            filteredFlag = filteredFlag, verbose = F),
                                     mc.cores = ncores, SIMPLIFY = F)
 
  print("Converting Frames") 
  ### number of ChIP sampels
  ndim = length(ChIPBamlist)

  ### convert bamsignals output to a list
  for(i in 1:ndim){ChIPBamCounts[[i]] = unlist(as.list(ChIPBamCounts[[i]]))}
  for(i in 1:ndim){InputBamCounts[[i]] = unlist(as.list(InputBamCounts[[i]]))}
  
  ### convert bamsignal lists to matrices
  ChIPBamMatrix = matrix(nrow = length(ChIPBamCounts[[1]]), ncol = ndim)
  for(i in 1:ndim){ChIPBamMatrix[,i] = ChIPBamCounts[[i]]}
  colnames(ChIPBamMatrix) = ChIPBamlist
  InputBamMatrix = matrix(nrow = length(InputBamCounts[[1]]), ncol = ndim)
  for(i in 1:ndim){InputBamMatrix[,i] = InputBamCounts[[i]]}
  colnames(InputBamMatrix) = InputBamList
  
  N = ChIPBamMatrix + InputBamMatrix
  
  NanEliminator = rowsums(N == 0) == 0
  
  ### call EM model
  print("Fitting model") 
  fit = BinomEMwrapperParallel(N[NanEliminator,], ChIPBamMatrix[NanEliminator,], k = nmodel, ncores = ncores, maxiter = 200)
  
  fit$Input = InputBamMatrix
  fit$InputFiles = InputBamList
  fit$ChIP = ChIPBamMatrix
  fit$ChIPFiles = ChIPBamlist
  fit$excluded = NanEliminator
  fit$Enrichment = EnrichCalc(ChIPBamMatrix[NanEliminator,], InputBamMatrix[NanEliminator,])
  fit$target = target
  
  reworkexport = unlist(tile(gr, width = binsize))
  reworkexport$name = "None"
  reworkexport$name[NanEliminator] = fit$Group
  
  fit$genome = reworkexport
  
  try(if(saveBed != "none"){
    export.bed(reworkexport, saveBed)
  })
  
  print("done")
  return(fit)
}

BinEMClusterFromBW = function(ChIPBamlist, InputBamList, gr, path, nmodel = 7, saveBed = F, binsize = F, ncores = parallel::detectCores()/2){
  # ChIPBamlist is a list of ChIP bw files to be used for processing
  # InputBamList is the list of matching input bw files (repeat input files where necessary)
  # gr is the GeneomicRanges object of the chromosome sizes. e.g.
      #genome = read.delim("/Users/streeck/Genomes/DmelBDGP6.91/chrNameLength.txt", header = F, stringsAsFactors = F)
      #genome = genome[1:7,]
      #gr = GRanges(genome[,1], IRanges(1, as.integer(genome[,2])))
  # path is the path to the bw file directory
  # saveBed is a path to a file for an output bed if desired
  
  require(rtracklayer)
  require(GenomicRanges)
  
  ### number of CHIP samples
  ndim = length(ChIPBamlist)
  
  ### read in bw files
  ChIPBamCounts = list()
  for (i in 1:ndim) {
    ChIPBamCounts[[i]] = import(paste(path, ChIPBamlist, sep = ""))
  }
  InputBamCounts = list()
  for (i in 1:ndim) {
    InputBamCounts[[i]] = import(paste(path, InputBamList, sep = ""))
  }
  
  ### determine binsize
  if(!as.logical(binsize)){
    binsize = round(median(ChIPBamCounts[1]@ranges@width))
  }
  
  ### convert gr genome into Tiled genome for detecting BW overlap
  TiledGenome = unlist(tile(gr, width = binsize))
  
  ### making fixed size matrixes out of bw files using tiled genome
  ChIPframe = matrix(ncol = ndim, nrow = length(TiledGenome))
  for (i in 1:ndim) {
    ChIPframe[,i] = ChIPBamCounts[[i]]$score[findOverlaps(TiledGenome, ChIPBamCounts[[i]], minoverlap = ceiling(binsize/2), select = "first")]
  }
  Inputframe = matrix(ncol = ndim, nrow = length(TiledGenome))
  for (i in 1:ndim) {
    Inputframe[,i] = InputBamCounts[[i]]$score[findOverlaps(TiledGenome, InputBamCounts[[i]], minoverlap = ceiling(binsize/2), select = "first")]
  }
  
  N = ChIPframe + Inputframe
  
  
  ### eliminate rows with 0
  NanEliminator = rowsums(N == 0) == 0
  
  ### fit model
  source("/Users/streeck/Desktop/EM-Project/EMcalc.R")
  fitBw = BinomEMwrapperParallel(N[NanEliminator,], ChIPframe[NanEliminator,], nmodel, ncores = ncores)
  
  ### write out resulting Bed files of clusters
  if(as.logical(saveBed)){
    reworkexport = TiledGenome
    reworkexport$name = "None"
    reworkexport$name[NanEliminator] = fitBw$Group
    export.bed(reworkexport, saveBed)
  }
  
  return(fitBw)
}
