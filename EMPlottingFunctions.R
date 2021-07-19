MeanSampleIntensityHeatmap = function(fit, returnData = FALSE){
  require(plyr)
  require(reshape2)
  require(ggplot2)
  E = as.data.frame(fit$Enrichment)
  colnames(E) = fit$ChIPFiles
  E$cluster = fit$Group
  E = ddply(E, .(cluster), numcolwise(mean))
  E = melt(E, id.vars = "cluster")
  E$target = fit$target[match(E$variable, fit$ChIPFiles)]
  E = ddply(E, .(cluster, target), summarize, Enrich = mean(value))
  Eclust = acast(E, target ~ cluster)
  rclust = hclust(dist(Eclust))
  E$target = factor(E$target, levels = rclust$labels[rclust$order])
  cclust = hclust(dist(t(Eclust)))
  E$cluster = factor(E$cluster, levels = cclust$labels[cclust$order])
  if(returnData){
    return(E)
  }else{
    p = ggplot(E, aes(cluster, target)) + geom_tile(aes(fill = Enrich)) + theme_bw() +
      ylab("Target Mark") + xlab("Cluster") +
      scale_fill_distiller(palette = "Spectral", direction = 1, name = "Enrichment") +
      theme(legend.position = "top", axis.title.x = element_text(size = 15), axis.title.y = element_blank(),
            axis.text = element_text(size = 10, face = "bold")) +
      guides(fill = guide_colorbar(barwidth = 15, barheight = 1, direction = "horizontal", title.vjust = .83)) +
      geom_text(aes(label = round(Enrich, 2)), alpha = .5)
    return(p)
  }
}

BinCount = function(fit, plot.order = NA){
  require(ggplot2)
  if(is.na(plot.order[1])){
    plot.order = unique(fit$Group)
  }
  E = data.frame(Cluster = factor(fit$Group, levels = plot.order))
  p = ggplot(E, aes(Cluster)) + geom_bar(aes(y= 100*..count../sum(..count..))) +
    theme_bw() + ylab("% of Genome") + xlab("Cluster") + ggtitle("Fraction of genome covered by cluster") + 
    theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10, face = "bold"))
  return(p)
}

GeneCounts = function(fit, gr, gtffile, returnData = FALSE, plot.order = NA){
  require(GenomicRanges)
  require(rtracklayer)
  require(ggplot2)
  rework = unlist(tile(gr, width = 200))
  rework$name = "None"
  rework$name[fit$excluded] = fit$Group
  GenomeGFF = import.gff(gtffile)
  GenomeGFF = GenomeGFF[GenomeGFF$type == "exon"]
  count = c()
  for(i in 1:dim(fit$Expectaion)[2]){
    count = c(count, length(unique(subsetByOverlaps(GenomeGFF, rework[rework$name == as.character(i)], minoverlap = 1)$gene_id)))
  }
  count = data.frame(Cluster = as.character(1:dim(fit$Expectaion)[2]), ngenes = count)
  if(is.na(plot.order[1])){
    count$Cluster = factor(count$Cluster, levels = plot.order)
  }
  if(returnData){
    return(count)
  }else{
    p = ggplot(count, aes(Cluster, ngenes)) + geom_bar(stat = "identity") +
      theme_bw() + ylab("Number of Genes") + xlab("Cluster") + ggtitle("Number of gene overlapping each cluster") +
      theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10, face = "bold"))
    return(p)
  }
}

GeneOverlap = function(fit, gr, gtffile, returnData = FALSE, plot.order = NA){
  require(GenomicRanges)
  require(rtracklayer)
  require(ggplot2)
  rework = unlist(tile(gr, width = 200))
  rework$name = "None"
  rework$name[fit$excluded] = fit$Group
  GenomeGFF = import.gff(gtffile)
  GenomeGFF = GenomeGFF[GenomeGFF$type == "exon"]
  Glist = list()
  for(i in 1:dim(fit$Expectaion)[2]){
    Glist[[i]] = subsetByOverlaps(rework[rework$name == as.character(i)], GenomeGFF, minoverlap = 1)$name
  }
  Glist = data.frame(Cluster = as.character(1:dim(fit$Expectaion)[2]),
                     GeneCount = unlist(lapply(Glist, length))/summary(as.factor(fit$Group)))
  if(is.na(plot.order[1])){
    Glist$Cluster = factor(Glist$Cluster, levels = plot.order)
  }
  if(returnData){
    return(count)
  }else{
    p = ggplot(Glist,aes(Cluster, GeneCount)) + geom_bar(stat = "identity") +
      theme_bw() + ylab("Exon density") + xlab("Cluster") + ggtitle("Exon Density in each cluster") +
      theme(axis.title = element_text(size = 15), axis.text = element_text(size = 10, face = "bold"))
    return(p)
  }
}

BinEMPCA = function(fit, returnData = FALSE){
  require(ggplot2)
  BinPCA = stats::prcomp(t(fit$Enrichment))
  sam = sample(1:dim(BinPCA$rotation)[1], 30000)
  BinPCArot = data.frame(BinPCA$rotation[sam,], group = fit$Group[sam])
  BinPCAAxisLabs = round(100*BinPCA$sdev^2/sum(BinPCA$sdev^2),1)
  p = ggplot(BinPCArot, aes(PC1, PC2, color = group)) + geom_point(alpha = .2, shape = 16) + geom_density2d(bins = 5) +
    scale_color_brewer(palette = "Dark2", name = "Cluster") +
    xlab(paste("PC1 (", BinPCAAxisLabs[1], "% of variance explained)", sep = "")) +
    ylab(paste("PC2 (", BinPCAAxisLabs[2], "% of variance explained)")) + theme_bw()
  return(p)
}

quickwrite = function(fit, file, gtf, chrsize, chrSet = "all"){
  require(rtracklayer)
  genome = read.delim(chrsize, header = F, stringsAsFactors = F)
  if(chrSet[1] != "all"){
    genome = genome[chrSet,]
  }
  gr = GRanges(genome[,1], IRanges(1, as.integer(genome[,2])))
  a = MeanSampleIntensityHeatmap(fit)
  ord = levels(a$data$cluster)
  pdf(file)
  print(a)
  print(BinCount(fit, plot.order = ord))
  print(GeneCounts(fit, gr, gtf, plot.order = ord))
  print(GeneOverlap(fit, gr, gtf, plot.order = ord))
  print(BinEMPCA(fit))
  dev.off()
}

Fingerprint = function(fit, samplesize = 10000, return.data = F, enrichment = T){
  require(reshape2)
  require(plyr)
  require(ggplot2)
  require(Rfast)
  require(gridExtra)
  if(enrichment){
    sam = sample(dim(fit$Enrichment)[1], samplesize)
    E = fit$Enrichment[sam,]
  }else{
    sam = sample(dim(fit$ChIP)[1], samplesize)
    E = fit$ChIP[sam,]
  }
  E = apply(E, 2, sort)
  E = apply(E, 2, cumsum)
  E = E/matrix(E[samplesize,], nrow = dim(E)[1], ncol = dim(E)[2], byrow = T)
  E = data.frame(E, x = 1:samplesize)
  E = melt(E, id.vars = "x")
  E$target = as.factor(fit$target[as.numeric(E$variable)])
  E = E[E$x %in% seq(1, samplesize, length.out = 100),]
  E = ddply(E, .(target, x), summarize, val = mean(value))
  if(return.data){
    if(enrichment){
      return(E)
    }else{
      f = rowSums(fit$Input[sam,])
      f = cumsum(f[order(f, decreasing = F)])
      f = f/f[samplesize]
      f = data.frame(target = "Input", x = 1:10000, val = f)
      f = f[f$x %in% seq(1, samplesize, length.out = 100),]
      E = rbind(E, f)
      return(E)
    }
  }else{
    E$fact = factor(rep(1:ceiling(length(levels(E$target))/5), each = 5))[as.numeric(E$target)]
    plots = list()
    for(i in levels(E$fact)){
      E = rbind(E, data.frame(target = "Input", x = 1:10000, val = f, fact = i))
    }
    for(i in levels(E$fact)){
      plots[[as.numeric(i)]] = ggplot(E[E$fact == i,], aes(x/samplesize, val, color = target)) + geom_line(size = 1.5) + theme_bw() +
        guides(color = guide_legend(name = "IP target")) + xlab("Fraction of bins") + ylab("Cumulative fraction of signal") +
        scale_color_brewer(palette = "Dark2")
    }
    return(grid.arrange(grobs = plots)) 
  }
}


GeneLevelEmPlotting = function(fit){
  require(ggplot2)
  ndim = dim(fit$Input)[2]
  if(ndim>1){
    Plot = as.data.frame(fit$Enrichment)
    colnames(Plot) = paste("Rep_", 1:ndim, sep = "")
    Plot$Clust = fit$Group
    
    PlotA = ggplot(Plot, aes(Rep_1, Rep_2, color = Clust)) + geom_point(alpha = .3, shape = 16) + theme_bw() +
      xlab("H3K27me3 Enrichment (Replicate A)") + ylab("H3K27me3 Enrichment (Replicate B)") +
      scale_color_brewer(palette = "Dark2", name = "Cluster") + guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
      theme(legend.text = element_text(size = 10))
    PlotB = ggplot(Plot, aes(Rep_1, fill = Clust)) + geom_histogram(aes(y = 100*..count../sum(..count..)), bins = 100) + theme_bw() +
      xlab("H3K27me3 Enrichment (Replicate A)") + ylab("% of Genes") + scale_y_continuous(breaks = c(0,3,6), labels = c("0.00", "3.00", "6.00")) +
      scale_fill_brewer(palette = "Dark2", name = "Cluster") + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                                                                     plot.margin = unit(c(5.5,62.5,5.5,5.5), "points")) + guides(fill = F)
    
    return(gridExtra::grid.arrange(PlotB, PlotA, ncol = 1, heights = c(.25, .75)))
  }
}