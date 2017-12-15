val2col <- function(x,gradientPalette=NULL,zlim=NULL,gradient.range.quantile=0.95) {
  if(all(sign(na.omit(x))>=0)) {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c('gray90','red'), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- as.numeric(quantile(x,p=c(1-gradient.range.quantile,gradient.range.quantile),na.rm=T))
      if(diff(zlim)==0) {
        zlim <- as.numeric(range(x))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])

  } else {
    if(is.null(gradientPalette)) {
      gradientPalette <- colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024)
    }
    if(is.null(zlim)) {
      zlim <- c(-1,1)*as.numeric(quantile(abs(x),p=gradient.range.quantile,na.rm=T))
      if(diff(zlim)==0) {
        zlim <- c(-1,1)*as.numeric(max(abs(x)))
      }
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    x <- (x-zlim[1])/(zlim[2]-zlim[1])

  }
  gradientPalette[x*(length(gradientPalette)-1)+1]
}


# note transpose is meant to speed up calculations when neither scaling nor centering is required
fast.pca <- function(m,nPcs=2,tol=1e-10,scale=F,center=F,transpose=F) {
  require(irlba)
  if(transpose) {
    if(center) { m <- m-Matrix::rowMeans(m)}; if(scale) { m <- m/sqrt(Matrix::rowSums(m*m)); }
    a <- irlba(tcrossprod(m)/(ncol(m)-1), nu=0, nv=nPcs,tol=tol);
    a$l <- t(t(a$v) %*% m)
  } else {
    if(scale||center) { m <- scale(m,scale=scale,center=center) }
    #a <- irlba((crossprod(m) - nrow(m) * tcrossprod(Matrix::colMeans(m)))/(nrow(m)-1), nu=0, nv=nPcs,tol=tol);
    a <- irlba(crossprod(m)/(nrow(m)-1), nu=0, nv=nPcs,tol=tol);
    a$l <- m %*% a$v
  }
  a
}

# quick utility to check if given character vector is colors
# thanks, Josh O'Brien: http://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
areColors <- function(x) {
  is.character(x) & sapply(x, function(X) {tryCatch(is.matrix(col2rgb(X)), error = function(e) FALSE)})
}

jw.disR <- function(x,y) { x <- x+1/length(x)/1e3; y <- y+1/length(y)/1e3; a <- x*log(x)  + y*log(y) - (x+y)*log((x+y)/2); sqrt(sum(a)/2)}

# translate multilevel segmentation into a dendrogram, with the lowest level of the dendrogram listing the cells
multi2dend <- function(cl,counts,deep=F) {
  if(deep) {
    clf <- as.integer(cl$memberships[1,]); # take the lowest level
  } else {
    clf <- as.integer(membership(cl));
  }
  names(clf) <- names(membership(cl))
  clf.size <- unlist(tapply(clf,factor(clf,levels=seq(1,max(clf))),length))
  rowFac <- rep(NA,nrow(counts));
  rowFac[match(names(clf),rownames(counts))] <- clf;
  lvec <- colSumByFac(counts,rowFac)[-1,];
  lvec.dist <- jsDist(t(lvec/pmax(1,Matrix::rowSums(lvec))));
  d <- as.dendrogram(hclust(as.dist(lvec.dist),method='ward.D'))
  # add cell info to the laves
  addinfo <- function(l,env) {
    v <- as.integer(mget("index",envir=env,ifnotfound=0)[[1]])+1;
    attr(l,'nodeId') <- v
    assign("index",v,envir=env)
    attr(l,'nCells') <- sum(clf.size[as.integer(unlist(l))]);
    if(is.leaf(l)) {
      attr(l,'cells') <- names(clf)[clf==attr(l,'label')];
    }
    attr(l,'root') <- FALSE;
    return(l);
  }
  d <- dendrapply(d,addinfo,env=environment())
  attr(d,'root') <- TRUE;
  d
}
# translate cell cluster dendrogram to an array, one row per node with 1/0 cluster membership
cldend2array <- function(d,cells=NULL) {
  if(is.null(cells)) { # figure out the total order of cells
    cells <- unlist(dendrapply(d,attr,'cells'))
  }
  getcellbin <- function(l) {
    if(is.leaf(l)) {
      vi <- match(attr(l,'cells'),cells)
      ra <- sparseMatrix(i=vi,p=c(0,length(vi)),x=rep(1,length(vi)),dims=c(length(cells),1),dimnames=list(NULL,attr(l,'nodeId')))
      return(ra);
    } else { # return rbind of the children arrays, plus your own
      ra <- do.call(cbind,lapply(l,getcellbin))
      ur <- unique(ra@i);
      ra <- cbind(sparseMatrix(ur+1,x=rep(1,length(ur)),p=c(0,length(ur)),dims=c(length(cells),1),dimnames=list(NULL,attr(l,'nodeId'))),ra);
      return(ra)
    }
  }
  a <- getcellbin(d)
  rownames(a) <- cells;
  return(t(a));
}


# peak annotation utils (based on Jean's code)

# hg38-specific utility function
# given a set of peaks ("chr:start-end", returns assignment to the ? gene using ChIPseeker)
hg38.peaks2Symbols <- function(peaks,includeDistal=FALSE, returnEntrezIds=FALSE,return.details=FALSE) {
  peak.df <- data.frame(do.call(rbind,strsplit(peaks,":|-")),stringsAsFactors=F); colnames(peak.df) <- c('chr','start','end');
  require(GenomicRanges)
  peaks.df <- with(peak.df, GRanges(chr, IRanges(as.numeric(start), as.numeric(end)), strand=NULL))

  ## Annotating peaks
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  require(ChIPseeker)
  spp.ann <- ChIPseeker::annotatePeak(peaks.df, TxDb=txdb, verbose = T)
  spp.ann.df <- as.data.frame(spp.ann)
  #rownames(spp.ann.df) <- paste(peaks)

  if(!includeDistal) {
    spp.ann.df$geneId[spp.ann.df$annotation=='Distal Intergenic'] <- NA
  }

  if(returnEntrezIds) {
    sym <- spp.ann.df$geneId;
  } else {
    ## get symbols
    library(org.Hs.eg.db)
    spp.ann.df$symbol <- sym <- select(org.Hs.eg.db,keys=spp.ann.df$geneId,columns=c('SYMBOL'))$SYMBOL
  }
  if(return.details) { return(spp.ann.df)}
  names(sym) <- peaks;
  sym
}

# group nearby peaks (returns a factor given peak names)
group.nearby.peaks <- function(peaks,window.size=1e3,verbose=T) {
  peaks.df <- data.frame(do.call(rbind,strsplit(peaks,":|-")),stringsAsFactors=F); colnames(peaks.df) <- c('chr','start','end');
  peaks.df$pos <- (as.numeric(peaks.df$start)+as.numeric(peaks.df$end))/2
  peak.f <- rep(NA,nrow(peaks.df));

  for(chr in unique(peaks.df$chr)) {
    ii <- which(peaks.df$chr==chr); ii <- ii[order(peaks.df$pos[ii],decreasing=F)];
    cf <- nearbyPointsGreedyCluster(peaks.df$pos[ii],window.size);
    peak.f[ii] <- cf+max(c(0,peak.f),na.rm=TRUE);
  }
  peak.f <- as.factor(peak.f);
  # translate into coordinate names
  peak.cfn <- unlist(tapply(1:nrow(peaks.df),peak.f,function(ii) paste(peaks.df$chr[ii[1]],':',min(peaks.df$pos[ii]),'-',max(peaks.df$pos[ii]),sep='')))
  peak.f <- as.factor(peak.cfn[peak.f]);
  names(peak.f) <- peaks;
  if(verbose) { cat(round(100-length(levels(peak.f))/length(peaks)*100,2),"% reduction\n")}
  return(peak.f);
}

hg38.closestGenes2Peaks <- function(peaks,n=10) {
  # parse out peak names
  peak.df <- data.frame(do.call(rbind,strsplit(peaks,":|-")),stringsAsFactors=F); colnames(peak.df) <- c('chr','start','end');
  peak.df$mid <- (as.numeric(peak.df$end)+as.numeric(peak.df$start))/2

  # get gene annotation
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(Homo.sapiens) # to look names homo sapiens can read
  # ... 20 packages and numberous overloaded functions later
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene;
  seqlevels(txdb) <- unique(peak.df$chr)
  genes <- genes(txdb)
  x <- tapply(1:nrow(peak.df),as.factor(peak.df$chr),function(ii) {
    ii <- ii[order(peak.df$mid[ii],decreasing=F)]
    chr <- peak.df$chr[ii[1]];
    gi <- seqnames(genes)==chr;
    gdf <- data.frame(start=start(genes[gi]),end=end(genes[gi]),entrezId=mcols(genes[gi])$gene_id,strand=as.integer(as.character(strand(genes[gi]))=="+"),stringsAsFactors=F)
    gdf <- gdf[order(gdf$start,decreasing=F),]
    x <- closestNSegmentsToPoints(gdf$start,gdf$end,peak.df$mid[ii],gdf$strand,n)
    x$i <- x$i+1; x$entrezId <- gdf$entrezId;
    rownames(x$i) <- rownames(x$d) <- rownames(x$s) <- peaks[ii];
    x
  })

  # combine joint index and distance matrices
  entrezIds <- unique(unlist(lapply(x,function(z) z$entrezId)))
  genes <- entrezIds;
  genes <- select(Homo.sapiens,keys=entrezIds,column=c("SYMBOL"),keytype="ENTREZID")$SYMBOL
  genes[is.na(genes)] <- entrezIds[is.na(genes)]

  # merge gene indices, translating them to the global index
  closestGenes <- do.call(rbind,lapply(x,function(z) {
    gids <- match(z$entrezId,entrezIds);
    gind <- apply(z$i,2,function(y) gids[y]);
    rownames(gind) <- rownames(z$i);
    return(gind)
  }))
  closestGenes <- closestGenes[peaks,]

  geneDistances <- do.call(rbind,lapply(x,function(z) {  z$d }))
  geneDistances <- geneDistances[peaks,]

  tssDistances <- do.call(rbind,lapply(x,function(z) {  z$s }))
  tssDistances <- tssDistances[peaks,]

  return(list(genes=genes,closestGenes=closestGenes,geneDistances=geneDistances,tssDistances=tssDistances))
}

hg38.closestPeaks2Genes <- function(peaks,n=10) {
  # parse out peak names
  peak.df <- data.frame(do.call(rbind,strsplit(peaks,":|-")),stringsAsFactors=F); colnames(peak.df) <- c('chr','start','end');
  peak.df$mid <- (as.numeric(peak.df$end)+as.numeric(peak.df$start))/2

  # get gene annotation
  require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  require(Homo.sapiens) # to look names homo sapiens can read
  # ... 20 packages and numberous overloaded functions later
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene;
  seqlevels(txdb) <- unique(peak.df$chr)
  genes <- genes(txdb)
  x <- tapply(1:nrow(peak.df),as.factor(peak.df$chr),function(ii) {
    ii <- ii[order(peak.df$mid[ii],decreasing=F)]
    chr <- peak.df$chr[ii[1]];
    gi <- seqnames(genes)==chr;
    gdf <- data.frame(start=start(genes[gi]),end=end(genes[gi]),entrezId=mcols(genes[gi])$gene_id,strand=as.integer(as.character(strand(genes[gi]))=="+"),stringsAsFactors=F)
    gdf <- gdf[order(gdf$start,decreasing=F),]
    x <- closestNPointsToSegments(gdf$start,gdf$end,peak.df$mid[ii],gdf$strand,n)
    x$i <- x$i+1;
    x$peaks <- peaks[ii];
    rownames(x$s) <- rownames(x$i) <- rownames(x$d) <- gdf$entrezId;
    x
  })

  # combine joint index and distance matrices
  entrezIds <- unique(unlist(lapply(x,function(z) rownames(z$i))))
  genes <- entrezIds;
  genes <- select(Homo.sapiens,keys=entrezIds,column=c("SYMBOL"),keytype="ENTREZID")$SYMBOL
  genes[is.na(genes)] <- entrezIds[is.na(genes)]

  # merge gene indices, translating them to the global index
  closestPeaks <- do.call(rbind,lapply(x,function(z) {
    gids <- match(z$peaks,peaks);
    gind <- apply(z$i,2,function(y) gids[y]);
    rownames(gind) <- rownames(z$i);
    return(gind)
  }))
  closestPeaks <- closestPeaks[entrezIds,]
  peakDistances <- do.call(rbind,lapply(x,function(z) {  z$d }))
  peakDistances <- peakDistances[entrezIds,]
  peakTssDistances <- do.call(rbind,lapply(x,function(z) {  z$s }))
  peakTssDistances <- peakTssDistances[entrezIds,]
  rownames(closestPeaks) <- rownames(peakDistances) <- rownames(peakTssDistances) <- genes;
  return(list(peaks=peaks,closestPeaks=closestPeaks,peakDistances=peakDistances,peakTssDistances=peakTssDistances))
}


# returns GO -> gene symbol map
hg38.getGo2Symbols <- function(symbols,entrezIds=FALSE,min.go.size=5,max.go.size=Inf,return.list=FALSE) {
  library(org.Hs.eg.db)
  # translate gene names to ids
  ids <- unlist(lapply(mget(unique(na.omit(symbols)),org.Hs.egALIAS2EG,ifnotfound=NA),function(x) x[1]))
  # reverse map
  rids <- names(ids); names(rids) <- ids;
  # list all the ids per GO category
  go.env <- eapply(org.Hs.egGO2ALLEGS,function(x) as.character(na.omit(rids[x])))
  sz <- unlist(lapply(go.env,length));
  go.env <- go.env[sz>=min.go.size & sz<=max.go.size];
  if(return.list) { go.env } else { list2env(go.env)}
}

# invert one-to-many string map (list)
invert.string.list <- function(map) {
  df <- cbind(rep(names(map),unlist(lapply(map,length))),unlist(map))
  tapply(df[,1],as.factor(df[,2]),I)
}

hg38.getSymbols2Go <- function(symbols,entrezIds=FALSE, return.list=FALSE, ...) {
  forward <- hg38.getGo2Symbols(na.omit(unique(symbols)),entrezIds=entrezIds,return.list = TRUE, ...)
  # reverse map
  gene2go <- invert.string.list(forward)
  if(return.list) { gene2go } else {return(list2env(gene2go))}
}

hg38.symbols2Peaks <- function( peaks, ... ) {
  symbols <- hg38.peaks2Symbols( peaks, ... );
  # build a reverse map
  tapply(peaks,as.factor(symbols),I)
}

#' Given mapping to peaks to symbols and symbols to gene sets, maps peaks to gene sets
peaks2GO <- function(peaks2Symbols, go.env, min.size=5, max.size=Inf) {
  go.env.peaks <- lapply(go.env, function(g) {
    na.omit(unique(unlist(symbols2Peaks[g])))
  })
  size <- unlist(lapply(go.env.peaks, length))
  go.env.peaks <- go.env.peaks[size > min.size & size < max.size]
  return(go.env.peaks)
}

hg38.peaks2GO <- function(peaks, includeDistal=F, ... ) {
  peaks2Symbols <- hg38.symbols2Peaks(peaks, includeDistal=includeDistal);
  library(liger)
  go.env <- org.Hs.GO2Symbol.list
  library(GO.db)
  desc <- select(GO.db, keys = names(go.env), columns = c("TERM"), multiVals = "CharacterList")
  stopifnot(all(names(go.env) == desc$GOID))
  names(go.env) <- paste(names(go.env), desc$TERM)
  go.env.peaks <- peaks2GO(peaks2Symbols, go.env, ...)
}




ths.collapse.aspect.clusters <- function(d,ct, scale = TRUE, pick.top = FALSE) {
  xvm <- do.call(rbind, tapply(seq_len(nrow(d)), factor(ct, levels = sort(unique(ct))), function(ii) {
    if(length(ii) == 1) return(d[ii, ])
    if(pick.top) {
      return(d[ii[which.max(apply(d[ii, ], 1, var))], ])
    }
    xp <- pcaMethods::pca(t(d[ii, ]), nPcs = 1, center = TRUE, scale = "none")
    xv <- pcaMethods::scores(xp)[, 1]
    if(sum(abs(diff(xv))) > 0 && cor(xv, Matrix::colMeans(d[ii, ]*abs(pcaMethods::loadings(xp)[, 1])))<0) { xv <- -1*xv }
    #set scale at top pathway?
    if(sum(abs(diff(xv))) > 0) {
      if(scale) {
        xv <- xv*sqrt(max(apply(d[ii, ], 1, var)))/sqrt(var(xv))
      }
      if(sum(abs(xv)) == 0) { xv <- abs(rnorm(length(xv), sd = 1e-6)) }
    } else {
      xv <- abs(rnorm(length(xv), sd = 1e-6))
    }
    #xv <- xv/sqrt(length(ii))
    xv
  }))
  rownames(xvm) <- unlist(tapply(seq_len(nrow(d)), factor(ct, levels = sort(unique(ct))), function(ii) {
    if(length(ii) == 1) return(rownames(d)[ii])
    return(rownames(d)[ii[which.max(apply(d[ii, ], 1, var))]])
  }))
  return(xvm);
}


draw.color.key <- function(x,y,w,h,col=colorRampPalette(c("green","black","red"),space="Lab")(100),n=300,labels=c("low","high"),lab="expression",labcol="white",labcex=1.2,vert=F,standalone=F,minor.ticks=NULL,useRaster=F,cex=1) {
  require(gridBase)
  require(grid)
  if(!standalone) {
    #opar <- par(no.readonly=TRUE)
    #grid.newpage();
    vp <- baseViewports()
    pushViewport(vp$inner,vp$figure,vp$plot)
    pushViewport(viewport(x=x,y=y,width=w,height=h,just=c("left","top")))
    par(plt=gridPLT(),new=T)
  }

  if(vert) {
    las=2; axd=4; srt=90;
    image(z=matrix(seq(-1,1,length=n),nrow=1),col=col,yaxt="n",xaxt="n",useRaster=useRaster);
  } else {
    las=0; axd=1; srt=0;
    image(z=matrix(seq(-1,1,length=n),ncol=1),col=col,yaxt="n",xaxt="n",useRaster=useRaster);
  }
  par(cex=cex)
  axis(axd,at=seq(0,(length(labels)-1))/(length(labels)-1),labels=labels,las=las);
  if(!is.null(minor.ticks)) {
    axis(axd,at=(minor.ticks-labels[1])/diff(range(labels)),labels=FALSE,tcl=-0.2)
  }
  text(0.5,0.1,labels=lab,col=labcol,adj=c(0.5,0.5),cex=labcex,srt=srt)
  if(!standalone) {
    popViewport(4)
    #par(opar);

  }
}

# a utility function to translate factor into colors
fac2col <- function(x,s=1,v=1,shuffle=FALSE,min.group.size=1,return.level.colors=F,unclassified.cell.color='gray50') {
  x <- as.factor(x);
  if(min.group.size>1) {
    x <- factor(x,exclude=levels(x)[unlist(tapply(rep(1,length(x)),x,length))<min.group.size])
  }
  col <- rainbow(length(levels(x)),s=s,v=v);
  if(shuffle) col <- sample(col);
  if(return.level.colors) { names(col) <- levels(x); return(col); }
  y <- col[as.integer(x)]; names(y) <- names(x);
  y[is.na(y)] <- unclassified.cell.color;
  y
}

