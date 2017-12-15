library(methods)
myPagoda2 <- setRefClass(
  "myPagoda2",
  contains = "Pagoda2",
  inheritPackage = T,
  methods = list(
    carefullyNormalizedReduction=function(type='counts',use.odgenes=FALSE, n.odgenes=NULL, odgenes=NULL, scale=F,center=T, name='normalized',n.depth.slices=20,depthScale=1e4,sparse=FALSE) {
      
      if(type=='counts') {
        #x <- counts;
        x <- misc[['rawCounts']]
      } else {
        if(!type %in% names(reductions)) { stop("reduction ",type,' not found')}
        x <- reductions[[type]]
      }
      if((use.odgenes || !is.null(n.odgenes)) && is.null(odgenes)) {
        if(is.null(misc[['devsites']] )) { stop("please run adjustVariance() first")}
        odgenes <- misc[['devsites']];
        if(!is.null(n.odgenes)) {
          if(n.odgenes>length(odgenes)) {
            #warning("number of specified odgenes is higher than the number of the statistically significant sites, will take top ",n.odgenes,' sites')
            odgenes <- rownames(misc[['varinfo']])[(order(misc[['varinfo']]$lp,decreasing=F)[1:min(ncol(counts),n.odgenes)])]
          } else {
            odgenes <- odgenes[1:n.odgenes]
          }
        }
      }
      if(!is.null(odgenes)) { x <- x[,odgenes] }
      
      if(is.null(misc[['cpois']])) {
        warning("using naive guesses for library size/abundance. getRefinedLibSizes() is recommended")
        dep <- Matrix::rowSums(x); names(dep) <- rownames(x);
        if(is.null(batch)) {
          p <- Matrix::colSums(x)/sum(dep); p <- p/sum(p); names(p) <- colnames(x)
        } else {
          # in this case, p is a matrix with the number of rows equal to the number of batches
          p <- colSumByFac(x,as.integer(batch))[-1,];
          p <- p/Matrix::rowSums(p);
          p <- p/sum(p);
          colnames(p) <- colnames(x)
          rownames(p) <- levels(batch)
        }
        rp <- list(p=p,depth=dep);
      } else {
        rp <- misc[['cpois']]
      }
      
      if(sparse) {
        
        rx <- x
        rx.gene <- rep(1:rx@Dim[2],diff(rx@p))
        if(is.null(batch)) {
          # -log prob of seeing more than 0 observations for all the non-zero observations, given site specific p and cell-specific depth
          rx@x <- -1*ppois(rx@x-1,rp$depth[rownames(rx)[rx@i+1]]*rp$p[colnames(rx)][rx.gene],lower.tail=FALSE,log.p=TRUE)
        } else {
          ps <- rp$p[,colnames(rx)]; gbi <- cbind(as.integer(batch[rownames(rx)[rx@i+1]]), rx.gene); # rows/column index of batch/gene elements in p
          rx@x <- -1*ppois(rx@x-1,rp$depth[rownames(rx)[rx@i+1]]*ps[gbi],lower.tail=FALSE,log.p=TRUE)
          
          # even out the residuals between batches for each gene
          tc <- colSumByFac(rx,as.integer(batch))[-1,]; # total residuals per batch (row) per gene (column)
          batch.sizes <- as.numeric(tapply(batch,batch,length));
          dataset.average <- Matrix::colSums(tc)/sum(batch.sizes); # dataset-wide averages
          batch.scaling <- t(t(tc/batch.sizes)/dataset.average);
          # find discrepant genes
          mins <- apply(batch.scaling,2,min); maxs <- apply(batch.scaling,2,max);
          batch.scaling[,mins<0.5 | maxs>2] <- 1e10; # take out the genes.
          
          rx@x <- rx@x/batch.scaling[gbi];
        }
      } else {
        if(!is.null(batch)) { stop("batch mode is not yet supported in the dense calculation")}
        rx <- do.call(cbind,mclapply(colnames(x), function(i) {
          nzi <- x[,i]>0;
          lp <- as.numeric(rp$depth*rp$p[i])
          lp[nzi] <- ppois(0,rp$depth[nzi]*rp$p[i],lower.tail=FALSE,log.p=TRUE)
          return(-1*lp)
        },mc.preschedule=T,mc.cores=n.cores))
      }
      
      # sdepth <- sort(depth,decreasing=T)
      # rx <- do.call(cbind,mclapply(colnames(x), function(i) {
      #   nzi <- x[,i]>0;
      #   plp <- lp <- as.numeric(rp$depth*rp$p[i])
      #   lp[nzi] <- ppois(0,rp$depth[nzi]*rp$p[i],lower.tail=FALSE,log.p=TRUE)
      #   names(plp) <- rownames(x);
      #   # scale relative to a perfect model where the occurrences are all within the top cells
      #   tn <- names(sdepth)[1:sum(nzi)];
      #   plp[tn] <- ppois(0,sdepth[1:sum(nzi)]*rp$p[i],lower.tail=FALSE,log.p=TRUE)
      #   #lp <- -1*lp; plp <- -1*plp;
      #   tl <- (sum(plp[plp<0]))/sum(lp<0)
      #   lp[lp<0] <- pmin(0,lp[lp<0]-tl);
      #   tl <- (sum(plp[plp>0]))/sum(lp>0)
      #   lp[lp>0] <- pmax(0,lp[lp>0]-tl);
      #   return(-1*lp)
      # },mc.preschedule=T,mc.cores=n.cores))
      
      rownames(rx) <- rownames(x)
      colnames(rx) <- colnames(x)
      
      reductions[[name]] <<- rx;
      
      return(invisible(rx))
    },
    
    # determine subpopulation-specific genes
    getDifferentialGenes=function(type='counts',clusterType=NULL,groups=NULL,testType='Fisher',name='customClustering', z.threshold=3,correct.z.scores=TRUE, upregulated.only=FALSE,long.form=FALSE) {
      if(is.null(groups)) {
        # look up the clustering based on a specified type
        if(is.null(clusterType)) {
          # take the first one
          cols <- clusters[[type]][[1]]
        } else {
          cols <- clusters[[type]][[clusterType]]
          if(is.null(cols)) { stop("clustering ",clusterType," for type ", type," doesn't exist")}
        }
        cols <- as.factor(cols);
      } else {
        # use clusters information
        if(!all(rownames(counts) %in% names(groups))) { warning("provided cluster vector doesn't list groups for all of the cells")}
        cols <- as.factor(groups);
      }
      cat("running differential expression with ",length(levels(cols))," clusters ... ")
      # use offsets based on the base model
      
      if(testType=='Fisher') {
        cm <- misc[['rawCounts']];
      } else {
        cm <- counts;
      }
      if(!all(rownames(cm) %in% names(cols))) { warning("cluster vector doesn't specify groups for all of the cells, dropping missing cells from comparison")}
      # determine a subset of cells that's in the cols and cols[cell]!=NA
      valid.cells <- rownames(cm) %in% names(cols)[!is.na(cols)];
      if(!all(valid.cells)) {
        # take a subset of the count matrix
        cm <- cm[valid.cells,]
      }
      # reorder cols
      cols <- as.factor(cols[match(rownames(cm),names(cols))]);
      
      if(testType=='Fisher') {
        # Fisher's exact test
        # total number of gene observations
        pickedN <- Matrix::colSums(cm)
        # number of gene observations within the group
        pickedRed <- colSumByFac(cm,as.integer(cols))[-1,]
        # total depth within each group
        redN <- Matrix::rowSums(pickedRed);
        # total depth within other cells
        whiteN <- sum(redN)-redN;
        lpv <- do.call(rbind,lapply(1:nrow(pickedRed),function(i) {
          phyper(pickedRed[i,],redN[i],whiteN[i],pickedN+1,lower.tail=FALSE,log.p=TRUE)
        }))
        lpv[pickedRed==0] <- 0;
        lower.lpv.limit <- -100;
        lpv[lpv<lower.lpv.limit] <- lower.lpv.limit;
        lpv <- -1*lpv;
        if(!upregulated.only) {
          llpv <- do.call(rbind,lapply(1:nrow(pickedRed),function(i) {
            phyper(pickedRed[i,],redN[i],whiteN[i],pickedN,lower.tail=TRUE,log.p=TRUE)
          }))
          llpv[llpv<lower.lpv.limit] <- lower.lpv.limit;
          di <- llpv < -1*lpv;
          lpv[di] <- llpv[di];
        }
        x <- t(lpv);
      } else {
        # Wilcoxon rank test
        # calculate rank per-column (per-gene) average rank matrix
        xr <- sparse_matrix_column_ranks(cm);
        # calculate rank sums per group
        grs <- colSumByFac(xr,as.integer(cols))[-1,]
        # calculate number of non-zero entries per group
        xr@x <- numeric(length(xr@x))+1
        gnzz <- colSumByFac(xr,as.integer(cols))[-1,]
        #group.size <- as.numeric(tapply(cols,cols,length));
        group.size <- as.numeric(tapply(cols,cols,length))[1:nrow(gnzz)]; group.size[is.na(group.size)]<-0; # trailing empty levels are cut off by colSumByFac
        # add contribution of zero entries to the grs
        gnz <- (group.size-gnzz)
        # rank of a 0 entry for each gene
        zero.ranks <- (nrow(xr)-diff(xr@p)+1)/2 # number of total zero entries per gene
        ustat <- t((t(gnz)*zero.ranks)) + grs - group.size*(group.size+1)/2
        # standardize
        n1n2 <- group.size*(nrow(cm)-group.size);
        # usigma <- sqrt(n1n2*(nrow(cm)+1)/12) # without tie correction
        # correcting for 0 ties, of which there are plenty
        usigma <- sqrt(n1n2*(nrow(cm)+1)/12)
        usigma <- sqrt((nrow(cm) +1 - (gnz^3 - gnz)/(nrow(cm)*(nrow(cm)-1)))*n1n2/12)
        x <- t((ustat - n1n2/2)/usigma); # standardized U value- z score
      }
      
      
      ## # run quick fisher test on the count matrix for each group
      ## lower.lpv.limit <- -100;
      ## x <- do.call(cbind,tapply(1:nrow(cm),cols,function(ii) {
      ##   # total number of matching, non-matching rows
      ##   redN <- length(ii); whiteN <- length(cols)-redN
      ##   pickedN <- diff(cm@p)
      ##   pickedRed <- diff(cm[ii,]@p)
      ##   lpv <- -1*pmax(lower.lpv.limit,phyper(pickedRed,redN,whiteN,pickedN+1,lower.tail=FALSE,log.p=TRUE))
      ##   if(!upregulated.only) {
      ##     lpvl <- pmax(lower.lpv.limit,phyper(pickedRed,redN,whiteN,pickedN,lower.tail=TRUE,log.p=TRUE))
      ##     di <- lpvl < -1*lpv;
      ##     lpv[di] <- lpvl[di];
      ##   }
      ##   lpv
      ## }))
      
      # TODO: batch correction
      
      
      # correct for multiple hypothesis
      cat("adjusting p-values ... ")
      if(correct.z.scores) {
        x <- matrix(qnorm(scde:::bh.adjust(pnorm(as.numeric(abs(x)), lower.tail = FALSE, log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE),ncol=ncol(x))*sign(x)
      }
      rownames(x) <- colnames(cm); colnames(x) <- levels(cols)[1:ncol(x)];
      cat("done.\n")
      if(upregulated.only) {
        ds <- apply(x,2,function(z) {vi <- which(z>=z.threshold); r <- z[vi]; names(r) <- rownames(x)[vi]; sort(r,decreasing=T)})
      } else {
        ds <- apply(x,2,function(z) {vi <- which(abs(z)>=z.threshold); r <- z[vi]; names(r) <- rownames(x)[vi]; sort(r,decreasing=T)})
      }
      if(is.null(groups)) {
        if(is.null(clusterType)) {
          diffgenes[[type]][[names(clusters[[type]])[1]]] <<- ds;
        } else {
          diffgenes[[type]][[clusterType]] <<- ds;
        }
      } else {
        diffgenes[[type]][[name]] <<- ds;
      }
      return(invisible(ds))
    }
  )
)
