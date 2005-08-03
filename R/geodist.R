piecewiselin <- function(distmatrix, maxdist=0.1*max(distmatrix)){
  ncd <- ncol(distmatrix)
  nmatrix <- distmatrix/maxdist  
  for (i in 1:(ncd-1))
    for (j in (i+1):ncd)
      if (nmatrix[i,j]>1)
        nmatrix[i,j] <- nmatrix[j,i] <- 1
  nmatrix
}
          
# Geographic Kulczynski distance, rows are regions
geco <- function(regmat,
                   geodist=as.dist(matrix(as.integer(!diag(nrow(regmat)))),
                             nrow=nrow(regmat)),transform="piece",
                   tf=0.1,
                   countmode=ncol(regmat)+1){
#  print(tf)
  nart <- ncol(regmat)
  ncell <- nrow(regmat)
  jdist <- rep(0, nart * nart)
  dim(jdist) <- c(nart, nart)
  if (transform=="none")
    geodist <- as.matrix(geodist)
  if (transform=="log")
    geodist <- log(as.matrix(tf*geodist)+1)
  if (transform=="sqrt")
    geodist <- sqrt(as.matrix(tf*geodist))
  if (transform=="piece"){
    mg <- max(geodist)
    geodist <- piecewiselin(as.matrix(geodist),tf*mg)
  }
  fi <- list()
  nci <- numeric(0)
  for (i in 1:nart){
    fi[[i]] <- (1:ncell)[as.logical(regmat[,i])]
    nci[[i]] <- sum(regmat[, i])
  }
  for (i in 1:(nart - 1)) {
    if (round(i/countmode)==(i/countmode))
      cat("Computing gb-distances for  species ",i,"\n")
    for (j in (i + 1):nart) {
      nsi <- numeric(0)
      gfi <- geodist[fi[[i]],fi[[j]],drop=FALSE]
#      print(gfi)
      nsi[1] <- sum(apply(gfi,1,min))
      nsi[2] <- sum(apply(gfi,2,min))
      jdist[i,j] <- jdist[j,i] <- (nsi[1]/nci[i]+nsi[2]/nci[j])/2
      if (is.na(jdist[i, j])) 
        cat("Warning! NA at i=", i, ", j=", j, "\n")
    }
  }
  jdist
}

geo2neighbor <- function(geodist,cut=0.1*max(geodist)){
  geodist <- as.matrix(geodist)
  n <- nrow(geodist)
  nblist <- list()
  for (i in 1:n) nblist[[i]] <- numeric(0)
  for (i in 1:(n-1))
    for(j in (i+1):n)
      if (geodist[i,j]<=cut){
        nblist[[i]] <- c(nblist[[i]],j)
        nblist[[j]] <- c(nblist[[j]],i)
      }
  nblist
}
        

# Quantitative Kulcynski distance (works also as Kulczynski distance)
qkulczynski <- function(regmat){
  nart <- ncol(regmat)
  jdist <- rep(0, nart*nart)
  dim(jdist) <- c(nart,nart)
  for (i in 1:(nart-1)){
#    cat("Row ",i,"\n")
    for (j in (i+1):nart){
      ri <- sum(regmat[,i])
      rj <- sum(regmat[,j])
      srij <- sum(pmin(regmat[,i],regmat[,j])) 
      jdist[j,i] <- jdist[i,j] <- 1 - 0.5* (srij/ri + srij/rj)
      if (is.na(jdist[i,j]))
        cat("Warning! NA at i=",i,", j=", j,"\n")
    }
  }
  jdist
}

  

hprabclust <- function (prabobj, cutdist=0.4, cutout=cutdist,
                        method="complete", nnout=2) 
{
#    require(mva)
#    require(MASS)
#    require(mclust)
    n <- prabobj$n.species
    nnd <- c()
    nout <- rep(TRUE,n)
    for (i in 1:n){
      nnd[i] <- sort(prabobj$distmat[i, ])[nnout + 1]
      if (nnd[i]>cutout) nout[i] <- FALSE
    }
    noisen <- n-sum(nout)
    dm <- as.dist(prabobj$distmat[nout,nout])
    cl1 <- hclust(dm, method=method)
    rclustering <- cl2 <- cutree(cl1, h=cutdist)
    nc <- max(cl2)
#    ncl <- max(cl2)
    csum <- function(nx, cv) {
        out <- c()
        for (i in 1:length(nx)) out[i] <- sum(cv == nx[i])
        out
    }
    cs <- csum(1:nc, cl2)
    ocs <- order(-cs)
    for (i in 1:nc) cl2[rclustering == ocs[i]] <- i
    clustering <- rep(nc+1,n)
    clustering[nout] <- cl2
    nmr <- max((1:nc)[cs[ocs]>nnout])
#    print(cs[ocs])
#    print(nmr)
    rclustering <- rep(nmr+1,n)
    rclustering[clustering<=nmr] <- clustering[clustering<=nmr]
    symbols <- c("N")
    if (nmr>0) symbols <- c(sapply(1:nmr, toString),"N")
    clsym <- symbols[rclustering]
    out <- list(clustering = clustering, rclustering=rclustering,
                cutdist=cutdist,
                cutout=cutout,nnout=nnout,noisen=noisen,
                symbols = clsym, hclustering=cl1)
    class(out) <- "comprabclust"
    out
}

