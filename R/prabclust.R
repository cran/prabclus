"prabclust" <-
function(prabobj, mdsmethod="classical", mdsdim=4,
                      nnk=ceiling(prabobj$n.species/40), nclus=0:9,
                      modelid="noVVV"){
  require(mva)
  require(MASS)
  require(mclust)
  dm <- prabobj$distmat
  if (mdsmethod!="classical"){
    mindm <- min(dm[dm>0])/10
    for (i in 1:(prabobj$n.species-1))
      for (j in (i+1):prabobj$n.species)
        if (dm[i,j]<mindm) dm[i,j] <- dm[j,i] <- mindm
  }
  mds <- switch(mdsmethod,
                classical = cmdscale(dm, k=mdsdim),
                kruskal = isoMDS(dm, k=mdsdim)$points,
                sammon = sammon(dm, k=mdsdim)$points)
  kn <- NNclean(mds, k=nnk)
  if (modelid=="all")
    kem <- EMclustN(mds, G=nclus, noise=1-kn$z)
  else{
    if (modelid=="noVVV")
      kem <- EMclustN(mds, G=nclus,
                           emModelNames=c("EII","VII","EEI","VEI",
                             "EVI", "VVI", "EEE","EEV", "VEV"),
                           noise=1-kn$z)
    else
      kem <- EMclustN(mds, G=nclus,
                     emModelNames=modelid, noise=1-kn$z)
  }
  skem <- summary(kem,mds)
  uclustering <- skem[[4]]
#  print(kem)
  nc <- max(uclustering)
  noisec <- ifelse(nc==ncol(skem$mu),FALSE,TRUE)
#   print(nc)
#   print(skem[[5]])
#   print(skem[[6]])
#   print(skem[[7]])
#   print(skem[[4]])
  clustering <- uclustering
  csum <- function(n, cv){
    out <- c()
    for (i in 1:length(n))
      out[i] <- sum(cv==n[i])
    out
  }
  ncl <- ifelse(noisec, nc-1, nc)
  cs <- csum(1:ncl,clustering)
  ocs <- order(-cs)
  for (i in 1:ncl)
    clustering[uclustering==ocs[i]] <- i
  if (noisec & nc==1) symbols <- c("N")
  else{
    if (noisec)
      symbols <- c(sapply(1:ncl, toString),"N")
    else
      symbols <- sapply(1:nc, toString)
  }
  clsym <- symbols[clustering]
  for (i in 1:ncl)
    if (sum(clustering==i)<2)
      clsym[clustering==i] <- "N"      
  plot(mds, pch=clsym)
  out <- list(clustering=clustering, clustsummary=skem, bicsummary=kem,
              points=mds, nnk=nnk, mdsdim=mdsdim, mdsmethod=mdsmethod,
              symbols=clsym)
  class(out) <- "prabclust"
  out
}
