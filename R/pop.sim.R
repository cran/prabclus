"pop.sim" <-
function(regmat, neighbors, h0c=1, times=200, dist="kulczynski",
                    teststat="isovertice",testc=NULL,geodist=NULL,gtf=0.1,
                    n.species=ncol(regmat), specperreg=NULL,
                    regperspec=NULL, species.fixed=FALSE, pdfnb=FALSE){
  statres <- rep(0,times)
  if (is.null(specperreg))
      nregions <- apply(regmat,1,sum)
  else
      nregions <- specperreg
  if (is.null(regperspec))
      nspecies <- apply(regmat,2,sum)
  else
      nspecies <- regperspec
  for (i in 1:times)
  {
    cat("Simulation run ",i)
    mat <- randpop.nb(neighbors,p.nb=h0c,n.species=n.species,
                        vector.species=nspecies, species.fixed=species.fixed,
                        pdf.regions=nregions/sum(nregions),count=FALSE,
                      pdfnb=pdfnb)
    if (teststat!="inclusions"){
      if (dist=="jaccard")
        distm <- jaccard(mat)
      if (dist=="kulczynski")
        distm <- kulczynski(mat)
      if (dist=="geco")
        distm <- geco(mat,geodist,tf=gtf)
    }
    else
      statres[i] <- incmatrix(mat)$ninc
    if (teststat=="isovertice")
    {
      test <- homogen.test(distm,ne=testc)
      statres[i] <- test$iv
    }
    if (teststat=="lcomponent")
      statres[i] <- lcomponent(distm,ne=testc)$lc
    if (teststat=="distratio")
      statres[i] <- distratio(distm,prop=testc)$dr
    if (teststat=="nn")
      statres[i] <- nn(distm,ne=testc)
    cat(" statistics value=",statres[i],"\n")
  }
  if (teststat!="inclusions"){
    if (dist=="jaccard")
      distm <- jaccard(regmat)
    if (dist=="kulczynski")
      distm <- kulczynski(regmat)
    if (dist=="geco")
      distm <- geco(regmat,geodist,tf=gtf)
  }
  else{
    test <- incmatrix(regmat)$ninc
    p.above <- (1+sum(statres>=test))/(1+times)
    p.below <- (1+sum(statres<=test))/(1+times)
    datac <- test
  }
  if (teststat=="isovertice")
  {
    test <- homogen.test(distm,ne=testc)
    p.above <- (1+sum(statres>=test$iv))/(1+times)
    p.below <- (1+sum(statres<=test$iv))/(1+times)
    pb <- min(p.above,p.below)*2
    p.above <- max(p.above,p.below)
    p.below <- pb
    datac <- test$iv
    testc <- test$ne
  }
  if (teststat=="lcomponent")
  {
    test <- lcomponent(distm,ne=testc)
    p.above <- (1+sum(statres>=test$lc))/(1+times)
    p.below <- (1+sum(statres<=test$lc))/(1+times)
    datac <- test$lc
    testc <- test$ne
  }
  if (teststat=="nn")
  {
    test <- nn(distm,ne=testc)
    p.above <- (1+sum(statres>=test))/(1+times)
    p.below <- (1+sum(statres<=test))/(1+times)
    datac <- test
  }
  if (teststat=="distratio")
  {
    test <- distratio(distm,prop=testc)
    p.above <- (1+sum(statres>=test$dr))/(1+times)
    p.below <- (1+sum(statres<=test$dr))/(1+times)
    datac <- test$dr
    testc <- test$prop
  }
  cat("Data value: ",datac,"\n")
  out <- list(results=statres,p.above=p.above,p.below=p.below,datac=datac,
              testc=testc)
  out
}
