"prabclust" <- function (prabobj, mdsmethod = "classical", mdsdim = 4,
                       nnk = ceiling(prabobj$n.species/40), 
    nclus = 0:9, modelid = "all", permutations=0) 
{
    require(MASS)
    require(mclust)
    # "Data-alphabetical" ordering 
    oregions <- order(prabobj$specperreg)
    prabo1 <- prabobj$prab[oregions,]
    ospecies <- do.call("order",as.data.frame(t(prabo1)))    
    dm <- prabobj$distmat[ospecies,ospecies]
    if (mdsmethod != "classical") {
        mindm <- min(dm[dm > 0])/10
        for (i in 1:(prabobj$n.species - 1))
          for (j in (i + 1):prabobj$n.species) if (dm[i, j] < mindm) 
            dm[i, j] <- dm[j, i] <- mindm
    }
    mdsbest <- mdsout <-
      switch(mdsmethod, classical = cmdscale(dm, k = mdsdim), 
        kruskal = isoMDS(dm, k = mdsdim), sammon = sammon(dm, 
            k = mdsdim))
    if (mdsmethod=="classical") mds <- mdsout
    else mds <- mdsout$points
    permchange=FALSE
    operm <- NULL
    if (permutations>0){
      if (mdsmethod!="classical"){
        for (i in 1:permutations){
          inumbers <- sample(1:prabobj$n.species,prabobj$n.species)
          dmperm <- dm[inumbers,inumbers]
          mdsout <- switch(mdsmethod,  
            kruskal = isoMDS(dmperm, k = mdsdim),
                        sammon = sammon(dm, k = mdsdim))
          if (mdsout$stress<mdsbest$stress){
            mdsbest <- mdsout
            operm <- inumbers
            permchange <- TRUE
          } # if mds$stress<beststress
        } # for i
        mdsout <- mdsbest
        if (!is.null(operm)){
          mdsout$points[operm,] <- mdsbest$points
          mds <- mdsout$points
        }
        operm <- NULL
      } # if mdsmethod!="classical"
    } # if permutations>0  
    if (nnk==0){
      if (0==nclus[1])
        nclus <- nclus[-1]
      if (modelid == "all") 
        kem <- mclustBIC(mds, G = nclus, modelNames=c("EII", 
                "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
                "VEV","VVV"))
      else {
        if (modelid == "noVVV") 
            kem <- mclustBIC(mds, G = nclus, modelNames = c("EII", 
                "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
                "VEV"))
        else kem <- mclustBIC(mds, G = nclus, modelNames = modelid)
      }
    }
    else{  
      kn <- NNclean(mds, k = nnk)
      if (modelid == "all") 
        kem <- mclustBIC(mds, G = nclus, modelNames=c("EII", 
                "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
                "VEV","VVV"),
                         initialization=list(noise = as.logical(1 - kn$z)))
      else {
        if (modelid == "noVVV") 
            kem <- mclustBIC(mds, G = nclus, modelNames = c("EII", 
                "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
                "VEV"), initialization=list(noise = as.logical(1 - kn$z)))
        else kem <- mclustBIC(mds, G = nclus, modelNames = modelid, 
            initialization=list(noise =as.logical( 1 - kn$z)))
      }
    }
    skembest <- skem <- summary(kem, mds)
    if (permutations>0){
      for (i in 1:permutations){
        inumbers <- sample(1:prabobj$n.species,prabobj$n.species)
        mdsperm <- mds[inumbers,]
        if (nnk==0){
          if (0==nclus[1])
            nclus <- nclus[-1]
          if (modelid == "all") 
            kem <- mclustBIC(mdsperm, G = nclus, modelNames=c("EII", 
                "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
                "VEV","VVV"))
          else {
            if (modelid == "noVVV") 
                kem <- mclustBIC(mdsperm, G = nclus, modelNames = c("EII", 
                    "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
                    "VEV"))
            else kem <- mclustBIC(mdsperm, G = nclus, modelNames = modelid)
          } # else (modelid!="all")
        } # if nnk==0
        else{  
          kn <- NNclean(mdsperm, k = nnk)
          if (modelid == "all") 
            kem <- mclustBIC(mdsperm, G = nclus, modelNames=c("EII", 
                    "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
                    "VEV","VVV"),initialization=
                             list(noise = as.logical(1 - kn$z)))
          else {
            if (modelid == "noVVV") 
              kem <- mclustBIC(mdsperm, G = nclus, modelNames = c("EII", 
                    "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", 
                    "VEV"), initialization=list(noise = as.logical(1 - kn$z)))
            else kem <- mclustBIC(mdsperm, G = nclus, modelNames = modelid, 
                initialization=list(noise =as.logical(1 - kn$z)))
          } # else - modelid!="all"
        } # else - nnk>0
        skem <- summary(kem, mdsperm)
        if (skem$bic>skembest$bic){
          skembest <- skem
          operm <- inumbers
          permchange <- TRUE
        }
      } # for i
    } # if permutations>0    
    skem <- skembest
    if (!is.null(operm)){
      skem$classification[operm] <- skembest$classification
      skem$z[operm,] <- skembest$z
      if (!is.null(attr(skembest,"initialization")$noise))
        attr(skem,"initialization")$noise[operm] <-
          attr(skembest,"initialization")$noise
      skembest <- skem
    }
    mdsr <- mds
    mds[ospecies,] <- mdsr
    skem$classification[ospecies] <- skembest$classification
    skem$z[ospecies,] <- skembest$z
    if (!is.null(attr(skembest,"initialization")$noise))
      attr(skem,"initialization")$noise[ospecies] <-
        attr(skembest,"initialization")$noise
    uclustering <- skem$classification
    ncl <- max(uclustering)
    nc <- ncl+1
    noisec <- min(uclustering)==0
    clustering <- uclustering
    csum <- function(n, cv) {
        out <- c()
        for (i in 1:length(n)) out[i] <- sum(cv == n[i])
        out
    }
#    print(skem$z)
#    str(skem)
#    print(ncl)
#    print(nc)
    cs <- csum(1:ncl, clustering)
    ocs <- order(-cs)
    for (i in 1:ncl) clustering[uclustering == ocs[i]] <- i
    if (noisec & nc == 1)
      clsym <- rep("N",prabobj$n.species)
    else {
      if (noisec){ 
         symbols <- c("N",sapply(1:ncl, toString))
         clsym <- symbols[clustering+1]
       }
       else{
         symbols <- sapply(1:ncl, toString)
         clsym <- symbols[clustering]
       }
    }
    for (i in 1:ncl) if (sum(clustering == i) < 2) 
        clsym[clustering == i] <- "N"
    plot(mds, pch = clsym)
    out <- list(clustering = clustering, clustsummary = skem, 
        bicsummary = kem, points = mds, nnk = nnk, mdsdim = mdsdim, 
        mdsmethod = mdsmethod, symbols = clsym, permutations=permutations,
                permchange=permchange)
    class(out) <- "prabclust"
    out
}
