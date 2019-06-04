# library(bootstrap)

# This is no good because for comparing regressions
# we don't have a single null model.
# compareregdist <- function(x,dmx,dmy,grouping){
#   lcl <- levels(as.factor(grouping)[x])
#   cld <- as.factor(grouping[x])==lcl[1]
#   n <- length(grouping[x])
#   distx <- dmx[x,x]
#   disty <- dmy[x,x]
#   nd1 <- length(as.vector(as.dist(distx[cld,cld])))
#   nd2 <- length(as.vector(as.dist(distx[!cld,!cld])))
#   ndummy <- c(rep(1,nd1),rep(0,nd2))
#   xv <- c(as.vector(as.dist(distx[cld,cld])),as.vector(as.dist(distx[!cld,!cld])))
#   yv <- c(as.vector(as.dist(disty[cld,cld])),as.vector(as.dist(disty[!cld,!cld])))
#   cxi <- xv*ndummy
#   lmfull <- lm(yv~ndummy+xv+cxi)
#   lmtogether <- lm(yv~xv)
#   logf <- log(anova(lmtogether,lmfull)$F[2])
#   logf
# }

# This defines communities from a given geographical distance matrix or dist
# object. cutoff is the geographical distance below which individuals are put
# together (by hclust/cutree method, default is single linkage).
# individuals in different groups in grouping
# cannot be in the same community. Probably for later testing 2 groups
# it is best to run this with a 2 groups grouping only.
communities <- function(geodist,grouping=NULL,
                        cutoff=1e-5,method="single"){
  geodist <- as.dist(geodist)
  hg <- hclust(geodist,method=method)
  hgc <- cutree(hg,h=cutoff)
  maxcl <- max(hgc)
  mc <- maxcl+1
  if (!is.null(grouping)){
    fgrouping <- as.factor(grouping)
    for (i in 1:maxcl){
      ilev <- as.vector(levels(droplevels(fgrouping[hgc==i])))
      q <- length(ilev)
      if (q>1)
        for (j in 2:q){
          hgc[(hgc==i) & (fgrouping==ilev[j])] <- mc
          mc <- mc+1
        }
    }   
  }
  hgc
}

phipt <- function(alleleobj,comvector,i,j){
    ssg <- numeric(0)
    ssg[i] <-sum(as.dist(alleleobj$distmat[comvector==i,comvector==i]))/sum(comvector==i)
    ssg[j] <-sum(as.dist(alleleobj$distmat[comvector==j,comvector==j]))/sum(comvector==j)
    sst <- sum(as.dist(alleleobj$distmat[comvector %in% c(i,j),comvector %in% c(i,j)]))/sum(comvector %in% c(i,j))
    msa <- ssa <- sst-ssg[i]-ssg[j]
    ssw <- ssg[i]+ssg[j]
    msw <- ssw/(sum(comvector %in% c(i,j))-2)
    if (is.na(msw)) msw <- 0
    n0 <- sum(comvector %in% c(i,j))-(sum(comvector==i)^2+sum(comvector==j)^2)/sum(comvector %in% c(i,j))
    vap <- (msa-msw)/n0
    phipt <- vap/(vap+msw)
    out <- list(phipt=phipt,vap=vap,n0=n0,sst=sst,ssg=ssg,msa=msa,msw=msw)
    out
}


# Cavalli-Sforza-Edwards Chord distance between communities
# from lists p1, p2. Every list entry
# corresponds to a locus and is a vector holding relative frequencies for
# every allele. A list entry can be a single NA in case all community members
# have NA on that locus.
cfchord <- function(p1,p2){
  n <- length(p1)
  n2 <- length(p2)
  if (n!=n2) stop("Lists p1 and p2 must have the same length.")
  naloci <- numeric(0)
  vd <- rep(NA,n)
  for (i in 1:n)
    if (!is.na(p1[[i]][1]) & !is.na(p2[[i]][1])){
      cosphi <- sum(sqrt(p1[[i]]*p2[[i]]))
      vd[i] <- 2*sqrt(2*(1-cosphi))/pi
    }
  validd <- n-sum(is.na(vd))
  out <- sqrt(n*sum(vd[!is.na(vd)]^2)/validd)
  out
}

shared.problist <- function(p1,p2){
  n <- length(p1)
  n2 <- length(p2)
  if (n!=n2) stop("Lists p1 and p2 must have the same length.")
  naloci <- numeric(0)
  vd <- rep(NA,n)
  for (i in 1:n)
    if (!is.na(p1[[i]][1]) & !is.na(p2[[i]][1])){
      vd[i] <- sum(pmin(p1[[i]],p2[[i]]))
    }
  out <- mean(vd,na.rm=TRUE)
  out
}
 

# Construct input lists for cfchord
diploidcomlist <- function(alleleobj,comvector,diploid=TRUE){
  ncom <- max(comvector)
  if (diploid)
    dcomvector <- rep(comvector,each=2)
  else
    dcomvector <- comvector
  out <- list()
#  print(str(alleleobj$charmatrix))
  for (i in 1:ncom){
    out[[i]]  <- list()
    for (j in 1:alleleobj$n.variables){
#      cat("i=",i,"j=",j,"\n")
#      print(dcomvector==i)
      dvec <- factor(
        as.vector(alleleobj$charmatrix[dcomvector==i,j]),
        levels=alleleobj$alevels)
      tdvec <- as.vector(table(dvec))
      ntdvec <- sum(tdvec)
      if (ntdvec>0)
        out[[i]][[j]] <- tdvec/ntdvec
      else
        out[[i]][[j]] <- NA
    }
  }
  out
}

# This constructs a chord-distance matrix between communities.
# comvector: indicates community for individuals; must be numbered from
# 1 to maximum, or "auto" in which case the communities function is
# used. ... are additional arguments to that function; if grouping=NA,
# communities will ignore grouping.
# distance can be "chord", "phipt", "shared.average", "shared.chakraborty",
# "shared.problist".
# communities requires a geographical distance geodist.
# compute.geodist: if TRUE, geographical distance
# between communities will be created, based on distance geodist between
# individuals.
# out.dist: Will dist objects be given out or matrices?
# out: comvector (community indicator), comgroup (group indicator for
# communities), dist, cgeodist
# phiptna: Value where phipt cannot be computed because of communities of size 1.
communitydist <- function(alleleobj,comvector="auto",distance="chord",
                          compute.geodist=TRUE,out.dist=FALSE,
                          grouping=NULL,geodist=NA,diploid=TRUE,
                          phiptna=NA,...){
  if (identical(comvector,"auto")){
    if(is.null(grouping))
      comvector <- communities(geodist,...)
    else
      comvector <- communities(geodist,grouping,...)
  }
  ncom <- max(comvector)
  dclist <- diploidcomlist(alleleobj,comvector,diploid)
  chorddist <- matrix(0,ncol=ncom,nrow=ncom)
  if (distance=="chord"){
    for(i in 1:(ncom-1))
      for (j in (i+1):ncom)
        chorddist[i,j] <- chorddist[j,i] <- cfchord(dclist[[i]],dclist[[j]])
  }
  if (distance=="shared.problist"){
    for(i in 1:(ncom-1))
      for (j in (i+1):ncom)
        chorddist[i,j] <- chorddist[j,i] <- 1-shared.problist(dclist[[i]],dclist[[j]])
  }
  if (distance=="phipt"){
    ssg <- numeric(0)
    for (i in 1:ncom)
      ssg[i] <-sum(as.dist(alleleobj$distmat[comvector==i,comvector==i]))/sum(comvector==i)
    for(i in 1:(ncom-1))
      for (j in (i+1):ncom){
        sst <- sum(as.dist(alleleobj$distmat[comvector %in% c(i,j),comvector %in% c(i,j)]))/sum(comvector %in% c(i,j))
        msa <- ssa <- sst-ssg[i]-ssg[j]
        ssw <- ssg[i]+ssg[j]
        msw <- ssw/(sum(comvector %in% c(i,j))-2)
#        if (is.na(msw)) msw <- phiptna
        n0 <- sum(comvector %in% c(i,j))-(sum(comvector==i)^2+sum(comvector==j)^2)/sum(comvector %in% c(i,j))
        vap <- (msa-msw)/n0
        chorddist[i,j] <- chorddist[j,i] <- vap/(vap+msw)
        if (is.na(chorddist[i,j])) chorddist[i,j] <- chorddist[j,i] <- phiptna
      }
  }
  if (distance=="shared.average"){
    for(i in 1:(ncom-1))
      for (j in (i+1):ncom)
        chorddist[i,j] <- chorddist[j,i] <- mean(alleleobj$distmat[comvector==i,comvector==j])
  }
  if (distance=="shared.chakraborty"){
    avewithin <- numeric(0)
    for (i in 1:ncom){
      imat <-  alleleobj$distmat[comvector==i,comvector==i]
      avewithin[i] <- mean(imat[upper.tri(imat,diag=TRUE)])
    }
    for(i in 1:(ncom-1))
      for (j in (i+1):ncom)
        chorddist[i,j] <- chorddist[j,i] <-
          1-(2-2*mean(alleleobj$distmat[comvector==i,comvector==j]))/
           (2-avewithin[i]-avewithin[j])
  }
  if (out.dist)
    chorddist <- as.dist(chorddist)
  out <- list(comvector=comvector,dist=chorddist)
  if (compute.geodist){
   cgeodist <- matrix(0,ncol=ncom,nrow=ncom)
   geodist=as.matrix(geodist)
   for(i in 1:(ncom-1))
     for (j in (i+1):ncom)
       cgeodist[i,j] <- cgeodist[j,i] <-
         mean(geodist[comvector==i,comvector==j])
     if (out.dist)
       cgeodist <- as.dist(cgeodist)
   out$cgeodist <- cgeodist 
  }
  if (is.null(grouping))
    out$comgroup  <- rep(1,ncom)
  else{
    out$comgroup <- numeric(0)
    for (i in 1:ncom)
      out$comgroup[i] <- grouping[comvector==i][1]
  }
  out
}


# This computes a regression between distance matrices extracted
# from dmx and dmy by indicator vector x
# param (output parameter) can be 1 or 2, 1 for intercept, 2 for slope.
# x centered by xcenter
regdist <- function(x,dmx,dmy,xcenter=0,param=NULL){
  dmx <- as.matrix(dmx)
  dmy <- as.matrix(dmy)
  distx <- dmx[x,x]
  disty <- dmy[x,x]
#   nd1 <- length(as.vector(as.dist(distx[cld,cld])))
#   nd2 <- length(as.vector(as.dist(distx[!cld,!cld])))
#   ndummy <- c(rep(1,nd1),rep(0,nd2))
  xv <- as.vector(as.dist(distx))-xcenter
  yv <- as.vector(as.dist(disty))
  lmdist <- lm(yv~xv)
  if (is.null(param)) out <- lmdist
  else out <- coef(lmdist)[param]
  out
}
 



# Jackknife-based test for equality of 2 regressions
# distance matrices dmx and dmy, grouping is vector with 1 and 2
# defining the data for the two regressions
regeqdist <- function(dmx,dmy,grouping,groups=levels(as.factor(grouping))[1:2]){ 
  dmx <- as.matrix(dmx)
  dmy <- as.matrix(dmy)
  grouping <- as.factor(grouping)
  n <- length(grouping)
  computable <- sum(grouping==groups[1])>2 & sum(grouping==groups[2])>2
#  cld <- as.factor(grouping[x])==groups[1]
  if (computable){
    dmxc <- dmyc <- jr <- lmfit <- xvi <- yvi <- list()
    nc <- sediff <- coefdiff <- pval <- condition <- numeric(0)
    for (i in 1:2){
      dmxc[[i]] <- dmx[grouping==groups[i],grouping==groups[i]]
      dmyc[[i]] <- dmy[grouping==groups[i],grouping==groups[i]]
      jr[[i]] <- list()
      nc[i] <- sum(grouping==groups[i])
      xvi[[i]] <- as.vector(as.dist(dmxc[[i]]))
      yvi[[i]] <- as.vector(as.dist(dmyc[[i]]))
    }
    xall <- c(xvi[[1]],xvi[[2]])
    xcenter <- mean(xall)
    clm <- jackpseudo <- jackestcl <- jackvarcl <- list()
    jackse <- jackest <- tstat <- tdf <- numeric(0)
    for (i in 1:2){
      jackpseudo[[i]] <- list()
      jackestcl[[i]] <- jackvarcl[[i]] <- numeric(0)
      for (j in 1:2)
        jackpseudo[[i]][[j]] <- numeric(0)
    }
    for (i in 1:2){
      xvi[[i]] <- xvi[[i]]-xcenter
      lmfit[[i]] <- lm(yvi[[i]]~xvi[[i]])
      mm <- model.matrix(~xvi[[i]])
      condition[i] <- kappa(mm)
      clm[[i]] <- coef(lmfit[[i]])
    }
    for (i in 1:2)
      for (j in 1:2)
        jr[[i]][[j]] <- bootstrap::jackknife(1:nc[i],regdist,dmx=dmxc[[i]],
                                  dmy=dmyc[[i]],xcenter=xcenter,param=j)
    for (j in 1:2){
      for (i in 1:2)
        if (is.na(jr[[i]][[j]]$jack.se)) computable <- FALSE
      coefdiff[j] <- clm[[1]][j]-clm[[2]][j]    
#      if(computable){
      for (i in 1:2){
        for (k in 1:nc[i])
          jackpseudo[[i]][[j]][k] <-
            nc[i]*clm[[i]][j]-(nc[i]-1)*jr[[i]][[j]]$jack.values[k]
        jackestcl[[i]][j] <- mean(jackpseudo[[i]][[j]])
        jackvarcl[[i]][j] <- var(jackpseudo[[i]][[j]])
      }
      jackest[j] <- jackestcl[[1]][j]-jackestcl[[2]][j]
      jackse[j] <- sqrt(jackvarcl[[1]][j]/nc[1]+jackvarcl[[2]][j]/nc[2])
      tstat[j] <- jackest[j]/jackse[j]
      # Run Welch's t-test, 2-sided
      tdf[j] <-
        jackse[j]^4/((jackvarcl[[1]][j]/nc[1])^2/(nc[1]-1)+
                     (jackvarcl[[2]][j]/nc[2])^2/(nc[2]-1))
      if (tstat[j]>0)
        pval[j] <- 2*(pt(tstat[j],tdf[j],lower.tail=FALSE))
      else
        pval[j] <- 2*(pt(tstat[j],tdf[j]))
 #     }
 #     else
 #       pval <- tstat <- tdf <- jackest <- jackse <- jackpseudo <- NA
    }
    out <- list(pval=pval,coefdiff=coefdiff,condition=condition,lmfit=lmfit,jr=jr,xcenter=xcenter,tstat=tstat,tdf=tdf,jackest=jackest,jackse=jackse,jackpseudo=jackpseudo,groups=groups)
  }
  else
    out <- NA
  class(out) <- "regeqdist"
  out
}

print.regeqdist <- function(x,...){
  if (identical(x[[1]], NA))
    print("Too few individuals in one group. Regression could not be computed.")
  else{
  cat("Testing equality for distance-based regressions in two groups\n")
  cat(x$groups[1]," and ",x$groups[2],"\n")
  if (all(is.na(x$pval)))
    cat("p-value could not be computed because of ill-conditioned regressions\n (probably there were too few or too geographically concentrated individuals\n in a group).\n")
  else{
    cat("Approx. p-values (intercept, slope): ", x$pval,"\n")
    cat("Approx. Bonferroni p-value: ", min(c(1,min(x$pval)*2)),"\n")
    cat("Difference between coefficients, plain (intercept,slope): ",x$coefdiff,"\n")
    cat("Difference between coefficients, jackknife (intercept,slope): ",x$jackest,"\n")
    cat("Standard error of difference, jackknife (intercept, slope): ",x$jackse,"\n")
    cat("Welch-Satterthwaite degrees of freedom: ",x$tdf,"\n")
  }
  }
#  cat("Condition numbers of regressions (clusters): ",x$condition,"\n")
}
  
# For testing whether distanzes between clusters are compatible with   
# joint regression within clusters, compute regression within clusters,
# compute jackknife se (both parameters?? check formula for confidence/
# prediction intervals), and then check whether between-cluster distances
# are outside prediction interval.


# Compute regression from within-cluster distances only; from all
# distances, and compute difference between the two
# where x is the center of between-cluster distances
# This should vary around 0 under H0 (equality) and se can be estimated.
regdistdiff <- function(x,dmx,dmy,grouping,xcenter=0,xcenterbetween=0){
  clx <- grouping[x]
  groups <- unique(clx)
  dmx <- as.matrix(dmx)
  dmy <- as.matrix(dmy)
  dmxx <- dmx[x,x]
  dmyx <- dmy[x,x]
  dmxc <- dmyc <- lmfit <- list()
  dmxc[[1]] <- c(as.vector(as.dist(dmxx[clx==groups[1],clx==groups[1]])),
                 as.vector(as.dist(dmxx[clx==groups[2],clx==groups[2]])))-xcenter
  dmxc[[2]] <- as.vector(as.dist(dmxx))-xcenter
  dmyc[[1]] <- c(as.vector(as.dist(dmyx[clx==groups[1],clx==groups[1]])),
                 as.vector(as.dist(dmyx[clx==groups[2],clx==groups[2]])))
  dmyc[[2]] <- as.vector(as.dist(dmyx))
  for (i in 1:2)
    lmfit[[i]] <- lm(dmyc[[i]]~dmxc[[i]])
  out <- coef(lmfit[[2]])[1]-coef(lmfit[[1]])[1]+(coef(lmfit[[2]])[2]-coef(lmfit[[1]])[2])*xcenterbetween
  out
}



# Tests whether  regression for between-cluster distances
# is compatible with regressions for within-cluster distances
# at between-cluster central value of x
# by calling regdistdiff and using jackknife 
regdistbetween <- function(dmx,dmy,grouping,groups=levels(as.factor(grouping))[1:2]){
  dmx <- as.matrix(dmx)
  dmy <- as.matrix(dmy)
  grouping <- as.factor(grouping)
  n <- sum(grouping %in% groups)
  computable <- (n>2)
#  cld <- as.factor(grouping[x])==groups[1]
  dmxc <- dmyc <- lmfit <- list()
  coefdiff <- pval <- condition <- numeric(0)
  dmxc[[1]] <- c(as.vector(as.dist(dmx[grouping==groups[1],grouping==groups[1]])),
                 as.vector(as.dist(dmx[grouping==groups[2],grouping==groups[2]])))
  xcenter <- mean(dmxc[[1]])
  dmxc[[1]] <- dmxc[[1]]-xcenter
  dmxc[[2]] <- as.vector(as.dist(dmx[grouping %in% groups, grouping %in% groups]))-xcenter
  dmxcbetween <- as.vector(dmx[grouping==groups[1],grouping==groups[2]])-xcenter
  xcenterbetween <- mean(dmxcbetween)
  dmyc[[1]] <- c(as.vector(as.dist(dmy[grouping==groups[1],grouping==groups[1]])),
                 as.vector(as.dist(dmy[grouping==groups[2],grouping==groups[2]])))
  dmyc[[2]] <- as.vector(as.dist(dmy[grouping %in% groups, grouping %in% groups]))
  jackpseudo <- numeric(0)
  for (i in 1:2){
    lmfit[[i]] <- lm(dmyc[[i]]~dmxc[[i]])
#    jr[[i]] <- list()
#    nc[i] <- sum(grouping==groups[i])
    mm <- model.matrix(~dmxc[[i]])
    condition[i] <- kappa(mm)
  }
  jr <- bootstrap::jackknife(1:n,regdistdiff,
                  dmx=dmx[grouping %in% groups, grouping %in% groups],
                  dmy=dmy[grouping %in% groups, grouping %in% groups],
                  grouping=grouping[grouping %in% groups],xcenter=xcenter,
                         xcenterbetween=xcenterbetween)
#    sediff[j] <- sqrt(jr[[1]][[j]]$jack.se^2+jr[[2]][[j]]$jack.se^2)
  if (is.na(jr$jack.se)) computable <- FALSE
  coefdiff <- (coef(lmfit[[2]])[1]-coef(lmfit[[1]])[1])+(coef(lmfit[[2]])[2]-coef(lmfit[[1]])[2])*xcenterbetween
  if (computable){
      for (k in 1:n)
        jackpseudo[k] <- n*coefdiff-(n-1)*jr$jack.values[k]
      jackest <- mean(jackpseudo)
      jackse <- sd(jackpseudo)/sqrt(n)
      tstat <- jackest/jackse
      # One-sided t-test
      tdf <- n-1
      pval <- pt(tstat,tdf,lower.tail=FALSE)
  }
  else
      pval <- tstat <- tdf <- jackest <- jackse <- jackpseudo <- NA
  testname <- "Testing whether regression for between-group distances\n is compatible with regressions for within-group distances"  
  out <- list(pval=pval,coefdiff=coefdiff,condition=condition,lmfit=lmfit,jr=jr,xcenter=xcenter,xcenterbetween=xcenterbetween,tstat=tstat,tdf=tdf,jackest=jackest,jackse=jackse,jackpseudo=jackpseudo,groups=groups,testname=testname)
  class(out) <- "regdistbetween"
  out
}
 
print.regdistbetween <- function(x,...){
  if (identical(x[[1]], NA))
    print("Too few individuals in one group. Regression could not be computed.")
  else{
    cat(x$testname,"\n")
    cat("Groups: ",x$groups[1]," and ",x$groups[2],"\n")
    if (all(is.na(x$pval)))
      cat("p-value could not be computed because of ill-conditioned regressions\n (probably there were too few or too geographically concentrated individuals\n in a species).\n")
    else{
      cat("Approx. p-value: ", x$pval,"\n")
      cat("Difference between coefficients at between groups center (plain): ",x$coefdiff,"\n")
      cat("Difference between coefficients at between groups center (jackknife): ",x$jackest,"\n")
      cat("Standard error of difference (jackknife): ",x$jackse,"\n")
      cat("Welch-Satterthwaite degrees of freedom: ",x$tdf,"\n")
    }
  }
#  cat("Condition numbers of regressions (within-cluster/all): ",x$condition,"\n")
}


# Compute regression from within-cluster distances of one species only
# (species gives number); from
# those and the distances to species 2,
# and compute difference between the two at center of between-cluster
# distances.
# This should vary around 0 under H0 (equality) and se can be estimated.
regdistdiffone <- function(x,dmx,dmy,grouping,xcenter=0,xcenterbetween=0,rgroup){
  clx <- grouping[x]
  groups <- unique(clx)
  dmxx <- dmx[x,x]
  dmyx <- dmy[x,x]  
#  print(x)
#  print(str(dmxx))
#  print(groups)
#  print(grouping)
#  print(clx)
  dmxc <- dmyc <- lmfit <- list()
  dmxc[[1]] <- as.vector(as.dist(dmxx[clx==rgroup,clx==rgroup]))-xcenter
  dmxc[[2]] <- c(dmxc[[1]],as.vector(dmxx[clx==groups[1],clx==groups[2]])-xcenter)
  dmyc[[1]] <- as.vector(as.dist(dmyx[clx==rgroup,clx==rgroup]))
  dmyc[[2]] <- c(dmyc[[1]],as.vector(dmyx[clx==groups[1],clx==groups[2]]))
  for (i in 1:2)
    lmfit[[i]] <- lm(dmyc[[i]]~dmxc[[i]])
#  print(str(lmfit))
  out <- coef(lmfit[[2]])[1]-coef(lmfit[[1]])[1]+(coef(lmfit[[2]])[2]-coef(lmfit[[1]])[2])*xcenterbetween
  out
}


# Tests whether  regression for between-cluster distances
# is compatible with regressions for within-cluster distances
# of a single rgroup (rgroup given number) 
# by calling regdistdiffone and using jackknife
regdistbetweenone <- function(dmx,dmy,grouping,groups=levels(as.factor(grouping))[1:2],rgroup){
  dmx <- as.matrix(dmx)
  dmy <- as.matrix(dmy)
  grouping <- as.factor(grouping)
  n <- sum(grouping %in% groups)
  n1 <- sum(grouping==rgroup)
  computable <- (n1>2)
  if (computable){
  n2 <- n-n1
#  cld <- as.factor(grouping[x])==groups[1]
  jackse <- jackest <- tstat <- tdf <- numeric(0)
  dmxc <- dmyc <- jr <- lmfit <- list()
  coefdiff <- pval <- condition <- sediff <- numeric(0)
  dmxc[[1]] <- as.vector(as.dist(dmx[grouping==rgroup,grouping==rgroup]))
  xcenter <- mean(dmxc[[1]])
  dmxc[[1]] <- dmxc[[1]]-xcenter
  dmxcbetween <- as.vector(dmx[grouping==groups[1],grouping==groups[2]])-xcenter
  xcenterbetween <- mean(dmxcbetween)
  dmxc[[2]] <- c(dmxc[[1]],dmxcbetween)
  dmyc[[1]] <- as.vector(as.dist(dmy[grouping==rgroup,grouping==rgroup]))
  dmyc[[2]] <- c(dmyc[[1]],as.vector(dmy[grouping==groups[1],grouping==groups[2]]))
  jackpseudo <- numeric(0)
  for (i in 1:2){
    lmfit[[i]] <- lm(dmyc[[i]]~dmxc[[i]])
#    jr[[i]] <- list()
#    nc[i] <- sum(grouping==groups[i])
    mm <- model.matrix(~dmxc[[i]])
    condition[i] <- kappa(mm)
  }
  jr <- bootstrap::jackknife(1:n,regdistdiffone,
                  dmx=dmx[grouping %in% groups, grouping %in% groups],
                  dmy=dmy[grouping %in% groups, grouping %in% groups],
                  grouping=grouping[grouping %in% groups],xcenter=xcenter,
                  xcenterbetween=xcenterbetween,rgroup=rgroup)
  if (is.na(jr$jack.se)) computable <- FALSE
  coefdiff <- (coef(lmfit[[2]])[1]-coef(lmfit[[1]])[1])+(coef(lmfit[[2]])[2]-coef(lmfit[[1]])[2])*xcenterbetween
#  if (computable){
  for (k in 1:n)
    jackpseudo[k] <- n*coefdiff-(n-1)*jr$jack.values[k]
  jackest <- mean(jackpseudo)
  jackvar1 <- var(jackpseudo[grouping[grouping %in% groups]==groups[1]])
  jackvar2 <- var(jackpseudo[grouping[grouping %in% groups]==groups[2]])
  jackvar <- (n1*jackvar1+n2*jackvar2)/(n1+n2)
  jackse <- sqrt(jackvar/n)
  tstat <- jackest/jackse
  # One-sided t-test, df according to Welch-Satterthwaite equation
  tdf <- jackvar^2/((n1*jackvar1/(n1+n2))^2/(n1-1)+(n2*jackvar2/(n1+n2))^2/(n2-1))
  pval <- pt(tstat,tdf,lower.tail=FALSE)
    
#  }
#  else
#      pval <- tstat <- tdf <- jackest <- jackse <- jackpseudo <- jackvar1 <- jackvar2 <- NA

  testname <- paste("Testing whether regression for between-group distances\n is compatible with regressions for within-group distances for species ",rgroup,"\n The maximum Bonferroni p-value over both species can be used to test\n the H0 that the between-group distances are compatible with within-group\n distances of at least one group.")
  out <- list(pval=pval,coefdiff=coefdiff,condition=condition,lmfit=lmfit,jr=jr,xcenter=xcenter,xcenterbetween=xcenterbetween,rgroup=rgroup,tstat=tstat,tdf=tdf,jackest=jackest,jackse=jackse,jackpseudo=jackpseudo,jackvar1=jackvar1,jackvar2=jackvar2,groups=groups,rgroup=rgroup,testname=testname)
  }
  else
    out <- NA
  class(out) <- "regdistbetween"
  out
}


plotdistreg <- function(dmx,dmy,grouping,
                        groups=levels(as.factor(grouping))[1:2],cols=c(1,2,3,4),
                        pchs=rep(1,3),
                        ltys=c(1,2,1,2),
                        individual=TRUE,jointwithin=TRUE,jointall=TRUE,
                        oneplusjoint=TRUE,jittering=TRUE,bcenterline=TRUE,
                        xlim=NULL,ylim=NULL,xlab="geographical distance",
                        ylab="genetic distance",...){
#  if (is.null(groups) & length(levels(as.factor(grouping))==2)){
#    if (is.numeric(grouping))
#      groups <- range(grouping)
#    else
#      groups <- levels(grouping)
  #  }
  grouping <- as.factor(grouping)
  dmx <- as.matrix(dmx)
  dmy <- as.matrix(dmy)
  selector <- grouping %in% groups
  selclust <- grouping[selector]
  ci <- list()
  ci[[1]] <- grouping==groups[1]
  ci[[2]] <- grouping==groups[2]
  betweencenter <- mean(dmx[ci[[1]],ci[[2]]])
  if (is.null(xlim))
    xlim=c(min(as.dist(dmx[ci[[1]]|ci[[2]],ci[[1]]|ci[[2]]])),max(as.dist(dmx[ci[[1]]|ci[[2]],ci[[1]]|ci[[2]]])))
  if (is.null(ylim))
    ylim=c(min(as.dist(dmy[ci[[1]]|ci[[2]],ci[[1]]|ci[[2]]])),max(as.dist(dmy[ci[[1]]|ci[[2]],ci[[1]]|ci[[2]]])))  
  if (jittering){
    plot(jitter(as.dist(dmx[ci[[1]],ci[[1]]])),
         jitter(as.dist(dmy[ci[[1]],ci[[1]]])),
         col=cols[1],pch=pchs[1],xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
    points(jitter(as.dist(dmx[ci[[2]],ci[[2]]])),
         jitter(as.dist(dmy[ci[[2]],ci[[2]]])),
         col=cols[2],pch=pchs[2])
    points(jitter(as.vector(dmx[ci[[1]],ci[[2]]])),
         jitter(as.vector(dmy[ci[[1]],ci[[2]]])),
         col=cols[3],pch=pchs[3])
  }
  else{
    plot(as.dist(dmx[ci[[1]],ci[[1]]]),
         as.dist(dmy[ci[[1]],ci[[1]]]),
         col=cols[1],pch=pchs[1],xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
    points(as.dist(dmx[ci[[2]],ci[[2]]]),
         as.dist(dmy[ci[[2]],ci[[2]]]),
         col=cols[2],pch=pchs[2])
    points(as.vector(dmx[ci[[1]],ci[[2]]]),
         as.vector(dmy[ci[[1]],ci[[2]]]),
         col=cols[3],pch=pchs[3])
  }    
  if(individual){
      lind <- list()
      for (i in 1:2){
        lind[[i]] <- lm(as.dist(dmy[ci[[i]],ci[[i]]])~as.dist(dmx[ci[[i]],ci[[i]]]))
        abline(coef(lind[[i]]),col=cols[i],lty=ltys[1])
      }
  }
  if (bcenterline)
    lines(c(betweencenter,betweencenter),c(-1,1.5*max(ylim)),col=cols[4])
  if (jointwithin){
      lmjw <- lm(c(as.dist(dmy[ci[[1]],ci[[1]]]),as.dist(dmy[ci[[2]],ci[[2]]]))~c(as.dist(dmx[ci[[1]],ci[[1]]]),as.dist(dmx[ci[[2]],ci[[2]]])))
      abline(coef(lmjw),col=cols[3],lty=ltys[2])
  }
  if (jointall){
    lmall <- lm(as.dist(dmy[ci[[1]]|ci[[2]],ci[[1]]|ci[[2]]])~as.dist(dmx[ci[[1]]|ci[[2]],ci[[1]]|ci[[2]]]))
    abline(coef(lmall),col=cols[3],lty=ltys[3])
  }
  if(oneplusjoint){
      lop <- list()
      for (i in 1:2){
        lop[[i]] <- lm(c(as.dist(dmy[ci[[i]],ci[[i]]]),as.vector(dmy[ci[[1]],ci[[2]]]))~c(as.dist(dmx[ci[[i]],ci[[i]]]),as.vector(dmx[ci[[1]],ci[[2]]])))
        abline(coef(lop[[i]]),col=cols[i],lty=ltys[4])
      }
  }
}
    
    
alleleconvert <- function (file = NULL, strmatrix = NULL, format.in = "genepop", 
    format.out = "prabclus", alength = 3, orig.nachar = "000", 
    new.nachar = "-", rows.are.individuals = TRUE, firstcolname = FALSE, 
    aletters = intToUtf8(c(65:90, 97:122), multiple = TRUE), 
    outfile = NULL, skip=0) 
{
    if (is.null(file)) 
        m1 <- strmatrix
    else m1 <- read.table(file, colClasses = "character", skip=skip)
    if (firstcolname) {
        if (format.in == "structure") 
            colnamesv <- m1[seq(1, nrow(m1) - 1, by = 2), 1]
        if (format.in %in% c("genepop","structureb")) 
            colnamesv <- m1[, 1]
        if (format.in == "structureb")
          m1 <- m1[, 3:ncol(m1)]
        else
          m1 <- m1[, 2:ncol(m1)] 
    }
    if (!rows.are.individuals) 
        m1 <- t(m1)
    if (format.in == "genepop") {
        n.individuals <- nrow(m1)
        n.variables <- ncol(m1)
    }
    if (format.in == "structure") {
        n.individuals <- nrow(m1)/2
        n.variables <- ncol(m1)
    }
    if (format.in == "structureb") {
        n.individuals <- nrow(m1)
        n.variables <- ncol(m1)/2
    }
    first.allele <- second.allele <- outmatrix <- matrix("", 
        ncol = n.variables, nrow = n.individuals)
    double.out <- matrix("", ncol = n.variables, nrow = 2 * n.individuals)
    if (format.in == "genepop") {
        for (i in 1:n.individuals) for (j in 1:n.variables) {
            first.allele[i, j] <- substring(m1[i, j], 1, alength)
            second.allele[i, j] <- substring(m1[i, j], alength + 
                1, 2 * alength)
        }
    }
    if (format.in == "structure") {
        for (i in 1:n.individuals) for (j in 1:n.variables) {
            first.allele[i, j] <- m1[2 * i - 1, j]
            second.allele[i, j] <- m1[2 * i, j]
        }
    }
    if (format.in == "structureb") {
        for (i in 1:n.individuals) for (j in 1:n.variables) {
            first.allele[i, j] <- m1[i, 2*j-1]
            second.allele[i, j] <- m1[i, 2*j]
        }
    }
    double.allele <- rbind(first.allele, second.allele)
    maxk <- 0
    for (j in 1:n.variables) {
        jlevels <- levels(as.factor(double.allele[, j]))
        jlevels <- jlevels[jlevels != orig.nachar]
        if (length(jlevels) > maxk) 
            maxk <- length(jlevels)
        double.out[double.allele[, j] == orig.nachar, j] <- new.nachar
        if (format.out == "prabclus") {
            for (k in 1:length(jlevels)) double.out[double.allele[, 
                j] == jlevels[k], j] <- aletters[k]
            outmatrix[, j] <- mapply(paste, double.out[1:n.individuals, 
                j], double.out[(n.individuals + 1):(2 * n.individuals), 
                j], sep = "")
            out <- structure(outmatrix, alevels = aletters[1:maxk])
        }
        if (format.out == "structurama") {
            for (k in 1:length(jlevels)) double.out[double.allele[, 
                j] == jlevels[k], j] <- jlevels[k]
            outmatrix[, j] <- mapply(paste, "(", double.out[1:n.individuals, 
                j], ",", double.out[(n.individuals + 1):(2 * 
                n.individuals), j], ")", sep = "")
            out <- structure(outmatrix, alevels = jlevels)
        }
    }
    if (!is.null(outfile)) {
        if (firstcolname) 
            outmatrix <- cbind(colnamesv, outmatrix)
        write.table(outmatrix, outfile, quote = FALSE, col.names = FALSE)
    }
    out
}
  
    
    
    
         

