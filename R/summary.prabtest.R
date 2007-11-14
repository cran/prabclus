"summary.prabtest" <-
function(object, above.p=object$teststat %in% c("groups","inclusions","mean"),
         group.outmean=FALSE,...){
  if (object$teststat!="groups"){
    rrange <- range(object$results)
    rmean <- mean(object$results)
    groupinfo <- NULL
  }
  else{
    rrange <- range(object$results$overall)
    rmean <- mean(object$results$overall)
    rrangem <- range(object$results$mean)
    rmeanm <- mean(object$results$mean)
    rrangeg <- matrix(0,ncol=2,nrow=object$groupinfo$ng)
    rmeang <- numeric(0)
    for (i in 1:object$groupinfo$ng){
      rrangeg[i,] <- range(object$results$gr[i,])
      rmeang[i] <- mean(object$results$gr[i,])
    }    
    groupinfo <- c(object$groupinfo,list(rrangeg=rrangeg,rmeang=rmeang,
                                         rrangem=rrangem,rmeanm=rmeanm))
  }
  p.value <- if (is.null(object$abund)) object$p.value
             else{
               if (above.p) object$p.above
               else object$p.below
             }
  out <- list(rrange=rrange, rmean=rmean, datac=object$datac,
              p.value=p.value,
              pd=object$pd, tuning=object$tuning,
              teststat=object$teststat, distance=object$distance,
              times=object$times,
              pdfnb=object$pdfnb,
              abund=object$abund, sarlambda=object$sarlambda,
              groupinfo=groupinfo, group.outmean=group.outmean)
  class(out) <- "summary.prabtest"
  out
}


"print.summary.prabtest" <-
function(x, ...){
  if (is.null(x$abund))
    cat("* Parametric bootstrap test for presence-absence data *\n\n")
  else
    cat("* Parametric bootstrap test for spatial abundance data *\n\n")
  cat("Test statistics: ",x$teststat,", Tuning constant=",x$tuning,"\n")
  cat("Distance: ",x$distance,"\n")
  cat("Simulation runs: ",x$times,"\n")
  if(is.null(x$pd))    
    cat("A model without spatial autocorrelation was used.\n")
  else{
    cat("Disjunction parameter for presence-absence pattern: ",x$pd,"\n")
    if (!is.null(x$abund))
      cat("Neighborhood parameter lambda for SAR-model: ",x$sarlambda,"\n")
  }
  if (!is.null(x$pdfnb))
    if (x$pdfnb)
      cat("Neighbor-based correction of region probabilities was used.\n")
  if(x$teststat=="groups"){
    cat("Mean within group distances for original data: ",x$datac$overall,"\n")
    cat("Mean of mean within group distances for null data: ",x$rmean,
        ",\n range: ",x$rrange,"\n")
    cat("p= ",x$p.value,"\n")
    if (x$group.outmean){
      cat("Overall mean for original data: ",x$groupinfo$testm,"\n")
      cat("Mean of overall mean for null data: ",x$groupinfo$rmeanm,
          ", range: ",x$groupinfo$rrangem,"\n")
      cat("p= ",x$groupinfo$pma,"\n")
    }
    for (i in 1:x$groupinfo$ng)
      if (x$groupinfo$nsg[i]>1){
        cat("  Group ",x$groupinfo$lg[i],
            " statistics value for original data: ",
          x$datac$gr[i],"\n")
        cat("  Group ",x$groupinfo$lg[i],
            " mean for null data: ",x$groupinfo$rmeang[i],
          ", range: ",x$groupinfo$rrangeg[i,],"\n")
        cat("  p= ",x$groupinfo$pa[i],"\n")
      }
      
  }
  else{
    cat("Statistics value for original data: ",x$datac,"\n")
    cat("Mean for null data: ",x$rmean,", range: ",x$rrange,"\n")
    cat("p= ",x$p.value,"\n")
  }
}
