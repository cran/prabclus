"prabtest" <-
function(x, teststat="distratio",tuning=switch(teststat,distratio=0.25,
                        lcomponent=floor(3*ncol(x$distmat)/4),
                        isovertice=ncol(x$distmat),nn=4,NA), times=1000,
                        pd=NULL,
                      prange=c(0,1), nperp=4, step=0.1, twostep=TRUE,
                      sf.sim=FALSE, sf.const=sf.sim,
         pdfnb=FALSE){
  if (is.null(pd) & x$spatial)
    ac <- autoconst(x, twostep=twostep, prange=prange, nperp=nperp,
                    step1=step, species.fixed=sf.const)
  else{
    if (is.null(pd))
      pd <- 1
    ac <- list(pd=pd, coef=NA)
  }
#  if (!is.logical(pdfnb))
#    pdfnb <- nbdiag(x,pd=ac$pd,plot=FALSE)$pdfnb
  psim <- pop.sim(x$prab,x$nb,teststat=teststat, h0c=ac$pd,
                  dist=x$distance,
                  times=times, testc=tuning, n.species=x$n.species,
                  specperreg=x$specperreg, regperspec=x$regperspec,
                  species.fixed=sf.sim, pdfnb=pdfnb)
  out <- list(results=psim$results, datac=psim$datac,
              p.value=ifelse(teststat=="inclusions", 
              psim$p.above,psim$p.below), tuning=tuning, pd=ac$pd,
              reg=ac$coef,
              teststat=teststat, distance=x$distance, times=times,
              pdfnb=pdfnb)
  class(out) <- "prabtest"
  out
}
