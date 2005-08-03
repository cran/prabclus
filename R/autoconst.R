"autoconst" <-
function(x, prange=c(0,1), twostep=TRUE,
                      step1=0.1, step2=0.01, plot=TRUE, nperp=4,
                      ejprob=NULL, species.fixed=TRUE, pdfnb=FALSE){
  probs <- prange[1]+step1*(0:round((prange[2]-prange[1])/step1))
  if (is.null(ejprob)){
    cat("  Calculating disjunction probability for original data ")
    cn <- con.regmat(x$prab,x$nb)
    ejumps <- sum(cn-1)
    den <- sum(x$regperspec-1)
    if (den==0 | ejumps==0)
      ejprob <- 0
    else
      ejprob <- ejumps/sum(x$regperspec-1)
    cat(ejprob,"\n")
  }
  out <- list()
  if (ejprob>0)
    out <-  autoreg(x, probs, ejprob, plot, nperp,
                    species.fixed=species.fixed, pdfnb=pdfnb)
  if (ejprob==0)
    out$pd <- 0
  if (out$pd<0)
    out$pd <- 0
  if (out$pd>1)
    out$pd <- 1
  if (twostep & out$pd>0){
    out1 <- out
    prange2 <- c(max(0,out$pd-5*step2),min(1,out$pd+5*step2))
    probs <- prange2[1]+step2*(0:round((prange2[2]-prange2[1])/step2))
    out <- autoreg(x, probs, ejprob, plot, nperp,
                     species.fixed=species.fixed, pdfnb=pdfnb)
    if (out$pd<0 | out$pd>1)
      out <- out1
  }
  out <- c(out,list(ejprob=ejprob))
  out
}











