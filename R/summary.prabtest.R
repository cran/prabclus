"summary.prabtest" <-
function(object, ...){
  rrange <- range(object$results)
  rmean <- mean(object$results)
  out <- list(rrange=rrange, rmean=rmean, datac=object$datac,
              p.value=object$p.value, pd=object$pd, tuning=object$tuning,
              teststat=object$teststat, distance=object$distance,
              times=object$times, pdfnb=object$pdfnb)
  class(out) <- "summary.prabtest"
  out
}
